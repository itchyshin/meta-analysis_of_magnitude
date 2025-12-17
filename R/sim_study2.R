# ===================================================================
# SAFE-only meta-analysis simulation (IV vs Multiplicative), K fixed
#   - K fixed per scenario; n1,n2 vary per study around n_mean (min 3)
#   - Per-study effect: lnM via SAFE; point = SAFE-BC (2*Delta1 - SAFE_mean)
#   - Per-study var:    SAFE variance only
#   - Estimators (rma.mv):
#       (A) Inverse-variance (IV): V = diag(vi)
#       (B) Multiplicative (n0-based): V = 0, R = diag(vtilde), vtilde ∝ 1/n0
#   - Coverage computed from confint() CIs (not 1.96*SE)
#   - Grid over K_fixed_vec, n_mean_vec, lnM_targets, tau_vec
#   - Run with 16 cores:
#       Sys.setenv(N_WORKERS = "16"); source("fixedK_varyN_SAFE.R")
# ===================================================================

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(metafor)
  library(MASS)
  library(tibble)
  library(dplyr)
  library(purrr)
  library(furrr)
  library(future)
  library(progressr)
})

# ------------------------- Progress handlers -------------------------
progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

# ------------------------- RNG + helpers -----------------------------
SEED_MAIN <- 202501
set.seed(SEED_MAIN)

safe_gap <- function(gap) ifelse(gap <= 0, NA_real_, gap)

# --------------------- Parallel plan helper -------------------------
set_parallel_plan <- function(n_workers = NULL, backend = c("multisession","multicore")) {
  backend <- match.arg(backend)
  if (is.null(n_workers)) n_workers <- as.integer(Sys.getenv("N_WORKERS", parallel::detectCores()))
  n_workers <- max(1L, n_workers)
  if (backend == "multicore" && .Platform$OS.type != "windows" && !identical(Sys.getenv("RSTUDIO"), "1")) {
    future::plan(future::multicore, workers = n_workers)
  } else {
    future::plan(future::multisession, workers = n_workers)
  }
  message(sprintf("[plan] backend=%s; workers=%d", backend, n_workers))
  invisible(n_workers)
}

.FOPTS <- furrr_options(
  seed = TRUE,
  packages = c("metafor","MASS","tibble","dplyr","purrr"),
  scheduling = 2
)

# --------------------- SAFE front-door (+ fallback) ------------------
SAFE_FILE <- "SAFE_fun.R"
SAFE_HAS_NEW <- FALSE
if (file.exists(SAFE_FILE)) {
  source(SAFE_FILE)
  SAFE_HAS_NEW <- exists("safe_lnM_indep")
  message(sprintf("[SAFE] Using %s API", if (SAFE_HAS_NEW) "NEW" else "fallback"))
} else {
  message("[SAFE] SAFE_fun.R not found; using fallback SAFE")
}

.safe_fallback_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                                 B = 2000, chunk = 4000, max_chunks = 50) {
  df1 <- n1 - 1L; df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)
  lnM_star <- numeric(0)
  kept <- total <- attempts <- 0L
  while (length(lnM_star) < B && attempts < max_chunks) {
    attempts <- attempts + 1L
    m1 <- rnorm(chunk, mean = x1bar, sd = s1/sqrt(n1))
    m2 <- rnorm(chunk, mean = x2bar, sd = s2/sqrt(n2))
    v1 <- s1^2 * stats::rchisq(chunk, df = df1) / df1
    v2 <- s2^2 * stats::rchisq(chunk, df = df2) / df2
    total <- total + chunk
    MSB   <- h * (m1 - m2)^2
    MSW   <- ((df1) * v1 + (df2) * v2) / (df1 + df2)
    good  <- MSB > MSW
    if (any(good)) {
      kept <- kept + sum(good)
      lnM_star <- c(lnM_star,
                    0.5 * (log((MSB[good] - MSW[good])/(2*h)) - log(MSW[good])))
    }
    if (attempts %% 5L == 0L) {
      message(sprintf("[SAFE-fallback] attempts=%d kept=%d total=%d", attempts, kept, total))
    }
  }
  lnM_star <- if (length(lnM_star)) lnM_star[seq_len(min(B, length(lnM_star)))] else lnM_star
  list(
    point    = if (length(lnM_star)) mean(lnM_star) else NA_real_,
    var      = if (length(lnM_star)) stats::var(lnM_star) else NA_real_,
    kept     = kept,
    total    = total,
    attempts = attempts,
    status   = if (length(lnM_star) >= B) "ok" else "stopped_early"
  )
}

safe_call_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                            min_kept = 2000, chunk_init = 4000,
                            chunk_max = 2e6, max_draws = Inf,
                            patience_noaccept = 5,
                            B_fallback = 2000, chunk_fallback = 4000,
                            max_chunks_fallback = 50) {
  if (SAFE_HAS_NEW) {
    res <- tryCatch(
      safe_lnM_indep(x1bar, x2bar, s1, s2, n1, n2,
                     min_kept = min_kept, chunk_init = chunk_init,
                     chunk_max = chunk_max, max_draws = max_draws,
                     patience_noaccept = patience_noaccept),
      error = function(e) NULL
    )
    if (!is.null(res)) return(res)
  }
  .safe_fallback_indep(x1bar, x2bar, s1, s2, n1, n2,
                       B = B_fallback, chunk = chunk_fallback,
                       max_chunks = max_chunks_fallback)
}

# ----------------- Delta-1 (for SAFE-BC point) ----------------------
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h    <- n1 * n2 / (n1 + n2)
  MSB  <- h * (x1bar - x2bar)^2
  MSW  <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA_real_))
  lnM <- 0.5 * (log(Delta / (2*h)) - log(MSW))
  c(point = lnM)
}

# -------------------- Summary draws (independent) --------------------
draw_summaries_indep <- function(mu1, mu2, sd1, sd2, n1, n2) {
  x1bar <- rnorm(1L, mu1, sd1 / sqrt(n1))
  x2bar <- rnorm(1L, mu2, sd2 / sqrt(n2))
  s1sq  <- sd1^2 * stats::rchisq(1L, df = n1 - 1) / (n1 - 1)
  s2sq  <- sd2^2 * stats::rchisq(1L, df = n2 - 1) / (n2 - 1)
  c(x1bar = x1bar, x2bar = x2bar, s1 = sqrt(s1sq), s2 = sqrt(s2sq))
}

# --------------- True lnM + inverse mapping helpers -----------------
delta_from_target_lnM <- function(lnM_target, n1_bar, n2_bar, sdW = 1) {
  # deq^2 = 2*exp(2*lnM) + 2/n0, with deq = |delta|/sW ; n0 = 2*h
  n0_bar <- 2 * n1_bar * n2_bar / (n1_bar + n2_bar)
  deq2   <- 2 * exp(2 * lnM_target) + 2 / n0_bar
  sdW * sqrt(pmax(deq2, 0))
}

true_lnM_indep <- function(mu1, mu2, sd1, sd2, n1, n2, tau_delta = 0) {
  h          <- n1 * n2 / (n1 + n2)
  n0         <- 2 * h
  delta_mean <- mu1 - mu2
  omega      <- ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  Delta_true <- h * (delta_mean^2 + tau_delta^2)
  0.5 * (log(Delta_true / n0) - log(omega))
}

# --------------------- n drawing utilities --------------------------
draw_n_vec <- function(K, n_mean, dist = c("poisson","nbinom"), nb_size = 10L, n_min = 3L) {
  dist <- match.arg(dist)
  if (dist == "poisson") {
    n <- rpois(K, lambda = n_mean)
  } else {
    n <- rnbinom(K, size = nb_size, mu = n_mean)
  }
  pmax(n_min, as.integer(n))
}

# ------------------- One meta-analysis replicate (with CIs) ---------
fit_meta_once_fixedK <- function(K_fixed,
                                 n1_mean, n2_mean,
                                 n_dist = c("poisson","nbinom"),
                                 n_nb_size = 10L,
                                 delta_mu, tau_delta,
                                 mu1 = 0, sd1 = 1, sd2 = 1,
                                 # SAFE tuning
                                 min_kept = 2000, chunk_init = 4000,
                                 chunk_max = 2e6, max_draws = Inf, patience_noaccept = 5,
                                 B_fallback = 2000, chunk_fallback = 4000, max_chunks_fallback = 50,
                                 scale_mult_to_vi_mean = FALSE,
                                 n_min = 3L) {
  n_dist <- match.arg(n_dist)
  
  # Per-study sample sizes (vary), floor at n_min
  n1_k <- draw_n_vec(K_fixed, n1_mean, dist = n_dist, nb_size = n_nb_size, n_min = n_min)
  n2_k <- draw_n_vec(K_fixed, n2_mean, dist = n_dist, nb_size = n_nb_size, n_min = n_min)
  
  # Study-level deltas
  delta_k <- if (tau_delta > 0) rnorm(K_fixed, mean = delta_mu, sd = tau_delta) else rep(delta_mu, K_fixed)
  mu2_k   <- mu1 + delta_k
  
  yi  <- numeric(K_fixed)
  vi  <- numeric(K_fixed)
  n0k <- numeric(K_fixed)
  
  for (k in seq_len(K_fixed)) {
    sm <- draw_summaries_indep(mu1, mu2_k[k], sd1, sd2, n1_k[k], n2_k[k])
    d1 <- lnM_delta1_indep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1_k[k], n2_k[k])
    
    SAFE <- safe_call_indep(
      sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1_k[k], n2_k[k],
      min_kept = min_kept, chunk_init = chunk_init, chunk_max = chunk_max,
      max_draws = max_draws, patience_noaccept = patience_noaccept,
      B_fallback = B_fallback, chunk_fallback = chunk_fallback, max_chunks_fallback = max_chunks_fallback
    )
    
    # SAFE-BC point; fall back to SAFE mean if needed
    yi[k]  <- if (is.na(d1["point"]) || is.na(SAFE$point)) SAFE$point else 2*d1["point"] - SAFE$point
    vi[k]  <- SAFE$var
    n0k[k] <- 2 * n1_k[k] * n2_k[k] / (n1_k[k] + n2_k[k])
  }
  
  good <- is.finite(yi) & is.finite(vi) & (vi > 0)
  if (!any(good) || sum(good) < 4) return(list(ok = FALSE, m_effects = sum(good)))
  
  yi <- yi[good]; vi <- pmax(vi[good], 1e-10); n0k <- n0k[good]
  m  <- length(yi)
  dat <- data.frame(yi = yi, vi = vi, n0 = n0k, ID = factor(seq_len(m)))
  
  # (A) IV
  fit_iv <- tryCatch(
    rma.mv(yi ~ 1, V = vi, random = ~ 1 | ID, data = dat,
           method = "REML", test = "t"),
    error = function(e) NULL
  )
  if (is.null(fit_iv)) return(list(ok = FALSE, m_effects = m))
  
  mu_iv   <- as.numeric(fit_iv$b[1])
  se_iv   <- sqrt(vcov(fit_iv))[1,1]
  tau2_iv <- sum(fit_iv$sigma2)
  
  ci_iv_lb <- ci_iv_ub <- NA_real_
  iv_ci <- tryCatch(confint(fit_iv, level = 0.95), error = function(e) NULL)
  if (!is.null(iv_ci)) {
    mat <- if (!is.null(iv_ci$beta)) iv_ci$beta else iv_ci
    if (!is.null(dim(mat))) {
      ci_iv_lb <- suppressWarnings(as.numeric(mat[1, "ci.lb"]))
      ci_iv_ub <- suppressWarnings(as.numeric(mat[1, "ci.ub"]))
    }
  }
  
  # (B) Multiplicative (weights ∝ n0 -> vtilde ∝ 1/n0)
  vtilde <- 1 / (dat$n0/2)
  if (isTRUE(scale_mult_to_vi_mean)) vtilde <- vtilde * (mean(dat$vi) / mean(vtilde))
  Vf <- diag(as.numeric(vtilde)); levs <- levels(dat$ID); rownames(Vf) <- levs; colnames(Vf) <- levs
  
  fit_mult <- tryCatch(
    rma.mv(yi ~ 1, V = 0, random = ~ 1 | ID, data = dat,
           R = list(ID = Vf), Rscale = FALSE,
           method = "REML", control = list(stepadj = 0.5, maxiter = 1000),
           test = "t"),
    error = function(e) NULL
  )
  if (is.null(fit_mult)) return(list(ok = FALSE, m_effects = m))
  
  mu_mult   <- as.numeric(fit_mult$b[1])
  se_mult   <- sqrt(vcov(fit_mult))[1,1]
  tau2_mult <- sum(fit_mult$sigma2)
  
  ci_mult_lb <- ci_mult_ub <- NA_real_
  mult_ci <- tryCatch(confint(fit_mult, level = 0.95), error = function(e) NULL)
  if (!is.null(mult_ci)) {
    mat <- if (!is.null(mult_ci$beta)) mult_ci$beta else mult_ci
    if (!is.null(dim(mat))) {
      ci_mult_lb <- suppressWarnings(as.numeric(mat[1, "ci.lb"]))
      ci_mult_ub <- suppressWarnings(as.numeric(mat[1, "ci.ub"]))
    }
  }
  
  list(ok = TRUE,
       m_effects = m,
       # IV
       mu_iv = mu_iv, se_iv = se_iv, tau2_iv = tau2_iv,
       ci_iv_lb = ci_iv_lb, ci_iv_ub = ci_iv_ub,
       # Multiplicative
       mu_mult = mu_mult, se_mult = se_mult, tau2_mult = tau2_mult,
       ci_mult_lb = ci_mult_lb, ci_mult_ub = ci_mult_ub)
}

# -------------------- One scenario driver (CI-based coverage) --------
run_scenario_fixedK <- function(R_meta,
                                K_fixed,
                                n1_mean, n2_mean,
                                n_dist = c("poisson","nbinom"),
                                n_nb_size = 10L,
                                lnM_target, tau_delta,
                                mu1 = 0, sd1 = 1, sd2 = 1,
                                min_kept = 2000, chunk_init = 4000,
                                chunk_max = 2e6, max_draws = Inf, patience_noaccept = 5,
                                B_fallback = 2000, chunk_fallback = 4000, max_chunks_fallback = 50,
                                scale_mult_to_vi_mean = FALSE,
                                parallel_replicates = TRUE,
                                n_min = 3L) {
  n_dist <- match.arg(n_dist)
  
  # Back-solve δ to hit lnM target using mean n
  delta_mu <- delta_from_target_lnM(lnM_target, n1_bar = n1_mean, n2_bar = n2_mean, sdW = 1)
  
  # "Truth" for summaries (large-sample analytic at mean n)
  lnM_true0 <- true_lnM_indep(mu1, mu1 + delta_mu, sd1, sd2, n1_mean, n2_mean, tau_delta = tau_delta)
  
  one_rep <- function(r) {
    fit_meta_once_fixedK(
      K_fixed = K_fixed,
      n1_mean = n1_mean, n2_mean = n2_mean,
      n_dist = n_dist, n_nb_size = n_nb_size,
      delta_mu = delta_mu, tau_delta = tau_delta,
      mu1 = mu1, sd1 = sd1, sd2 = sd2,
      min_kept = min_kept, chunk_init = chunk_init,
      chunk_max = chunk_max, max_draws = max_draws, patience_noaccept = patience_noaccept,
      B_fallback = B_fallback, chunk_fallback = chunk_fallback, max_chunks_fallback = max_chunks_fallback,
      scale_mult_to_vi_mean = scale_mult_to_vi_mean,
      n_min = n_min
    )
  }
  
  reps <- seq_len(R_meta)
  out_list <- if (parallel_replicates) {
    progressr::with_progress({
      future_map(reps, one_rep, .options = .FOPTS, .progress = TRUE)
    })
  } else {
    pb <- txtProgressBar(min = 0, max = length(reps), style = 3)
    on.exit(close(pb), add = TRUE)
    map(reps, ~{res <- one_rep(.x); setTxtProgressBar(pb, .x); res})
  }
  
  ok <- vapply(out_list, function(z) isTRUE(z$ok), logical(1))
  invalid_rate <- 1 - mean(ok)
  
  if (!any(ok)) {
    return(tibble(
      K_fixed = K_fixed, n1_mean = n1_mean, n2_mean = n2_mean,
      lnM_target = lnM_target, tau_delta = tau_delta,
      lnM_true0 = lnM_true0,
      bias_iv = NA_real_, bias_mult = NA_real_,
      rmse_iv = NA_real_, rmse_mult = NA_real_,
      cover_iv = NA_real_, cover_mult = NA_real_,
      tau2_iv_mean = NA_real_, tau2_mult_mean = NA_real_,
      mean_m_effects = NA_real_,
      invalid_meta_rate = 1
    ))
  }
  
  ext <- function(name) vapply(out_list[ok], `[[`, numeric(1), name)
  
  mu_iv     <- ext("mu_iv")
  mu_mult   <- ext("mu_mult")
  tau_iv    <- ext("tau2_iv")
  tau_mult  <- ext("tau2_mult")
  m_eff     <- ext("m_effects")
  
  ci_iv_lb   <- ext("ci_iv_lb");   ci_iv_ub   <- ext("ci_iv_ub")
  ci_mult_lb <- ext("ci_mult_lb"); ci_mult_ub <- ext("ci_mult_ub")
  
  cover_iv_vec   <- (ci_iv_lb   <= lnM_true0) & (lnM_true0 <= ci_iv_ub)
  cover_mult_vec <- (ci_mult_lb <= lnM_true0) & (lnM_true0 <= ci_mult_ub)
  
  tibble(
    K_fixed = K_fixed, n1_mean = n1_mean, n2_mean = n2_mean,
    lnM_target = lnM_target, tau_delta = tau_delta,
    lnM_true0 = lnM_true0,
    bias_iv   = mean(mu_iv   - lnM_true0),
    bias_mult = mean(mu_mult - lnM_true0),
    rmse_iv   = sqrt(mean((mu_iv   - lnM_true0)^2)),
    rmse_mult = sqrt(mean((mu_mult - lnM_true0)^2)),
    cover_iv   = mean(cover_iv_vec,   na.rm = TRUE),
    cover_mult = mean(cover_mult_vec, na.rm = TRUE),
    tau2_iv_mean   = mean(tau_iv),
    tau2_mult_mean = mean(tau_mult),
    mean_m_effects = mean(m_eff),
    invalid_meta_rate = invalid_rate
  )
}

# ---------------------- Scenario grid + runner -----------------------
seq_progress_map <- function(idx, fun) {
  pb <- txtProgressBar(min = 0, max = length(idx), style = 3)
  on.exit(close(pb), add = TRUE)
  res <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    res[[i]] <- fun(idx[i])
    setTxtProgressBar(pb, i)
  }
  res
}

run_grid_fixedK <- function(
    n_workers = NULL,                      # plan set outside; cosmetic only
    parallel_replicates = TRUE,
    parallel_scenarios  = FALSE,
    # Grid
    R_meta        = 400,
    K_fixed_vec   = c(20, 50, 100),
    n_mean_vec    = c(10, 20, 30, 50),
    lnM_targets   = c(-2.0, -1.5, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0),
    tau_vec       = c(0.2, 0.4),
    # n distribution & overdispersion
    n_dist        = c("poisson","nbinom"),
    n_nb_size     = 10L,
    n_min         = 3L,
    # SAFE tuning
    min_kept = 2000, chunk_init = 4000,
    chunk_max = 2e6, max_draws = Inf, patience_noaccept = 5,
    B_fallback = 2000, chunk_fallback = 4000, max_chunks_fallback = 50,
    scale_mult_to_vi_mean = FALSE
) {
  n_dist <- match.arg(n_dist)
  
  grid <- expand.grid(
    K_fixed = K_fixed_vec,
    n_mean  = n_mean_vec,
    lnM     = lnM_targets,
    tau     = tau_vec,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid <- tibble::as_tibble(grid)
  message(sprintf("[grid] Scenarios: %d (K={%s}; n̄={%s}; lnM=%d; tau={%s}; ndist=%s)",
                  nrow(grid),
                  paste(K_fixed_vec, collapse=","),
                  paste(n_mean_vec, collapse=","),
                  length(lnM_targets),
                  paste(tau_vec, collapse=","),
                  n_dist))
  
  run_one <- function(irow) {
    par <- grid[irow,]
    message(sprintf("[scenario %d/%d] K=%d; n̄=%d; lnM=%.3f; tau=%.2f  @%s",
                    irow, nrow(grid), par$K_fixed, par$n_mean, par$lnM, par$tau, format(Sys.time(), "%H:%M:%S")))
    res <- run_scenario_fixedK(
      R_meta = R_meta,
      K_fixed = par$K_fixed,
      n1_mean = par$n_mean, n2_mean = par$n_mean,
      n_dist = n_dist, n_nb_size = n_nb_size,
      lnM_target = par$lnM, tau_delta = par$tau,
      min_kept = min_kept, chunk_init = chunk_init,
      chunk_max = chunk_max, max_draws = max_draws, patience_noaccept = patience_noaccept,
      B_fallback = B_fallback, chunk_fallback = chunk_fallback, max_chunks_fallback = max_chunks_fallback,
      scale_mult_to_vi_mean = scale_mult_to_vi_mean,
      parallel_replicates = parallel_replicates,
      n_min = n_min
    ) %>% mutate(scenario_id = irow)
    message(sprintf("[scenario %d] done @%s", irow, format(Sys.time(), "%H:%M:%S")))
    res
  }
  
  idx <- seq_len(nrow(grid))
  res_list <- if (parallel_scenarios) {
    progressr::with_progress({
      future_map(idx, run_one, .options = .FOPTS, .progress = TRUE)
    })
  } else {
    seq_progress_map(idx, run_one)
  }
  
  bind_rows(res_list) %>%
    relocate(scenario_id) %>%
    arrange(scenario_id, K_fixed, n1_mean, lnM_target, tau_delta)
}

# --------------------------- Demo run --------------------------------
if (sys.nframe() == 0) {
  message("=== Quick demo (K fixed; n varies; progress enabled) ===")
  Sys.setenv(N_WORKERS = Sys.getenv("N_WORKERS", "16"))
  set_parallel_plan(n_workers = as.integer(Sys.getenv("N_WORKERS")), backend = "multisession")
  message(sprintf("[start] %s", format(Sys.time(), "%H:%M:%S")))
  
  demo <- run_grid_fixedK(
    n_workers = as.integer(Sys.getenv("N_WORKERS")),
    R_meta = 200,
    K_fixed_vec = c(20, 50, 100),
    n_mean_vec = c(10, 20, 30, 50),
    lnM_targets = c(-1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0),
    tau_vec = c(0.2, 0.4),
    n_dist = "nbinom", n_nb_size = 8L,  # overdispersed n; use "poisson" if preferred
    n_min = 3L,
    parallel_replicates = TRUE,
    parallel_scenarios = FALSE,
    min_kept = 1000, chunk_init = 3000,
    B_fallback = 1000, chunk_fallback = 3000
  )
  
  message(sprintf("[end]   %s", format(Sys.time(), "%H:%M:%S")))
  print(demo)
}

# ========================= Plotting section ==========================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

# Map columns for plotting
res <- demo %>%
  rename(K     = K_fixed,
         nbar  = n1_mean,
         lnM_x = lnM_true0,
         tau   = tau_delta) %>%
  mutate(
    setting = paste0("K=", K, ", n=", nbar),
    tau_f   = factor(tau, levels = sort(unique(tau)),
                     labels = paste0("tau = ", sort(unique(tau))))
  )

# Bias
df_bias <- res %>%
  select(setting, tau_f, lnM_x, bias_iv, bias_mult) %>%
  pivot_longer(starts_with("bias_"), names_to = "method", values_to = "bias") %>%
  mutate(method = recode(method, bias_iv = "IV", bias_mult = "Multiplicative"))

p_bias <- ggplot(df_bias, aes(x = lnM_x, y = bias, color = method, group = method)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + geom_point(size = 1.8) +
  facet_grid(tau_f ~ setting, scales = "free_x") +
  labs(x = "lnM (true)", y = "Bias", color = "Estimator",
       title = "Bias: IV vs Multiplicative") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

# RMSE
df_rmse <- res %>%
  select(setting, tau_f, lnM_x, rmse_iv, rmse_mult) %>%
  pivot_longer(starts_with("rmse_"), names_to = "method", values_to = "rmse") %>%
  mutate(method = recode(method, rmse_iv = "IV", rmse_mult = "Multiplicative"))

p_rmse <- ggplot(df_rmse, aes(x = lnM_x, y = rmse, color = method, group = method)) +
  geom_line() + geom_point(size = 1.8) +
  facet_grid(tau_f ~ setting, scales = "free_x") +
  labs(x = "lnM (true)", y = "RMSE", color = "Estimator",
       title = "RMSE: IV vs Multiplicative") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

# Coverage from CIs
if (all(c("cover_iv","cover_mult") %in% names(res))) {
  df_cov <- res %>%
    select(setting, tau_f, lnM_x, cover_iv, cover_mult) %>%
    pivot_longer(starts_with("cover_"), names_to = "method", values_to = "coverage") %>%
    mutate(method = recode(method, cover_iv = "IV", cover_mult = "Multiplicative"))
  
  p_cover <- ggplot(df_cov, aes(x = lnM_x, y = coverage, color = method, group = method)) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_line() + geom_point(size = 1.8) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    facet_grid(tau_f ~ setting, scales = "free_x") +
    labs(x = "lnM (true)", y = "Coverage (95% nominal)", color = "Estimator",
         title = "Interval coverage: IV vs Multiplicative (CI-based)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
}

# Estimated tau^2 vs truth
if (all(c("tau2_iv_mean","tau2_mult_mean") %in% names(res))) {
  df_tau <- res %>%
    select(setting, tau_f, tau, lnM_x, tau2_iv_mean, tau2_mult_mean) %>%
    pivot_longer(ends_with("_mean"), names_to = "method", values_to = "tau2_est") %>%
    mutate(method = recode(method, "tau2_iv_mean"="IV","tau2_mult_mean"="Multiplicative"))
  
  tau_lines <- res %>% distinct(tau_f, tau)
  
  p_tau <- ggplot(df_tau, aes(x = lnM_x, y = tau2_est, color = method, group = method)) +
    geom_line() + geom_point(size = 1.8) +
    geom_hline(data = tau_lines, aes(yintercept = tau), linetype = 3,
               inherit.aes = FALSE, colour = "grey40") +
    facet_grid(tau_f ~ setting, scales = "free_x") +
    labs(x = "lnM (true)", y = expression(hat(tau)^2~"(mean across reps)"),
         color = "Estimator",
         title = expression(paste("Between-study variance ", tau^2, ": estimated vs truth"))) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
}

# Invalid replicate rate
p_invalid <- ggplot(res, aes(x = lnM_x, y = setting, fill = invalid_meta_rate)) +
  geom_tile(color = "white") +
  facet_wrap(~ tau_f, nrow = 1) +
  scale_fill_gradient(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(x = "lnM (true)", y = "Scenario", fill = "Invalid rate",
       title = "Proportion of failed meta-analysis replicates") +
  theme_bw(base_size = 12) + theme(legend.position = "right")

# RMSE ratio (Mult / IV)
df_ratio <- res %>%
  transmute(setting, tau_f, x = lnM_x,
            rmse_ratio = rmse_mult / pmax(rmse_iv, .Machine$double.eps))

p_ratio <- ggplot(df_ratio, aes(x = x, y = rmse_ratio, group = 1)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line() + geom_point(size = 1.8) +
  facet_grid(tau_f ~ setting, scales = "free_x") +
  scale_y_log10() +
  labs(x = "lnM (true)", y = "RMSE ratio (Mult / IV, log scale)",
       title = "Relative performance (<1 favours Multiplicative)") +
  theme_bw(base_size = 12)

# Print
p_bias
p_rmse
if (exists("p_cover")) p_cover
if (exists("p_tau"))   p_tau
p_invalid
p_ratio

# # Optional saves:
# ggsave("bias_iv_vs_mult.pdf", p_bias, width = 10, height = 6)
# ggsave("rmse_iv_vs_mult.pdf", p_rmse, width = 10, height = 6)
# if (exists("p_cover")) ggsave("coverage_iv_vs_mult.pdf", p_cover, width = 10, height = 6)
# if (exists("p_tau"))   ggsave("tau2_iv_vs_mult.pdf", p_tau, width = 10, height = 6)
# ggsave("invalid_rate.pdf", p_invalid, width = 10, height = 4.5)
# ggsave("rmse_ratio_mult_over_iv.pdf", p_ratio, width = 10, height = 6)