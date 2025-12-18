# ===================================================================
# SAFE-only meta-analysis simulation (IV Only), K fixed
#   - K fixed per scenario; n1,n2 vary per study around n_mean (min 3)
#   - Per-study effect: lnM via SAFE; point = SAFE-BC (2*Delta1 - SAFE_mean)
#   - Per-study var:    SAFE variance only
#   - Estimator (metafor::rma.mv):
#       (A) Inverse-variance (IV):         V = diag(vi)
#   - δμ is calibrated so that true lnM (given τδ) ≈ lnM_target at mean n
#   - Coverage computed from these CIs
#   - Run with 16 cores:
#       Sys.setenv(N_WORKERS = "16"); source("fixedK_varyN_SAFE_IV_ONLY.R")
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
  library(tidyr)
  library(ggplot2)
  library(scales)
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

# --------------------- SAFE sampler (stand-alone) -------------------
# Accept-reject over summary stats for two-group independent design.
# Returns mean and variance of lnM across accepted draws (SAFE).
.safe_sampler_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                                min_kept = 2000, chunk_init = 4000,
                                chunk_max = 2e6, max_draws = Inf,
                                patience_noaccept = 5) {
  df1 <- n1 - 1L; df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)
  
  kept <- 0L
  lnM_star <- numeric(0)
  draws <- 0L
  chunk <- chunk_init
  noacc_streak <- 0L
  
  while (kept < min_kept && draws < max_draws && chunk <= chunk_max) {
    m1 <- rnorm(chunk, mean = x1bar, sd = s1/sqrt(n1))
    m2 <- rnorm(chunk, mean = x2bar, sd = s2/sqrt(n2))
    v1 <- s1^2 * stats::rchisq(chunk, df = df1) / df1
    v2 <- s2^2 * stats::rchisq(chunk, df = df2) / df2
    
    MSB <- h * (m1 - m2)^2
    MSW <- ((df1) * v1 + (df2) * v2) / (df1 + df2)
    good <- MSB > MSW
    
    if (any(good)) {
      noacc_streak <- 0L
      gidx <- which(good)
      lnM_star <- c(lnM_star,
                    0.5 * (log((MSB[gidx] - MSW[gidx])/(2*h)) - log(MSW[gidx])))
      kept <- length(lnM_star)
    } else {
      noacc_streak <- noacc_streak + 1L
      if (noacc_streak >= patience_noaccept) {
        chunk <- min(as.integer(chunk * 2), chunk_max)
        noacc_streak <- 0L
      }
    }
    draws <- draws + chunk
  }
  
  if (!length(lnM_star)) {
    return(list(point = NA_real_, var = NA_real_,
                kept = 0L, total = draws, status = "no_accept"))
  }
  list(point = mean(lnM_star), var = stats::var(lnM_star),
       kept = length(lnM_star), total = draws,
       status = if (kept >= min_kept) "ok" else "stopped_early")
}

# Wrapper used below
safe_call_indep <- .safe_sampler_indep

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

# --------------- Truth & calibration helpers ------------------------
true_lnM_indep <- function(mu1, mu2, sd1, sd2, n1, n2, tau_delta = 0) {
  # Large-sample expectation of lnM given mean difference (mu1 - mu2)
  # and between-study SD on the delta-scale (tau_delta).
  h           <- n1 * n2 / (n1 + n2)
  n0          <- 2 * h
  delta_mean  <- mu1 - mu2
  omega       <- ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  Delta_true  <- h * (delta_mean^2 + tau_delta^2)
  0.5 * (log(Delta_true / n0) - log(omega))
}

# Calibrate delta_mu so that true lnM (given tau_delta) hits lnM_target
delta_mu_from_target <- function(lnM_target, tau_delta, n1_bar, n2_bar, sd1=1, sd2=1,
                                 upper_init = 10, max_upper = 1e4) {
  f <- function(dmu) {
    true_lnM_indep(mu1 = 0, mu2 = dmu, sd1 = sd1, sd2 = sd2,
                   n1 = n1_bar, n2 = n2_bar, tau_delta = tau_delta) - lnM_target
  }
  
  # FIX: Start at a small epsilon instead of 0 to avoid log(0) when tau=0
  eps <- 1e-8
  f0 <- f(eps)
  
  if (!is.finite(f0)) stop("Non-finite baseline in delta_mu_from_target() even with epsilon.")
  if (f0 >= 0) return(eps)  # target ≤ attainable minimum
  
  up <- upper_init
  while (f(up) < 0 && up < max_upper) up <- up * 2
  if (f(up) < 0) stop("Target lnM too large; increase max_upper.")
  uniroot(f, c(eps, up))$root
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

# ----------------------- CI helper (robust) -------------------------
ci_95_from_fit <- function(fit) {
  b  <- as.numeric(fit$b[1])
  se <- as.numeric(sqrt(vcov(fit)[1,1]))
  df <- NA_real_
  if (!is.null(fit$ddf)) {
    df <- suppressWarnings(as.numeric(fit$ddf[1]))
  }
  crit <- if (is.finite(df)) qt(0.975, df) else qnorm(0.975)
  c(lb = b - crit*se, ub = b + crit*se)
}

# ------------------- One meta-analysis replicate --------------------
fit_meta_once_fixedK <- function(K_fixed,
                                 n1_mean, n2_mean,
                                 n_dist = c("poisson","nbinom"),
                                 n_nb_size = 10L,
                                 delta_mu, tau_delta,
                                 mu1 = 0, sd1 = 1, sd2 = 1,
                                 # SAFE tuning
                                 min_kept = 2000, chunk_init = 4000,
                                 chunk_max = 2e6, max_draws = Inf, patience_noaccept = 5,
                                 n_min = 3L) {
  n_dist <- match.arg(n_dist)
  
  # Per-study sample sizes (vary)
  n1_k <- draw_n_vec(K_fixed, n1_mean, dist = n_dist, nb_size = n_nb_size, n_min = n_min)
  n2_k <- draw_n_vec(K_fixed, n2_mean, dist = n_dist, nb_size = n_nb_size, n_min = n_min)
  
  # Study-level deltas (between-study SD = tau_delta) on delta-scale
  delta_k <- if (tau_delta > 0) rnorm(K_fixed, mean = delta_mu, sd = tau_delta) else rep(delta_mu, K_fixed)
  mu2_k   <- mu1 + delta_k
  
  yi  <- numeric(K_fixed)
  vi  <- numeric(K_fixed)
  n0k <- 2 * n1_k * n2_k / (n1_k + n2_k)
  
  for (k in seq_len(K_fixed)) {
    sm <- draw_summaries_indep(mu1, mu2_k[k], sd1, sd2, n1_k[k], n2_k[k])
    d1 <- lnM_delta1_indep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1_k[k], n2_k[k])
    SAFE <- safe_call_indep(
      sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1_k[k], n2_k[k],
      min_kept = min_kept, chunk_init = chunk_init,
      chunk_max = chunk_max, max_draws = max_draws,
      patience_noaccept = patience_noaccept
    )
    # SAFE-BC point; fall back to SAFE mean if needed
    yi[k] <- if (is.na(d1["point"]) || is.na(SAFE$point)) SAFE$point else 2*d1["point"] - SAFE$point
    vi[k] <- SAFE$var
  }
  
  good <- is.finite(yi) & is.finite(vi) & (vi > 0)
  if (!any(good) || sum(good) < 4) return(list(ok = FALSE, m_effects = sum(good)))
  
  yi  <- yi[good]
  vi  <- pmax(vi[good], 1e-12)
  m   <- length(yi)
  dat <- data.frame(yi = yi, vi = vi, ID = factor(seq_len(m)))
  
  # (A) IV (standard)
  fit_iv <- tryCatch(
    rma.mv(yi ~ 1, V = vi, random = ~ 1 | ID, data = dat,
           method = "REML", test = "t"),
    error = function(e) NULL
  )
  if (is.null(fit_iv)) return(list(ok = FALSE, m_effects = m))
  
  mu_iv   <- as.numeric(fit_iv$b[1])
  se_iv   <- as.numeric(sqrt(vcov(fit_iv)[1,1]))
  tau2_iv <- sum(fit_iv$sigma2)
  ci_iv   <- ci_95_from_fit(fit_iv)
  
  # Return ONLY IV results (Multiplicative removed)
  list(ok = TRUE,
       m_effects = m,
       # IV
       mu_iv = mu_iv, se_iv = se_iv, tau2_iv = tau2_iv,
       ci_iv_lb = ci_iv["lb"], ci_iv_ub = ci_iv["ub"])
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
                                parallel_replicates = TRUE,
                                n_min = 3L) {
  n_dist <- match.arg(n_dist)
  
  # Calibrate delta_mu so that true lnM (given tau_delta) matches lnM_target
  delta_mu <- delta_mu_from_target(lnM_target, tau_delta,
                                   n1_bar = n1_mean, n2_bar = n2_mean,
                                   sd1 = sd1, sd2 = sd2)
  
  # Analytic "truth" at mean n (should ≈ lnM_target)
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
      bias_iv = NA_real_,
      rmse_iv = NA_real_,
      cover_iv = NA_real_,
      tau2_iv_mean = NA_real_,
      mean_m_effects = NA_real_,
      invalid_meta_rate = 1
    ))
  }
  
  ext <- function(name) vapply(out_list[ok], `[[`, numeric(1), name)
  
  mu_iv     <- ext("mu_iv")
  tau_iv    <- ext("tau2_iv")
  m_eff     <- ext("m_effects")
  
  ci_iv_lb   <- ext("ci_iv_lb");   ci_iv_ub   <- ext("ci_iv_ub")
  
  cover_iv_vec   <- (ci_iv_lb   <= lnM_true0) & (lnM_true0 <= ci_iv_ub)
  
  tibble(
    K_fixed = K_fixed, n1_mean = n1_mean, n2_mean = n2_mean,
    lnM_target = lnM_target, tau_delta = tau_delta,
    lnM_true0 = lnM_true0,
    bias_iv    = mean(mu_iv    - lnM_true0),
    rmse_iv    = sqrt(mean((mu_iv    - lnM_true0)^2)),
    cover_iv   = mean(cover_iv_vec,    na.rm = TRUE),
    tau2_iv_mean     = mean(tau_iv),
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
    chunk_max = 2e6, max_draws = Inf, patience_noaccept = 5
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
  message("=== Quick demo (K fixed; n varies; IV Only; robust CIs; calibrated δμ) ===")
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
    n_dist = "nbinom", n_nb_size = 8L,
    n_min = 3L,
    parallel_replicates = TRUE,
    parallel_scenarios = FALSE,
    min_kept = 1000, chunk_init = 3000
  )
  
  message(sprintf("[end]    %s", format(Sys.time(), "%H:%M:%S")))
  print(demo)
}

# ========================= Plotting section ==========================
res <- demo %>%
  rename(K     = K_fixed,
         nbar  = n1_mean,
         lnM_x = lnM_target,
         tau   = tau_delta) %>%
  mutate(
    setting = paste0("K=", K, ", n=", nbar),
    tau_f   = factor(tau, levels = sort(unique(tau)),
                     labels = paste0("tau = ", sort(unique(tau))))
  )

# Bias: IV Only
p_bias <- ggplot(res, aes(x = lnM_x, y = bias_iv)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(color = "blue") + geom_point(size = 1.8, color = "blue") +
  facet_grid(tau_f ~ setting, scales = "free_x") +
  labs(x = "lnM (target ≈ truth)", y = "Bias",
       title = "Bias: IV Estimator") +
  theme_bw(base_size = 12)

# RMSE: IV Only
p_rmse <- ggplot(res, aes(x = lnM_x, y = rmse_iv)) +
  geom_line(color = "red") + geom_point(size = 1.8, color = "red") +
  facet_grid(tau_f ~ setting, scales = "free_x") +
  labs(x = "lnM (target ≈ truth)", y = "RMSE",
       title = "RMSE: IV Estimator") +
  theme_bw(base_size = 12)

# Coverage: IV Only
if ("cover_iv" %in% names(res)) {
  p_cover <- ggplot(res, aes(x = lnM_x, y = cover_iv)) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_line(color = "darkgreen") + geom_point(size = 1.8, color = "darkgreen") +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    facet_grid(tau_f ~ setting, scales = "free_x") +
    labs(x = "lnM (target ≈ truth)", y = "Coverage (95% nominal)",
         title = "Interval coverage: IV Estimator") +
    theme_bw(base_size = 12)
}

# Estimated tau^2: IV Only
if ("tau2_iv_mean" %in% names(res)) {
  p_tau <- ggplot(res, aes(x = lnM_x, y = tau2_iv_mean)) +
    geom_line(color = "purple") + geom_point(size = 1.8, color = "purple") +
    facet_grid(tau_f ~ setting, scales = "free_x") +
    labs(x = "lnM (target ≈ truth)", y = expression(mean(hat(tau)^2)~"(REML)"),
         title = expression(paste("Between-study variance ", tau^2, ": estimated (IV)"))) +
    theme_bw(base_size = 12)
}

# Invalid replicate rate heatmap
p_invalid <- ggplot(res, aes(x = lnM_x, y = setting, fill = invalid_meta_rate)) +
  geom_tile(color = "white") +
  facet_wrap(~ tau_f, nrow = 1) +
  scale_fill_gradient(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(x = "lnM (target ≈ truth)", y = "Scenario", fill = "Invalid rate",
       title = "Proportion of failed meta-analysis replicates") +
  theme_bw(base_size = 12) + theme(legend.position = "right")

# Print
p_bias
p_rmse
if (exists("p_cover")) p_cover
if (exists("p_tau"))   p_tau
p_invalid