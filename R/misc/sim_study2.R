# ===================================================================
# SAFE-only meta-analysis simulation (Inverse-Variance only), K fixed
#   • K fixed per scenario; n1, n2 vary per study around n_mean (min 3)
#   • Per-study effect: lnM via SAFE; point = SAFE-BC (2*Delta1 − SAFE_mean)
#   • Per-study var:    SAFE variance only
#   • Estimator (metafor::rma.mv): IV with V = diag(vi), random = ~1|ID
#   • δμ calibrated to hit lnM_target GIVEN τδ (at mean n)
#   • 95% CIs via predict(fit) (fallback: confint() → ±1.96*SE)
#   • Coverage uses those CIs
#   • MC uncertainty: MCSE for all reported metrics + 95% MC intervals
#   • Reduced grid demo: K ∈ {20, 50}; lnM_target ∈ {−1.0, −0.7, −0.4, 0.0, 0.4};
#                         τδ ∈ {0, 0.2}; n̄ ∈ {10, 20, 30, 50, 100}; n-dist = nbinom(size=8)
#   • Run:
#       Sys.setenv(N_WORKERS = "16"); source("SAFE_IV_only_MC.R")
# ===================================================================

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(metafor)
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
  packages = c("metafor","tibble","dplyr","purrr"),
  scheduling = 2
)

# --------------------- SAFE sampler (stand-alone) -------------------
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
  h          <- n1 * n2 / (n1 + n2)
  n0         <- 2 * h
  delta_mean <- mu1 - mu2
  omega      <- ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  Delta_true <- h * (delta_mean^2 + tau_delta^2)
  0.5 * (log(Delta_true / n0) - log(omega))
}

# ---- Calibrate delta_mu to hit lnM_target given tau_delta (mean n) --
# • Closed form when tau_delta == 0 (avoids log(0) baseline)
# • Robust root-finding for tau_delta > 0 without evaluating at 0
delta_mu_from_target <- function(lnM_target, tau_delta,
                                 n1_bar, n2_bar, sd1 = 1, sd2 = 1,
                                 upper_init = 10, max_upper = 1e4) {
  # pooled within-group variance at the mean n
  omega <- ((n1_bar - 1) * sd1^2 + (n2_bar - 1) * sd2^2) / (n1_bar + n2_bar - 2)
  # harmonic sample size factor
  h <- n1_bar * n2_bar / (n1_bar + n2_bar)
  n0 <- 2 * h
  
  if (isTRUE(all.equal(tau_delta, 0))) {
    # Solve lnM = 0.5*(log((h*d^2)/n0) - log(omega))  for d
    # ⇒ d = sqrt(2*omega) * exp(lnM)
    return(sqrt(2 * omega) * exp(lnM_target))
  }
  
  f <- function(dmu) 0.5 * (log((h * (dmu^2 + tau_delta^2)) / n0) - log(omega)) - lnM_target
  
  low <- 1e-12
  f_low <- f(low)
  if (!is.finite(f_low)) stop("Calibration failed: non-finite objective at lower bound.")
  
  if (f_low >= 0) return(0)  # target ≤ minimum achievable ⇒ dmu = 0
  
  up <- upper_init
  while (f(up) < 0 && up < max_upper) up <- up * 2
  if (f(up) < 0) stop("Target lnM too large at given tau_delta; increase max_upper.")
  
  uniroot(f, c(low, up))$root
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

# ------------------- One meta-analysis replicate (IV) ----------------
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
  
  for (k in seq_len(K_fixed)) {
    sm <- draw_summaries_indep(mu1, mu2_k[k], sd1, sd2, n1_k[k], n2_k[k])
    d1 <- lnM_delta1_indep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1_k[k], n2_k[k])
    SAFE <- safe_call_indep(
      sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1_k[k], n2_k[k],
      min_kept = min_kept, chunk_init = chunk_init,
      chunk_max = chunk_max, max_draws = max_draws,
      patience_noaccept = patience_noaccept
    )
    yi[k] <- if (is.na(d1["point"]) || is.na(SAFE$point)) SAFE$point else 2*d1["point"] - SAFE$point
    vi[k] <- SAFE$var
  }
  
  good <- is.finite(yi) & is.finite(vi) & (vi > 0)
  if (!any(good) || sum(good) < 4) return(list(ok = FALSE, m_effects = sum(good)))
  
  yi <- yi[good]; vi <- pmax(vi[good], 1e-12)
  m  <- length(yi)
  dat <- data.frame(yi = yi, vi = vi, ID = factor(seq_len(m)))
  
  # IV model
  fit_iv <- tryCatch(
    rma.mv(yi ~ 1, V = vi, random = ~ 1 | ID, data = dat,
           method = "REML", test = "t"),
    error = function(e) NULL
  )
  if (is.null(fit_iv)) return(list(ok = FALSE, m_effects = m))
  
  mu_iv   <- as.numeric(fit_iv$b[1])
  se_iv   <- as.numeric(sqrt(vcov(fit_iv)[1,1]))
  tau2_iv <- sum(fit_iv$sigma2)
  
  # 95% CI via predict(); fallback to confint(); then ±1.96*SE
  ci_iv_lb <- ci_iv_ub <- NA_real_
  pr <- tryCatch(predict(fit_iv), error = function(e) NULL)
  if (!is.null(pr) && all(c("ci.lb","ci.ub") %in% names(pr))) {
    ci_iv_lb <- as.numeric(pr$ci.lb)
    ci_iv_ub <- as.numeric(pr$ci.ub)
  } else {
    iv_ci <- tryCatch(confint(fit_iv, level = 0.95), error = function(e) NULL)
    if (!is.null(iv_ci)) {
      mat <- if (!is.null(iv_ci$beta)) iv_ci$beta else iv_ci
      if (!is.null(dim(mat)) && all(c("ci.lb","ci.ub") %in% colnames(mat))) {
        r <- 1L
        rn <- rownames(mat)
        if (!is.null(rn)) {
          hit <- which(grepl("intrc|Intercept", rn, ignore.case = TRUE))
          if (length(hit)) r <- hit[1]
        }
        ci_iv_lb <- suppressWarnings(as.numeric(mat[r, "ci.lb"]))
        ci_iv_ub <- suppressWarnings(as.numeric(mat[r, "ci.ub"]))
      }
    }
    if (!is.finite(ci_iv_lb) || !is.finite(ci_iv_ub)) {
      ci_iv_lb <- mu_iv - 1.96 * se_iv
      ci_iv_ub <- mu_iv + 1.96 * se_iv
    }
  }
  
  list(ok = TRUE,
       m_effects = m,
       mu_iv = mu_iv, se_iv = se_iv, tau2_iv = tau2_iv,
       ci_iv_lb = ci_iv_lb, ci_iv_ub = ci_iv_ub)
}

# -------------------- One scenario driver + MC error -----------------
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
  
  # Calibrate delta_mu to hit lnM_target GIVEN tau_delta (at mean n)
  delta_mu <- delta_mu_from_target(lnM_target, tau_delta,
                                   n1_bar = n1_mean, n2_bar = n2_mean,
                                   sd1 = sd1, sd2 = sd2)
  
  # Analytic "truth" at mean n (≈ lnM_target after calibration)
  lnM_true0 <- true_lnM_indep(mu1, mu1 + delta_mu, sd1, sd2,
                              n1_mean, n2_mean, tau_delta = tau_delta)
  
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
  R_ok <- sum(ok)
  invalid_rate <- 1 - mean(ok)
  
  # If nothing converged
  if (!any(ok)) {
    return(tibble(
      K_fixed = K_fixed, n1_mean = n1_mean, n2_mean = n2_mean,
      lnM_target = lnM_target, tau_delta = tau_delta,
      lnM_true0 = lnM_true0,
      bias_iv = NA_real_, bias_iv_mcse = NA_real_, bias_iv_lb = NA_real_, bias_iv_ub = NA_real_,
      rmse_iv = NA_real_, rmse_iv_mcse = NA_real_, rmse_iv_lb = NA_real_, rmse_iv_ub = NA_real_,
      cover_iv = NA_real_, cover_iv_mcse = NA_real_, cover_iv_lb = NA_real_, cover_iv_ub = NA_real_,
      tau2_iv_mean = NA_real_, tau2_iv_mean_mcse = NA_real_, tau2_iv_lb = NA_real_, tau2_iv_ub = NA_real_,
      mean_m_effects = NA_real_, mean_m_effects_mcse = NA_real_,
      invalid_meta_rate = 1, invalid_meta_rate_mcse = 0,
      R_ok = 0L, R_meta = R_meta
    ))
  }
  
  ext <- function(name) vapply(out_list[ok], `[[`, numeric(1), name)
  mu_iv   <- ext("mu_iv")
  tau_iv  <- ext("tau2_iv")
  m_eff   <- ext("m_effects")
  ci_lb   <- ext("ci_iv_lb")
  ci_ub   <- ext("ci_iv_ub")
  
  diffs   <- mu_iv - lnM_true0
  d2      <- diffs^2
  cover   <- (ci_lb <= lnM_true0) & (lnM_true0 <= ci_ub)
  p_cover <- mean(cover)
  
  # --- MCSEs ---
  se_bias    <- if (R_ok > 1) sd(diffs) / sqrt(R_ok) else NA_real_
  mean_d2    <- mean(d2)
  var_d2     <- if (R_ok > 1) var(d2) else NA_real_
  se_rmse    <- if (is.finite(var_d2) && mean_d2 > 0) sqrt(var_d2 / (4 * R_ok * mean_d2)) else NA_real_
  se_cover   <- if (R_ok > 0) sqrt(p_cover * (1 - p_cover) / R_ok) else NA_real_
  se_tau2    <- if (R_ok > 1) sd(tau_iv) / sqrt(R_ok) else NA_real_
  se_meff    <- if (R_ok > 1) sd(m_eff)  / sqrt(R_ok) else NA_real_
  se_invalid <- sqrt(invalid_rate * (1 - invalid_rate) / R_meta)
  
  bias_hat  <- mean(diffs)
  rmse_hat  <- sqrt(mean_d2)
  tau2_hat  <- mean(tau_iv)
  meff_hat  <- mean(m_eff)
  
  # 95% MC intervals (Normal approx); clamp coverage to [0,1]
  ci95 <- function(est, se) c(lb = est - 1.96 * se, ub = est + 1.96 * se)
  bias_ci   <- if (is.finite(se_bias))  ci95(bias_hat, se_bias) else c(lb = NA_real_, ub = NA_real_)
  rmse_ci   <- if (is.finite(se_rmse))  ci95(rmse_hat, se_rmse) else c(lb = NA_real_, ub = NA_real_)
  cover_ci  <- if (is.finite(se_cover)) ci95(p_cover, se_cover) else c(lb = NA_real_, ub = NA_real_)
  tau2_ci   <- if (is.finite(se_tau2))  ci95(tau2_hat, se_tau2) else c(lb = NA_real_, ub = NA_real_)
  
  cover_ci[1] <- max(0, cover_ci[1]); cover_ci[2] <- min(1, cover_ci[2])
  
  tibble(
    K_fixed = K_fixed, n1_mean = n1_mean, n2_mean = n2_mean,
    lnM_target = lnM_target, tau_delta = tau_delta,
    lnM_true0 = lnM_true0,
    
    bias_iv  = bias_hat,
    bias_iv_mcse = se_bias,
    bias_iv_lb = bias_ci["lb"], bias_iv_ub = bias_ci["ub"] ,
    
    rmse_iv  = rmse_hat,
    rmse_iv_mcse = se_rmse,
    rmse_iv_lb = rmse_ci["lb"], rmse_iv_ub = rmse_ci["ub"],
    
    cover_iv = p_cover,
    cover_iv_mcse = se_cover,
    cover_iv_lb = cover_ci["lb"], cover_iv_ub = cover_ci["ub"],
    
    tau2_iv_mean = tau2_hat,
    tau2_iv_mean_mcse = se_tau2,
    tau2_iv_lb = tau2_ci["lb"], tau2_iv_ub = tau2_ci["ub"],
    
    mean_m_effects = meff_hat,
    mean_m_effects_mcse = se_meff,
    
    invalid_meta_rate = invalid_rate,
    invalid_meta_rate_mcse = se_invalid,
    
    R_ok = R_ok, R_meta = R_meta
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
    n_workers = NULL,
    parallel_replicates = TRUE,
    parallel_scenarios  = FALSE,
    # Grid
    R_meta        = 400,
    K_fixed_vec   = c(20, 50),
    n_mean_vec    = c(10, 20, 30, 50, 100),
    lnM_targets   = c(-1.0, -0.7, -0.4, 0.0, 0.4),
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
  message(sprintf("[grid] Scenarios: %d (K={%s}; n̄={%s}; lnM={%s}; tau={%s}; ndist=%s)",
                  nrow(grid),
                  paste(K_fixed_vec, collapse=","),
                  paste(n_mean_vec, collapse=","),
                  paste(lnM_targets, collapse=","),
                  paste(tau_vec, collapse=","),
                  n_dist))
  
  run_one <- function(irow) {
    par <- grid[irow,]
    message(sprintf("[scenario %d/%d] K=%d; n̄=%d; lnM=%.3f; tau=%.2f  @%s",
                    irow, nrow(grid), par$K_fixed, par$n_mean, par$lnM, par$tau,
                    format(Sys.time(), "%H:%M:%S")))
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
  Sys.setenv(N_WORKERS = Sys.getenv("N_WORKERS", "16"))
  set_parallel_plan(n_workers = as.integer(Sys.getenv("N_WORKERS")), backend = "multisession")
  message(sprintf("[start] %s", format(Sys.time(), "%H:%M:%S")))
  
  demo <- run_grid_fixedK(
    n_workers = as.integer(Sys.getenv("N_WORKERS")),
    R_meta = 250,
    K_fixed_vec = c(20, 50),
    n_mean_vec = c(10, 20, 30, 50, 100),
    lnM_targets = c(-1.0, -0.7, -0.4, 0.0, 0.4),
    tau_vec = c(0, 0.2),
    n_dist = "nbinom", n_nb_size = 8L,
    n_min = 3L,
    parallel_replicates = TRUE,
    parallel_scenarios = FALSE,
    min_kept = 1000, chunk_init = 3000
  )
  
  message(sprintf("[end]   %s", format(Sys.time(), "%H:%M:%S")))
  print(demo)
}

# ========================= Plotting section ==========================
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

# Bias (IV) with MC ribbon
p_bias <- ggplot(res, aes(x = lnM_x, y = bias_iv, group = 1)) +
  geom_ribbon(aes(ymin = bias_iv_lb, ymax = bias_iv_ub), alpha = 0.18, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + geom_point(size = 1.8) +
  facet_grid(tau_f ~ setting, scales = "free_x") +
  labs(x = "lnM (true ≈ target)", y = "Bias (IV)",
       title = "Bias: inverse-variance (lines) with 95% Monte Carlo intervals (ribbons)") +
  theme_bw(base_size = 12)

# RMSE (IV) with MC ribbon
p_rmse <- ggplot(res, aes(x = lnM_x, y = rmse_iv, group = 1)) +
  geom_ribbon(aes(ymin = pmax(0, rmse_iv_lb), ymax = rmse_iv_ub), alpha = 0.18, na.rm = TRUE) +
  geom_line() + geom_point(size = 1.8) +
  facet_grid(tau_f ~ setting, scales = "free_x") +
  labs(x = "lnM (true ≈ target)", y = "RMSE (IV)",
       title = "RMSE: inverse-variance with 95% Monte Carlo intervals") +
  theme_bw(base_size = 12)

# Coverage (IV) with MC ribbon
p_cover <- ggplot(res, aes(x = lnM_x, y = cover_iv, group = 1)) +
  geom_ribbon(aes(ymin = pmax(0, cover_iv_lb), ymax = pmin(1, cover_iv_ub)), alpha = 0.18, na.rm = TRUE) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  geom_line() + geom_point(size = 1.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  facet_grid(tau_f ~ setting, scales = "free_x") +
  labs(x = "lnM (true ≈ target)", y = "Coverage (95% nominal)",
       title = "CI coverage: inverse-variance with 95% Monte Carlo intervals") +
  theme_bw(base_size = 12)

# Estimated tau^2 (IV) with MC ribbon
p_tau <- ggplot(res, aes(x = lnM_x, y = tau2_iv_mean, group = 1)) +
  geom_ribbon(aes(ymin = tau2_iv_lb, ymax = tau2_iv_ub), alpha = 0.18, na.rm = TRUE) +
  geom_line() + geom_point(size = 1.8) +
  facet_grid(tau_f ~ setting, scales = "free_x") +
  labs(x = "lnM (true ≈ target)", y = expression(mean(hat(tau)^2)~"(REML, IV)"),
       title = expression(paste("Between-study variance ", tau^2, ": mean with 95% Monte Carlo intervals"))) +
  theme_bw(base_size = 12)

# Invalid replicate rate heatmap
p_invalid <- ggplot(res, aes(x = lnM_x, y = setting, fill = invalid_meta_rate)) +
  geom_tile(color = "white") +
  facet_wrap(~ tau_f, nrow = 1) +
  scale_fill_gradient(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(x = "lnM (true ≈ target)", y = "Scenario", fill = "Invalid rate",
       title = "Proportion of failed meta-analysis replicates") +
  theme_bw(base_size = 12)

# Print
p_bias
p_rmse
p_cover
# p_tau
# p_invalid