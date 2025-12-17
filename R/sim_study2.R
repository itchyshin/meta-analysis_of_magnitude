# ===================================================================
# SAFE-only meta-analysis simulation (parallel, 16+ cores ready)
#   - Per-study effect: lnM via SAFE; point = SAFE-BC (2*Delta1 - SAFE_mean)
#   - Per-study var:    SAFE variance only (no delta-method anywhere)
#   - Compare meta-analytic estimators (all via rma.mv):
#       (A) Inverse-variance (IV): V = diag(vi)
#       (B) Multiplicative (n0-based): V = 0, R = diag(vtilde) with vtilde ∝ 1/n0
#     (also report two references: unweighted mean; n0-weighted mean)
#   - Scenario grid includes positives: 0.4, 0.7, 1.0 (and more if desired)
#   - Parallelization via future/furrr; plan set once via set_parallel_plan()
#   - To run on 16 cores locally:
#       Sys.setenv(N_WORKERS = "16"); source("this_script.R")
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
})

# ------------------------- RNG + helpers ----------------------------

SEED_MAIN <- 202501
set.seed(SEED_MAIN)

posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(gap) ifelse(gap <= 0, NA_real_, gap)

# --------------------- Parallel plan helper -------------------------
# Use multisession by default (safe everywhere). You can switch to multicore
# on non-Windows and non-RStudio by passing backend="multicore".
set_parallel_plan <- function(n_workers = NULL, backend = c("multisession","multicore")) {
  backend <- match.arg(backend)
  if (is.null(n_workers)) n_workers <- as.integer(Sys.getenv("N_WORKERS", parallel::detectCores()))
  n_workers <- max(1L, n_workers)
  if (backend == "multicore" && .Platform$OS.type != "windows" && !identical(Sys.getenv("RSTUDIO"), "1")) {
    future::plan(future::multicore, workers = n_workers)
  } else {
    future::plan(future::multisession, workers = n_workers)
  }
  invisible(n_workers)
}

# A single furrr options object to ensure workers load needed packages
.FOPTS <- furrr_options(
  seed = TRUE,
  packages = c("metafor","MASS","tibble","dplyr","purrr")
)

# --------------------- Load SAFE (new API if present) ----------------
# If SAFE_fun.R with new API is present, we'll call it; otherwise fall back.

SAFE_FILE <- "SAFE_fun.R"
SAFE_HAS_NEW <- FALSE
if (file.exists(SAFE_FILE)) {
  source(SAFE_FILE)
  SAFE_HAS_NEW <- exists("safe_lnM_indep")
}

# Fallback SAFE (B/chunk) if new API not available
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

# A single front-door SAFE call (new API if available; fallback otherwise)
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

# ----------------- Delta-1 (for bias-corrected point) ----------------
# We use Delta-1 only for the plug-in point to compute SAFE-BC = 2*Delta1 - SAFE_mean.

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

# --------------- "True" lnM and inverse mapping helpers --------------

# Target lnM -> delta using the n0 term:
#   deq^2 = 2*exp(2*lnM) + 2/n0, with deq = |delta|/sW (take sW ≈ 1 via sd1=sd2=1)
delta_from_target_lnM <- function(lnM_target, n1_bar, n2_bar, sdW = 1) {
  n0_bar <- 2 * n1_bar * n2_bar / (n1_bar + n2_bar)
  deq2   <- 2 * exp(2 * lnM_target) + 2 / n0_bar
  sdW * sqrt(pmax(deq2, 0))
}

# Large-sample analytic "truth" for lnM under random-effects in delta:
true_lnM_indep <- function(mu1, mu2, sd1, sd2, n1, n2, tau_delta = 0) {
  h          <- n1 * n2 / (n1 + n2)
  n0         <- 2 * h
  delta_mean <- mu1 - mu2
  omega      <- ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  Delta_true <- h * (delta_mean^2 + tau_delta^2)
  0.5 * (log(Delta_true / n0) - log(omega))
}

# ------------------- One meta-analysis replicate ---------------------

fit_meta_once <- function(K_study,
                          n1_mean, n2_mean,
                          delta_mu, tau_delta,
                          mu1 = 0, sd1 = 1, sd2 = 1,
                          # SAFE tuning
                          min_kept = 2000, chunk_init = 4000,
                          chunk_max = 2e6, max_draws = Inf, patience_noaccept = 5,
                          B_fallback = 2000, chunk_fallback = 4000, max_chunks_fallback = 50,
                          # multiplicative scaling option
                          scale_mult_to_vi_mean = FALSE) {
  # Per-study sample sizes (Poisson around means, min 3)
  n1_k <- pmax(3L, rpois(K_study, lambda = n1_mean))
  n2_k <- pmax(3L, rpois(K_study, lambda = n2_mean))
  # True study-level deltas
  delta_k <- if (tau_delta > 0) rnorm(K_study, mean = delta_mu, sd = tau_delta) else rep(delta_mu, K_study)
  mu2_k   <- mu1 + delta_k
  
  yi  <- numeric(K_study)
  vi  <- numeric(K_study)
  n0k <- numeric(K_study)
  
  for (k in seq_len(K_study)) {
    sm <- draw_summaries_indep(mu1, mu2_k[k], sd1, sd2, n1_k[k], n2_k[k])
    d1 <- lnM_delta1_indep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1_k[k], n2_k[k])
    SAFE <- safe_call_indep(
      sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1_k[k], n2_k[k],
      min_kept = min_kept, chunk_init = chunk_init, chunk_max = chunk_max,
      max_draws = max_draws, patience_noaccept = patience_noaccept,
      B_fallback = B_fallback, chunk_fallback = chunk_fallback, max_chunks_fallback = max_chunks_fallback
    )
    # SAFE-BC point: 2*Delta1 - SAFE_mean; if Delta1 missing, fall back to SAFE mean
    yi[k] <- if (is.na(d1["point"]) || is.na(SAFE$point)) SAFE$point else 2*d1["point"] - SAFE$point
    vi[k] <- SAFE$var
    # harmonic n0 (for multiplicative)
    n0k[k] <- 2 * n1_k[k] * n2_k[k] / (n1_k[k] + n2_k[k])
  }
  
  good <- is.finite(yi) & is.finite(vi) & (vi > 0)
  if (!any(good) || sum(good) < 4) return(list(ok = FALSE))
  
  yi <- yi[good]; vi <- pmax(vi[good], 1e-10); n0k <- n0k[good]
  m  <- length(yi)
  dat <- data.frame(yi = yi, vi = vi, n0 = n0k, ID = factor(seq_len(m)))
  
  # (A) IV meta-analysis (rma.mv with V = diag(vi))
  fit_iv <- tryCatch(
    rma.mv(yi ~ 1, V = vi, random = ~ 1 | ID, data = dat,
           method = "REML", control = list(stepadj = 0.5, maxiter = 1000)),
    error = function(e) NULL
  )
  if (is.null(fit_iv)) return(list(ok = FALSE))
  mu_iv   <- as.numeric(fit_iv$b[1]); se_iv <- sqrt(vcov(fit_iv))[1,1]; tau2_iv <- sum(fit_iv$sigma2)
  
  # (B) Multiplicative (n0-based) via rma.mv using R = diag(vtilde) with vtilde ∝ 1/n0
  inv_w_n0 <- 1 / dat$n0  # ∝ 1 / n0
  vtilde   <- inv_w_n0
  if (isTRUE(scale_mult_to_vi_mean)) vtilde <- vtilde * (mean(dat$vi) / mean(vtilde))
  Vf <- diag(as.numeric(vtilde)); levs <- levels(dat$ID); rownames(Vf) <- levs; colnames(Vf) <- levs
  
  fit_mult <- tryCatch(
    rma.mv(yi ~ 1, V = 0, random = ~ 1 | ID, data = dat,
           R = list(ID = Vf), Rscale = FALSE,
           method = "REML", control = list(stepadj = 0.5, maxiter = 1000)),
    error = function(e) NULL
  )
  if (is.null(fit_mult)) return(list(ok = FALSE))
  mu_mult <- as.numeric(fit_mult$b[1]); se_mult <- sqrt(vcov(fit_mult))[1,1]; tau2_mult <- sum(fit_mult$sigma2)
  
  # References
  mu_uw  <- mean(dat$yi); se_uw  <- sqrt(var(dat$yi) / m)
  mu_wn0 <- weighted.mean(dat$yi, w = dat$n0)
  se_wn0 <- sqrt( sum(dat$n0^2 * (dat$yi - mu_wn0)^2) / (sum(dat$n0)^2 * (m - 1)) )
  
  list(ok = TRUE,
       mu_iv = mu_iv,   se_iv = se_iv,   tau2_iv = tau2_iv,
       mu_mult = mu_mult, se_mult = se_mult, tau2_mult = tau2_mult,
       mu_uw = mu_uw, se_uw = se_uw, mu_wn0 = mu_wn0, se_wn0 = se_wn0,
       m_effects = m)
}

# -------------------- One scenario driver (parallel) -----------------

run_scenario <- function(R_meta,
                         K_study,
                         n1_mean, n2_mean,
                         lnM_target, tau_delta,
                         mu1 = 0, sd1 = 1, sd2 = 1,
                         min_kept = 2000, chunk_init = 4000,
                         chunk_max = 2e6, max_draws = Inf, patience_noaccept = 5,
                         B_fallback = 2000, chunk_fallback = 4000, max_chunks_fallback = 50,
                         scale_mult_to_vi_mean = FALSE,
                         parallel_replicates = TRUE) {
  # Back-solve delta_mu for the requested lnM_target using n0 identity
  delta_mu <- delta_from_target_lnM(lnM_target, n1_bar = n1_mean, n2_bar = n2_mean, sdW = 1)
  # Large-sample analytic "truth" for reporting
  lnM_true0 <- true_lnM_indep(mu1, mu1 + delta_mu, sd1, sd2, n1_mean, n2_mean, tau_delta = tau_delta)
  
  one_rep <- function(r) {
    fit_meta_once(
      K_study = K_study,
      n1_mean = n1_mean, n2_mean = n2_mean,
      delta_mu = delta_mu, tau_delta = tau_delta,
      mu1 = mu1, sd1 = sd1, sd2 = sd2,
      min_kept = min_kept, chunk_init = chunk_init,
      chunk_max = chunk_max, max_draws = max_draws, patience_noaccept = patience_noaccept,
      B_fallback = B_fallback, chunk_fallback = chunk_fallback, max_chunks_fallback = max_chunks_fallback,
      scale_mult_to_vi_mean = scale_mult_to_vi_mean
    )
  }
  
  reps <- seq_len(R_meta)
  out_list <- if (parallel_replicates) {
    future_map(reps, one_rep, .options = .FOPTS)
  } else {
    map(reps, one_rep)
  }
  
  ok <- vapply(out_list, function(z) isTRUE(z$ok), logical(1))
  invalid_rate <- 1 - mean(ok)
  
  if (!any(ok)) {
    return(tibble(
      K_study = K_study, n1 = n1_mean, n2 = n2_mean,
      lnM_target = lnM_target, tau_delta = tau_delta,
      lnM_true0 = lnM_true0,
      bias_iv = NA_real_, bias_mult = NA_real_, bias_uw = NA_real_, bias_wn0 = NA_real_,
      rmse_iv = NA_real_, rmse_mult = NA_real_, rmse_uw = NA_real_, rmse_wn0 = NA_real_,
      cover_iv = NA_real_, cover_mult = NA_real_,
      tau2_iv_mean = NA_real_, tau2_mult_mean = NA_real_,
      invalid_meta_rate = 1
    ))
  }
  
  ext <- function(name) vapply(out_list[ok], `[[`, numeric(1), name)
  mu_iv   <- ext("mu_iv");   se_iv   <- ext("se_iv");   tau_iv   <- ext("tau2_iv")
  mu_mult <- ext("mu_mult"); se_mult <- ext("se_mult"); tau_mult <- ext("tau2_mult")
  mu_uw   <- ext("mu_uw");   mu_wn0  <- ext("mu_wn0")
  
  # Performance summaries vs lnM_true0
  bias_iv   <- mean(mu_iv   - lnM_true0)
  bias_mult <- mean(mu_mult - lnM_true0)
  bias_uw   <- mean(mu_uw   - lnM_true0)
  bias_wn0  <- mean(mu_wn0  - lnM_true0)
  
  rmse_iv   <- sqrt(mean((mu_iv   - lnM_true0)^2))
  rmse_mult <- sqrt(mean((mu_mult - lnM_true0)^2))
  rmse_uw   <- sqrt(mean((mu_uw   - lnM_true0)^2))
  rmse_wn0  <- sqrt(mean((mu_wn0  - lnM_true0)^2))
  
  cover_iv   <- mean(abs(mu_iv   - lnM_true0) <= 1.96 * se_iv)
  cover_mult <- mean(abs(mu_mult - lnM_true0) <= 1.96 * se_mult)
  
  tibble(
    K_study = K_study, n1 = n1_mean, n2 = n2_mean,
    lnM_target = lnM_target, tau_delta = tau_delta,
    lnM_true0 = lnM_true0,
    bias_iv = bias_iv, bias_mult = bias_mult, bias_uw = bias_uw, bias_wn0 = bias_wn0,
    rmse_iv = rmse_iv, rmse_mult = rmse_mult, rmse_uw = rmse_uw, rmse_wn0 = rmse_wn0,
    cover_iv = cover_iv, cover_mult = cover_mult,
    tau2_iv_mean = mean(tau_iv), tau2_mult_mean = mean(tau_mult),
    invalid_meta_rate = invalid_rate
  )
}

# ---------------------- Master grid (parallel) -----------------------

run_grid <- function(
    # Parallel config (plan is set outside via set_parallel_plan(); we don't touch it here)
  n_workers = NULL,                      # kept for interface symmetry; not used internally
  parallel_replicates = TRUE,
  parallel_scenarios  = FALSE,
  # Scenario grid
  R_meta      = 400,
  K_study_vec = c(20, 50),
  n_vec       = c(10, 20),
  lnM_targets = c(-2.0, -1.5, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0),
  tau_vec     = c(0.2, 0.4),
  # SAFE tuning
  min_kept = 2000, chunk_init = 4000,
  chunk_max = 2e6, max_draws = Inf, patience_noaccept = 5,
  B_fallback = 2000, chunk_fallback = 4000, max_chunks_fallback = 50,
  scale_mult_to_vi_mean = FALSE
) {
  grid <- expand.grid(
    K_study = K_study_vec,
    n       = n_vec,
    lnM     = lnM_targets,
    tau     = tau_vec,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid <- tibble::as_tibble(grid)
  
  run_one <- function(irow) {
    par <- grid[irow,]
    run_scenario(
      R_meta = R_meta,
      K_study = par$K_study,
      n1_mean = par$n, n2_mean = par$n,
      lnM_target = par$lnM,
      tau_delta = par$tau,
      # SAFE tuning
      min_kept = min_kept, chunk_init = chunk_init,
      chunk_max = chunk_max, max_draws = max_draws, patience_noaccept = patience_noaccept,
      B_fallback = B_fallback, chunk_fallback = chunk_fallback, max_chunks_fallback = max_chunks_fallback,
      scale_mult_to_vi_mean = scale_mult_to_vi_mean,
      parallel_replicates = parallel_replicates
    ) %>% mutate(scenario_id = irow)
  }
  
  idx <- seq_len(nrow(grid))
  res_list <- if (parallel_scenarios) {
    future_map(idx, run_one, .options = .FOPTS)
  } else {
    map(idx, run_one)
  }
  
  bind_rows(res_list) %>%
    relocate(scenario_id) %>%
    arrange(scenario_id, K_study, n1, lnM_target, tau_delta)
}

# --------------------------- Demo run --------------------------------
# Use 16 cores locally:
#   Sys.setenv(N_WORKERS = "16")
#   source("this_script.R")
# Or run manually:
#   set_parallel_plan(16, backend="multisession")
#   results <- run_grid(R_meta=400, parallel_replicates=TRUE, parallel_scenarios=FALSE)

if (sys.nframe() == 0) {
  message("=== Quick demo (16 cores, multisession) ===")
  Sys.setenv(N_WORKERS = Sys.getenv("N_WORKERS", "16"))
  set_parallel_plan(n_workers = as.integer(Sys.getenv("N_WORKERS")), backend = "multisession")
  
  demo <- run_grid(
    n_workers = as.integer(Sys.getenv("N_WORKERS")),  # not used internally, just for clarity
    R_meta = 200,
    K_study_vec = c(20, 50),
    n_vec = c(10, 20),
    lnM_targets = c(-1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0),
    tau_vec = c(0.2, 0.4),
    parallel_replicates = TRUE,
    parallel_scenarios = FALSE,
    min_kept = 1000, chunk_init = 3000,
    B_fallback = 1000, chunk_fallback = 3000
  )
  print(demo)
}