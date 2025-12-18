# ================================================================
# simulation16_parallel_SAFE_with_plots_FIXED.R
# lnM simulations (PI vs SAFE-BC), parallel, with plots
#
# FIXES vs your current sim16:
#  1) PI is no longer NA when MSB-MSW <= 0:
#     we truncate the gap: Delta <- max(MSB - MSW, EPS_GAP)
#     (sim15-style "make it defined", and we track truncation rate)
#  2) Delta-method variance is capped at DELTA_VAR_CAP = 20 (paper/sim15)
#     and we track cap-hit rate/n
#  3) Relative-bias is simplified to ONE baseline:
#     true variance = MC variance of PI point estimator (Var_MC_PI)
# ================================================================

suppressPackageStartupMessages({
  library(MASS)      # mvrnorm, rWishart
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(pbapply)
  library(parallel)
  library(here)
})

theme_set(theme_bw(11))

# ------------------- Load SAFE ----------------------------------
SAFE_FILE <- here("R","SAFE_fun.R")
if (!file.exists(SAFE_FILE)) {
  stop("SAFE_fun.R not found. Please place it at: ", SAFE_FILE)
}
source(SAFE_FILE)
if (!exists("safe_lnM_indep") || !exists("safe_lnM_dep")) {
  stop("SAFE_fun.R must define safe_lnM_indep() and safe_lnM_dep().")
}

# Optional: prevent oversubscription (OpenBLAS/MKL)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

# ------------------- Global knobs --------------------------------
RNGkind("L'Ecuyer-CMRG")
SEED_MAIN <- 12345
set.seed(SEED_MAIN)

K_default <- 2000

MIN_KEPT_SAFE    <- 2000
CHUNK_INIT_SAFE  <- 4000
CHUNK_MAX_SAFE   <- 2e6
MAX_DRAWS_SAFE   <- Inf
PATIENCE_SAFE    <- 5

# Paper/sim15 convention
DELTA_VAR_CAP <- 20

# Critical sim15-style “define PI even if MSB-MSW <= 0”
EPS_GAP <- 1e-12

r_default <- 0.8

theta_grid <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,4,5)

designs_indep <- list(
  c(5,5), c(10,10), c(20,20), c(100,100),
  c(3,7), c(6,14), c(12,28), c(40,160)
)
designs_paired_n <- c(5,10,20,100)

# ------------------- Helpers ------------------------------------
posify <- function(x, eps = 1e-12) pmax(x, eps)

safe_mean <- function(x) {
  y <- x[is.finite(x)]
  if (length(y)) mean(y) else NA_real_
}
safe_var <- function(x) {
  y <- x[is.finite(x)]
  if (length(y) >= 2) var(y) else NA_real_
}
safe_prop <- function(x, empty_value = NA_real_) {
  y <- x[is.finite(x)]
  if (length(y)) mean(y) else empty_value
}
safe_rmse <- function(err) {
  y <- err[is.finite(err)]
  if (length(y)) sqrt(mean(y^2)) else NA_real_
}
relbias_pct <- function(est, target) {
  if (is.finite(est) && is.finite(target) && target > 0) 100 * (est - target) / target else NA_real_
}

# ------------------- lnM Delta-1 (independent) ------------------
# PI point + delta-method variance with:
#  - gap truncation to EPS_GAP (so PI is always defined)
#  - variance cap at DELTA_VAR_CAP
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                             cap = DELTA_VAR_CAP, eps_gap = EPS_GAP) {
  h    <- n1 * n2 / (n1 + n2)
  MSB  <- h * (x1bar - x2bar)^2
  MSW  <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  
  Delta_raw <- MSB - MSW
  trunc_flag <- as.numeric(is.finite(Delta_raw) && (Delta_raw <= eps_gap))
  Delta <- if (is.finite(Delta_raw)) max(Delta_raw, eps_gap) else NA_real_
  if (!is.finite(Delta) || !is.finite(MSW) || MSW <= 0) {
    return(c(point = NA_real_, var = NA_real_, capped = NA_real_, trunc = NA_real_))
  }
  
  lnM <- 0.5 * (log(Delta / (2*h)) - log(MSW))
  
  sigmaD2 <- s1^2 / n1 + s2^2 / n2
  delta   <- x1bar - x2bar
  
  Var_B   <- h^2 * (2 * sigmaD2^2 + 4 * sigmaD2 * delta^2)
  Var_W   <- 2 * MSW^2 / (n1 + n2 - 2)
  
  gB   <- 0.5 / Delta
  gW   <- -0.5 * MSB / (Delta * MSW)
  
  Var1_raw <- posify(gB^2 * Var_B + gW^2 * Var_W)
  cap_flag <- as.numeric(is.finite(Var1_raw) && (Var1_raw > cap))
  Var1 <- pmin(Var1_raw, cap)
  
  c(point = lnM, var = Var1, capped = cap_flag, trunc = trunc_flag)
}

# ------------------- lnM Delta-1 (paired) -----------------------
lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r,
                           cap = DELTA_VAR_CAP, eps_gap = EPS_GAP) {
  h     <- n / 2
  MSB   <- h * (x1bar - x2bar)^2
  MSW   <- (s1^2 + s2^2) / 2
  
  Delta_raw <- MSB - MSW
  trunc_flag <- as.numeric(is.finite(Delta_raw) && (Delta_raw <= eps_gap))
  Delta <- if (is.finite(Delta_raw)) max(Delta_raw, eps_gap) else NA_real_
  if (!is.finite(Delta) || !is.finite(MSW) || MSW <= 0) {
    return(c(point = NA_real_, var = NA_real_, capped = NA_real_, trunc = NA_real_))
  }
  
  lnM <- 0.5 * (log(Delta / n) - log(MSW))
  
  sigmaD2 <- s1^2 + s2^2 - 2 * r * s1 * s2
  delta   <- x1bar - x2bar
  
  Var_B   <- h^2 * (2 * sigmaD2^2 / n^2 + 4 * delta^2 * sigmaD2 / n)
  Var_W   <- (s1^4 + s2^4 + 2 * r^2 * s1^2 * s2^2) / (2 * (n - 1))
  
  gB   <- 0.5 / Delta
  gW   <- -0.5 * MSB / (Delta * MSW)
  
  Var1_raw <- posify(gB^2 * Var_B + gW^2 * Var_W)
  cap_flag <- as.numeric(is.finite(Var1_raw) && (Var1_raw > cap))
  Var1 <- pmin(Var1_raw, cap)
  
  c(point = lnM, var = Var1, capped = cap_flag, trunc = trunc_flag)
}

# ------------------- "True" lnM (sim16 definition) --------------
true_lnM_indep <- function(mu1, mu2, sd1, sd2, n1, n2) {
  h     <- n1*n2/(n1+n2)
  n0    <- 2*h
  delta <- mu1 - mu2
  omega <- ((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2)
  Delta_true <- h * delta^2
  0.5 * (log(Delta_true / n0) - log(omega))
}
true_lnM_dep <- function(mu1, mu2, sd1, sd2, n, r) {
  delta <- mu1 - mu2
  sigmaD2 <- sd1^2 + sd2^2 - 2*r*sd1*sd2
  E_MSW   <- (sd1^2 + sd2^2)/2
  Delta_true <- (n/2)*delta^2 + (1/2)*sigmaD2 - E_MSW
  n0 <- n
  0.5 * (log(Delta_true / n0) - log(E_MSW))
}

# ------------------- Summary draws (Normal theory) --------------
draw_summaries_indep <- function(mu1, mu2, sd1, sd2, n1, n2) {
  x1bar <- rnorm(1L, mu1, sd1/sqrt(n1))
  x2bar <- rnorm(1L, mu2, sd2/sqrt(n2))
  s1sq  <- sd1^2 * rchisq(1L, df = n1-1) / (n1-1)
  s2sq  <- sd2^2 * rchisq(1L, df = n2-1) / (n2-1)
  c(x1bar = x1bar, x2bar = x2bar, s1 = sqrt(s1sq), s2 = sqrt(s2sq))
}

draw_summaries_dep <- function(mu1, mu2, sd1, sd2, n, r) {
  Sig_pop <- matrix(c(sd1^2, r*sd1*sd2,
                      r*sd1*sd2, sd2^2), 2, 2)
  mu  <- MASS::mvrnorm(1L, mu = c(mu1, mu2), Sigma = Sig_pop / n)
  W   <- stats::rWishart(1L, df = n-1, Sigma = Sig_pop)
  S11 <- W[1,1,1] / (n-1)
  S22 <- W[2,2,1] / (n-1)
  c(x1bar = mu[1], x2bar = mu[2], s1 = sqrt(S11), s2 = sqrt(S22))
}

# ------------------- One parameter set: independent -------------
run_param_indep <- function(theta, n1, n2,
                            K = K_default,
                            min_kept = MIN_KEPT_SAFE,
                            chunk_init = CHUNK_INIT_SAFE,
                            chunk_max = CHUNK_MAX_SAFE,
                            max_draws = MAX_DRAWS_SAFE,
                            patience_noaccept = PATIENCE_SAFE,
                            mu1 = 0, sd1 = 1, sd2 = 1) {
  
  mu2 <- theta
  lnM_true <- true_lnM_indep(mu1, mu2, sd1, sd2, n1, n2)
  
  lnM_PI    <- rep(NA_real_, K)
  Var_delta <- rep(NA_real_, K)
  cap_hit   <- rep(NA_real_, K)
  trunc_hit <- rep(NA_real_, K)
  
  lnM_BC   <- rep(NA_real_, K)
  Var_SAFE <- rep(NA_real_, K)
  kept_vec   <- integer(K)
  total_vec  <- integer(K)
  status_vec <- character(K)
  
  for (k in seq_len(K)) {
    sm <- draw_summaries_indep(mu1, mu2, sd1, sd2, n1, n2)
    
    d1 <- lnM_delta1_indep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1, n2)
    
    SAFE <- safe_lnM_indep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1, n2,
                           min_kept = min_kept,
                           chunk_init = chunk_init,
                           chunk_max = chunk_max,
                           max_draws = max_draws,
                           patience_noaccept = patience_noaccept)
    
    lnM_PI[k]    <- unname(d1["point"])
    Var_delta[k] <- unname(d1["var"])
    cap_hit[k]   <- unname(d1["capped"])
    trunc_hit[k] <- unname(d1["trunc"])
    
    kept_vec[k]   <- SAFE$kept
    total_vec[k]  <- SAFE$total
    status_vec[k] <- SAFE$status
    
    # SAFE-BC: 2*PI - SAFE(point); if PI missing, fall back to SAFE(point)
    lnM_BC[k]   <- if (!is.finite(lnM_PI[k])) SAFE$point else 2*lnM_PI[k] - SAFE$point
    Var_SAFE[k] <- SAFE$var
  }
  
  pi_ok <- is.finite(lnM_PI) & is.finite(Var_delta) & (Var_delta > 0)
  bc_ok <- is.finite(lnM_BC) & is.finite(Var_SAFE)  & (Var_SAFE > 0)
  
  # "True variance" baseline = MC variance of PI point estimator
  Var_MC_PI <- safe_var(lnM_PI[pi_ok])
  
  Var_delta_bar <- safe_mean(Var_delta[pi_ok])
  Var_SAFE_bar  <- safe_mean(Var_SAFE[bc_ok])
  
  # ONE relative-bias set (baseline = Var_MC_PI)
  rb_delta <- relbias_pct(Var_delta_bar, Var_MC_PI)
  rb_SAFE  <- relbias_pct(Var_SAFE_bar,  Var_MC_PI)
  
  mu_PI <- safe_mean(lnM_PI[pi_ok])
  mu_BC <- safe_mean(lnM_BC[bc_ok])
  
  bias_PI <- if (is.finite(mu_PI)) mu_PI - lnM_true else NA_real_
  bias_BC <- if (is.finite(mu_BC)) mu_BC - lnM_true else NA_real_
  
  cover_PI <- safe_prop(abs(lnM_PI[pi_ok] - lnM_true) <= 1.96 * sqrt(Var_delta[pi_ok]))
  cover_BC <- safe_prop(abs(lnM_BC[bc_ok] - lnM_true) <= 1.96 * sqrt(Var_SAFE[bc_ok]))
  
  rmse_PI <- safe_rmse(lnM_PI[pi_ok] - lnM_true)
  rmse_BC <- safe_rmse(lnM_BC[bc_ok] - lnM_true)
  
  # Diagnostics
  PI_valid_rate        <- safe_prop(pi_ok, empty_value = 0)
  PI_gap_trunc_rate    <- safe_prop(trunc_hit[pi_ok], empty_value = 0)
  delta_var_capped_rate<- safe_prop(cap_hit[pi_ok], empty_value = 0)
  delta_var_capped_n   <- sum(cap_hit[pi_ok] == 1, na.rm = TRUE)
  
  SAFE_kept_rate <- safe_mean(kept_vec / pmax(1L, total_vec))
  SAFE_status_ok <- safe_prop(status_vec == "ok")
  
  # Old-style aliases
  delta_mean <- mu_PI
  safe_mean_ <- mu_BC
  delta_bias <- bias_PI
  safe_bias  <- bias_BC
  mean_var_delta <- Var_delta_bar
  mean_var_safe  <- Var_SAFE_bar
  relbias_delta  <- rb_delta
  relbias_safe   <- rb_SAFE
  rmse_delta     <- rmse_PI
  rmse_safe      <- rmse_BC
  cover_delta    <- cover_PI
  cover_safe     <- cover_BC
  
  tibble(
    design = "indep",
    n1 = n1, n2 = n2, r = NA_real_,
    theta = theta,
    lnM_true = lnM_true,
    
    delta_var_cap = DELTA_VAR_CAP,
    EPS_GAP = EPS_GAP,
    PI_gap_trunc_rate = PI_gap_trunc_rate,
    delta_var_capped_rate = delta_var_capped_rate,
    delta_var_capped_n    = delta_var_capped_n,
    
    bias_PI = bias_PI,
    bias_BC = bias_BC,
    Var_MC_PI = Var_MC_PI,
    Var_delta_bar = Var_delta_bar,
    Var_SAFE_bar  = Var_SAFE_bar,
    
    rb_delta = rb_delta,
    rb_SAFE  = rb_SAFE,
    
    cover_PI = cover_PI,
    cover_BC = cover_BC,
    rmse_PI  = rmse_PI,
    rmse_BC  = rmse_BC,
    
    PI_valid_rate = PI_valid_rate,
    SAFE_kept_rate = SAFE_kept_rate,
    SAFE_status_ok = SAFE_status_ok,
    
    # old-style aliases
    delta_mean = delta_mean,
    safe_mean  = safe_mean_,
    delta_bias = delta_bias,
    safe_bias  = safe_bias,
    mean_var_delta = mean_var_delta,
    mean_var_safe  = mean_var_safe,
    relbias_delta  = relbias_delta,
    relbias_safe   = relbias_safe,
    rmse_delta     = rmse_delta,
    rmse_safe      = rmse_safe,
    cover_delta    = cover_delta,
    cover_safe     = cover_safe
  )
}

# ------------------- One parameter set: paired ------------------
run_param_dep <- function(theta, n,
                          r = r_default,
                          K = K_default,
                          min_kept = MIN_KEPT_SAFE,
                          chunk_init = CHUNK_INIT_SAFE,
                          chunk_max = CHUNK_MAX_SAFE,
                          max_draws = MAX_DRAWS_SAFE,
                          patience_noaccept = PATIENCE_SAFE,
                          mu1 = 0, sd1 = 1, sd2 = 1) {
  
  mu2 <- theta
  lnM_true <- true_lnM_dep(mu1, mu2, sd1, sd2, n, r)
  
  lnM_PI    <- rep(NA_real_, K)
  Var_delta <- rep(NA_real_, K)
  cap_hit   <- rep(NA_real_, K)
  trunc_hit <- rep(NA_real_, K)
  
  lnM_BC   <- rep(NA_real_, K)
  Var_SAFE <- rep(NA_real_, K)
  kept_vec   <- integer(K)
  total_vec  <- integer(K)
  status_vec <- character(K)
  
  for (k in seq_len(K)) {
    sm <- draw_summaries_dep(mu1, mu2, sd1, sd2, n, r)
    
    d1 <- lnM_delta1_dep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n, r)
    
    SAFE <- safe_lnM_dep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n, r,
                         min_kept = min_kept,
                         chunk_init = chunk_init,
                         chunk_max = chunk_max,
                         max_draws = max_draws,
                         patience_noaccept = patience_noaccept)
    
    lnM_PI[k]    <- unname(d1["point"])
    Var_delta[k] <- unname(d1["var"])
    cap_hit[k]   <- unname(d1["capped"])
    trunc_hit[k] <- unname(d1["trunc"])
    
    kept_vec[k]   <- SAFE$kept
    total_vec[k]  <- SAFE$total
    status_vec[k] <- SAFE$status
    
    lnM_BC[k]   <- if (!is.finite(lnM_PI[k])) SAFE$point else 2*lnM_PI[k] - SAFE$point
    Var_SAFE[k] <- SAFE$var
  }
  
  pi_ok <- is.finite(lnM_PI) & is.finite(Var_delta) & (Var_delta > 0)
  bc_ok <- is.finite(lnM_BC) & is.finite(Var_SAFE)  & (Var_SAFE > 0)
  
  Var_MC_PI <- safe_var(lnM_PI[pi_ok])
  
  Var_delta_bar <- safe_mean(Var_delta[pi_ok])
  Var_SAFE_bar  <- safe_mean(Var_SAFE[bc_ok])
  
  rb_delta <- relbias_pct(Var_delta_bar, Var_MC_PI)
  rb_SAFE  <- relbias_pct(Var_SAFE_bar,  Var_MC_PI)
  
  mu_PI <- safe_mean(lnM_PI[pi_ok])
  mu_BC <- safe_mean(lnM_BC[bc_ok])
  
  bias_PI <- if (is.finite(mu_PI)) mu_PI - lnM_true else NA_real_
  bias_BC <- if (is.finite(mu_BC)) mu_BC - lnM_true else NA_real_
  
  cover_PI <- safe_prop(abs(lnM_PI[pi_ok] - lnM_true) <= 1.96 * sqrt(Var_delta[pi_ok]))
  cover_BC <- safe_prop(abs(lnM_BC[bc_ok] - lnM_true) <= 1.96 * sqrt(Var_SAFE[bc_ok]))
  
  rmse_PI <- safe_rmse(lnM_PI[pi_ok] - lnM_true)
  rmse_BC <- safe_rmse(lnM_BC[bc_ok] - lnM_true)
  
  PI_valid_rate         <- safe_prop(pi_ok, empty_value = 0)
  PI_gap_trunc_rate     <- safe_prop(trunc_hit[pi_ok], empty_value = 0)
  delta_var_capped_rate <- safe_prop(cap_hit[pi_ok], empty_value = 0)
  delta_var_capped_n    <- sum(cap_hit[pi_ok] == 1, na.rm = TRUE)
  
  SAFE_kept_rate <- safe_mean(kept_vec / pmax(1L, total_vec))
  SAFE_status_ok <- safe_prop(status_vec == "ok")
  
  # old-style
  delta_mean <- mu_PI
  safe_mean_ <- mu_BC
  delta_bias <- bias_PI
  safe_bias  <- bias_BC
  mean_var_delta <- Var_delta_bar
  mean_var_safe  <- Var_SAFE_bar
  relbias_delta  <- rb_delta
  relbias_safe   <- rb_SAFE
  rmse_delta     <- rmse_PI
  rmse_safe      <- rmse_BC
  cover_delta    <- cover_PI
  cover_safe     <- cover_BC
  
  tibble(
    design = "paired",
    n1 = n, n2 = n, r = r,
    theta = theta,
    lnM_true = lnM_true,
    
    delta_var_cap = DELTA_VAR_CAP,
    EPS_GAP = EPS_GAP,
    PI_gap_trunc_rate = PI_gap_trunc_rate,
    delta_var_capped_rate = delta_var_capped_rate,
    delta_var_capped_n    = delta_var_capped_n,
    
    bias_PI = bias_PI,
    bias_BC = bias_BC,
    Var_MC_PI = Var_MC_PI,
    Var_delta_bar = Var_delta_bar,
    Var_SAFE_bar  = Var_SAFE_bar,
    
    rb_delta = rb_delta,
    rb_SAFE  = rb_SAFE,
    
    cover_PI = cover_PI,
    cover_BC = cover_BC,
    rmse_PI  = rmse_PI,
    rmse_BC  = rmse_BC,
    
    PI_valid_rate = PI_valid_rate,
    SAFE_kept_rate = SAFE_kept_rate,
    SAFE_status_ok = SAFE_status_ok,
    
    # old-style aliases
    delta_mean = delta_mean,
    safe_mean  = safe_mean_,
    delta_bias = delta_bias,
    safe_bias  = safe_bias,
    mean_var_delta = mean_var_delta,
    mean_var_safe  = mean_var_safe,
    relbias_delta  = relbias_delta,
    relbias_safe   = relbias_safe,
    rmse_delta     = rmse_delta,
    rmse_safe      = rmse_safe,
    cover_delta    = cover_delta,
    cover_safe     = cover_safe
  )
}

# ------------------- Parameter grid -----------------------------
build_param_grid <- function() {
  indep_tbl <- tidyr::expand_grid(theta = theta_grid, dn = designs_indep) %>%
    mutate(n1 = purrr::map_int(dn, 1L), n2 = purrr::map_int(dn, 2L)) %>%
    select(-dn) %>%
    mutate(design = "indep")
  
  paired_tbl <- tidyr::expand_grid(theta = theta_grid, n = designs_paired_n) %>%
    transmute(theta = theta, n1 = n, n2 = n, design = "paired")
  
  bind_rows(indep_tbl, paired_tbl)
}

# ------------------- Parallel driver ----------------------------
run_full_simulation_parallel <- function(
    K = K_default,
    min_kept = MIN_KEPT_SAFE,
    chunk_init = CHUNK_INIT_SAFE,
    chunk_max = CHUNK_MAX_SAFE,
    max_draws = MAX_DRAWS_SAFE,
    patience_noaccept = PATIENCE_SAFE,
    r = r_default,
    n_cores = NULL,
    cluster_type = c("auto","psock","fork"),
    export_SAFE_fun = TRUE,
    safe_fun_path = SAFE_FILE
) {
  pg <- build_param_grid()
  n_tasks <- nrow(pg)
  
  if (is.null(n_cores)) {
    env_cores <- suppressWarnings(as.integer(Sys.getenv("N_CORES", "")))
    n_cores <- if (!is.na(env_cores) && env_cores > 0) env_cores else parallel::detectCores()
  }
  n_cores <- min(n_cores, n_tasks)
  cluster_type <- match.arg(cluster_type)
  
  .run_one_param <- function(row) {
    if (row$design == "indep") {
      run_param_indep(theta = row$theta, n1 = row$n1, n2 = row$n2,
                      K = K,
                      min_kept = min_kept,
                      chunk_init = chunk_init,
                      chunk_max = chunk_max,
                      max_draws = max_draws,
                      patience_noaccept = patience_noaccept)
    } else {
      run_param_dep(theta = row$theta, n = row$n1, r = r,
                    K = K,
                    min_kept = min_kept,
                    chunk_init = chunk_init,
                    chunk_max = chunk_max,
                    max_draws = max_draws,
                    patience_noaccept = patience_noaccept)
    }
  }
  
  if (n_cores <= 1) {
    message("Running serially...")
    res_list <- pbapply::pblapply(seq_len(n_tasks), function(i) .run_one_param(pg[i, ]))
    return(dplyr::bind_rows(res_list))
  }
  
  on_windows <- (.Platform$OS.type == "windows")
  use_fork   <- (!on_windows) && (cluster_type %in% c("fork","auto"))
  use_psock  <- on_windows || cluster_type == "psock" || (!use_fork && cluster_type == "auto")
  
  if (use_fork) {
    options(mc.cores = n_cores)
    if (export_SAFE_fun && file.exists(safe_fun_path)) source(safe_fun_path)
    res_list <- pbapply::pblapply(seq_len(n_tasks), function(i) .run_one_param(pg[i, ]), cl = n_cores)
  } else if (use_psock) {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, iseed = SEED_MAIN)
    
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(MASS); library(dplyr); library(tidyr); library(purrr); library(tibble)
      })
    })
    
    if (export_SAFE_fun && file.exists(safe_fun_path)) {
      parallel::clusterExport(cl, varlist = "safe_fun_path", envir = environment())
      parallel::clusterEvalQ(cl, source(safe_fun_path))
    }
    
    parallel::clusterExport(
      cl,
      varlist = c(
        "K","min_kept","chunk_init","chunk_max","max_draws","patience_noaccept","r",
        "DELTA_VAR_CAP","EPS_GAP",
        "pg","theta_grid","designs_indep","designs_paired_n",
        "posify",
        "safe_mean","safe_var","safe_prop","safe_rmse","relbias_pct",
        "lnM_delta1_indep","lnM_delta1_dep",
        "true_lnM_indep","true_lnM_dep",
        "draw_summaries_indep","draw_summaries_dep",
        "run_param_indep","run_param_dep",
        "SEED_MAIN"
      ),
      envir = environment()
    )
    
    res_list <- pbapply::pblapply(seq_len(n_tasks), function(i) .run_one_param(pg[i, ]), cl = cl)
  } else {
    stop("Unknown cluster_type")
  }
  
  dplyr::bind_rows(res_list)
}

# ------------------- Plot helpers -------------------------------
facet_ordering <- function(results) {
  results %>%
    mutate(facet_label = paste0(design, " n1=", n1, " n2=", n2)) %>%
    {
      facet_info <- distinct(., facet_label, design, n1, n2)
      balanced_ind <- facet_info %>%
        filter(design == "indep", n1 == n2) %>% arrange(n1) %>% pull(facet_label)
      unbalanced_ind <- facet_info %>%
        filter(design == "indep", n1 != n2) %>% arrange(n1) %>% pull(facet_label)
      paired_balanced <- facet_info %>%
        filter(design == "paired") %>% arrange(n1) %>% pull(facet_label)
      new_levels <- c(balanced_ind, unbalanced_ind, paired_balanced)
      mutate(., facet_label = factor(facet_label, levels = new_levels))
    }
}

plot_bias <- function(results) {
  df <- bind_rows(
    results %>% transmute(theta, facet_label, estimator = "PI",      bias = bias_PI),
    results %>% transmute(theta, facet_label, estimator = "SAFE-BC", bias = bias_BC)
  )
  ggplot(df, aes(theta, bias, colour = estimator, group = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = "theta", y = "Bias (estimate - true lnM)", colour = NULL) +
    scale_colour_manual(values = c("PI" = "firebrick", "SAFE-BC" = "steelblue"))
}

plot_relbias_simple <- function(results) {
  df <- bind_rows(
    results %>% transmute(theta, facet_label, estimator = "delta-var", relbias = rb_delta),
    results %>% transmute(theta, facet_label, estimator = "SAFE-var",  relbias = rb_SAFE)
  )
  ggplot(df, aes(theta, relbias, colour = estimator, group = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = "theta", y = "Relative bias of variance (%) vs MC(PI)", colour = NULL) +
    scale_colour_manual(values = c("delta-var" = "firebrick", "SAFE-var" = "steelblue"))
}

plot_coverage <- function(results) {
  df <- bind_rows(
    results %>% transmute(theta, facet_label, estimator = "PI",      cover = cover_PI),
    results %>% transmute(theta, facet_label, estimator = "SAFE-BC", cover = cover_BC)
  )
  ggplot(df, aes(theta, cover, colour = estimator, group = estimator)) +
    geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50") +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = "theta", y = "Empirical coverage", colour = NULL) +
    scale_colour_manual(values = c("PI" = "firebrick", "SAFE-BC" = "steelblue"))
}

plot_rmse <- function(results) {
  df <- bind_rows(
    results %>% transmute(theta, facet_label, estimator = "PI",      rmse = rmse_PI),
    results %>% transmute(theta, facet_label, estimator = "SAFE-BC", rmse = rmse_BC)
  )
  ggplot(df, aes(theta, rmse, colour = estimator, group = estimator)) +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = "theta", y = "RMSE (ln Mhat - ln M)", colour = NULL) +
    scale_colour_manual(values = c("PI" = "firebrick", "SAFE-BC" = "steelblue"))
}

plot_gap_trunc <- function(results) {
  ggplot(results, aes(theta, PI_gap_trunc_rate)) +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = "theta", y = sprintf("Gap truncation rate (Delta<=%g)", EPS_GAP))
}

plot_delta_cap <- function(results) {
  ggplot(results, aes(theta, delta_var_capped_rate)) +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = "theta", y = sprintf("Delta-var cap-hit rate (cap=%g)", DELTA_VAR_CAP))
}

plot_safe_accept <- function(results) {
  ggplot(results, aes(theta, SAFE_kept_rate)) +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = "theta", y = "SAFE acceptance rate (kept/total)")
}

# ------------------- Main ---------------------------------------
if (sys.nframe() == 0) {
  message("Running a parallel demo; increase K and MIN_KEPT_SAFE for HPC-scale runs.")
  
  K_demo   <- as.integer(Sys.getenv("K_DEMO",    "200"))
  MIN_demo <- as.integer(Sys.getenv("MIN_KEPT",  "2000"))
  CH_demo  <- as.integer(Sys.getenv("CHUNK_INIT","2000"))
  N_CORES  <- suppressWarnings(as.integer(Sys.getenv("N_CORES", "")))
  if (is.na(N_CORES) || N_CORES <= 0) N_CORES <- min(parallel::detectCores(), 216L)
  
  results <- run_full_simulation_parallel(
    K = K_demo,
    min_kept = MIN_demo,
    chunk_init = CH_demo,
    chunk_max = CHUNK_MAX_SAFE,
    max_draws = MAX_DRAWS_SAFE,
    patience_noaccept = PATIENCE_SAFE,
    r = r_default,
    n_cores = N_CORES,
    cluster_type = "auto",
    export_SAFE_fun = TRUE,
    safe_fun_path = SAFE_FILE
  )
  
  results <- facet_ordering(results)
  
  stamp <- format(Sys.Date(), "%Y-%m-%d")
  saveRDS(results, file = sprintf("lnM_sim16_SAFE_FIXED_%s.rds", stamp))
  write.csv(results, file = sprintf("lnM_sim16_SAFE_FIXED_%s.csv", stamp), row.names = FALSE)
  
  dir.create("plots", showWarnings = FALSE)
  
  p_bias    <- plot_bias(results)
  p_relb    <- plot_relbias_simple(results)
  p_cov     <- plot_coverage(results)
  p_rmse    <- plot_rmse(results)
  p_trunc   <- plot_gap_trunc(results)
  p_cap     <- plot_delta_cap(results)
  p_accept  <- plot_safe_accept(results)
  
  print(p_bias); print(p_relb); print(p_cov); print(p_rmse)
  print(p_trunc); print(p_cap); print(p_accept)
  
  ggsave(sprintf("plots/lnM_bias_%s.pdf", stamp),        p_bias,   width = 11, height = 8.5)
  ggsave(sprintf("plots/lnM_relbias_%s.pdf", stamp),     p_relb,   width = 11, height = 8.5)
  ggsave(sprintf("plots/lnM_cover_%s.pdf", stamp),       p_cov,    width = 11, height = 8.5)
  ggsave(sprintf("plots/lnM_rmse_%s.pdf", stamp),        p_rmse,   width = 11, height = 8.5)
  ggsave(sprintf("plots/lnM_gaptrunc_%s.pdf", stamp),    p_trunc,  width = 11, height = 8.5)
  ggsave(sprintf("plots/lnM_deltacap_%s.pdf", stamp),    p_cap,    width = 11, height = 8.5)
  ggsave(sprintf("plots/lnM_safeaccept_%s.pdf", stamp),  p_accept, width = 11, height = 8.5)
  
  message(sprintf("Delta-var cap used: %g", DELTA_VAR_CAP))
  message(sprintf("Gap truncation EPS_GAP: %g", EPS_GAP))
  message(sprintf("Mean |bias| PI: %.4f",      mean(abs(results$bias_PI), na.rm = TRUE)))
  message(sprintf("Mean |bias| SAFE-BC: %.4f", mean(abs(results$bias_BC), na.rm = TRUE)))
  message(sprintf("Mean |rb_delta| (vs MC(PI)): %.2f%%", mean(abs(results$rb_delta), na.rm = TRUE)))
  message(sprintf("Mean |rb_SAFE|  (vs MC(PI)): %.2f%%", mean(abs(results$rb_SAFE),  na.rm = TRUE)))
  message(sprintf("Mean PI_valid_rate: %.3f", mean(results$PI_valid_rate, na.rm = TRUE)))
  message(sprintf("Mean PI_gap_trunc_rate: %.3f", mean(results$PI_gap_trunc_rate, na.rm = TRUE)))
  message(sprintf("Mean delta_var_capped_rate: %.3f", mean(results$delta_var_capped_rate, na.rm = TRUE)))
  message(sprintf("Mean SAFE_kept_rate: %.3f; SAFE ok: %.3f",
                  mean(results$SAFE_kept_rate, na.rm = TRUE),
                  mean(results$SAFE_status_ok, na.rm = TRUE)))
}

# ------------------- How to run (examples) -----------------------
# 1) Linux/macOS (fork):  N_CORES=48 K_DEMO=2000 MIN_KEPT=2000 CHUNK_INIT=4000 Rscript simulation16_parallel_SAFE_with_plots_FIXED.R
# 2) Windows (PSOCK):     N_CORES=48 K_DEMO=2000 MIN_KEPT=2000 CHUNK_INIT=4000 Rscript simulation16_parallel_SAFE_with_plots_FIXED.R
# 3) Accuracy/HPC:        N_CORES=200 K_DEMO=10000 MIN_KEPT=100000 CHUNK_INIT=4000 Rscript simulation16_parallel_SAFE_with_plots_FIXED.R