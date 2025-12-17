## ============================================================
## lnM weighting simulation with SAFE (Delta-1 + BC) per study
## Varying sample sizes per study, compare:
##   (1) IV meta-analysis: rma.uni with vi = SAFE var
##   (2) Multiplicative method (harmonic n0-based):
##       rma.mv with V = 0, R = list(ID2 = Vf) from n0
##   (3) Unweighted meta-analysis: lm(yi ~ 1)
##   (4) n0–weighted meta-analysis: lm(yi ~ 1, weights = n0)
## ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(metafor)
  library(tibble)
  library(MASS)
})

## ------------------- Core lnM + SAFE (independent) -----------------

posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(gap) ifelse(gap <= 0, NA_real_, gap)

## Delta-1 lnM (independent)
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2)
{
  h    <- n1 * n2 / (n1 + n2)                         # harmonic component
  MSB  <- h * (x1bar - x2bar)^2                       # Eq. (1)
  MSW  <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)   # Eq. (2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA, var = NA, se = NA))
  
  ## lnM: Eq. (4); n0 = 2h for independent case
  lnM <- 0.5 * (log(Delta / (2 * h)) - log(MSW))
  
  ## Delta-1 variance (for completeness; we mostly use SAFE var)
  sigmaD2 <- s1^2 / n1 + s2^2 / n2
  delta   <- x1bar - x2bar
  Var_B   <- h^2 * (2 * sigmaD2^2 + 4 * sigmaD2 * delta^2)
  Var_W   <- 2 * MSW^2 / (n1 + n2 - 2)
  
  gB   <- 0.5 / Delta
  gW   <- -0.5 * MSB / (Delta * MSW)
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

## SAFE bootstrap for lnM (independent)
safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           B = 1e4, chunk = 4e3, max_chunks = 50)
{
  df1 <- n1 - 1L
  df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)
  
  lnM_star <- numeric(0)
  total    <- 0L
  kept     <- 0L
  attempts <- 0L
  
  while (length(lnM_star) < B && attempts < max_chunks) {
    attempts <- attempts + 1L
    
    m1 <- rnorm(chunk, mean = x1bar, sd = s1 / sqrt(n1))
    m2 <- rnorm(chunk, mean = x2bar, sd = s2 / sqrt(n2))
    v1 <- s1^2 * stats::rchisq(chunk, df = df1) / df1
    v2 <- s2^2 * stats::rchisq(chunk, df = df2) / df2
    
    total <- total + chunk
    
    MSB  <- h * (m1 - m2)^2
    MSW  <- ((df1) * v1 + (df2) * v2) / (df1 + df2)
    good <- MSB > MSW
    
    if (any(good)) {
      kept <- kept + sum(good)
      lnM_star <- c(
        lnM_star,
        0.5 * (log((MSB[good] - MSW[good]) / (2 * h)) -
                 log(MSW[good]))
      )
    }
  }
  
  if (length(lnM_star) > 0L) {
    lnM_star <- lnM_star[seq_len(min(B, length(lnM_star)))]
  }
  
  list(
    point    = if (length(lnM_star)) mean(lnM_star) else NA_real_,
    var      = if (length(lnM_star)) stats::var(lnM_star) else NA_real_,
    kept     = kept,
    total    = total,
    attempts = attempts,
    status   = if (length(lnM_star) >= B) "ok" else "stopped_early"
  )
}

## ------------------- Sampling & "truth" ----------------------------

draw_summaries_indep <- function(mu1, mu2, sd1, sd2, n1, n2) {
  x1bar <- rnorm(1L, mu1, sd1 / sqrt(n1))
  x2bar <- rnorm(1L, mu2, sd2 / sqrt(n2))
  s1sq  <- sd1^2 * stats::rchisq(1L, df = n1 - 1) / (n1 - 1)
  s2sq  <- sd2^2 * stats::rchisq(1L, df = n2 - 1) / (n2 - 1)
  c(x1bar = x1bar, x2bar = x2bar, s1 = sqrt(s1sq), s2 = sqrt(s2sq))
}

true_lnM_indep <- function(mu1, mu2, sd1, sd2, n1, n2, tau_delta = 0) {
  h          <- n1 * n2 / (n1 + n2)
  n0         <- 2 * h
  delta_mean <- mu1 - mu2
  omega      <- ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  Delta_true <- h * (delta_mean^2 + tau_delta^2)
  0.5 * (log(Delta_true / n0) - log(omega))
}

estimate_lnM_true_MC_SAFE <- function(B_true   = 20000,
                                      n1_mean, n2_mean,
                                      delta_mu, tau_delta,
                                      mu1 = 0, sd1 = 1, sd2 = 1,
                                      B_SAFE = 2000, chunk_SAFE = 4000,
                                      max_chunks = 50) {
  lnM_vals <- numeric(B_true)
  for (b in seq_len(B_true)) {
    n1_b <- pmax(3L, rpois(1L, lambda = n1_mean))
    n2_b <- pmax(3L, rpois(1L, lambda = n2_mean))
    delta_b <- if (tau_delta > 0) {
      rnorm(1L, mean = delta_mu, sd = tau_delta)
    } else {
      delta_mu
    }
    mu2_b <- mu1 + delta_b
    sm_b  <- draw_summaries_indep(mu1, mu2_b, sd1, sd2, n1_b, n2_b)
    SAFE  <- safe_lnM_indep(sm_b["x1bar"], sm_b["x2bar"],
                            sm_b["s1"], sm_b["s2"],
                            n1_b, n2_b,
                            B = B_SAFE, chunk = chunk_SAFE,
                            max_chunks = max_chunks)
    lnM_vals[b] <- SAFE$point
  }
  mean(lnM_vals, na.rm = TRUE)
}

## ------------------- One meta-analysis replicate -----------------

simulate_meta_once_mv_SAFE <- function(K_study,
                                       n1, n2,
                                       delta_mu,
                                       tau_delta,
                                       mu1 = 0,
                                       sd1 = 1, sd2 = 1,
                                       B_SAFE = 2000,
                                       chunk_SAFE = 4000,
                                       max_chunks = 50) {
  
  ## Per-study sample sizes (Poisson around means, min 3)
  n1_k <- pmax(3L, rpois(K_study, lambda = n1))
  n2_k <- pmax(3L, rpois(K_study, lambda = n2))
  
  ## True study-level deltas
  delta_k <- if (tau_delta > 0) {
    rnorm(K_study, mean = delta_mu, sd = tau_delta)
  } else {
    rep(delta_mu, K_study)
  }
  mu2_k <- mu1 + delta_k
  
  lnM_BC_vec <- numeric(K_study)
  vi_vec     <- numeric(K_study)
  n0_vec     <- numeric(K_study)
  
  for (k in seq_len(K_study)) {
    sm <- draw_summaries_indep(mu1, mu2_k[k], sd1, sd2, n1_k[k], n2_k[k])
    
    d1 <- lnM_delta1_indep(
      x1bar = sm["x1bar"],
      x2bar = sm["x2bar"],
      s1    = sm["s1"],
      s2    = sm["s2"],
      n1    = n1_k[k],
      n2    = n2_k[k]
    )
    
    SAFE <- safe_lnM_indep(
      x1bar = sm["x1bar"],
      x2bar = sm["x2bar"],
      s1    = sm["s1"],
      s2    = sm["s2"],
      n1    = n1_k[k],
      n2    = n2_k[k],
      B     = B_SAFE,
      chunk = chunk_SAFE,
      max_chunks = max_chunks
    )
    
    ## Bias-corrected lnM:
    ## lnM_BC = 2 * lnM_delta1 - SAFE_mean (if both defined),
    ## otherwise fall back to SAFE mean.
    if (is.na(d1["point"]) || is.na(SAFE$point)) {
      lnM_BC_vec[k] <- SAFE$point
    } else {
      lnM_BC_vec[k] <- 2 * d1["point"] - SAFE$point
    }
    vi_vec[k] <- SAFE$var
    
    h_k       <- n1_k[k] * n2_k[k] / (n1_k[k] + n2_k[k])
    n0_vec[k] <- 2 * h_k
  }
  
  ## Keep only usable studies
  good <- is.finite(lnM_BC_vec) & is.finite(vi_vec) & (vi_vec > 0)
  yi   <- lnM_BC_vec[good]
  vi   <- vi_vec[good]
  n0   <- n0_vec[good]
  
  m <- length(yi)
  if (m < 4) {
    return(list(ok = FALSE))
  }
  
  vi <- pmax(vi, 1e-8)
  
  dat <- data.frame(
    yi  = yi,
    vi  = vi,
    n0  = n0,
    ID2 = factor(seq_len(m))
  )
  
  ## -------- (1) IV meta-analysis: rma.uni with vi = SAFE var ----
  fit_iv <- tryCatch(
    metafor::rma.uni(
      yi     = dat$yi,
      vi     = dat$vi,
      method = "REML",
      control = list(stepadj = 0.5, maxiter = 1000)
    ),
    error = function(e) NULL
  )
  if (is.null(fit_iv)) return(list(ok = FALSE))
  
  ## -------- (2) Multiplicative, harmonic-n0-based ---------------
  ## We want weights ∝ n0; in rma.mv + R, weight ∝ 1 / diag(Vf),
  ## so set vtilde ∝ 1 / n0, then rescale to match mean(vi).
  inv_w_n0 <- mean(dat$n0) / dat$n0          # ∝ 1 / n0
  vtilde   <- inv_w_n0
  #vtilde   <- vtilde * (mean(dat$vi) / mean(vtilde))  # scale ~ vi
  
  Vf   <- diag(as.numeric(vtilde))
  levs <- levels(dat$ID2)
  rownames(Vf) <- levs
  colnames(Vf) <- levs
  
  fit_mult <- tryCatch(
    metafor::rma.mv(
      yi     = yi,
      V      = 0,
      random = list(~1 | ID2),
      data   = dat,
      R      = list(ID2 = Vf),
      Rscale = FALSE,
      method = "REML",
      control = list(stepadj = 0.5, maxiter = 1000)
    ),
    error = function(e) NULL
  )
  if (is.null(fit_mult)) return(list(ok = FALSE))
  
  ## -------- (3) Unweighted meta-analysis (simple mean) ----------
  fit_uw <- tryCatch(
    lm(yi ~ 1, data = dat),
    error = function(e) NULL
  )
  if (is.null(fit_uw)) return(list(ok = FALSE))
  
  ## -------- (4) n0–weighted meta-analysis ----------------------
  fit_wn0 <- tryCatch(
    lm(yi ~ 1, weights = n0, data = dat),
    error = function(e) NULL
  )
  if (is.null(fit_wn0)) return(list(ok = FALSE))
  
  list(
    ok           = TRUE,
    mu_hat_iv    = as.numeric(fit_iv$b[1]),
    mu_hat_n0    = as.numeric(fit_mult$b[1]),
    mu_hat_uw    = as.numeric(coef(fit_uw)[1]),
    mu_hat_wn0   = as.numeric(coef(fit_wn0)[1]),
    tau2_iv      = as.numeric(fit_iv$tau2),
    tau2_n0      = sum(fit_mult$sigma2),
    m_effects    = m
  )
}

## ------------------- Simulation driver ----------------------------

run_weighting_sim_mv_SAFE <- function(R_meta    = 200,
                                      K_study   = 20,
                                      n1        = 10,
                                      n2        = 10,
                                      delta_mu  = 0.8,
                                      tau_delta = 0.4,
                                      mu1       = 0,
                                      sd1       = 1,
                                      sd2       = 1,
                                      use_MC_truth = FALSE,
                                      B_true       = 20000,
                                      B_SAFE       = 2000,
                                      chunk_SAFE   = 4000,
                                      max_chunks   = 50) {
  
  ## "True" lnM target: analytic or MC-SAFE
  if (use_MC_truth) {
    lnM_true0 <- estimate_lnM_true_MC_SAFE(
      B_true    = B_true,
      n1_mean   = n1,
      n2_mean   = n2,
      delta_mu  = delta_mu,
      tau_delta = tau_delta,
      mu1       = mu1,
      sd1       = sd1,
      sd2       = sd2,
      B_SAFE    = B_SAFE,
      chunk_SAFE = chunk_SAFE,
      max_chunks = max_chunks
    )
  } else {
    lnM_true0 <- true_lnM_indep(
      mu1, mu1 + delta_mu, sd1, sd2,
      n1, n2, tau_delta = tau_delta
    )
  }
  
  out_list <- vector("list", R_meta)
  for (r in seq_len(R_meta)) {
    out_list[[r]] <- simulate_meta_once_mv_SAFE(
      K_study   = K_study,
      n1        = n1,
      n2        = n2,
      delta_mu  = delta_mu,
      tau_delta = tau_delta,
      mu1       = mu1,
      sd1       = sd1,
      sd2       = sd2,
      B_SAFE    = B_SAFE,
      chunk_SAFE = chunk_SAFE,
      max_chunks = max_chunks
    )
  }
  
  ok <- vapply(out_list, function(z) isTRUE(z$ok), logical(1))
  invalid_rate <- 1 - mean(ok)
  
  if (!any(ok)) {
    warning("All meta-analysis replicates were invalid.")
    return(
      tibble(
        K_study           = K_study,
        n1                = n1,
        n2                = n2,
        delta_mu          = delta_mu,
        tau_delta         = tau_delta,
        lnM_true0         = lnM_true0,
        bias_iv           = NA_real_,
        bias_n0           = NA_real_,
        bias_uw           = NA_real_,
        bias_wn0          = NA_real_,
        rmse_iv           = NA_real_,
        rmse_n0           = NA_real_,
        rmse_uw           = NA_real_,
        rmse_wn0          = NA_real_,
        tau2_iv_mean      = NA_real_,
        tau2_n0_mean      = NA_real_,
        invalid_meta_rate = 1
      )
    )
  }
  
  mu_iv   <- vapply(out_list[ok], `[[`, numeric(1), "mu_hat_iv")
  mu_n0   <- vapply(out_list[ok], `[[`, numeric(1), "mu_hat_n0")
  mu_uw   <- vapply(out_list[ok], `[[`, numeric(1), "mu_hat_uw")
  mu_wn0  <- vapply(out_list[ok], `[[`, numeric(1), "mu_hat_wn0")
  tau_iv  <- vapply(out_list[ok], `[[`, numeric(1), "tau2_iv")
  tau_n0  <- vapply(out_list[ok], `[[`, numeric(1), "tau2_n0")
  
  bias_iv  <- mean(mu_iv  - lnM_true0)
  bias_n0  <- mean(mu_n0  - lnM_true0)
  bias_uw  <- mean(mu_uw  - lnM_true0)
  bias_wn0 <- mean(mu_wn0 - lnM_true0)
  
  rmse_iv  <- sqrt(mean((mu_iv  - lnM_true0)^2))
  rmse_n0  <- sqrt(mean((mu_n0  - lnM_true0)^2))
  rmse_uw  <- sqrt(mean((mu_uw  - lnM_true0)^2))
  rmse_wn0 <- sqrt(mean((mu_wn0 - lnM_true0)^2))
  
  tibble(
    K_study           = K_study,
    n1                = n1,
    n2                = n2,
    delta_mu          = delta_mu,
    tau_delta         = tau_delta,
    lnM_true0         = lnM_true0,
    bias_iv           = bias_iv,
    bias_n0           = bias_n0,
    bias_uw           = bias_uw,
    bias_wn0          = bias_wn0,
    rmse_iv           = rmse_iv,
    rmse_n0           = rmse_n0,
    rmse_uw           = rmse_uw,
    rmse_wn0          = rmse_wn0,
    tau2_iv_mean      = mean(tau_iv),
    tau2_n0_mean      = mean(tau_n0),
    invalid_meta_rate = invalid_rate
  )
}

## ------------------- Example pilot runs ----------------------

set.seed(123)

## Example 1: lnM_true ≈ -0.5 region
res_mid <- run_weighting_sim_mv_SAFE(
  R_meta    = 50,
  K_study   = 50,
  n1        = 20,
  n2        = 20,
  delta_mu  = 0.76,
  tau_delta = 0.4,
  use_MC_truth = FALSE,
  B_SAFE       = 1000
)
print(res_mid)

## Example 2: lnM_true ≈ -0.1 region
res_shallow <- run_weighting_sim_mv_SAFE(
  R_meta    = 50,
  K_study   = 50,
  n1        = 20,
  n2        = 20,
  delta_mu  = 2,
  tau_delta = 0.4,
  use_MC_truth = FALSE,
  B_SAFE       = 1000
)
print(res_shallow)


####################
####################

## ============================================================
## Multi-scenario grid for SAFE-based lnM weighting simulation
## ============================================================

run_weighting_sim_mv_SAFE_grid <- function(
    R_meta      = 200,
    K_study_vec = c(20, 50),
    n_vec       = c(10, 20),
    delta_vec   = c(0.7, 1.0, 1.3),
    tau_vec     = c(0.2, 0.4),
    mu1         = 0,
    sd1         = 1,
    sd2         = 1,
    use_MC_truth = FALSE,
    B_true       = 20000,
    B_SAFE       = 1000,
    chunk_SAFE   = 4000,
    max_chunks   = 50
) {
  grid <- expand.grid(
    K_study   = K_study_vec,
    n         = n_vec,
    delta_mu  = delta_vec,
    tau_delta = tau_vec,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  res_list <- vector("list", nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    par_i <- grid[i, ]
    cat("\n-------------------------------------------------\n")
    cat("Scenario", i, "of", nrow(grid), "\n")
    print(par_i)
    
    res_list[[i]] <- run_weighting_sim_mv_SAFE(
      R_meta      = R_meta,
      K_study     = par_i$K_study,
      n1          = par_i$n,
      n2          = par_i$n,
      delta_mu    = par_i$delta_mu,
      tau_delta   = par_i$tau_delta,
      mu1         = mu1,
      sd1         = sd1,
      sd2         = sd2,
      use_MC_truth = use_MC_truth,
      B_true       = B_true,
      B_SAFE       = B_SAFE,
      chunk_SAFE   = chunk_SAFE,
      max_chunks   = max_chunks
    )
    
    res_list[[i]]$scenario_id <- i
  }
  
  out <- do.call(rbind, res_list)
  
  out <- out[, c("scenario_id",
                 "K_study", "n1", "n2",
                 "delta_mu", "tau_delta",
                 "lnM_true0",
                 "bias_iv", "bias_n0", "bias_uw", "bias_wn0",
                 "rmse_iv", "rmse_n0", "rmse_uw", "rmse_wn0",
                 "tau2_iv_mean", "tau2_n0_mean",
                 "invalid_meta_rate")]
  
  rownames(out) <- NULL
  out
}

## ------------------- Run the 24-scenario grid ----------------

set.seed(123)

grid_res <- run_weighting_sim_mv_SAFE_grid(
  R_meta      = 200,
  K_study_vec = c(20, 50),
  n_vec       = c(10, 20),
  delta_vec   = c(0.3, 1.2, 2.5),
  tau_vec     = c(0.2, 0.4),
  use_MC_truth = FALSE,
  B_SAFE       = 1000
)

print(grid_res)

