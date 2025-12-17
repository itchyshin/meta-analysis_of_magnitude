## ------------------------------------------------------------
## Packages
## ------------------------------------------------------------
suppressPackageStartupMessages({
  library(MASS)      # mvrnorm not strictly needed now, but useful later
})

## ------------------------------------------------------------
## lnM (delta-1, independent groups) from *summary* stats
##  - Point estimate + delta-method sampling variance
## ------------------------------------------------------------
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h   <- n1 * n2 / (n1 + n2)
  MSB <- h * (x1bar - x2bar)^2
  MSW <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  
  Delta <- MSB - MSW
  if (!is.finite(Delta) || Delta <= 0) {
    return(c(point = NA_real_, var = NA_real_))
  }
  
  ## lnM (Eq. 4; n0 = 2h for independent designs)
  lnM <- 0.5 * (log(Delta / (2 * h)) - log(MSW))
  
  ## delta-method variance (Eq. 6)
  sigmaD2 <- s1^2 / n1 + s2^2 / n2
  delta   <- x1bar - x2bar
  Var_B   <- h^2 * (2 * sigmaD2^2 + 4 * sigmaD2 * delta^2)
  Var_W   <- 2 * MSW^2 / (n1 + n2 - 2)
  
  gB   <- 0.5 / Delta
  gW   <- -0.5 * MSB / (Delta * MSW)
  Var1 <- gB^2 * Var_B + gW^2 * Var_W
  
  if (!is.finite(Var1) || Var1 <= 0) Var1 <- NA_real_
  
  c(point = lnM, var = Var1)
}

## ------------------------------------------------------------
## "True" lnM for a given delta (no among-study heterogeneity)
##  (Normal-theory expected MSs, independent design)
## ------------------------------------------------------------
true_lnM_indep <- function(mu1, mu2, sd1, sd2, n1, n2) {
  h    <- n1 * n2 / (n1 + n2)
  n0   <- 2 * h
  delta <- mu1 - mu2
  
  omega <- ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  Delta_true <- h * delta^2
  
  0.5 * (log(Delta_true / n0) - log(omega))
}

## ------------------------------------------------------------
## One meta-analysis replicate: generate raw data and pool lnM
##  - K_study studies in one meta-analysis
##  - delta_i ~ Normal(mu_delta, tau_delta^2)
##  - Each study: simulate raw x1, x2; compute lnM + var
##  - Fit meta-analysis with:
##      (i) inverse-variance weights
##     (ii) harmonic-n0 weights (multiplicative method)
## ------------------------------------------------------------
simulate_meta_once <- function(K_study,
                               n1, n2,
                               mu_delta, tau_delta,
                               sigma   = 1,
                               mu_base = 0) {
  
  ## Study-specific true deltas
  delta_i <- rnorm(K_study, mean = mu_delta, sd = tau_delta)
  
  y  <- vi <- n0 <- rep(NA_real_, K_study)
  
  for (k in seq_len(K_study)) {
    mu1 <- mu_base
    mu2 <- mu_base + delta_i[k]
    
    ## Raw data for this study
    x1 <- rnorm(n1, mean = mu1, sd = sigma)
    x2 <- rnorm(n2, mean = mu2, sd = sigma)
    
    m1 <- mean(x1)
    m2 <- mean(x2)
    s1 <- stats::sd(x1)
    s2 <- stats::sd(x2)
    
    est <- lnM_delta1_indep(m1, m2, s1, s2, n1, n2)
    
    y[k]  <- est["point"]
    vi[k] <- est["var"]
    
    h      <- n1 * n2 / (n1 + n2)
    n0[k]  <- 2 * h        # harmonic n0 used in lnM definition
  }
  
  ## Keep only valid lnM estimates
  keep <- is.finite(y) & is.finite(vi) & (vi > 0)
  if (!any(keep)) {
    return(list(
      mu_hat_iv    = NA_real_,
      mu_hat_n0    = NA_real_,
      tau2_iv      = NA_real_,
      tau2_n0      = NA_real_,
      invalid_meta = TRUE
    ))
  }
  
  y  <- y[keep]
  vi <- vi[keep]
  n0 <- n0[keep]
  
  ## Optional: mild trimming of extreme vi to avoid insane weight ratios
  ## (keeps simulation stable but you can comment this out if you prefer)
  # q_lo <- stats::quantile(vi, 0.01)
  # q_hi <- stats::quantile(vi, 0.99)
  # vi   <- pmin(pmax(vi, q_lo), q_hi)
  
  ## Fixed-effect pooling (safe baseline)
  w_iv <- 1 / vi
  w_n0 <- n0
  
  mu_hat_iv <- sum(w_iv * y) / sum(w_iv)
  mu_hat_n0 <- sum(w_n0 * y) / sum(w_n0)
  tau2_iv   <- NA_real_
  tau2_n0   <- NA_real_
  
  ## Optional: REML random-effects via metafor, but robust to failures
  if (requireNamespace("metafor", quietly = TRUE)) {
    ## IV-REML
    re_iv <- try(
      suppressWarnings(
        metafor::rma.uni(yi = y, vi = vi, method = "REML")
      ),
      silent = TRUE
    )
    if (!inherits(re_iv, "try-error")) {
      mu_hat_iv <- as.numeric(re_iv$b[1])
      tau2_iv   <- re_iv$tau2
    }
    
    ## Multiplicative weights (n0) in metafor
    re_n0 <- try(
      suppressWarnings(
        metafor::rma.uni(yi = y, vi = vi, method = "REML",
                         weights = n0)
      ),
      silent = TRUE
    )
    if (!inherits(re_n0, "try-error")) {
      mu_hat_n0 <- as.numeric(re_n0$b[1])
      tau2_n0   <- re_n0$tau2
    }
  }
  
  list(
    mu_hat_iv    = mu_hat_iv,
    mu_hat_n0    = mu_hat_n0,
    tau2_iv      = tau2_iv,
    tau2_n0      = tau2_n0,
    invalid_meta = FALSE
  )
}

## ------------------------------------------------------------
## Top-level simulation over many meta-analyses
##  - Compares IV vs n0 (multiplicative) weighting
## ------------------------------------------------------------
run_weighting_sim <- function(R_meta   = 1000,
                              K_study  = 20,
                              n1       = 10,
                              n2       = 10,
                              delta_mu = 0.8,
                              tau_delta = 0.0,   # >0 for heterogeneity
                              sigma    = 1,
                              mu_base  = 0) {
  
  ## Approximate target: lnM at the mean delta
  lnM_true0 <- true_lnM_indep(mu1 = mu_base,
                              mu2 = mu_base + delta_mu,
                              sd1 = sigma, sd2 = sigma,
                              n1 = n1, n2 = n2)
  
  mu_iv   <- mu_n0 <- tau2_iv <- tau2_n0 <- rep(NA_real_, R_meta)
  invalid <- logical(R_meta)
  
  for (r in seq_len(R_meta)) {
    sim <- simulate_meta_once(K_study  = K_study,
                              n1       = n1,
                              n2       = n2,
                              mu_delta = delta_mu,
                              tau_delta = tau_delta,
                              sigma    = sigma,
                              mu_base  = mu_base)
    mu_iv[r]   <- sim$mu_hat_iv
    mu_n0[r]   <- sim$mu_hat_n0
    tau2_iv[r] <- sim$tau2_iv
    tau2_n0[r] <- sim$tau2_n0
    invalid[r] <- sim$invalid_meta
  }
  
  ok <- !invalid & is.finite(mu_iv) & is.finite(mu_n0)
  if (!any(ok)) {
    warning("All meta-analysis replicates were invalid (no usable lnM).")
    return(data.frame(
      lnM_true0         = lnM_true0,
      bias_iv           = NA_real_,
      bias_n0           = NA_real_,
      rmse_iv           = NA_real_,
      rmse_n0           = NA_real_,
      tau2_iv_mean      = NA_real_,
      tau2_n0_mean      = NA_real_,
      invalid_meta_rate = 1
    ))
  }
  
  mu_iv   <- mu_iv[ok]
  mu_n0   <- mu_n0[ok]
  tau2_iv <- tau2_iv[ok]
  tau2_n0 <- tau2_n0[ok]
  
  bias_iv <- mean(mu_iv) - lnM_true0
  bias_n0 <- mean(mu_n0) - lnM_true0
  rmse_iv <- sqrt(mean((mu_iv  - lnM_true0)^2))
  rmse_n0 <- sqrt(mean((mu_n0  - lnM_true0)^2))
  
  data.frame(
    lnM_true0         = lnM_true0,
    bias_iv           = bias_iv,
    bias_n0           = bias_n0,
    rmse_iv           = rmse_iv,
    rmse_n0           = rmse_n0,
    tau2_iv_mean      = mean(tau2_iv,  na.rm = TRUE),
    tau2_n0_mean      = mean(tau2_n0,  na.rm = TRUE),
    invalid_meta_rate = mean(invalid)
  )
}

## ------------------------------------------------------------
## Small pilot run
## ------------------------------------------------------------
set.seed(123)

demo_res <- run_weighting_sim(
  R_meta    = 200,   # number of meta-analyses
  K_study   = 20,    # studies per meta-analysis
  n1        = 10,
  n2        = 10,
  delta_mu  = 0.8,
  tau_delta = 0.4    # try 0 first, then add heterogeneity
)

print(demo_res)