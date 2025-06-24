########################################################################
##  lnM  – simulation, summary & graphics (delta vs SAFE)
##  24 independent + 12 paired parameter sets  (36 total)
##  ASCII-only code, works on Linux / macOS / Windows
##  last tested: 25-Jun-2025  (R ≥ 4.3, MASS 7.3-60, ggplot2 3.5-0)
########################################################################

library(MASS)       # mvrnorm()
library(parallel)   # mclapply / parLapply
library(ggplot2)    # plots

## ---------- 0. helpers ------------------------------------------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)           # force > 0
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)         # NA if MSB≤MSW

## ---------- 0a. true-lnM functions -----------------------------------
lnM_true_ind <- function(theta, n1, n2, sigma = 1) {
  msb <- (n1*n2)/(n1+n2) * theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  0.5*( log(msb-msw) - log(2*n1*n2/(n1+n2)) - log(msw) )
}

# the same as ind but less computation
lnM_true_dep <- function(theta, n, sigma = 1) {
  msb <- (n/2)*theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  0.5*( log(msb-msw) - log(n) - log(msw) )
}

## ---------- 1. delta-method plug-in (indep / paired) ------------------
## … (unchanged – truncated for brevity in this snippet)
## keep your lnM_delta_ind() and lnM_delta_dep() definitions here

## ---------- 1. delta-method plug-in  (independent / paired) -----------

lnM_delta_ind <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h      <- n1 * n2 / (n1 + n2)
  MSB    <- h * (x1bar - x2bar)^2
  MSW    <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  n0     <- 2 * n1 * n2 / (n1 + n2)
  Delta  <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA, var = NA))
  
  point  <- 0.5 * (log(Delta) - log(n0) - log(MSW))
  
  sD2    <- s1^2 / n1 + s2^2 / n2
  difM   <- x1bar - x2bar
  var_MSB <- h^2 * (2 * sD2^2 + 4 * sD2 * difM^2)
  var_MSW <- 2 * MSW^2 / (n1 + n2 - 2)
  
  g1 <- 0.5 / Delta
  g2 <- -0.5 * MSB / (Delta * MSW)
  var_est <- posify(g1^2 * var_MSB + g2^2 * var_MSW)
  
  c(pt = point, var = var_est)
}

lnM_delta_dep <- function(x1bar, x2bar, s1, s2, n, rho) {
  MSB    <- (n / 2) * (x1bar - x2bar)^2
  MSW    <- (s1^2 + s2^2) / 2
  Delta  <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA, var = NA))
  
  point  <- 0.5 * (log(Delta) - log(n) - log(MSW))
  
  sD2    <- s1^2 + s2^2 - 2 * rho * s1 * s2
  difM   <- x1bar - x2bar
  var_MSB <- (n / 2)^2 * (2 * sD2^2 / n^2 + 4 * difM^2 * sD2 / n)
  var_MSW <- (s1^4 + s2^4 + 2 * rho^2 * s1^2 * s2^2) / (2 * (n - 1))
  
  g1 <- 0.5 / Delta
  g2 <- -0.5 * MSB / (Delta * MSW)
  var_est <- posify(g1^2 * var_MSB + g2^2 * var_MSW)
  
  c(pt = point, var = var_est)
}


## ---------- 2. SAFE bootstrap  (indep / paired) ----------------------
safe_ind <- function(x1bar, x2bar, s1, s2, n1, n2,
                     B = 1e4, chunk = 5e3) {
  
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- diag(c(s1^2/n1, s2^2/n2,
                2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  
  cloud <- numeric(); kept <- 0L
  while (length(cloud) < B) {
    d   <- mvrnorm(chunk, mu, Sig)
    okV <- d[,3] > 0 & d[,4] > 0
    d   <- d[okV,,drop = FALSE]
    if (!nrow(d)) next
    
    m1  <- d[,1];  m2 <- d[,2];  v1 <- d[,3];  v2 <- d[,4]
    MSB <- (n1*n2)/(n1+n2) * (m1 - m2)^2
    MSW <- ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2)
    ok  <- MSB > MSW
    kept <- kept + sum(ok)
    if (!any(ok)) next
    
    cloud <- c(cloud,
               0.5*(log(MSB[ok]-MSW[ok]) -
                      log(2*n1*n2/(n1+n2)) - log(MSW[ok])))
  }
  cloud <- cloud[seq_len(B)]
  list(pt = mean(cloud),
       var= var(cloud),
       lo = quantile(cloud, 0.025),          ### CHG
       hi = quantile(cloud, 0.975),          ### CHG
       kept = kept,
       total= B)
}

safe_dep <- function(x1bar, x2bar, s1, s2, n, rho,
                     B = 1e4, chunk = 5e3) {
  
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- matrix(0, 4, 4)
  Sig[1,1] <- s1^2/n;  Sig[2,2] <- s2^2/n
  Sig[1,2] <- Sig[2,1] <- rho*s1*s2/n
  Sig[3,3] <- 2*s1^4/(n-1);  Sig[4,4] <- 2*s2^4/(n-1)
  Sig[3,4] <- Sig[4,3] <- 2*rho^2*s1^2*s2^2/(n-1)
  
  cloud <- numeric(); kept <- 0L
  while (length(cloud) < B) {
    d   <- mvrnorm(chunk, mu, Sig)
    okV <- d[,3] > 0 & d[,4] > 0
    d   <- d[okV,,drop = FALSE]
    if (!nrow(d)) next
    
    m1  <- d[,1];  m2 <- d[,2];  v1 <- d[,3];  v2 <- d[,4]
    MSB <- (n/2) * (m1 - m2)^2
    MSW <- (v1 + v2)/2
    ok  <- MSB > MSW
    kept <- kept + sum(ok)
    if (!any(ok)) next
    
    cloud <- c(cloud,
               0.5*(log(MSB[ok]-MSW[ok]) - log(n) - log(MSW[ok])))
  }
  cloud <- cloud[seq_len(B)]
  list(pt = mean(cloud),
       var= var(cloud),
       lo = quantile(cloud, 0.025),          ### CHG
       hi = quantile(cloud, 0.975),          ### CHG
       kept = kept,
       total= B)
}

## ---------- 3. one replicate  (now stores lo / hi) --------------------
one_rep <- function(mu1, mu2, sd1, sd2,
                    n1, n2 = NULL, rho = 0, B = 1e4) {
  
  if (is.null(n2)) {
    Sig <- matrix(c(sd1^2, rho*sd1*sd2,
                    rho*sd1*sd2, sd2^2), 2, 2)
    xy  <- mvrnorm(n1, c(mu1, mu2), Sig)
    x1 <- xy[,1]; x2 <- xy[,2]
    rho_hat <- cor(x1, x2)
  } else {
    x1 <- rnorm(n1, mu1, sd1)
    x2 <- rnorm(n2, mu2, sd2)
    rho_hat <- 0
  }
  
  x1bar <- mean(x1); s1 <- sd(x1)
  x2bar <- mean(x2); s2 <- sd(x2)
  
  delta <- if (is.null(n2))
    lnM_delta_dep(x1bar, x2bar, s1, s2, n1, rho_hat)
  else
    lnM_delta_ind(x1bar, x2bar, s1, s2, n1, n2)
  
  safe  <- if (is.null(n2))
    safe_dep(x1bar, x2bar, s1, s2, n1, rho_hat, B)
  else
    safe_ind(x1bar, x2bar, s1, s2, n1, n2, B)
  
  c(delta_pt  = delta["pt"],
    delta_var = delta["var"],
    safe_pt   = safe$pt,
    safe_var  = safe$var,
    safe_lo   = safe$lo,      ### CHG
    safe_hi   = safe$hi,      ### CHG
    safe_kept = safe$kept,
    safe_total= safe$total)
}

## ---------- 4. parameter grid  (same 36 sets as before) ---------------
theta_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 3, 4, 5)
pairs_ind  <- data.frame(n1 = c(5, 10, 20, 100, 3, 6, 12, 60),
                         n2 = c(5, 10, 20, 100, 7, 14, 28, 140))

grid_ind <- expand.grid(theta = theta_vals,
                        idx    = seq_len(nrow(pairs_ind)),
                        KEEP.OUT.ATTRS = FALSE)
grid_ind$n1 <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2 <- pairs_ind$n2[grid_ind$idx]
grid_ind$rho <- 0
grid_ind$design <- "indep"
grid_ind$idx <- NULL

grid_dep <- expand.grid(theta = theta_vals,
                        n     = c(5, 10, 20, 100),
                        KEEP.OUT.ATTRS = FALSE)
grid_dep$rho <- 0.8
grid_dep$design <- "paired"
grid_dep$n1 <- grid_dep$n; grid_dep$n2 <- grid_dep$n; grid_dep$n <- NULL
grid_dep <- grid_dep[, names(grid_ind)]

param_grid <- rbind(grid_ind, grid_dep)
stopifnot(nrow(param_grid) == 168)

## ---------- 5. simulation driver (z-Wald + percentile) ---------------
set.seed(20250625)

K_repl  <- 1000
B_boot  <- 1e4
n_cores <- max(1, detectCores()-1)
use_fork <- (.Platform$OS.type != "windows")

####################################################################
##  --  replace everything from here down to the parallel call  -- ##
####################################################################

runner <- function(i) {
  p   <- param_grid[i, ]
  sd0 <- 1
  
  ## -------- true lnM ----------------------------------------------------
  true_lnM <- if (p$design == "indep")
    lnM_true_ind(p$theta, p$n1, p$n2, sd0) else
      lnM_true_dep(p$theta, p$n1,          sd0)
  
  ## -------- generate K_repl replicates ----------------------------------
  rep_list <- replicate(
    K_repl,
    if (p$design == "indep")
      one_rep(0, p$theta, sd0, sd0,
              n1 = p$n1, n2 = p$n2, B = B_boot) else
                one_rep(0, p$theta, sd0, sd0,
                        n1 = p$n1, rho = p$rho, B = B_boot),
    simplify = FALSE)
  
  M <- do.call(cbind, rep_list)
  rownames(M) <- c("delta_pt","delta_var",
                   "safe_pt","safe_var",
                   "safe_lo","safe_hi",
                   "safe_kept","safe_total")
  
  ok  <- !is.na(M["delta_pt", ])
  Mok <- M[, ok, drop = FALSE]
  
  ## ---------- NEW-1 : Monte-Carlo “true” variances ----------------------
  true_var_delta <- var(Mok["delta_pt", ])
  true_var_safe  <- var(Mok["safe_pt",  ])
  
  ## ---------- coverage (both use z-Wald now) ---------------------------
  cover_delta_z <- mean(abs(Mok["delta_pt", ] - true_lnM) <=
                          1.96 * sqrt(Mok["delta_var", ]), na.rm = TRUE)
  cover_safe_z  <- mean(abs(Mok["safe_pt",  ] - true_lnM) <=
                          1.96 * sqrt(Mok["safe_var", ]),  na.rm = TRUE)
  ## (safe_lo / safe_hi remain in the record but are no longer used) ### OUT
  
  ## ---------- NEW-2 : relative bias of each variance estimator ---------
  relbias_delta <- 100 * (mean(Mok["delta_var", ]) / true_var_delta - 1)
  relbias_safe  <- 100 * (mean(Mok["safe_var",  ]) / true_var_safe  - 1)
  
  ## ---------- summary row ----------------------------------------------
  data.frame(
    theta       = p$theta,
    design      = p$design,
    n1          = p$n1,
    n2          = p$n2,
    rho         = p$rho,
    
    ## true values
    true_lnM        = true_lnM,
    true_var_delta  = true_var_delta,
    true_var_safe   = true_var_safe,
    
    replicates  = ncol(Mok),
    
    ## means & bias of point estimates
    delta_mean  = mean(Mok["delta_pt", ]),
    safe_mean   = mean(Mok["safe_pt",  ]),
    delta_bias  = mean(Mok["delta_pt", ]) - true_lnM,
    safe_bias   = mean(Mok["safe_pt",  ]) - true_lnM,
    
    ## variance estimators
    mean_var_delta = mean(Mok["delta_var", ]),
    mean_var_safe  = mean(Mok["safe_var",  ]),
    relbias_delta  = relbias_delta,    ### NEW-2
    relbias_safe   = relbias_safe,     ### NEW-2
    
    ## RMSE
    rmse_delta = sqrt(mean((Mok["delta_pt", ] - true_lnM)^2)),
    rmse_safe  = sqrt(mean((Mok["safe_pt",  ] - true_lnM)^2)),
    
    ## coverage (z-Wald for both)
    cover_delta = cover_delta_z,       ### NEW-3
    cover_safe  = cover_safe_z,        ### NEW-3
    
    ## bootstrap diagnostics
    boot_keep   = sum(Mok["safe_kept", ]),
    boot_total  = sum(Mok["safe_total", ])
  )
}

## ---------- run in parallel exactly as before --------------------------
if (use_fork) {
  res_list <- mclapply(seq_len(nrow(param_grid)),
                       runner, mc.cores = n_cores)
} else {
  res_list <- lapply(seq_len(nrow(param_grid)), runner)
}

results <- do.call(rbind, res_list)

## save the richer table
write.csv(results, "lnM_sim_results_full.csv", row.names = FALSE)

read.csv("lnM_sim_results_full.csv", stringsAsFactors = FALSE) -> results


## ---------- 6. graphics ----------------------------------------------## =====================================================================
##  A.  Bias of the point estimators  (Δ-method vs SAFE)
## =====================================================================

bias_df <- rbind(
  data.frame(results[, c("theta","design","n1","n2")],
             estimator = "delta",
             bias      = results$delta_bias),
  data.frame(results[, c("theta","design","n1","n2")],
             estimator = "SAFE",
             bias      = results$safe_bias))

p_bias <- ggplot(bias_df,
                 aes(theta, bias,
                     colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line(size = 0.8) +
  facet_wrap(~ paste0(design, "  n1=", n1, "  n2=", n2), ncol = 4) +
  scale_colour_manual(values = c(delta = "firebrick",
                                 SAFE  = "steelblue")) +
  labs(x     = expression(theta),
       y     = "Bias (estimate − true lnM)",
       colour = NULL) +
  theme_bw(11)

p_bias
# ggsave("fig_bias_point.pdf", p_bias, width = 9, height = 5)


## =====================================================================
##  B.  Relative bias of the variance estimators
##      (each plotted separately)
## =====================================================================

relbias_df <- rbind(
  data.frame(results[, c("theta","design","n1","n2")],
             estimator = "delta",
             relbias   = results$relbias_delta),
  data.frame(results[, c("theta","design","n1","n2")],
             estimator = "SAFE",
             relbias   = results$relbias_safe))

p_var <- ggplot(relbias_df,
                aes(theta, relbias,
                    colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line(size = 0.8) +
  facet_wrap(~ paste0(design, "  n1=", n1, "  n2=", n2), ncol = 4) +
  scale_colour_manual(values = c(delta = "firebrick",
                                 SAFE  = "steelblue")) +
  labs(x     = expression(theta),
       y     = "Relative bias of Var  (%)",
       colour = NULL) +
  theme_bw(11)

p_var
# ggsave("fig_bias_variance.pdf", p_var, width = 9, height = 5)


## =====================================================================
##  C.  Coverage of 95 % z-Wald CI (both estimators)
## =====================================================================

cover_df <- rbind(
  data.frame(results[, c("theta","design","n1","n2")],
             estimator = "delta",
             cover     = results$cover_delta),
  data.frame(results[, c("theta","design","n1","n2")],
             estimator = "SAFE",
             cover     = results$cover_safe))

p_cov <- ggplot(cover_df,
                aes(theta, cover,
                    colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50") +
  geom_line(size = 0.8) +
  facet_wrap(~ paste0(design, "  n1=", n1, "  n2=", n2), ncol = 4) +
  scale_colour_manual(values = c(delta = "firebrick",
                                 SAFE  = "steelblue")) +
  labs(x     = expression(theta),
       y     = "Empirical coverage",
       colour = NULL) +
  theme_bw(11)

p_cov
# ggsave("fig_coverage.pdf", p_cov, width = 9, height = 5)


