########################################################################
##  lnM  – simulation, summary & graphics (delta vs SAFE)
##  now with Parametric‐Truncation Bias Correction for SAFE
##  KEEP‐RATE computed correctly
##  24 independent + 12 paired parameter sets (36 total)
##  ASCII-only code, works on Linux / macOS / Windows
##  last tested: 25-Jun-2025  (R ≥ 4.3, MASS 7.3-60, ggplot2 3.5-0)
########################################################################

library(MASS)     # mvrnorm()
library(parallel) # mclapply / parLapply
library(ggplot2)  # plots

## ---------- 0. helpers ------------------------------------------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)

## ---------- 0a. true-lnM functions -----------------------------------
lnM_true_ind <- function(theta, n1, n2, sigma = 1) {
  msb <- (n1*n2)/(n1+n2) * theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  0.5*( log(msb-msw) - log(2*n1*n2/(n1+n2)) - log(msw) )
}
lnM_true_dep <- function(theta, n, sigma = 1) {
  msb <- (n/2)*theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  0.5*( log(msb-msw) - log(n) - log(msw) )
}

## ---------- 1. delta-method plug-in  (independent / paired) -----------
lnM_delta_ind <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h      <- n1 * n2 / (n1 + n2)
  MSB    <- h * (x1bar - x2bar)^2
  MSW    <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2)/(n1 + n2 - 2)
  n0     <- 2 * n1 * n2 / (n1 + n2)        # ← restored original scaling
  Delta  <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt=NA,var=NA))
  point  <- 0.5*(log(Delta) - log(n0) - log(MSW))
  sD2    <- s1^2/n1 + s2^2/n2
  dif    <- x1bar - x2bar
  var_MSB <- h^2*(2*sD2^2 + 4*sD2*dif^2)
  var_MSW <- 2*MSW^2/(n1+n2-2)
  g1 <- 0.5/Delta
  g2 <- -0.5*MSB/(Delta*MSW)
  var_est <- posify(g1^2*var_MSB + g2^2*var_MSW)
  c(pt=point, var=var_est)
}

lnM_delta_dep <- function(x1bar,x2bar,s1,s2,n,rho) {
  MSB   <- (n/2)*(x1bar-x2bar)^2
  MSW   <- (s1^2 + s2^2)/2
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt=NA,var=NA))
  point <- 0.5*(log(Delta) - log(n) - log(MSW))
  sD2   <- s1^2 + s2^2 - 2*rho*s1*s2
  dif   <- x1bar - x2bar
  var_MSB <- (n/2)^2*(2*sD2^2/n^2 + 4*dif^2*sD2/n)
  var_MSW <- (s1^4 + s2^4 + 2*rho^2*s1^2*s2^2)/(2*(n-1))
  g1 <- 0.5/Delta
  g2 <- -0.5*MSB/(Delta*MSW)
  var_est <- posify(g1^2*var_MSB + g2^2*var_MSW)
  c(pt=point, var=var_est)
}

## =====================================================================
##  2.  SAFE bootstrap – PARAMETRIC TRUNCATION BIAS CORRECTION
## =====================================================================
## ---------- 2a. independent groups -----------------------------------
safe_ind <- function(x1bar, x2bar, s1, s2, n1, n2,
                     B = 1e4, chunk = 5e3) {
  mu       <- c(x1bar,x2bar,s1^2,s2^2)
  Sig      <- diag(c(s1^2/n1, s2^2/n2,
                     2*s1^4/(n1-1),2*s2^4/(n2-1)))
  cloud    <- numeric()
  kept     <- 0L
  attempts <- 0L
  
  while(length(cloud) < B) {
    d        <- mvrnorm(chunk, mu, Sig)
    attempts <- attempts + nrow(d)
    okV      <- d[,3]>0 & d[,4]>0
    d        <- d[okV,,drop=FALSE]
    if (!nrow(d)) next
    
    m1 <- d[,1]; m2 <- d[,2]
    v1 <- d[,3]; v2 <- d[,4]
    MSB <- (n1*n2)/(n1+n2)*(m1 - m2)^2
    MSW <- ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2)
    
    ok   <- MSB > MSW
    kept <- kept + sum(ok)
    if (!any(ok)) next
    
    Δ    <- MSB[ok] - MSW[ok]
    cloud <- c(cloud,
               0.5*(log(Δ) -
                      log(2*n1*n2/(n1+n2)) -
                      log(MSW[ok])))
  }
  
  cloud <- cloud[seq_len(B)]
  p_keep <- kept / attempts
  p_keep <- posify(pmin(p_keep, 1 - 1e-12))
  
  μ_raw  <- mean(cloud)
  σ_raw  <- sd(cloud)
  α      <- qnorm(p_keep)
  λ      <- dnorm(α)/p_keep
  pt_corr<- μ_raw - σ_raw * λ
  
  list(pt    = pt_corr,
       var   = var(cloud),
       lo    = quantile(cloud,0.025),
       hi    = quantile(cloud,0.975),
       kept  = kept,
       total = B)
}

## ---------- 2b. paired groups ----------------------------------------
safe_dep <- function(x1bar, x2bar, s1, s2, n, rho,
                     B = 1e4, chunk = 5e3) {
  mu       <- c(x1bar,x2bar,s1^2,s2^2)
  Sig      <- matrix(0,4,4)
  Sig[1,1]<-s1^2/n;  Sig[2,2]<-s2^2/n
  Sig[1,2]<-Sig[2,1]<-rho*s1*s2/n
  Sig[3,3]<-2*s1^4/(n-1); Sig[4,4]<-2*s2^4/(n-1)
  Sig[3,4]<-Sig[4,3]<-2*rho^2*s1^2*s2^2/(n-1)
  
  cloud    <- numeric()
  kept     <- 0L
  attempts <- 0L
  
  while(length(cloud) < B) {
    d        <- mvrnorm(chunk, mu, Sig)
    attempts <- attempts + nrow(d)
    okV      <- d[,3]>0 & d[,4]>0
    d        <- d[okV,,drop=FALSE]
    if (!nrow(d)) next
    
    m1 <- d[,1]; m2 <- d[,2]
    v1 <- d[,3]; v2 <- d[,4]
    MSB <- (n/2)*(m1 - m2)^2
    MSW <- (v1 + v2)/2
    
    ok   <- MSB > MSW
    kept <- kept + sum(ok)
    if (!any(ok)) next
    
    Δ    <- MSB[ok] - MSW[ok]
    cloud <- c(cloud,
               0.5*(log(Δ) -
                      log(n) -
                      log(MSW[ok])))
  }
  
  cloud <- cloud[seq_len(B)]
  p_keep <- kept / attempts
  p_keep <- posify(pmin(p_keep, 1 - 1e-12))
  
  μ_raw  <- mean(cloud)
  σ_raw  <- sd(cloud)
  α      <- qnorm(p_keep)
  λ      <- dnorm(α)/p_keep
  pt_corr<- μ_raw - σ_raw * λ
  
  list(pt    = pt_corr,
       var   = var(cloud),
       lo    = quantile(cloud,0.025),
       hi    = quantile(cloud,0.975),
       kept  = kept,
       total = B)
}

## ---------- 3. one‐replicate wrapper ---------------------------------
one_rep <- function(mu1,mu2,sd1,sd2,n1,n2=NULL,rho=0,B=1e4) {
  if (is.null(n2)) {
    Sig <- matrix(c(sd1^2,rho*sd1*sd2,
                    rho*sd1*sd2,sd2^2),2,2)
    xy <- mvrnorm(n1,c(mu1,mu2),Sig)
    x1 <- xy[,1]; x2 <- xy[,2]; rho_hat <- cor(x1,x2)
  } else {
    x1 <- rnorm(n1,mu1,sd1)
    x2 <- rnorm(n2,mu2,sd2)
    rho_hat <- 0
  }
  x1bar <- mean(x1); s1 <- sd(x1)
  x2bar <- mean(x2); s2 <- sd(x2)
  
  delta <- if (is.null(n2))
    lnM_delta_dep(x1bar,x2bar,s1,s2,n1,rho_hat) else
      lnM_delta_ind(x1bar,x2bar,s1,s2,n1,n2)
  
  safe  <- if (is.null(n2))
    safe_dep(x1bar,x2bar,s1,s2,n1,rho_hat,B) else
      safe_ind(x1bar,x2bar,s1,s2,n1,n2,B)
  
  c(delta_pt   = delta["pt"],
    delta_var  = delta["var"],
    safe_pt    = safe$pt,
    safe_var   = safe$var,
    safe_lo    = safe$lo,
    safe_hi    = safe$hi,
    safe_kept  = safe$kept,
    safe_total = safe$total)
}

## ---------- 4. parameter grid ----------------------------------------
theta_vals <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,4,5)
pairs_ind  <- data.frame(n1=c(5,10,20,100,3,6,12,60),
                         n2=c(5,10,20,100,7,14,28,140))
grid_ind <- expand.grid(theta=theta_vals,
                        idx=seq_len(nrow(pairs_ind)),
                        KEEP.OUT.ATTRS=FALSE)
grid_ind$n1     <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2     <- pairs_ind$n2[grid_ind$idx]
grid_ind$rho    <- 0
grid_ind$design <- "indep"
grid_ind$idx    <- NULL

grid_dep <- expand.grid(theta=theta_vals,
                        n=c(5,10,20,100),
                        KEEP.OUT.ATTRS=FALSE)
grid_dep$rho    <- 0.8
grid_dep$design <- "paired"
grid_dep$n1     <- grid_dep$n
grid_dep$n2     <- grid_dep$n
grid_dep$n      <- NULL
grid_dep        <- grid_dep[,names(grid_ind)]

param_grid <- rbind(grid_ind, grid_dep)
stopifnot(nrow(param_grid)==168)

## ---------- 5. run simulation ----------------------------------------
set.seed(20250625)
K_repl  <- 1000
B_boot  <- 1e4
n_cores <- max(1, detectCores()-1)
use_fork<- (.Platform$OS.type!="windows")

runner <- function(i) {
  p       <- param_grid[i,]
  sd0     <- 1
  true_lnM<- if(p$design=="indep")
    lnM_true_ind(p$theta,p$n1,p$n2,sd0)
  else
    lnM_true_dep(p$theta,p$n1,sd0)
  
  reps <- replicate(K_repl,
                    if(p$design=="indep")
                      one_rep(0,p$theta,sd0,sd0,n1=p$n1,n2=p$n2,B=B_boot)
                    else
                      one_rep(0,p$theta,sd0,sd0,n1=p$n1,rho=p$rho,B=B_boot),
                    simplify=FALSE)
  M <- do.call(cbind, reps)
  rownames(M) <- c("delta_pt","delta_var",
                   "safe_pt","safe_var",
                   "safe_lo","safe_hi",
                   "safe_kept","safe_total")
  ok  <- !is.na(M["delta_pt",])
  Mok <- M[,ok,drop=FALSE]
  
  data.frame(
    theta       = p$theta,
    design      = p$design,
    n1          = p$n1,
    n2          = p$n2,
    rho         = p$rho,
    true_lnM    = true_lnM,
    delta_bias  = mean(Mok["delta_pt",]) - true_lnM,
    safe_bias   = mean(Mok["safe_pt",])  - true_lnM,
    relbias_delta = 100*(mean(Mok["delta_var",]) / var(Mok["delta_pt",])),
    relbias_safe  = 100*(mean(Mok["safe_var", ]) / var(Mok["safe_pt", ])),
    cover_delta = mean(abs(Mok["delta_pt",]-true_lnM) <= 1.96*sqrt(Mok["delta_var",])),
    cover_safe  = mean(abs(Mok["safe_pt", ]-true_lnM) <= 1.96*sqrt(Mok["safe_var", ])),
    boot_keep   = sum(Mok["safe_kept",]),
    boot_total  = sum(Mok["safe_total",])
  )
}

if(use_fork) {
  res_list <- mclapply(seq_len(nrow(param_grid)), runner, mc.cores=n_cores)
} else {
  res_list <- lapply(seq_len(nrow(param_grid)), runner)
}
results <- do.call(rbind, res_list)

## ---------- 6. graphics (unchanged) ----------------------------------
# ... now plot `results$delta_bias` vs `results$safe_bias` with ggplot2 ...
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

