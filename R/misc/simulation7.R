########################################################################
##  lnM  – simulation, summary & graphics (delta vs SAFE)
##  now with SAFE-median & SAFE-mode point estimators
########################################################################
library(MASS)       # mvrnorm()
library(parallel)   # mclapply / parLapply
library(ggplot2)    # plots

## ---------- 0. helpers ------------------------------------------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)      # force > 0
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)    # NA if MSB≤MSW
d_mode   <- function(x) {                              # quick KDE mode
  d <- density(x, adjust = 1)
  d$x[ which.max(d$y) ]
}

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

## ---------- 1. delta-method plug-in  (indep / paired) -----------------
## … unchanged … (keep your lnM_delta_ind / lnM_delta_dep definitions)

## =====================================================================
##  2.  SAFE bootstrap – returns mean / median / mode               <––
## =====================================================================
## -------- 2a. independent groups -------------------------------------
safe_ind <- function(x1bar, x2bar, s1, s2, n1, n2,
                     B = 1e4, chunk = 5e3)
{
  mu  <- c(x1bar, x2bar,  s1^2, s2^2)
  Sig <- diag(c(s1^2/n1,  s2^2/n2,
                2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  
  cloud <- numeric()
  kept  <- 0L
  while (length(cloud) < B) {
    d   <- mvrnorm(chunk, mu, Sig)
    okV <- d[,3] > 0 & d[,4] > 0
    d   <- d[okV, , drop = FALSE]
    if (!nrow(d)) next
    
    m1  <- d[,1];  m2 <- d[,2];  v1 <- d[,3];  v2 <- d[,4]
    MSB <- (n1*n2)/(n1+n2) * (m1 - m2)^2
    MSW <- ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2)
    
    ok   <- MSB > MSW
    kept <- kept + sum(ok)
    if (!any(ok)) next
    
    Delta <- MSB[ok] - MSW[ok]
    cloud <- c(cloud,
               0.5*(log(Delta) -
                      log(2*n1*n2/(n1+n2)) - log(MSW[ok])))
  }
  cloud <- cloud[seq_len(B)]
  list(pt    = mean(cloud),
       med   = median(cloud),          ### NEW-MED
       mode  = d_mode(cloud),          ### NEW-MODE
       var   = var(cloud),
       lo    = quantile(cloud, 0.025),
       hi    = quantile(cloud, 0.975),
       kept  = kept,
       total = B)
}

## -------- 2b. paired groups ------------------------------------------
safe_dep <- function(x1bar, x2bar, s1, s2, n, rho,
                     B = 1e4, chunk = 5e3)
{
  mu  <- c(x1bar, x2bar,  s1^2, s2^2)
  Sig <- matrix(0, 4, 4)
  Sig[1,1] <- s1^2/n;  Sig[2,2] <- s2^2/n
  Sig[1,2] <- Sig[2,1] <- rho*s1*s2/n
  Sig[3,3] <- 2*s1^4/(n-1);  Sig[4,4] <- 2*s2^4/(n-1)
  Sig[3,4] <- Sig[4,3] <- 2*rho^2*s1^2*s2^2/(n-1)
  
  cloud <- numeric(); kept <- 0L
  while (length(cloud) < B) {
    d   <- mvrnorm(chunk, mu, Sig)
    okV <- d[,3] > 0 & d[,4] > 0
    d   <- d[okV, , drop = FALSE]
    if (!nrow(d)) next
    
    m1  <- d[,1];  m2 <- d[,2];  v1 <- d[,3];  v2 <- d[,4]
    MSB <- (n/2) * (m1 - m2)^2
    MSW <- (v1 + v2)/2
    
    ok   <- MSB > MSW
    kept <- kept + sum(ok)
    if (!any(ok)) next
    
    Delta <- MSB[ok] - MSW[ok]
    cloud <- c(cloud,
               0.5*(log(Delta) - log(n) - log(MSW[ok])))
  }
  cloud <- cloud[seq_len(B)]
  list(pt    = mean(cloud),
       med   = median(cloud),          ### NEW-MED
       mode  = d_mode(cloud),          ### NEW-MODE
       var   = var(cloud),
       lo    = quantile(cloud, 0.025),
       hi    = quantile(cloud, 0.975),
       kept  = kept,
       total = B)
}

## ---------- 3. one replicate  (collect mean / med / mode) ------------
one_rep <- function(mu1, mu2, sd1, sd2,
                    n1, n2 = NULL, rho = 0, B = 1e4)
{
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
  
  c(delta_pt   = delta["pt"],
    delta_var  = delta["var"],
    safe_mean  = safe$pt,
    safe_med   = safe$med,            ### NEW-MED
    safe_mode  = safe$mode,           ### NEW-MODE
    safe_var   = safe$var,
    safe_lo    = safe$lo,
    safe_hi    = safe$hi,
    safe_kept  = safe$kept,
    safe_total = safe$total)
}

## ---------- 4. parameter grid (unchanged) -----------------------------

## ---------- 4. parameter grid  ----------------------------------------
theta_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                1, 1.5, 2, 3, 4, 5)

## 8 independent designs: four balanced + four un-balanced
pairs_ind  <- data.frame(
  n1 = c(5, 10, 20, 100,   3,  6, 12, 60),
  n2 = c(5, 10, 20, 100,   7, 14, 28,140)
)

## build grid for independent designs
grid_ind <- expand.grid(theta = theta_vals,
                        idx    = seq_len(nrow(pairs_ind)),
                        KEEP.OUT.ATTRS = FALSE)

grid_ind$n1     <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2     <- pairs_ind$n2[grid_ind$idx]
grid_ind$rho    <- 0
grid_ind$design <- "indep"
grid_ind$idx    <- NULL     # drop helper column

## four paired (balanced) designs
grid_dep <- expand.grid(theta = theta_vals,
                        n     = c(5, 10, 20, 100),
                        KEEP.OUT.ATTRS = FALSE)

grid_dep$rho    <- 0.8
grid_dep$design <- "paired"
grid_dep$n1     <- grid_dep$n
grid_dep$n2     <- grid_dep$n
grid_dep$n      <- NULL      # drop helper column
grid_dep        <- grid_dep[, names(grid_ind)]  # reorder cols

## final 168-row grid
param_grid <- rbind(grid_ind, grid_dep)
stopifnot(nrow(param_grid) == 168)

## ---------- 5. simulation driver (robust + extra estimators) ---------
set.seed(20250625)
K_repl  <- 1000
B_boot  <- 1e4
n_cores <- max(1, detectCores() - 1)
use_fork <- (.Platform$OS.type != "windows")

runner <- function(i) {
  p   <- param_grid[i, ]
  sd0 <- 1
  
  true_lnM <- if (p$design == "indep")
    lnM_true_ind(p$theta, p$n1, p$n2, sd0)
  else
    lnM_true_dep(p$theta, p$n1, sd0)
  
  reps <- replicate(
    K_repl,
    if (p$design == "indep")
      one_rep(0, p$theta, sd0, sd0,
              n1 = p$n1, n2 = p$n2, B = B_boot)
    else
      one_rep(0, p$theta, sd0, sd0,
              n1 = p$n1, rho = p$rho, B = B_boot),
    simplify = FALSE
  )
  M <- do.call(cbind, reps)
  rownames(M) <- c("delta_pt","delta_var",
                   "safe_mean","safe_med","safe_mode",
                   "safe_var","safe_lo","safe_hi",
                   "safe_kept","safe_total")
  
  true_var_delta <- var(M["delta_pt", ], na.rm = TRUE)
  true_var_safe  <- var(M["safe_mean",], na.rm = TRUE) # same cloud
  
  ok   <- !is.na(M["delta_pt", ])
  Mok  <- M[, ok, drop = FALSE]
  
  ## ------ biases / RMSE for 3 SAFE point estimators -------------------
  safe_stats <- function(rowname) {
    x <- Mok[rowname, ]
    c(bias = mean(x) - true_lnM,
      rmse = sqrt(mean((x - true_lnM)^2)))
  }
  b_mean <- safe_stats("safe_mean")
  b_med  <- safe_stats("safe_med")
  b_mode <- safe_stats("safe_mode")
  
  ## ----- variance summaries (unchanged) ------------------------------
  trim_var_delta <- mean(Mok["delta_var", ], trim = 0.10)
  trim_var_safe  <- mean(Mok["safe_var",  ], trim = 0.10)
  relbias_delta  <- 100*(trim_var_delta/true_var_delta - 1)
  relbias_safe   <- 100*(trim_var_safe /true_var_safe  - 1)
  
  cover_delta <- mean(abs(Mok["delta_pt", ] - true_lnM) <=
                        1.96*sqrt(Mok["delta_var", ]), na.rm = TRUE)
  cover_safe  <- mean(abs(Mok["safe_mean",] - true_lnM) <=
                        1.96*sqrt(Mok["safe_var", ]),  na.rm = TRUE)
  
  data.frame(
    theta = p$theta, design = p$design, n1 = p$n1, n2 = p$n2, rho = p$rho,
    true_lnM = true_lnM,
    
    ## Δ-method
    delta_bias = mean(Mok["delta_pt", ]) - true_lnM,
    rmse_delta = sqrt(mean((Mok["delta_pt", ] - true_lnM)^2)),
    
    ## SAFE – mean / median / mode
    safe_bias_mean = b_mean["bias"],  rmse_mean = b_mean["rmse"],
    safe_bias_med  = b_med["bias"],   rmse_med  = b_med["rmse"],
    safe_bias_mode = b_mode["bias"],  rmse_mode = b_mode["rmse"],
    
    ## variance diagnostics
    relbias_delta = relbias_delta,
    relbias_safe  = relbias_safe,
    
    cover_delta = cover_delta,
    cover_safe  = cover_safe,
    
    boot_keep = sum(Mok["safe_kept", ]),
    boot_total= sum(Mok["safe_total", ])
  )
}

## ---------- run in parallel ------------------------------------------
if (use_fork) {
  res_list <- mclapply(seq_len(nrow(param_grid)),
                       runner, mc.cores = n_cores)
} else {
  res_list <- lapply(seq_len(nrow(param_grid)), runner)
}
results <- do.call(rbind, res_list)
write.csv(results, "lnM_sim_results_SAFEmean-med-mode.csv", row.names = FALSE)

## ---------- 6. graphics ----------------------------------------------
# results should contain at least:
#   theta, design, n1, n2,
#   delta_bias,
#   safe_bias_mean,
#   safe_bias_med,
#   safe_bias_mode

library(dplyr)
library(ggplot2)
library(tidyverse)
## 1. reshape into long form -------------------------------------------------

bias_df <- results %>% 
  select(theta, design, n1, n2,
         delta_bias,
         safe_bias_mean,
         safe_bias_med,
         safe_bias_mode) %>%
  pivot_longer(
    cols = c(delta_bias, safe_bias_mean, safe_bias_med, safe_bias_mode),
    names_to  = "estimator",
    values_to = "bias"
  ) %>%
  mutate(
    estimator = factor(estimator,
                       levels = c("delta_bias",
                                  "safe_bias_mean",
                                  "safe_bias_med",
                                  "safe_bias_mode"),
                       labels = c("Δ-method",
                                  "SAFE (mean)",
                                  "SAFE (median)",
                                  "SAFE (mode)")
    )
  )

## 2. plot --------------------------------------------------------------------

p_bias_all <- ggplot(bias_df,
                     aes(x = theta, y = bias,
                         colour   = estimator,
                         linetype = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ paste0(design, " (n1=", n1, ", n2=", n2, ")"),
             ncol = 4, scales = "free_y") +
  scale_colour_manual(values = c(
    "Δ-method"      = "firebrick",
    "SAFE (mean)"   = "steelblue",
    "SAFE (median)" = "darkgreen",
    "SAFE (mode)"   = "purple"
  )) +
  scale_linetype_manual(values = c(
    "Δ-method"      = "solid",
    "SAFE (mean)"   = "dashed",
    "SAFE (median)" = "dotdash",
    "SAFE (mode)"   = "twodash"
  )) +
  labs(x        = expression(theta),
       y        = "Bias (estimate − true lnM)",
       colour   = NULL,
       linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_bias_all
# ggsave("fig_bias_all_estimators.pdf", p_bias_all, width = 9, height = 6)
