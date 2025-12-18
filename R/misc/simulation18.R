########################################################################
## lnM simulation: PI vs SAFE_raw vs SAFE_BC vs SAFE_mix
## - SAFE uses Normal for means + Chi-square/Wishart for variances
## - SAFE cloud is conditional on Delta*>0 (only usable draws)
## - SAFE_raw = mean(cloud), SAFE_var = var(cloud)
## - SAFE_BC  = 2*PI - SAFE_raw   (when PI exists)
## - SAFE_mix = smooth blend of raw and BC in small-n & small-signal regime
## - Optional sensitivity analysis over SAFE_mix tuning parameters
##
## Run examples (macOS/Linux):
##   N_CORES=24 K_REPL=300 MIN_KEPT=5000 Rscript sim_lnM_compare_allinone.R
##
## Switch on sensitivity:
##   RUN_SENS=1 N_CORES=24 K_REPL=200 MIN_KEPT=2000 Rscript sim_lnM_compare_allinone.R
########################################################################

suppressPackageStartupMessages({
  library(MASS)
  library(ggplot2)
  library(dplyr)
  library(parallel)
  library(pbapply)
  library(tidyr)
})

theme_set(theme_bw(11))

## -------- 0. globals & helpers ---------------------------------------
maxVar   <- 20
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)

lnM_core <- function(Delta, MSW, n0) 0.5 * (log(Delta) - log(n0) - log(MSW))

## -------- 0a. true lnM (for bias/RMSE) -------------------------------
lnM_true_ind <- function(theta, n1, n2, sigma = 1) {
  msb <- (n1 * n2) / (n1 + n2) * theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  lnM_core(msb - msw, msw, 2 * n1 * n2 / (n1 + n2))
}
lnM_true_dep <- function(theta, n, sigma = 1) {
  msb <- (n / 2) * theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  lnM_core(msb - msw, msw, n)
}

## -------- 1. Delta-method plug-in (PI) -------------------------------
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h    <- n1 * n2 / (n1 + n2)
  MSB  <- h * (x1bar - x2bar)^2
  MSW  <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA_real_, var = NA_real_, capped = NA_integer_))
  
  pt <- lnM_core(Delta, MSW, 2 * h)
  
  sigmaD2 <- s1^2 / n1 + s2^2 / n2
  dif <- x1bar - x2bar
  vB  <- h^2 * (2 * sigmaD2^2 + 4 * sigmaD2 * dif^2)
  vW  <- 2 * MSW^2 / (n1 + n2 - 2)
  
  g1 <- 0.5 / Delta
  g2 <- -0.5 * MSB / (Delta * MSW)
  var_raw <- posify(g1^2 * vB + g2^2 * vW)
  
  capped <- as.integer(is.finite(var_raw) && var_raw > maxVar)
  c(pt = pt, var = pmin(var_raw, maxVar), capped = capped)
}

lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r) {
  h    <- n / 2
  MSB  <- h * (x1bar - x2bar)^2
  MSW  <- (s1^2 + s2^2) / 2
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA_real_, var = NA_real_, capped = NA_integer_))
  
  pt <- lnM_core(Delta, MSW, n)
  
  sigmaD2 <- s1^2 + s2^2 - 2 * r * s1 * s2
  dif <- x1bar - x2bar
  vB  <- h^2 * (2 * sigmaD2^2 / n^2 + 4 * dif^2 * sigmaD2 / n)
  vW  <- (s1^4 + s2^4 + 2 * r^2 * s1^2 * s2^2) / (2 * (n - 1))
  
  g1 <- 0.5 / Delta
  g2 <- -0.5 * MSB / (Delta * MSW)
  var_raw <- posify(g1^2 * vB + g2^2 * vW)
  
  capped <- as.integer(is.finite(var_raw) && var_raw > maxVar)
  c(pt = pt, var = pmin(var_raw, maxVar), capped = capped)
}

## -------- 2. SAFE bootstrap (raw mean) -------------------------------
## NOTE: cloud is conditional on MSB*>MSW* (usable draws only)

safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           min_kept = 5000,
                           chunk_init = 4000,
                           chunk_max  = 2e6,
                           max_draws  = Inf,
                           patience_noaccept = 5) {
  
  df1 <- n1 - 1L
  df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)
  n_eff <- 2 * h
  
  MSW0 <- ((df1) * s1^2 + (df2) * s2^2) / (df1 + df2)
  MSB0 <- h * (x1bar - x2bar)^2
  Delta0 <- MSB0 - MSW0
  d_hat <- abs(x1bar - x2bar) / sqrt(MSW0)
  
  lnM_star <- numeric(0L)
  total <- 0L; kept <- 0L; attempts <- 0L
  zero_streak <- 0L
  chunk <- as.integer(chunk_init)
  status <- "ok"
  
  while (kept < min_kept && total < max_draws) {
    attempts <- attempts + 1L
    
    m1 <- rnorm(chunk, mean = x1bar, sd = s1 / sqrt(n1))
    m2 <- rnorm(chunk, mean = x2bar, sd = s2 / sqrt(n2))
    v1 <- s1^2 * rchisq(chunk, df = df1) / df1
    v2 <- s2^2 * rchisq(chunk, df = df2) / df2
    
    total <- total + chunk
    
    MSB  <- h * (m1 - m2)^2
    MSW  <- ((df1) * v1 + (df2) * v2) / (df1 + df2)
    good <- MSB > MSW
    n_good <- sum(good)
    
    if (n_good > 0L) {
      zero_streak <- 0L
      vals <- 0.5 * (log((MSB[good] - MSW[good]) / (2 * h)) - log(MSW[good]))
      lnM_star <- c(lnM_star, vals)
      kept <- length(lnM_star)
    } else {
      zero_streak <- zero_streak + 1L
      if (zero_streak >= patience_noaccept) {
        status <- "no_usable_draws"
        break
      }
    }
    
    acc <- if (total > 0L) kept / total else 0
    remaining <- max(0L, min_kept - kept)
    next_needed <- if (acc > 0) ceiling(remaining / acc) else chunk * 2L
    chunk <- as.integer(max(chunk_init, min(chunk_max, next_needed)))
  }
  
  if (kept < min_kept && total >= max_draws) status <- "hit_max_draws"
  if (kept > min_kept) {
    lnM_star <- lnM_star[seq_len(min_kept)]
    kept <- min_kept
  }
  
  list(
    pt_raw = if (kept) mean(lnM_star) else NA_real_,
    var    = if (kept) var(lnM_star)  else NA_real_,
    kept   = kept,
    total  = total,
    status = status,
    n_eff  = n_eff,
    d_hat  = d_hat,
    Delta0 = Delta0
  )
}

safe_lnM_dep <- function(x1bar, x2bar, s1, s2, n, r,
                         min_kept = 5000,
                         chunk_init = 4000,
                         chunk_max  = 2e6,
                         max_draws  = Inf,
                         patience_noaccept = 5) {
  
  df <- n - 1L
  h  <- n / 2
  n_eff <- n
  
  Sig <- matrix(c(s1^2, r*s1*s2,
                  r*s1*s2, s2^2), 2, 2)
  
  MSW0 <- (s1^2 + s2^2) / 2
  MSB0 <- h * (x1bar - x2bar)^2
  Delta0 <- MSB0 - MSW0
  d_hat <- abs(x1bar - x2bar) / sqrt(MSW0)
  
  lnM_star <- numeric(0L)
  total <- 0L; kept <- 0L; attempts <- 0L
  zero_streak <- 0L
  chunk <- as.integer(chunk_init)
  status <- "ok"
  
  while (kept < min_kept && total < max_draws) {
    attempts <- attempts + 1L
    
    Mu <- MASS::mvrnorm(n = chunk, mu = c(x1bar, x2bar), Sigma = Sig / n)
    W  <- stats::rWishart(n = chunk, df = df, Sigma = Sig)
    S11 <- W[1,1,] / df
    S22 <- W[2,2,] / df
    
    total <- total + chunk
    
    MSB  <- h * (Mu[,1] - Mu[,2])^2
    MSW  <- (S11 + S22) / 2
    good <- MSB > MSW
    n_good <- sum(good)
    
    if (n_good > 0L) {
      zero_streak <- 0L
      vals <- 0.5 * (log((MSB[good] - MSW[good]) / n) - log(MSW[good]))
      lnM_star <- c(lnM_star, vals)
      kept <- length(lnM_star)
    } else {
      zero_streak <- zero_streak + 1L
      if (zero_streak >= patience_noaccept) {
        status <- "no_usable_draws"
        break
      }
    }
    
    acc <- if (total > 0L) kept / total else 0
    remaining <- max(0L, min_kept - kept)
    next_needed <- if (acc > 0) ceiling(remaining / acc) else chunk * 2L
    chunk <- as.integer(max(chunk_init, min(chunk_max, next_needed)))
  }
  
  if (kept < min_kept && total >= max_draws) status <- "hit_max_draws"
  if (kept > min_kept) {
    lnM_star <- lnM_star[seq_len(min_kept)]
    kept <- min_kept
  }
  
  list(
    pt_raw = if (kept) mean(lnM_star) else NA_real_,
    var    = if (kept) var(lnM_star)  else NA_real_,
    kept   = kept,
    total  = total,
    status = status,
    n_eff  = n_eff,
    d_hat  = d_hat,
    Delta0 = Delta0
  )
}

## -------- 2b. SAFE_mix weight (smooth trigger) -------------------------
mix_weight <- function(n_eff, d_hat, n_small = 20, d_small = 0.75,
                       s_n = 3, s_d = 0.15) {
  # w ~ 1 only when (n_eff << n_small) AND (d_hat << d_small)
  w_n <- plogis((n_small - n_eff) / s_n)
  w_d <- plogis((d_small - d_hat) / s_d)
  w_n * w_d
}

## -------- 3. one replicate -------------------------------------------
one_rep <- function(theta, design = c("indep","paired"),
                    n1, n2 = NULL, rho = 0.8,
                    min_kept = 5000,
                    chunk_init = 4000,
                    chunk_max  = 2e6,
                    max_draws  = Inf,
                    patience_noaccept = 5,
                    # mixing controls
                    n_small = 20, d_small = 0.75,
                    s_n = 3, s_d = 0.15) {
  
  design <- match.arg(design)
  
  if (design == "indep") {
    x1 <- rnorm(n1, 0, 1)
    x2 <- rnorm(n2, theta, 1)
    x1bar <- mean(x1); x2bar <- mean(x2)
    s1 <- sd(x1); s2 <- sd(x2)
    
    pi <- lnM_delta1_indep(x1bar, x2bar, s1, s2, n1, n2)
    safe <- safe_lnM_indep(x1bar, x2bar, s1, s2, n1, n2,
                           min_kept, chunk_init, chunk_max, max_draws, patience_noaccept)
  } else {
    Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    xy <- MASS::mvrnorm(n1, mu = c(0, theta), Sigma = Sigma)
    x1 <- xy[,1]; x2 <- xy[,2]
    x1bar <- mean(x1); x2bar <- mean(x2)
    s1 <- sd(x1); s2 <- sd(x2)
    rho_hat <- suppressWarnings(cor(x1, x2))
    if (!is.finite(rho_hat)) rho_hat <- rho
    
    pi <- lnM_delta1_dep(x1bar, x2bar, s1, s2, n1, rho_hat)
    safe <- safe_lnM_dep(x1bar, x2bar, s1, s2, n1, rho_hat,
                         min_kept, chunk_init, chunk_max, max_draws, patience_noaccept)
  }
  
  pt_pi <- unname(pi["pt"])
  pt_safe_raw <- safe$pt_raw
  pt_safe_bc <- if (is.finite(pt_pi) && is.finite(pt_safe_raw)) 2 * pt_pi - pt_safe_raw else NA_real_
  
  w <- if (is.finite(safe$n_eff) && is.finite(safe$d_hat) &&
           is.finite(pt_safe_raw) && is.finite(pt_safe_bc)) {
    mix_weight(safe$n_eff, safe$d_hat, n_small, d_small, s_n, s_d)
  } else NA_real_
  
  pt_safe_mix <- if (is.finite(w)) (1 - w) * pt_safe_raw + w * pt_safe_bc else NA_real_
  
  c(
    PI_pt = pt_pi,
    PI_var = unname(pi["var"]),
    SAFE_raw_pt = pt_safe_raw,
    SAFE_BC_pt  = pt_safe_bc,
    SAFE_mix_pt = pt_safe_mix,
    SAFE_var    = safe$var,
    w_mix       = w,
    kept        = safe$kept,
    tried       = safe$total,
    safe_ok     = as.integer(isTRUE(safe$status == "ok"))
  )
}

## -------- 4. parameter grid ------------------------------------------
theta_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4, 5)

pairs_ind <- data.frame(n1 = c(5, 10, 20, 100, 3, 6, 12, 40),
                        n2 = c(5, 10, 20, 100, 7, 14, 28, 160))

grid_ind <- expand.grid(theta = theta_vals, idx = seq_len(nrow(pairs_ind)))
grid_ind$n1 <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2 <- pairs_ind$n2[grid_ind$idx]
grid_ind$design <- "indep"
grid_ind$idx <- NULL

grid_dep <- expand.grid(theta = theta_vals, n = c(5, 10, 20, 100))
grid_dep$n1 <- grid_dep$n
grid_dep$n2 <- grid_dep$n
grid_dep$design <- "paired"
grid_dep$n <- NULL

param_grid <- rbind(grid_ind, grid_dep)

## -------- 5. simulation driver ---------------------------------------
set.seed(20250625)

K_REPL     <- as.integer(Sys.getenv("K_REPL", "300"))
MIN_KEPT   <- as.integer(Sys.getenv("MIN_KEPT", "5000"))
CHUNK_INIT <- as.integer(Sys.getenv("CHUNK_INIT", "4000"))
CHUNK_MAX  <- as.numeric(Sys.getenv("CHUNK_MAX", "2000000"))
MAX_DRAWS  <- as.numeric(Sys.getenv("MAX_DRAWS", "Inf"))
PATIENCE   <- as.integer(Sys.getenv("PATIENCE", "5"))
N_CORES    <- suppressWarnings(as.integer(Sys.getenv("N_CORES", "")))

N_SMALL <- as.numeric(Sys.getenv("N_SMALL", "20"))
D_SMALL <- as.numeric(Sys.getenv("D_SMALL", "0.75"))
S_N     <- as.numeric(Sys.getenv("S_N", "3"))
S_D     <- as.numeric(Sys.getenv("S_D", "0.15"))

mcse_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sqrt(var(x) / length(x))
}

runner <- function(i) {
  p <- param_grid[i, ]
  
  true_ln <- if (p$design == "indep") lnM_true_ind(p$theta, p$n1, p$n2, 1) else lnM_true_dep(p$theta, p$n1, 1)
  
  M <- matrix(NA_real_, 10, K_REPL,
              dimnames = list(
                c("PI_pt","PI_var",
                  "SAFE_raw_pt","SAFE_BC_pt","SAFE_mix_pt","SAFE_var",
                  "w_mix","kept","tried","safe_ok"),
                NULL))
  
  for (k in seq_len(K_REPL)) {
    M[, k] <- one_rep(theta = p$theta,
                      design = p$design,
                      n1 = p$n1,
                      n2 = if (p$design == "indep") p$n2 else NULL,
                      rho = 0.8,
                      min_kept = MIN_KEPT,
                      chunk_init = CHUNK_INIT,
                      chunk_max = CHUNK_MAX,
                      max_draws = MAX_DRAWS,
                      patience_noaccept = PATIENCE,
                      n_small = N_SMALL, d_small = D_SMALL,
                      s_n = S_N, s_d = S_D)
  }
  
  # baseline for Var relative bias: MC Var of SAFE_raw point estimator
  ok_raw <- is.finite(M["SAFE_raw_pt", ])
  Var_MC_SAFEraw <- if (sum(ok_raw) >= 2) var(M["SAFE_raw_pt", ok_raw]) else NA_real_
  
  mean_var_PI   <- mean(M["PI_var", ], na.rm = TRUE)
  mean_var_SAFE <- mean(M["SAFE_var", ], na.rm = TRUE)
  
  relbias_var_PI <- if (is.finite(Var_MC_SAFEraw) && Var_MC_SAFEraw > 0)
    100 * (mean_var_PI / Var_MC_SAFEraw - 1) else NA_real_
  
  relbias_var_SAFE <- if (is.finite(Var_MC_SAFEraw) && Var_MC_SAFEraw > 0)
    100 * (mean_var_SAFE / Var_MC_SAFEraw - 1) else NA_real_
  
  bias <- function(x) if (is.finite(true_ln)) mean(x, na.rm = TRUE) - true_ln else NA_real_
  rmse <- function(x) if (is.finite(true_ln)) sqrt(mean((x - true_ln)^2, na.rm = TRUE)) else NA_real_
  
  data.frame(
    theta = p$theta,
    design = p$design,
    n1 = p$n1,
    n2 = ifelse(p$design == "indep", p$n2, p$n1),
    true_lnM = true_ln,
    
    bias_PI       = bias(M["PI_pt", ]),
    bias_SAFE_raw = bias(M["SAFE_raw_pt", ]),
    bias_SAFE_BC  = bias(M["SAFE_BC_pt", ]),
    bias_SAFE_mix = bias(M["SAFE_mix_pt", ]),
    
    rmse_PI       = rmse(M["PI_pt", ]),
    rmse_SAFE_raw = rmse(M["SAFE_raw_pt", ]),
    rmse_SAFE_BC  = rmse(M["SAFE_BC_pt", ]),
    rmse_SAFE_mix = rmse(M["SAFE_mix_pt", ]),
    
    Var_MC_SAFEraw   = Var_MC_SAFEraw,
    mean_var_PI      = mean_var_PI,
    mean_var_SAFE    = mean_var_SAFE,
    relbias_var_PI   = relbias_var_PI,
    relbias_var_SAFE = relbias_var_SAFE,
    
    accept_prop = sum(M["kept", ], na.rm = TRUE) / max(1, sum(M["tried", ], na.rm = TRUE)),
    SAFE_ok_rate = mean(M["safe_ok", ], na.rm = TRUE),
    mean_w_mix   = mean(M["w_mix", ], na.rm = TRUE),
    
    mcse_bias_PI       = if (is.finite(true_ln)) mcse_mean(M["PI_pt", ]       - true_ln) else NA_real_,
    mcse_bias_SAFE_raw = if (is.finite(true_ln)) mcse_mean(M["SAFE_raw_pt", ] - true_ln) else NA_real_,
    mcse_bias_SAFE_BC  = if (is.finite(true_ln)) mcse_mean(M["SAFE_BC_pt", ]  - true_ln) else NA_real_,
    mcse_bias_SAFE_mix = if (is.finite(true_ln)) mcse_mean(M["SAFE_mix_pt", ] - true_ln) else NA_real_
  )
}

# parallel setup
if (is.na(N_CORES) || N_CORES <= 0) N_CORES <- max(1L, detectCores() - 2)
pbop <- pbapply::pboptions(type = "txt")

if (.Platform$OS.type == "windows") {
  cl <- makeCluster(N_CORES)
  clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)
  res_list <- pbapply::pblapply(seq_len(nrow(param_grid)), runner, cl = cl)
  stopCluster(cl)
} else {
  res_list <- pbapply::pblapply(seq_len(nrow(param_grid)), runner, cl = N_CORES)
}
pbapply::pboptions(pbop)

results <- do.call(rbind, res_list)

## -------- 6. plots ----------------------------------------------------
results <- results %>% mutate(facet_label = paste0(design, " n1=", n1, " n2=", n2))

# facet ordering: balanced indep, unbalanced indep, paired
facet_info <- results %>% distinct(facet_label, design, n1, n2)
balanced_ind <- facet_info %>% filter(design == "indep", n1 == n2) %>% arrange(n1) %>% pull(facet_label)
unbalanced_ind <- facet_info %>% filter(design == "indep", n1 != n2) %>% arrange(n1) %>% pull(facet_label)
paired_balanced <- facet_info %>% filter(design == "paired") %>% arrange(n1) %>% pull(facet_label)
new_levels <- c(balanced_ind, unbalanced_ind, paired_balanced)
results <- results %>% mutate(facet_label = factor(facet_label, levels = new_levels))

# Bias curves
bias_df <- results %>%
  select(theta, facet_label, starts_with("bias_")) %>%
  pivot_longer(cols = starts_with("bias_"),
               names_to = "estimator", values_to = "bias") %>%
  mutate(estimator = recode(estimator,
                            bias_PI = "PI",
                            bias_SAFE_raw = "SAFE_raw",
                            bias_SAFE_BC = "SAFE_BC",
                            bias_SAFE_mix = "SAFE_mix"))

p_bias <- ggplot(bias_df, aes(theta, bias, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
  labs(x = expression(theta), y = "Bias (estimate - true lnM)", colour = "Estimator")
print(p_bias)

# RMSE curves
rmse_df <- results %>%
  select(theta, facet_label, starts_with("rmse_")) %>%
  pivot_longer(cols = starts_with("rmse_"),
               names_to = "estimator", values_to = "rmse") %>%
  mutate(estimator = recode(estimator,
                            rmse_PI = "PI",
                            rmse_SAFE_raw = "SAFE_raw",
                            rmse_SAFE_BC = "SAFE_BC",
                            rmse_SAFE_mix = "SAFE_mix"))

p_rmse <- ggplot(rmse_df, aes(theta, rmse, colour = estimator, group = estimator)) +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
  labs(x = expression(theta), y = "RMSE", colour = "Estimator")
print(p_rmse)

# Relative bias of variance estimates, baseline = MC Var(SAFE_raw)
rb_df <- results %>%
  select(theta, facet_label, relbias_var_PI, relbias_var_SAFE) %>%
  pivot_longer(cols = c(relbias_var_PI, relbias_var_SAFE),
               names_to = "estimator", values_to = "relbias") %>%
  mutate(estimator = recode(estimator,
                            relbias_var_PI = "PI_var",
                            relbias_var_SAFE = "SAFE_var"))

p_relbias <- ggplot(rb_df, aes(theta, relbias, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
  labs(x = expression(theta),
       y = "Relative bias of Var (%) vs MC Var(SAFE_raw)",
       colour = "Variance estimator")
print(p_relbias)

# Quick diagnostics
message("Mean accept_prop (kept/tried): ", mean(results$accept_prop, na.rm = TRUE))
message("Mean SAFE_ok_rate: ", mean(results$SAFE_ok_rate, na.rm = TRUE))
message("Mean mixing weight (SAFE_mix): ", mean(results$mean_w_mix, na.rm = TRUE))

## -------- 7. OPTIONAL sensitivity analysis for SAFE_mix tuning ---------
RUN_SENS <- as.integer(Sys.getenv("RUN_SENS", "0"))

if (RUN_SENS == 1) {
  
  sens_grid <- expand.grid(
    N_SMALL = c(15, 20, 25),
    D_SMALL = c(0.6, 0.75, 0.9),
    S_N     = c(2.5, 3.0, 4.0),
    S_D     = c(0.10, 0.15, 0.20)
  )
  
  theta_weight <- function(theta) 1 / pmax(theta, 0.4)
  
  run_one_setting <- function(j) {
    Ns <- sens_grid$N_SMALL[j]
    Ds <- sens_grid$D_SMALL[j]
    sn <- sens_grid$S_N[j]
    sd <- sens_grid$S_D[j]
    
    runner_sens <- function(i) {
      p <- param_grid[i, ]
      true_ln <- if (p$design == "indep") lnM_true_ind(p$theta, p$n1, p$n2, 1) else lnM_true_dep(p$theta, p$n1, 1)
      
      M <- matrix(NA_real_, 10, K_REPL,
                  dimnames = list(
                    c("PI_pt","PI_var",
                      "SAFE_raw_pt","SAFE_BC_pt","SAFE_mix_pt","SAFE_var",
                      "w_mix","kept","tried","safe_ok"),
                    NULL))
      
      for (k in seq_len(K_REPL)) {
        M[, k] <- one_rep(theta = p$theta,
                          design = p$design,
                          n1 = p$n1,
                          n2 = if (p$design == "indep") p$n2 else NULL,
                          rho = 0.8,
                          min_kept = MIN_KEPT,
                          chunk_init = CHUNK_INIT,
                          chunk_max = CHUNK_MAX,
                          max_draws = MAX_DRAWS,
                          patience_noaccept = PATIENCE,
                          n_small = Ns, d_small = Ds,
                          s_n = sn, s_d = sd)
      }
      
      bias <- function(x) if (is.finite(true_ln)) mean(x, na.rm = TRUE) - true_ln else NA_real_
      rmse <- function(x) if (is.finite(true_ln)) sqrt(mean((x - true_ln)^2, na.rm = TRUE)) else NA_real_
      
      data.frame(
        theta = p$theta, design = p$design, n1 = p$n1,
        n2 = ifelse(p$design == "indep", p$n2, p$n1),
        bias_PI       = bias(M["PI_pt", ]),
        bias_SAFE_raw = bias(M["SAFE_raw_pt", ]),
        bias_SAFE_BC  = bias(M["SAFE_BC_pt", ]),
        bias_SAFE_mix = bias(M["SAFE_mix_pt", ]),
        rmse_PI       = rmse(M["PI_pt", ]),
        rmse_SAFE_raw = rmse(M["SAFE_raw_pt", ]),
        rmse_SAFE_BC  = rmse(M["SAFE_BC_pt", ]),
        rmse_SAFE_mix = rmse(M["SAFE_mix_pt", ]),
        mean_w_mix    = mean(M["w_mix", ], na.rm = TRUE),
        accept_prop   = sum(M["kept", ], na.rm = TRUE) / max(1, sum(M["tried", ], na.rm = TRUE))
      )
    }
    
    res_list_j <- pbapply::pblapply(seq_len(nrow(param_grid)), runner_sens, cl = N_CORES)
    Rj <- do.call(rbind, res_list_j)
    Rj$w_theta <- theta_weight(Rj$theta)
    
    pooled_rmse <- function(col) {
      ok <- is.finite(Rj[[col]]) & is.finite(Rj$w_theta)
      sqrt(sum(Rj$w_theta[ok] * (Rj[[col]][ok]^2)) / sum(Rj$w_theta[ok]))
    }
    
    worst_abs_bias <- function(col) {
      tmp <- Rj %>%
        mutate(facet = paste0(design," n1=",n1," n2=",n2),
               ab = abs(.data[[col]])) %>%
        group_by(facet) %>%
        summarise(worst = max(ab, na.rm = TRUE), .groups = "drop")
      mean(tmp$worst, na.rm = TRUE)
    }
    
    wins <- Rj %>%
      mutate(
        best = pmap_chr(
          list(rmse_PI, rmse_SAFE_raw, rmse_SAFE_BC, rmse_SAFE_mix),
          ~{
            v <- c(PI=..1, SAFE_raw=..2, SAFE_BC=..3, SAFE_mix=..4)
            if (all(!is.finite(v))) return(NA_character_)
            names(which.min(v))
          })
      ) %>%
      count(best) %>%
      tidyr::complete(best = c("PI","SAFE_raw","SAFE_BC","SAFE_mix"), fill = list(n = 0)) %>%
      mutate(prop = n / sum(n))
    
    out <- data.frame(
      N_SMALL=Ns, D_SMALL=Ds, S_N=sn, S_D=sd,
      pooled_RMSE_PI       = pooled_rmse("rmse_PI"),
      pooled_RMSE_SAFE_raw = pooled_rmse("rmse_SAFE_raw"),
      pooled_RMSE_SAFE_BC  = pooled_rmse("rmse_SAFE_BC"),
      pooled_RMSE_SAFE_mix = pooled_rmse("rmse_SAFE_mix"),
      worstAbsBias_PI       = worst_abs_bias("bias_PI"),
      worstAbsBias_SAFE_raw = worst_abs_bias("bias_SAFE_raw"),
      worstAbsBias_SAFE_BC  = worst_abs_bias("bias_SAFE_BC"),
      worstAbsBias_SAFE_mix = worst_abs_bias("bias_SAFE_mix"),
      mean_w_mix = mean(Rj$mean_w_mix, na.rm = TRUE),
      mean_accept = mean(Rj$accept_prop, na.rm = TRUE),
      winProp_PI       = wins$prop[wins$best=="PI"],
      winProp_SAFE_raw = wins$prop[wins$best=="SAFE_raw"],
      winProp_SAFE_BC  = wins$prop[wins$best=="SAFE_BC"],
      winProp_SAFE_mix = wins$prop[wins$best=="SAFE_mix"]
    )
    
    out
  }
  
  sens_out <- do.call(rbind, lapply(seq_len(nrow(sens_grid)), run_one_setting))
  sens_out <- sens_out %>% arrange(pooled_RMSE_SAFE_mix)
  
  print(head(sens_out, 15), row.names = FALSE)
}

########################################################################
## Why SAFE_mix is often the best (in this problem)
##
## 1) SAFE_raw is a conditional mean E[lnM* | Delta*>0]. When theta is small
##    and n is small, Delta*>0 is a “rare-event” filter that selects draws
##    with unusually large between-group signal, causing upward bias.
##
## 2) SAFE_BC (2*PI - SAFE_raw) can reduce that bias near the boundary,
##    but it can overshoot (become negatively biased), because PI itself is
##    noisy and the SAFE cloud is truncated/conditional.
##
## 3) SAFE_mix uses a smooth weight w(n_eff, d_hat) to apply BC only where
##    the mechanism in (1) is strongest (small n_eff AND small signal d_hat),
##    while reverting to SAFE_raw elsewhere. This typically improves the
##    bias/RMSE trade-off: you get boundary de-biasing without paying an
##    RMSE penalty in the moderate/large-signal regime.
########################################################################