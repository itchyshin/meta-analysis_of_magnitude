# ================================================================
# simulation16_using_SAFE.R — lnM simulations with NEW SAFE (direct calls)
# - Assumes SAFE_fun.R defines:
#     safe_lnM_indep(x1bar,x2bar,s1,s2,n1,n2,
#                    min_kept, chunk_init, chunk_max, max_draws, patience_noaccept)
#     safe_lnM_dep  (x1bar,x2bar,s1,s2,n,r,
#                    min_kept, chunk_init, chunk_max, max_draws, patience_noaccept)
# - Point estimators: PI (Δ plug-in) and SAFE-BC (2*PI - SAFE$point; fallback to SAFE$point if PI NA)
# - Variance estimators: Δ-method vs SAFE bootstrap
# - Computes Eq. (16)–(20): bias, MC variances, variance relative bias (4 combos), coverage, RMSE
# - Restores “simulation15” plots + adds SAFE acceptance and four-way rel. bias grid
# ================================================================

suppressPackageStartupMessages({
  library(MASS)      # mvrnorm, rWishart
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

# ------------------- Load the latest SAFE -----------------------
SAFE_FILE <- "SAFE_fun.R"
if (!file.exists(SAFE_FILE)) {
  stop("SAFE_fun.R not found in working directory. Please place it next to this script.")
}
source(SAFE_FILE)

if (!exists("safe_lnM_indep") || !exists("safe_lnM_dep")) {
  stop("SAFE_fun.R must define safe_lnM_indep() and safe_lnM_dep().")
}

# ------------------- Global controls / knobs --------------------
SEED_MAIN <- 12345
set.seed(SEED_MAIN)

# Monte Carlo replicates per parameter set (increase on HPC)
K_default <- 2000

# SAFE controls (new API)
MIN_KEPT_SAFE    <- 2000     # target usable lnM* draws per replicate
CHUNK_INIT_SAFE  <- 4000     # starting chunk size
CHUNK_MAX_SAFE   <- 2e6      # cap for adaptive chunk
MAX_DRAWS_SAFE   <- Inf      # safety cap
PATIENCE_SAFE    <- 5        # consecutive zero-accept chunks before early stop

# Cap for Δ-method variance
DELTA_VAR_CAP <- 20

# Paired correlation
r_default <- 0.8

# θ grid and designs
theta_grid <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,4,5)

designs_indep <- list(
  c(5,5), c(10,10), c(20,20), c(100,100),
  c(3,7), c(6,14), c(12,28), c(40,160)
)
designs_paired_n <- c(5,10,20,100)

# ------------------- Helpers -----------------------------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(gap) ifelse(gap <= 0, NA_real_, gap)

# ------------------- lnM Delta-1 (independent) ------------------
# MS Eq. (1)-(4); variance Eq. (6)
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h    <- n1 * n2 / (n1 + n2)
  MSB  <- h * (x1bar - x2bar)^2
  MSW  <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA, var = NA, se = NA))
  
  lnM <- 0.5 * (log(Delta / (2 * h)) - log(MSW))
  
  sigmaD2 <- s1^2 / n1 + s2^2 / n2
  delta   <- x1bar - x2bar
  Var_B   <- h^2 * (2 * sigmaD2^2 + 4 * sigmaD2 * delta^2)
  Var_W   <- 2 * MSW^2 / (n1 + n2 - 2)
  
  gB   <- 0.5 / Delta
  gW   <- -0.5 * MSB / (Delta * MSW)
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  Var1 <- pmin(Var1, DELTA_VAR_CAP)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ------------------- lnM Delta-1 (paired) -----------------------
# Paired forms; variance Eq. (7)
lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r) {
  h     <- n / 2
  MSB   <- h * (x1bar - x2bar)^2
  MSW   <- (s1^2 + s2^2) / 2
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA, var = NA, se = NA))
  
  lnM <- 0.5 * (log(Delta / n) - log(MSW))
  
  sigmaD2 <- s1^2 + s2^2 - 2 * r * s1 * s2
  delta   <- x1bar - x2bar
  Var_B   <- h^2 * (2 * sigmaD2^2 / n^2 + 4 * delta^2 * sigmaD2 / n)
  Var_W   <- (s1^4 + s2^4 + 2 * r^2 * s1^2 * s2^2) / (2 * (n - 1))
  
  gB   <- 0.5 / Delta
  gW   <- -0.5 * MSB / (Delta * MSW)
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  Var1 <- pmin(Var1, DELTA_VAR_CAP)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ------------------- “True” lnM (expected MS) -------------------
true_lnM_indep <- function(mu1, mu2, sd1, sd2, n1, n2) {
  h     <- n1*n2/(n1+n2)
  n0    <- 2*h
  delta <- mu1 - mu2
  omega <- ((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2)
  Delta_true <- h * delta^2
  0.5 * (log(Delta_true / n0) - log(omega))
}

# Paired: Δ_true = (n/2)δ^2 - r σ1 σ2 ; E[MSW] = (σ1^2+σ2^2)/2
true_lnM_dep <- function(mu1, mu2, sd1, sd2, n, r) {
  delta <- mu1 - mu2
  sigmaD2 <- sd1^2 + sd2^2 - 2*r*sd1*sd2
  E_MSW   <- (sd1^2 + sd2^2)/2
  Delta_true <- (n/2)*delta^2 + (1/2)*sigmaD2 - E_MSW  # = (n/2)δ^2 - r σ1 σ2
  n0 <- n
  0.5 * (log(Delta_true / n0) - log(E_MSW))
}

# ------------------- Summary draws (fast, Normal theory) --------
draw_summaries_indep <- function(mu1, mu2, sd1, sd2, n1, n2) {
  x1bar <- rnorm(1L, mu1, sd1/sqrt(n1))
  x2bar <- rnorm(1L, mu2, sd2/sqrt(n2))
  s1sq  <- sd1^2 * rchisq(1L, df = n1-1) / (n1-1)
  s2sq  <- sd2^2 * rchisq(1L, df = n2-1) / (n2-1)
  c(x1bar = x1bar, x2bar = x2bar, s1 = sqrt(s1sq), s2 = sqrt(s2sq))
}

draw_summaries_dep <- function(mu1, mu2, sd1, sd2, n, r) {
  Sig <- matrix(c(sd1^2, r*sd1*sd2,
                  r*sd1*sd2, sd2^2), 2, 2)
  mu  <- MASS::mvrnorm(1L, mu = c(mu1, mu2), Sigma = Sig / n)
  W   <- stats::rWishart(1L, df = n-1, Sigma = Sig)
  S11 <- W[1,1,1] / (n-1)
  S22 <- W[2,2,1] / (n-1)
  c(x1bar = mu[1], x2bar = mu[2], s1 = sqrt(S11), s2 = sqrt(S22))
}

# ------------------- One parameter set: independent --------------
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
  
  lnM_PI     <- numeric(K)
  Var_delta  <- numeric(K)
  lnM_BC     <- numeric(K)
  Var_SAFE   <- numeric(K)
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
    
    lnM_PI[k]     <- unname(d1["point"])
    Var_delta[k]  <- unname(d1["var"])
    kept_vec[k]   <- SAFE$kept
    total_vec[k]  <- SAFE$total
    status_vec[k] <- SAFE$status
    
    # SAFE-BC: 2*PI - SAFE point; fallback to SAFE point if PI undefined
    lnM_BC[k] <- if (is.na(lnM_PI[k])) SAFE$point else 2*lnM_PI[k] - SAFE$point
    Var_SAFE[k] <- SAFE$var
  }
  
  Var_MC_PI <- var(lnM_PI[is.finite(lnM_PI)], na.rm = TRUE)
  Var_MC_BC <- var(lnM_BC[is.finite(lnM_BC)], na.rm = TRUE)
  
  Var_delta_bar <- mean(Var_delta[is.finite(Var_delta)], na.rm = TRUE)
  Var_SAFE_bar  <- mean(Var_SAFE[is.finite(Var_SAFE)],   na.rm = TRUE)
  
  rb_delta_PI <- 100 * (Var_delta_bar - Var_MC_PI) / Var_MC_PI
  rb_delta_BC <- 100 * (Var_delta_bar - Var_MC_BC) / Var_MC_BC
  rb_SAFE_PI  <- 100 * (Var_SAFE_bar  - Var_MC_PI) / Var_MC_PI
  rb_SAFE_BC  <- 100 * (Var_SAFE_bar  - Var_MC_BC) / Var_MC_BC
  
  bias_PI <- mean(lnM_PI, na.rm = TRUE) - lnM_true
  bias_BC <- mean(lnM_BC, na.rm = TRUE) - lnM_true
  
  cover_PI  <- mean(abs(lnM_PI - lnM_true) <= 1.96*sqrt(Var_delta), na.rm = TRUE)
  cover_BC  <- mean(abs(lnM_BC - lnM_true) <= 1.96*sqrt(Var_SAFE),  na.rm = TRUE)
  rmse_PI   <- sqrt(mean((lnM_PI - lnM_true)^2, na.rm = TRUE))
  rmse_BC   <- sqrt(mean((lnM_BC - lnM_true)^2, na.rm = TRUE))
  
  tibble(
    design   = "indep",
    n1 = n1, n2 = n2, r = NA_real_,
    theta = theta,
    lnM_true = lnM_true,
    bias_PI = bias_PI,
    bias_BC = bias_BC,
    Var_MC_PI = Var_MC_PI,
    Var_MC_BC = Var_MC_BC,
    Var_delta_bar = Var_delta_bar,
    Var_SAFE_bar  = Var_SAFE_bar,
    rb_delta_PI = rb_delta_PI,
    rb_delta_BC = rb_delta_BC,
    rb_SAFE_PI  = rb_SAFE_PI,
    rb_SAFE_BC  = rb_SAFE_BC,
    cover_PI = cover_PI,
    cover_BC = cover_BC,
    rmse_PI  = rmse_PI,
    rmse_BC  = rmse_BC,
    SAFE_kept_rate = mean(kept_vec / pmax(1L, total_vec), na.rm = TRUE),
    SAFE_status_ok = mean(status_vec == "ok", na.rm = TRUE)
  )
}

# ------------------- One parameter set: paired -------------------
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
  
  lnM_PI     <- numeric(K)
  Var_delta  <- numeric(K)
  lnM_BC     <- numeric(K)
  Var_SAFE   <- numeric(K)
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
    
    lnM_PI[k]     <- unname(d1["point"])
    Var_delta[k]  <- unname(d1["var"])
    kept_vec[k]   <- SAFE$kept
    total_vec[k]  <- SAFE$total
    status_vec[k] <- SAFE$status
    
    lnM_BC[k] <- if (is.na(lnM_PI[k])) SAFE$point else 2*lnM_PI[k] - SAFE$point
    Var_SAFE[k] <- SAFE$var
  }
  
  Var_MC_PI <- var(lnM_PI[is.finite(lnM_PI)], na.rm = TRUE)
  Var_MC_BC <- var(lnM_BC[is.finite(lnM_BC)], na.rm = TRUE)
  
  Var_delta_bar <- mean(Var_delta[is.finite(Var_delta)], na.rm = TRUE)
  Var_SAFE_bar  <- mean(Var_SAFE[is.finite(Var_SAFE)],   na.rm = TRUE)
  
  rb_delta_PI <- 100 * (Var_delta_bar - Var_MC_PI) / Var_MC_PI
  rb_delta_BC <- 100 * (Var_delta_bar - Var_MC_BC) / Var_MC_BC
  rb_SAFE_PI  <- 100 * (Var_SAFE_bar  - Var_MC_PI) / Var_MC_PI
  rb_SAFE_BC  <- 100 * (Var_SAFE_bar  - Var_MC_BC) / Var_MC_BC
  
  bias_PI <- mean(lnM_PI, na.rm = TRUE) - lnM_true
  bias_BC <- mean(lnM_BC, na.rm = TRUE) - lnM_true
  
  cover_PI  <- mean(abs(lnM_PI - lnM_true) <= 1.96*sqrt(Var_delta), na.rm = TRUE)
  cover_BC  <- mean(abs(lnM_BC - lnM_true) <= 1.96*sqrt(Var_SAFE),  na.rm = TRUE)
  rmse_PI   <- sqrt(mean((lnM_PI - lnM_true)^2, na.rm = TRUE))
  rmse_BC   <- sqrt(mean((lnM_BC - lnM_true)^2, na.rm = TRUE))
  
  tibble(
    design   = "paired",
    n1 = n, n2 = n, r = r,
    theta = theta,
    lnM_true = lnM_true,
    bias_PI = bias_PI,
    bias_BC = bias_BC,
    Var_MC_PI = Var_MC_PI,
    Var_MC_BC = Var_MC_BC,
    Var_delta_bar = Var_delta_bar,
    Var_SAFE_bar  = Var_SAFE_bar,
    rb_delta_PI = rb_delta_PI,
    rb_delta_BC = rb_delta_BC,
    rb_SAFE_PI  = rb_SAFE_PI,
    rb_SAFE_BC  = rb_SAFE_BC,
    cover_PI = cover_PI,
    cover_BC = cover_BC,
    rmse_PI  = rmse_PI,
    rmse_BC  = rmse_BC,
    SAFE_kept_rate = mean(kept_vec / pmax(1L, total_vec), na.rm = TRUE),
    SAFE_status_ok = mean(status_vec == "ok", na.rm = TRUE)
  )
}

# ------------------- Run full grid -------------------------------
run_full_simulation <- function(K = K_default,
                                min_kept = MIN_KEPT_SAFE,
                                chunk_init = CHUNK_INIT_SAFE,
                                chunk_max = CHUNK_MAX_SAFE,
                                max_draws = MAX_DRAWS_SAFE,
                                patience_noaccept = PATIENCE_SAFE,
                                r = r_default) {
  indep_tbl <- expand_grid(theta = theta_grid, dn = designs_indep) %>%
    mutate(n1 = map_int(dn, 1L), n2 = map_int(dn, 2L)) %>%
    select(-dn)
  
  paired_tbl <- expand_grid(theta = theta_grid, n = designs_paired_n)
  
  res_indep <- pmap_dfr(indep_tbl,
                        ~ run_param_indep(theta = ..1, n1 = ..2, n2 = ..3,
                                          K = K,
                                          min_kept = min_kept,
                                          chunk_init = chunk_init,
                                          chunk_max = chunk_max,
                                          max_draws = max_draws,
                                          patience_noaccept = patience_noaccept))
  res_paired <- pmap_dfr(paired_tbl,
                         ~ run_param_dep(theta = ..1, n = ..2, r = r,
                                         K = K,
                                         min_kept = min_kept,
                                         chunk_init = chunk_init,
                                         chunk_max = chunk_max,
                                         max_draws = max_draws,
                                         patience_noaccept = patience_noaccept))
  bind_rows(res_indep, res_paired)
}

# ------------------- Plot helpers (15 + extras) ------------------
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
    results %>% select(theta, facet_label, bias_PI) %>%
      rename(bias = bias_PI) %>% mutate(estimator = "PI"),
    results %>% select(theta, facet_label, bias_BC) %>%
      rename(bias = bias_BC) %>% mutate(estimator = "SAFE-BC")
  )
  ggplot(df, aes(theta, bias, colour = estimator, group = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = expression(theta), y = "Bias (estimate \u2212 true lnM)", colour = NULL) +
    scale_colour_manual(values = c("PI" = "firebrick", "SAFE-BC" = "steelblue")) +
    theme_bw(11)
}

plot_relbias_15 <- function(results) {
  df <- bind_rows(
    results %>% select(theta, facet_label, rb_delta_BC) %>%
      rename(relbias = rb_delta_BC) %>% mutate(estimator = "Δ-var vs MC(BC)"),
    results %>% select(theta, facet_label, rb_SAFE_BC) %>%
      rename(relbias = rb_SAFE_BC)  %>% mutate(estimator = "SAFE-var vs MC(BC)")
  )
  ggplot(df, aes(theta, relbias, colour = estimator, group = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = expression(theta), y = "Relative bias of Var (%)", colour = NULL) +
    scale_colour_manual(values = c("Δ-var vs MC(BC)" = "firebrick",
                                   "SAFE-var vs MC(BC)" = "steelblue")) +
    theme_bw(11)
}

plot_coverage <- function(results) {
  df <- bind_rows(
    results %>% transmute(theta, facet_label, estimator = "PI", cover = cover_PI),
    results %>% transmute(theta, facet_label, estimator = "SAFE-BC", cover = cover_BC)
  )
  ggplot(df, aes(theta, cover, colour = estimator, group = estimator)) +
    geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50") +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = expression(theta), y = "Empirical coverage", colour = NULL) +
    scale_colour_manual(values = c("PI" = "firebrick", "SAFE-BC" = "steelblue")) +
    theme_bw(11)
}

plot_rmse <- function(results) {
  df <- bind_rows(
    results %>% transmute(theta, facet_label, estimator = "PI",     rmse = rmse_PI),
    results %>% transmute(theta, facet_label, estimator = "SAFE-BC", rmse = rmse_BC)
  )
  ggplot(df, aes(theta, rmse, colour = estimator, group = estimator)) +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    scale_colour_manual(values = c("PI" = "firebrick", "SAFE-BC" = "steelblue")) +
    labs(x = expression(theta), y = "RMSE ( ln M̂ − ln M )", colour = NULL) +
    theme_bw(11)
}

plot_safe_accept <- function(results) {
  ggplot(results, aes(theta, SAFE_kept_rate)) +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = expression(theta), y = "SAFE acceptance rate (kept/total)") +
    theme_bw(11)
}

plot_relbias_grid <- function(results) {
  df <- bind_rows(
    results %>% transmute(theta, facet_label, estimator = "Δ-var",    baseline = "MC(PI)", relbias = rb_delta_PI),
    results %>% transmute(theta, facet_label, estimator = "Δ-var",    baseline = "MC(BC)", relbias = rb_delta_BC),
    results %>% transmute(theta, facet_label, estimator = "SAFE-var", baseline = "MC(PI)", relbias = rb_SAFE_PI),
    results %>% transmute(theta, facet_label, estimator = "SAFE-var", baseline = "MC(BC)", relbias = rb_SAFE_BC)
  )
  ggplot(df, aes(theta, relbias, group = interaction(estimator, baseline),
                 colour = estimator, linetype = baseline)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line() +
    facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
    labs(x = expression(theta), y = "Relative bias of Var (%)",
         colour = "Estimator", linetype = "Baseline") +
    scale_colour_manual(values = c("Δ-var" = "firebrick", "SAFE-var" = "steelblue")) +
    theme_bw(11)
}

# ------------------- Run (demo) + Save + Plot --------------------
if (sys.nframe() == 0) {
  message("Running a small demo; increase K and MIN_KEPT_SAFE for HPC-scale runs.")
  K_demo  <- 200
  MIN_demo <- 2000
  CH_demo <- 2000
  
  results <- run_full_simulation(K = K_demo,
                                 min_kept = MIN_demo,
                                 chunk_init = CH_demo,
                                 chunk_max = CHUNK_MAX_SAFE,
                                 max_draws = MAX_DRAWS_SAFE,
                                 patience_noaccept = PATIENCE_SAFE)
  
  # Save
  stamp <- format(Sys.Date(), "%Y-%m-%d")
  saveRDS(results, file = sprintf("lnM_sim16_SAFE_%s.rds", stamp))
  write.csv(results, file = sprintf("lnM_sim16_SAFE_%s.csv", stamp), row.names = FALSE)
  
  # Facet ordering like simulation15
  results <- facet_ordering(results)
  
  # Plots from simulaiton15
  p_bias     <- plot_bias(results);       print(p_bias)
  p_relbias  <- plot_relbias_15(results); print(p_relbias)
  p_cover    <- plot_coverage(results);   print(p_cover)
  p_rmse     <- plot_rmse(results);       print(p_rmse)
  
  # Extras
  p_accept   <- plot_safe_accept(results); print(p_accept)
  p_grid4    <- plot_relbias_grid(results); print(p_grid4)
  
  # Quick roll-ups
  message(sprintf("Mean |bias| PI: %.4f",  mean(abs(results$bias_PI), na.rm = TRUE)))
  message(sprintf("Mean |bias| SAFE-BC: %.4f", mean(abs(results$bias_BC), na.rm = TRUE)))
  message(sprintf("Mean |relbias| Δ-var vs MC(BC): %.2f%%",
                  mean(abs(results$rb_delta_BC), na.rm = TRUE)))
  message(sprintf("Mean |relbias| SAFE-var vs MC(BC): %.2f%%",
                  mean(abs(results$rb_SAFE_BC),  na.rm = TRUE)))
}