########################################################################
## lnM simulation (INDEPENDENT ONLY): n=5 and n=10 (balanced)
## Estimators: PI, SAFE_raw, SAFE_BC, Qmid_05/10, Trim_05/10
## SAFE uses Normal means + Chi-square variances; cloud is conditional on Delta*>0
##
## Run:
##   K_REPL=500 MIN_KEPT=5000 THETA_MAX=5 Rscript sim_lnM_indep_qmid.R
########################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

theme_set(theme_bw(11))

## ---------------- helpers ----------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)

lnM_core <- function(Delta, MSW, n0) 0.5 * (log(Delta) - log(n0) - log(MSW))

## "true" lnM under the generative model: sigma=1, indep groups
## (this is the population analogue using MSB = h * theta^2, MSW = 1)
lnM_true_ind <- function(theta, n1, n2, sigma = 1) {
  h   <- n1 * n2 / (n1 + n2)
  MSB <- h * theta^2
  MSW <- sigma^2
  if (MSB <= MSW) return(NA_real_)  # lnM undefined when Delta<=0 at population level
  lnM_core(MSB - MSW, MSW, 2 * h)
}

## quantile-midpoint + trimmed mean from the SAFE cloud
qmid <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NA_real_)
  qs <- stats::quantile(x, probs = c(p, 1 - p), type = 8, names = FALSE)
  mean(qs)
}
trim_mean <- function(x, trim) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NA_real_)
  mean(x, trim = trim)
}

## ---------------- PI (Delta-1 plug-in), point only ----------------
lnM_PI_indep <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h    <- n1 * n2 / (n1 + n2)
  MSB  <- h * (x1bar - x2bar)^2
  MSW  <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(NA_real_)
  lnM_core(Delta, MSW, 2 * h)
}

## ---------------- SAFE (independent): return cloud + point summaries ----------------
safe_lnM_indep_cloud <- function(x1bar, x2bar, s1, s2, n1, n2,
                                 min_kept = 5000,
                                 chunk_init = 4000,
                                 chunk_max  = 2e6,
                                 max_draws  = Inf,
                                 patience_noaccept = 5) {
  df1 <- n1 - 1L
  df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)
  
  lnM_star <- numeric(0L)
  total <- 0L; kept <- 0L; attempts <- 0L
  zero_streak <- 0L
  chunk <- as.integer(chunk_init)
  status <- "ok"
  
  while (kept < min_kept && total < max_draws) {
    attempts <- attempts + 1L
    
    # Means ~ Normal
    m1 <- rnorm(chunk, mean = x1bar, sd = s1 / sqrt(n1))
    m2 <- rnorm(chunk, mean = x2bar, sd = s2 / sqrt(n2))
    
    # Variances ~ scaled Chi-square
    v1 <- s1^2 * stats::rchisq(chunk, df = df1) / df1
    v2 <- s2^2 * stats::rchisq(chunk, df = df2) / df2
    
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
    cloud    = lnM_star,
    kept     = kept,
    total    = total,
    attempts = attempts,
    status   = status
  )
}

## ---------------- one replicate ----------------
one_rep <- function(theta, n, min_kept) {
  x1 <- rnorm(n, 0, 1)
  x2 <- rnorm(n, theta, 1)
  x1bar <- mean(x1); x2bar <- mean(x2)
  s1 <- sd(x1); s2 <- sd(x2)
  
  PI <- lnM_PI_indep(x1bar, x2bar, s1, s2, n, n)
  
  safe <- safe_lnM_indep_cloud(x1bar, x2bar, s1, s2, n, n, min_kept = min_kept)
  
  cloud <- safe$cloud
  SAFE_raw <- if (length(cloud)) mean(cloud) else NA_real_
  
  # BC uses PI if finite (classic bootstrap bias correction form, but here using PI as "anchor")
  SAFE_BC <- if (is.finite(PI) && is.finite(SAFE_raw)) 2 * PI - SAFE_raw else NA_real_
  
  Qmid_05 <- if (length(cloud)) qmid(cloud, 0.05) else NA_real_
  Qmid_10 <- if (length(cloud)) qmid(cloud, 0.10) else NA_real_
  
  Trim_05 <- if (length(cloud)) trim_mean(cloud, 0.05) else NA_real_
  Trim_10 <- if (length(cloud)) trim_mean(cloud, 0.10) else NA_real_
  
  c(
    PI = PI,
    SAFE_raw = SAFE_raw,
    SAFE_BC  = SAFE_BC,
    Qmid_05  = Qmid_05,
    Qmid_10  = Qmid_10,
    Trim_05  = Trim_05,
    Trim_10  = Trim_10,
    kept     = safe$kept,
    tried    = safe$total,
    safe_ok  = as.integer(isTRUE(safe$status == "ok"))
  )
}

## ---------------- simulation settings ----------------
set.seed(20251218)

K_REPL   <- as.integer(Sys.getenv("K_REPL", "400"))
MIN_KEPT <- as.integer(Sys.getenv("MIN_KEPT", "5000"))
THETA_MAX <- as.numeric(Sys.getenv("THETA_MAX", "5"))

theta_vals <- c(0.25, 0.35, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5)
theta_vals <- theta_vals[theta_vals <= THETA_MAX]

n_vals <- c(5, 10)

param_grid <- expand.grid(theta = theta_vals, n = n_vals) %>%
  arrange(n, theta)

## ---------------- run ----------------
results_list <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {
  th <- param_grid$theta[i]
  n  <- param_grid$n[i]
  
  true_ln <- lnM_true_ind(th, n, n, sigma = 1)
  
  M <- matrix(NA_real_, nrow = 10, ncol = K_REPL,
              dimnames = list(c("PI","SAFE_raw","SAFE_BC","Qmid_05","Qmid_10","Trim_05","Trim_10","kept","tried","safe_ok"), NULL))
  
  for (k in seq_len(K_REPL)) {
    M[, k] <- one_rep(theta = th, n = n, min_kept = MIN_KEPT)
  }
  
  bias_fun <- function(x) if (is.finite(true_ln)) mean(x, na.rm = TRUE) - true_ln else NA_real_
  rmse_fun <- function(x) if (is.finite(true_ln)) sqrt(mean((x - true_ln)^2, na.rm = TRUE)) else NA_real_
  
  out <- data.frame(
    theta = th, n = n, true_lnM = true_ln,
    
    bias_PI       = bias_fun(M["PI",]),
    bias_SAFE_raw = bias_fun(M["SAFE_raw",]),
    bias_SAFE_BC  = bias_fun(M["SAFE_BC",]),
    bias_Qmid_05  = bias_fun(M["Qmid_05",]),
    bias_Qmid_10  = bias_fun(M["Qmid_10",]),
    bias_Trim_05  = bias_fun(M["Trim_05",]),
    bias_Trim_10  = bias_fun(M["Trim_10",]),
    
    rmse_PI       = rmse_fun(M["PI",]),
    rmse_SAFE_raw = rmse_fun(M["SAFE_raw",]),
    rmse_SAFE_BC  = rmse_fun(M["SAFE_BC",]),
    rmse_Qmid_05  = rmse_fun(M["Qmid_05",]),
    rmse_Qmid_10  = rmse_fun(M["Qmid_10",]),
    rmse_Trim_05  = rmse_fun(M["Trim_05",]),
    rmse_Trim_10  = rmse_fun(M["Trim_10",]),
    
    SAFE_ok_rate  = mean(M["safe_ok",], na.rm = TRUE),
    accept_prop   = sum(M["kept",], na.rm = TRUE) / max(1, sum(M["tried",], na.rm = TRUE))
  )
  
  results_list[[i]] <- out
}

results <- do.call(rbind, results_list)

## ---------------- plots ----------------
results <- results %>% mutate(facet_label = paste0("n=", n))

bias_df <- results %>%
  select(theta, facet_label, starts_with("bias_")) %>%
  pivot_longer(cols = starts_with("bias_"),
               names_to = "Estimator", values_to = "Bias") %>%
  mutate(Estimator = gsub("^bias_", "", Estimator),
         Estimator = recode(Estimator,
                            PI="PI",
                            SAFE_raw="SAFE_raw",
                            SAFE_BC="SAFE_BC",
                            Qmid_05="Qmid_05",
                            Qmid_10="Qmid_10",
                            Trim_05="Trim_05",
                            Trim_10="Trim_10"))

p_bias <- ggplot(bias_df, aes(theta, Bias, colour = Estimator, group = Estimator)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ facet_label, nrow = 1) +
  labs(x = expression(theta),
       y = "Bias (estimate - true lnM)",
       title = "Independent lnM: Bias curves",
       colour = "Estimator")

rmse_df <- results %>%
  select(theta, facet_label, starts_with("rmse_")) %>%
  pivot_longer(cols = starts_with("rmse_"),
               names_to = "Estimator", values_to = "RMSE") %>%
  mutate(Estimator = gsub("^rmse_", "", Estimator),
         Estimator = recode(Estimator,
                            PI="PI",
                            SAFE_raw="SAFE_raw",
                            SAFE_BC="SAFE_BC",
                            Qmid_05="Qmid_05",
                            Qmid_10="Qmid_10",
                            Trim_05="Trim_05",
                            Trim_10="Trim_10"))

p_rmse <- ggplot(rmse_df, aes(theta, RMSE, colour = Estimator, group = Estimator)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ facet_label, nrow = 1) +
  labs(x = expression(theta),
       y = "RMSE",
       title = "Independent lnM: RMSE curves",
       colour = "Estimator")

print(p_bias)
print(p_rmse)

## diagnostics
message("Mean SAFE_ok_rate:  ", mean(results$SAFE_ok_rate, na.rm = TRUE))
message("Mean accept_prop:   ", mean(results$accept_prop,  na.rm = TRUE))

## save
out_rds <- sprintf("lnM_indep_qmid_summary_%s.rds", Sys.Date())
out_csv <- sprintf("lnM_indep_qmid_summary_%s.csv", Sys.Date())
saveRDS(results, out_rds)
write.csv(results, out_csv, row.names = FALSE)
message("Saved: ", out_rds)
message("Saved: ", out_csv)