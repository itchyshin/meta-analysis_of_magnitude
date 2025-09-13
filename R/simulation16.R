# ================================================================
#  simulation15.R  —  lnM simulations with Chi-square / Wishart SAFE
#  - Single delta-1 only (indep/paired) per MS Eq. (6),(7)
#  - SAFE uses Normal + Chi-square (indep) / Wishart (paired)
#  - Computes all combinations for Eq. (18)-(20):
#       Var_delta vs Var_SAFE   ×   Var_MC from lnM_PI vs lnM_BC
#  - Closely follows MS §§2.2–2.4 (Eqs. (1)–(7), (11)–(14), (15)–(20))
#  - Design grid per MS: 144 indep + 72 paired = 216 parameter sets
# ================================================================

suppressPackageStartupMessages({
  library(MASS)     # mvrnorm
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# ------------------- Global controls / knobs --------------------
SEED_MAIN      <- 12345
set.seed(SEED_MAIN)

# Monte Carlo replicates per parameter set (K in MS §2.4)
K_default      <- 2000     # <-- set higher (e.g., 1e5) on HPC

# SAFE bootstrap size per replicate
B_SAFE         <- 2000     # <-- increase (e.g., 1e5) for final runs
chunk_SAFE     <- 4000     # vectorized chunk size inside SAFE
max_chunks     <- 50       # hard stop to avoid stuck while-loops

# Cap for delta-method variance (per MS note)
DELTA_VAR_CAP  <- 20

# Paired within-pair correlation
r_default      <- 0.8

# θ grid (MS §2.4)
theta_grid <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,4,5)

# Independent designs (n1,n2) — 8 cells (MS)
designs_indep <- list(
  c(5,5), c(10,10), c(20,20), c(100,100),
  c(3,7), c(6,14), c(12,28), c(40,160)
)

# Paired designs (balanced) — 4 cells (MS)
designs_paired_n <- c(5,10,20,100)

# ----------------------------------------------------------------
# Helpers (safe numerics)
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(gap) ifelse(gap <= 0, NA_real_, gap)

# ----------------------------------------------------------------
# lnM delta-1 (independent) — MS Eq. (1)-(4), variance Eq. (6)
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h    <- n1 * n2 / (n1 + n2)                           # harmonic component
  MSB  <- h * (x1bar - x2bar)^2                         # Eq. (1)
  MSW  <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)  # Eq. (2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA, var = NA, se = NA))
  
  # lnM — Eq. (4); n0 = 2h (indep)
  lnM <- 0.5 * (log(Delta / (2 * h)) - log(MSW))
  
  # Delta-1 variance — Eq. (6)
  sigmaD2 <- s1^2 / n1 + s2^2 / n2
  delta   <- x1bar - x2bar
  Var_B   <- h^2 * (2 * sigmaD2^2 + 4 * sigmaD2 * delta^2)
  Var_W   <- 2 * MSW^2 / (n1 + n2 - 2)
  
  gB   <- 0.5 / Delta
  gW   <- -0.5 * MSB / (Delta * MSW)
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  Var1 <- pmin(Var1, DELTA_VAR_CAP)   # cap per MS note
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# lnM delta-1 (paired) — MS paired forms, variance Eq. (7)
lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r) {
  h     <- n / 2
  MSB   <- h * (x1bar - x2bar)^2
  MSW   <- (s1^2 + s2^2) / 2
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA, var = NA, se = NA))
  
  # lnM — Eq. (4); n0 = n (paired)
  lnM <- 0.5 * (log(Delta / n) - log(MSW))
  
  # Delta-1 variance — Eq. (7)
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

# ----------------------------------------------------------------
# SAFE bootstrap — independent (Chi-square for variances)
safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           B = 1e4, chunk = 5e3, max_chunks = 50) {
  df1 <- n1 - 1L
  df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)
  
  lnM_star <- numeric(0)
  total    <- 0L
  kept     <- 0L
  attempts <- 0L
  
  while (length(lnM_star) < B && attempts < max_chunks) {
    attempts <- attempts + 1L
    
    # Means ~ Normal; variances ~ scaled Chi-square (strictly > 0)
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
      lnM_star <- c(lnM_star,
                    0.5 * (log((MSB[good] - MSW[good]) / (2 * h)) -
                             log(MSW[good])))
    }
  }
  
  status <- if (length(lnM_star) >= B) "ok" else "stopped_early"
  if (length(lnM_star) > 0L) lnM_star <- lnM_star[seq_len(min(B, length(lnM_star)))]
  
  list(mean_star = if (length(lnM_star)) mean(lnM_star) else NA_real_,
       var_star  = if (length(lnM_star)) var(lnM_star)  else NA_real_,
       kept      = kept,
       total     = total,
       attempts  = attempts,
       status    = status)
}

# SAFE bootstrap — paired (Wishart for variances)
safe_lnM_dep <- function(x1bar, x2bar, s1, s2, n, r,
                         B = 1e4, chunk = 5e3, max_chunks = 50) {
  df <- n - 1L
  h  <- n / 2
  Sig <- matrix(c(s1^2, r*s1*s2,
                  r*s1*s2, s2^2), 2, 2)
  
  lnM_star <- numeric(0)
  total    <- 0L
  kept     <- 0L
  attempts <- 0L
  
  while (length(lnM_star) < B && attempts < max_chunks) {
    attempts <- attempts + 1L
    
    Mu <- MASS::mvrnorm(n = chunk, mu = c(x1bar, x2bar), Sigma = Sig / n)
    W  <- stats::rWishart(n = chunk, df = df, Sigma = Sig)   # 2x2xchunk array
    S11 <- W[1,1,] / df
    S22 <- W[2,2,] / df
    
    total <- total + chunk
    
    MSB  <- h * (Mu[,1] - Mu[,2])^2
    MSW  <- (S11 + S22) / 2
    good <- MSB > MSW
    
    if (any(good)) {
      kept <- kept + sum(good)
      lnM_star <- c(lnM_star,
                    0.5 * (log((MSB[good] - MSW[good]) / n) -
                             log(MSW[good])))
    }
  }
  
  status <- if (length(lnM_star) >= B) "ok" else "stopped_early"
  if (length(lnM_star) > 0L) lnM_star <- lnM_star[seq_len(min(B, length(lnM_star)))]
  
  list(mean_star = if (length(lnM_star)) mean(lnM_star) else NA_real_,
       var_star  = if (length(lnM_star)) var(lnM_star)  else NA_real_,
       kept      = kept,
       total     = total,
       attempts  = attempts,
       status    = status)
}

# ----------------------------------------------------------------
# "True" lnM under Normal-theory expected mean squares.
# This avoids undefined cases by using E[MS_B] and E[MS_W].
# Indep: E[MSW] = ω; E[MSB] = ω + h*δ^2   => Δ_true = h*δ^2
true_lnM_indep <- function(mu1, mu2, sd1, sd2, n1, n2) {
  h    <- n1*n2/(n1+n2)
  n0   <- 2*h
  delta <- mu1 - mu2
  omega <- ((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2)
  Delta_true <- h * delta^2
  0.5 * (log(Delta_true / n0) - log(omega))
}

# Paired: E[MSW] = (σ1^2+σ2^2)/2; E[MSB] = (n/2)*(δ^2 + σ_D^2/n)
# with σ_D^2 = σ1^2 + σ2^2 - 2ρσ1σ2  => Δ_true = (n/2)δ^2 - ρσ1σ2
true_lnM_dep <- function(mu1, mu2, sd1, sd2, n, r) {
  delta <- mu1 - mu2
  sigmaD2 <- sd1^2 + sd2^2 - 2*r*sd1*sd2
  E_MSW   <- (sd1^2 + sd2^2)/2
  Delta_true <- (n/2)*delta^2 + (1/2)*sigmaD2 - E_MSW  # simplifies to (n/2)δ^2 - r σ1 σ2
  n0 <- n
  0.5 * (log(Delta_true / n0) - log(E_MSW))
}

# ----------------------------------------------------------------
# One Monte-Carlo replicate of sample summaries using Normal theory (fast)
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

# ----------------------------------------------------------------
# Run one parameter set (independent)
run_param_indep <- function(theta, n1, n2,
                            K = K_default, B = B_SAFE,
                            chunk = chunk_SAFE, max_chunks = max_chunks,
                            mu1 = 0, sd1 = 1, sd2 = 1) {
  
  mu2 <- theta
  lnM_true <- true_lnM_indep(mu1, mu2, sd1, sd2, n1, n2)
  
  # storage
  lnM_PI     <- numeric(K)
  Var_delta  <- numeric(K)
  lnM_BC     <- numeric(K)    # SAFE bias-corrected (or SAFE-only if PI undefined)
  Var_SAFE   <- numeric(K)
  kept_vec   <- integer(K)
  total_vec  <- integer(K)
  status_vec <- character(K)
  
  for (k in seq_len(K)) {
    sm <- draw_summaries_indep(mu1, mu2, sd1, sd2, n1, n2)
    d1 <- lnM_delta1_indep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"], n1, n2)
    SAFE <- safe_lnM_indep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"],
                           n1, n2, B = B, chunk = chunk, max_chunks = max_chunks)
    
    lnM_PI[k]    <- unname(d1["point"])
    Var_delta[k] <- unname(d1["var"])
    kept_vec[k]  <- SAFE$kept
    total_vec[k] <- SAFE$total
    status_vec[k] <- SAFE$status
    
    # Bias-corrected point (Eq. (13)); when PI undefined, use SAFE mean_star (Eq. (14))
    if (is.na(lnM_PI[k])) {
      lnM_BC[k] <- SAFE$mean_star
    } else {
      lnM_BC[k] <- 2*lnM_PI[k] - SAFE$mean_star
    }
    Var_SAFE[k] <- SAFE$var_star
  }
  
  # Monte-Carlo variances for both point estimators (Eq. (20))
  Var_MC_PI <- var(lnM_PI[is.finite(lnM_PI)], na.rm = TRUE)
  Var_MC_BC <- var(lnM_BC[is.finite(lnM_BC)], na.rm = TRUE)
  
  # Averages of variance estimators (Eq. (19))
  Var_delta_bar <- mean(Var_delta[is.finite(Var_delta)], na.rm = TRUE)
  Var_SAFE_bar  <- mean(Var_SAFE[is.finite(Var_SAFE)],   na.rm = TRUE)
  
  # Relative biases (Eq. (18)) — all 4 combos
  rb_delta_PI <- 100 * (Var_delta_bar - Var_MC_PI) / Var_MC_PI
  rb_delta_BC <- 100 * (Var_delta_bar - Var_MC_BC) / Var_MC_BC
  rb_SAFE_PI  <- 100 * (Var_SAFE_bar  - Var_MC_PI) / Var_MC_PI
  rb_SAFE_BC  <- 100 * (Var_SAFE_bar  - Var_MC_BC) / Var_MC_BC
  
  # Biases of point estimators (Eq. (16)-(17))
  bias_PI <- mean(lnM_PI, na.rm = TRUE) - lnM_true
  bias_BC <- mean(lnM_BC, na.rm = TRUE) - lnM_true
  
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
    SAFE_kept_rate = mean(kept_vec / pmax(1L, total_vec)),
    SAFE_status_ok = mean(status_vec == "ok")
  )
}

# Run one parameter set (paired)
run_param_dep <- function(theta, n,
                          r = r_default,
                          K = K_default, B = B_SAFE,
                          chunk = chunk_SAFE, max_chunks = max_chunks,
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
    SAFE <- safe_lnM_dep(sm["x1bar"], sm["x2bar"], sm["s1"], sm["s2"],
                         n, r, B = B, chunk = chunk, max_chunks = max_chunks)
    
    lnM_PI[k]    <- unname(d1["point"])
    Var_delta[k] <- unname(d1["var"])
    kept_vec[k]  <- SAFE$kept
    total_vec[k] <- SAFE$total
    status_vec[k] <- SAFE$status
    
    if (is.na(lnM_PI[k])) {
      lnM_BC[k] <- SAFE$mean_star
    } else {
      lnM_BC[k] <- 2*lnM_PI[k] - SAFE$mean_star
    }
    Var_SAFE[k] <- SAFE$var_star
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
    SAFE_kept_rate = mean(kept_vec / pmax(1L, total_vec)),
    SAFE_status_ok = mean(status_vec == "ok")
  )
}

# ----------------------------------------------------------------
# Run full grid
run_full_simulation <- function(K = K_default, B = B_SAFE,
                                chunk = chunk_SAFE, max_chunks = max_chunks,
                                r = r_default) {
  
  # Independent grid: 18 θ × 8 designs = 144 sets
  indep_tbl <- expand_grid(theta = theta_grid,
                           dn = designs_indep) %>%
    mutate(n1 = map_int(dn, 1L), n2 = map_int(dn, 2L)) %>%
    select(-dn)
  
  # Paired grid: 18 θ × 4 designs = 72 sets
  paired_tbl <- expand_grid(theta = theta_grid,
                            n = designs_paired_n)
  
  # Run
  res_indep <- pmap_dfr(indep_tbl,
                        ~ run_param_indep(theta = ..1, n1 = ..2, n2 = ..3,
                                          K = K, B = B, chunk = chunk, max_chunks = max_chunks))
  
  res_paired <- pmap_dfr(paired_tbl,
                         ~ run_param_dep(theta = ..1, n = ..2, r = r,
                                         K = K, B = B, chunk = chunk, max_chunks = max_chunks))
  
  bind_rows(res_indep, res_paired)
}

# --------------------------- Go! --------------------------------
# (You can source this file and call run_full_simulation with HPC-scale K,B)

if (sys.nframe() == 0) {
  message("Running a small demo (reduce K/B for speed; increase for final runs).")
  demo_res <- run_full_simulation(K = 200, B = 1000, chunk = 2000, max_chunks = 30, r = r_default)
  print(demo_res %>% slice_head(n = 10))
  # Write to disk if desired:
  # write.csv(demo_res, "simulation15_results.csv", row.names = FALSE)
}