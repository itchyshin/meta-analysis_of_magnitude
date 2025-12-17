# ================================================================
#  lnM: Delta-1 vs SAFE bootstrap (independent & paired)
#  - SAFE uses Chi-square / Wishart (no MVN for variances)
#  - Single delta version (delta-1) only
#  - Loop bounded by min_kept + max_draws + patience_noaccept
# ================================================================

library(MASS)    # for mvrnorm (paired means)
library(dplyr)

# ---------- tiny helpers ----------------------------------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(gap) ifelse(gap <= 0, NA_real_, gap)

# ---------- Delta-1 (independent) --------------------------------
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2)
{
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
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ---------- Delta-1 (paired) ------------------------------------
lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r)
{
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
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ---------- SAFE bootstrap : independent -------------------------
safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           min_kept   = 1e4,   # target usable lnM* draws
                           B          = NULL,  # alias for min_kept
                           chunk_init = 5e3,   # starting chunk size
                           chunk_max  = 2e6,   # cap for adaptive chunk
                           max_draws  = Inf,   # cap on total attempted draws
                           patience_noaccept = 5)
{
  if (!is.null(B) && missing(min_kept)) min_kept <- B
  
  df1 <- n1 - 1L
  df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)
  
  lnM_star <- numeric(0L)
  total    <- 0L
  kept     <- 0L
  attempts <- 0L
  zero_streak <- 0L
  chunk <- as.integer(chunk_init)
  status <- "ok"
  
  while (kept < min_kept && total < max_draws) {
    attempts <- attempts + 1L
    
    # Step 2: means ~ Normal; variances ~ Chi-square
    m1 <- rnorm(chunk, mean = x1bar, sd = s1 / sqrt(n1))
    m2 <- rnorm(chunk, mean = x2bar, sd = s2 / sqrt(n2))
    v1 <- s1^2 * stats::rchisq(chunk, df = df1) / df1
    v2 <- s2^2 * stats::rchisq(chunk, df = df2) / df2
    
    total <- total + chunk
    
    # Step 3: transform
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
    
    # Adaptive next chunk
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
  
  list(point    = if (kept) mean(lnM_star) else NA_real_,
       var      = if (kept) var(lnM_star)  else NA_real_,
       kept     = kept,
       total    = total,
       attempts = attempts,
       status   = status)
}

# ---------- SAFE bootstrap : paired --------------------------------
safe_lnM_dep <- function(x1bar, x2bar, s1, s2, n, r,
                         min_kept   = 1e4,
                         B          = NULL,
                         chunk_init = 5e3,
                         chunk_max  = 2e6,
                         max_draws  = Inf,
                         patience_noaccept = 5)
{
  if (!is.null(B) && missing(min_kept)) min_kept <- B
  
  df <- n - 1L
  h  <- n / 2
  Sig <- matrix(c(s1^2, r*s1*s2,
                  r*s1*s2, s2^2), 2, 2)
  
  lnM_star <- numeric(0L)
  total    <- 0L
  kept     <- 0L
  attempts <- 0L
  zero_streak <- 0L
  chunk <- as.integer(chunk_init)
  status <- "ok"
  
  while (kept < min_kept && total < max_draws) {
    attempts <- attempts + 1L
    
    # Step 2: means (μ*) and Wishart draws
    Mu <- MASS::mvrnorm(n = chunk, mu = c(x1bar, x2bar), Sigma = Sig / n)
    W  <- stats::rWishart(n = chunk, df = df, Sigma = Sig)
    S11 <- W[1,1,] / df
    S22 <- W[2,2,] / df
    
    total <- total + chunk
    
    # Step 3: transform
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
    
    # Adaptive next chunk
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
  
  list(point    = if (kept) mean(lnM_star) else NA_real_,
       var      = if (kept) var(lnM_star)  else NA_real_,
       kept     = kept,
       total    = total,
       attempts = attempts,
       status   = status)
}

# ---------- comparison wrapper -----------------------------------
compare_methods <- function(x1, x2, s1, s2,
                            n1, n2 = NULL, r = NULL,
                            min_kept   = 1e4,
                            B          = NULL,
                            chunk_init = 4e3,
                            chunk_max  = 2e6,
                            max_draws  = Inf)
{
  if (is.null(n2)) {
    if (is.null(r)) stop("For paired data supply 'r' (within-pair correlation).")
    d1   <- lnM_delta1_dep(x1, x2, s1, s2, n1, r)
    safe <- safe_lnM_dep  (x1, x2, s1, s2, n1, r,
                           min_kept, B, chunk_init, chunk_max, max_draws)
  } else {
    d1   <- lnM_delta1_indep(x1, x2, s1, s2, n1, n2)
    safe <- safe_lnM_indep(x1, x2, s1, s2, n1, n2,
                           min_kept, B, chunk_init, chunk_max, max_draws)
  }
  
  out <- dplyr::bind_rows(
    data.frame(method = "Delta-1",        lnM = d1["point"], Var = d1["var"]),
    data.frame(method = "SAFE bootstrap", lnM = safe$point,  Var = safe$var)
  )
  
  print(out, row.names = FALSE, digits = 4)
  pct <- 100 * safe$kept / max(1L, safe$total)
  cat(sprintf("\nSAFE kept %d of %d draws (%.2f%% usable) across %d chunk(s) [%s]\n",
              safe$kept, safe$total, pct, safe$attempts, safe$status))
  invisible(out)
}

# ---------- quick demo ------------------------------------------
#set.seed(123)

# # Independent examples (target ≥ 100,000 usable SAFE draws)
# compare_methods(x1 = 5, x2 = 5,
#                 s1 = 1,  s2 = 1,
#                 n1 = 40, n2 = 40,
#                 min_kept = 1e5, chunk_init = 4000, max_draws = 5e6)
# 
# compare_methods(x1 = 7, x2 = 5,
#                 s1 = 1,  s2 = 1,
#                 n1 = 40, n2 = 40,
#                 min_kept = 1e5, chunk_init = 4000, max_draws = 5e6)
# 
# compare_methods(x1 = 8, x2 = 5,
#                 s1 = 2,  s2 = 2,
#                 n1 = 40, n2 = 40,
#                 min_kept = 1e5, chunk_init = 4000, max_draws = 5e6)
# 
# # Paired examples (r = 0.8; target ≥ 100,000 usable SAFE draws)
# compare_methods(x1 = 5, x2 = 5,
#                 s1 = 1,  s2 = 1,
#                 n1 = 40, r  = 0.8,
#                 min_kept = 1e5, chunk_init = 4000, max_draws = 5e6)
# 
# compare_methods(x1 = 7, x2 = 5,
#                 s1 = 1,  s2 = 1,
#                 n1 = 40, r  = 0.8,
#                 min_kept = 1e5, chunk_init = 4000, max_draws = 5e6)
# 
# compare_methods(x1 = 8, x2 = 5,
#                 s1 = 2,  s2 = 2,
#                 n1 = 40, r  = 0.8,
#                 min_kept = 1e5, chunk_init = 4000, max_draws = 5e6)
# 
# compare_methods(x1 = 5, x2 = 5,
#                 s1 = 2,  s2 = 2.2,
#                 n1 = 40, r  = 0.8,
#                 min_kept = 1e5, chunk_init = 4000, max_draws = 5e6)