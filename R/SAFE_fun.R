# ================================================================
#  Functions for Delta-method (Delta-1) and SAFE bootstrap for lnM
# ================================================================
#
# PURPOSE
#   Compute lnM (log magnitude / “unsigned” mean-difference signal-to-noise,
#   derived from a between-groups ANOVA decomposition) and its uncertainty
#   using either:
#     (i)  Delta-1 approximation (fast, analytic), or
#     (ii) SAFE bootstrap (simulation-based, robust near boundaries)
#
# CONTEXT (informal)
#   lnM is defined (for both independent and paired designs) from a positive
#   “gap” Delta = MSB - MSW, where:
#     MSB = between-group mean square (signal)
#     MSW = within-group mean square (noise)
#   lnM is only defined when Delta > 0. When MSB <= MSW, lnM is undefined
#   (log of a non-positive quantity); we return NA (and SAFE discards such draws).
#
#   The SAFE bootstrap here is parametric: it draws group means and variances
#   from their sampling distributions under the usual normal-model assumptions.
#   It then computes lnM on each draw, keeping only draws with MSB > MSW.
#
# NOTES ON NUMERICAL/ALGORITHMIC CHOICES
#   - We guard against negative/zero “gap” values (Delta <= 0) explicitly.
#   - We “posify” variances in Delta-method outputs to avoid tiny negative values
#     from floating-point arithmetic.
#   - SAFE uses adaptive chunking to target a minimum number of *usable* draws
#     (i.e., those with MSB > MSW), rather than a fixed number of total draws.
#   - patience_noaccept stops SAFE if we repeatedly fail to obtain *any* usable
#     draws, which typically indicates the acceptance probability is extremely
#     small for the supplied parameters.
#
# DEPENDENCIES
#   - MASS: for mvrnorm() and rWishart() in the paired SAFE bootstrap.
#   - dplyr: only used for bind_rows() in compare_methods().
#
# ================================================================

library(MASS)
library(dplyr)

# ---------- tiny helpers ----------------------------------------
# Ensure positive numeric output (useful when small negative values appear due
# to floating-point roundoff, e.g., a computed variance of -1e-16).
posify <- function(x, eps = 1e-12) pmax(x, eps)

# SAFE "gap": return NA if gap <= 0, because lnM requires gap > 0.
safe_gap <- function(gap) ifelse(gap <= 0, NA_real_, gap)

# A small “is effectively positive?” check used for acceptance rates in SAFE.
# (Acceptance rates can be extremely small; we treat values <= 1e-10 as 0.)
is_positive <- function(x) { !is.na(x) && x > 1e-10 }


# ================================================================
#  Delta-1 (independent groups)
# ================================================================
#
# INPUTS
#   x1bar, x2bar : sample means of groups 1 and 2
#   s1, s2       : sample SDs of groups 1 and 2
#   n1, n2       : sample sizes
#
# OUTPUT
#   named numeric vector: c(point = lnM, var = Var(lnM), se = sqrt(var))
#
# DETAILS
#   - Define h = n1*n2/(n1+n2) (the harmonic-like factor used in MSB).
#   - MSB  = h*(x1bar - x2bar)^2
#   - MSW  = pooled within-group variance (ANOVA within mean square)
#   - Delta = MSB - MSW must be > 0 to define lnM
#   - lnM  = 0.5 * [ log(Delta/(2h)) - log(MSW) ]
#
#   Delta-1 variance uses a first-order delta-method approximation, treating
#   MSB and MSW as approximately independent and using analytic approximations
#   to Var(MSB) and Var(MSW). (This is a pragmatic “Delta-1” approximation.)
#
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2) {
  
  # h is the standard two-sample ANOVA factor in MSB
  h   <- n1 * n2 / (n1 + n2)
  
  # Between-group mean square (signal)
  MSB <- h * (x1bar - x2bar)^2
  
  # Within-group mean square (pooled variance; noise)
  MSW <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  
  # Gap must be strictly positive
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA, var = NA, se = NA))
  
  # Point estimate
  lnM <- 0.5 * (log(Delta / (2 * h)) - log(MSW))
  
  # Approximate components for delta-method variance
  # sigmaD2: variance of (x1bar - x2bar) under independence
  sigmaD2 <- s1^2 / n1 + s2^2 / n2
  delta   <- x1bar - x2bar
  
  # Var_B: approximate variance of MSB = h*(delta)^2, using moment formula
  Var_B <- h^2 * (2 * sigmaD2^2 + 4 * sigmaD2 * delta^2)
  
  # Var_W: approximate variance of MSW (scaled chi-square variance approximation)
  Var_W <- 2 * MSW^2 / (n1 + n2 - 2)
  
  # Delta-method gradients of lnM wrt MSB and MSW
  # lnM = 0.5*(log(MSB - MSW) - log(2h) - log(MSW))
  # d/dMSB: 0.5 * 1/(MSB - MSW) = 0.5/Delta
  gB <- 0.5 / Delta
  
  # d/dMSW: 0.5 * [ -1/(MSB - MSW) - 1/MSW ] = -0.5*(MSB)/(Delta*MSW)
  # (algebraic simplification used in original code)
  gW <- -0.5 * MSB / (Delta * MSW)
  
  # Combine; enforce positivity
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}


# ================================================================
#  Delta-1 (paired / dependent)
# ================================================================
#
# INPUTS
#   x1bar, x2bar : paired-condition means
#   s1, s2       : paired-condition SDs
#   n            : number of pairs
#   r            : within-pair correlation between conditions
#
# OUTPUT
#   named numeric vector: c(point = lnM, var = Var(lnM), se = sqrt(var))
#
# DETAILS
#   - For paired data, we use h = n/2 and MSW = (s1^2 + s2^2)/2 (mean of variances).
#     This matches the two-condition repeated-measures ANOVA decomposition in a
#     simplified form.
#   - sigmaD2 is Var(mean difference) accounting for correlation:
#       sigmaD2 = s1^2 + s2^2 - 2*r*s1*s2
#   - Var_B and Var_W are approximations suitable for the paired setting.
#
lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r) {
  
  h   <- n / 2
  
  # Between-condition signal
  MSB <- h * (x1bar - x2bar)^2
  
  # Within-condition “noise” (average of the two condition variances)
  MSW <- (s1^2 + s2^2) / 2
  
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA, var = NA, se = NA))
  
  # Point estimate
  lnM <- 0.5 * (log(Delta / n) - log(MSW))
  
  # Var of the paired mean difference accounting for correlation r
  sigmaD2 <- s1^2 + s2^2 - 2 * r * s1 * s2
  delta   <- x1bar - x2bar
  
  # Approx Var(MSB) in paired setting (scaled for mean difference sampling var)
  Var_B <- h^2 * (2 * sigmaD2^2 / n^2 + 4 * delta^2 * sigmaD2 / n)
  
  # Approx Var(MSW) allowing correlation in the two sample variances
  Var_W <- (s1^4 + s2^4 + 2 * r^2 * s1^2 * s2^2) / (2 * (n - 1))
  
  # Gradients
  gB <- 0.5 / Delta
  gW <- -0.5 * MSB / (Delta * MSW)
  
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}


# ================================================================
#  SAFE bootstrap (independent groups)
# ================================================================
#
# INPUTS
#   x1bar, x2bar, s1, s2, n1, n2 : observed summary stats
#   min_kept   : target number of *usable* draws (MSB > MSW) to keep
#   chunk_init : starting chunk size (number of draws per iteration)
#   chunk_max  : cap on chunk size (protects memory/time)
#   max_draws  : cap on total draws attempted (protects time)
#   patience_noaccept :
#       if we see 'patience_noaccept' consecutive iterations with 0 usable
#       draws, we stop early and return status = "no_usable_draws".
#
# OUTPUT (list)
#   point, var : mean and variance of kept lnM* draws (if kept >= 2; else NA)
#   kept       : number of usable draws retained
#   total      : total draws attempted (including unusable)
#   attempts   : number of while-loop iterations (chunks processed)
#   status     : "ok" or "no_usable_draws"
#
# SAMPLING MODEL (parametric)
#   - m1, m2 : draws of sample means, Normal(mean = xbar, sd = s/sqrt(n))
#   - v1, v2 : draws of sample variances, scaled chi-square:
#       v ~ s^2 * ChiSq(df)/df
#     so that E[v] = s^2 under normality.
#
safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           min_kept = 2000,
                           chunk_init = 4000,
                           chunk_max = 2e6,
                           max_draws = Inf,
                           patience_noaccept = 5) {
  
  df1 <- n1 - 1L
  df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)
  
  # Storage for usable draws (kept lnM* values).
  # We grow the vector dynamically; for big jobs you might later consider
  # preallocation, but dynamic growth is usually fine at these sizes.
  lnM_star <- numeric(0L)
  
  # Counters
  total <- 0L     # total draws attempted
  kept  <- 0L     # usable draws retained
  attempts <- 0L  # number of chunk iterations
  
  # Early-stop tracker: how many chunks in a row had zero usable draws?
  zero_streak <- 0L
  
  # Adaptive chunk size starts at chunk_init
  chunk <- as.integer(chunk_init)
  
  # Status flag
  status <- "ok"
  
  while (kept < min_kept && total < max_draws) {
    
    attempts <- attempts + 1L
    
    # ---- draw sampling distributions of means and variances ----
    # Sample means
    m1 <- rnorm(chunk, mean = x1bar, sd = s1 / sqrt(n1))
    m2 <- rnorm(chunk, mean = x2bar, sd = s2 / sqrt(n2))
    
    # Sample variances (chi-square sampling distribution under normality)
    v1 <- s1^2 * rchisq(chunk, df = df1) / df1
    v2 <- s2^2 * rchisq(chunk, df = df2) / df2
    
    # ---- compute MSB, MSW and accept draws with MSB > MSW ----
    MSB <- h * (m1 - m2)^2
    MSW <- (df1 * v1 + df2 * v2) / (df1 + df2)
    
    good <- which(MSB > MSW)
    n_good <- length(good)
    
    if (n_good > 0L) {
      zero_streak <- 0L
      
      # lnM* for accepted draws:
      # lnM* = 0.5 * [ log((MSB - MSW)/(2h)) - log(MSW) ]
      vals <- 0.5 * (log((MSB[good] - MSW[good]) / (2 * h)) - log(MSW[good]))
      
      lnM_star <- c(lnM_star, vals)
      kept <- length(lnM_star)
      
    } else {
      # No usable draws this chunk: track a streak and possibly stop
      zero_streak <- zero_streak + 1L
      if (zero_streak >= patience_noaccept) {
        status <- "no_usable_draws"
        break
      }
    }
    
    total <- total + chunk
    
    # ---- adaptive chunking ----
    # Acceptance rate acc = kept/total
    # If acceptance is non-negligible, we choose next chunk to aim for the
    # remaining usable draws; otherwise, we grow chunk size geometrically.
    acc <- if (total > 0) kept / total else 0
    
    if (is_positive(acc)) {
      # Remaining usable draws needed: (min_kept - kept)
      # Estimated total draws needed: remaining / acc
      # We cap to chunk_max, and also guard with min() to avoid overflow.
      next_needed <- min(chunk_max, ceiling(max(0, min_kept - kept) / acc))
    } else {
      # Practically zero acceptance so far: double chunk to probe feasibility
      next_needed <- chunk * 2L
    }
    
    # Keep chunk within [chunk_init, chunk_max]
    chunk <- as.integer(max(chunk_init, min(chunk_max, next_needed)))
  }
  
  list(
    point = if (kept >= 2) mean(lnM_star) else NA_real_,
    var   = if (kept >= 2) var(lnM_star)  else NA_real_,
    kept  = kept,
    total = total,
    attempts = attempts,
    status = status
  )
}

# ================================================================
#  SAFE bootstrap (paired / dependent)
# ================================================================
#
# Paired SAFE differs only in the parametric sampling model:
#   - Sample mean vector (x1bar, x2bar) has covariance Sig/n.
#   - Sample covariance matrix S has Wishart(df=n-1, Sigma=Sig),
#     so that the diagonal elements yield sampled condition variances.
#
# This corresponds to a multivariate normal model for paired observations.
#
safe_lnM_dep <- function(x1bar, x2bar, s1, s2, n, r,
                         min_kept = 2000,
                         chunk_init = 4000,
                         chunk_max = 2e6,
                         max_draws = Inf,
                         patience_noaccept = 5) {
  
  df <- n - 1L
  h  <- n / 2
  
  # Population covariance matrix implied by s1, s2, r
  Sig <- matrix(c(s1^2, r * s1 * s2,
                  r * s1 * s2, s2^2), 2, 2)
  
  lnM_star <- numeric(0L)
  total <- 0L
  kept  <- 0L
  attempts <- 0L
  zero_streak <- 0L
  chunk <- as.integer(chunk_init)
  status <- "ok"
  
  while (kept < min_kept && total < max_draws) {
    
    attempts <- attempts + 1L
    
    # Draw sample means (vector) for each bootstrap draw:
    # Mean of paired sample means: MVN(mu, Sigma/n)
    Mu <- mvrnorm(n = chunk, mu = c(x1bar, x2bar), Sigma = Sig / n)
    
    # Draw sample covariance matrices for each bootstrap draw:
    # W ~ Wishart(df, Sigma)
    W  <- rWishart(n = chunk, df = df, Sigma = Sig)
    
    # Extract diagonal (sample variances); divide by df to get unbiased-ish scale
    S11 <- W[1, 1, ] / df
    S22 <- W[2, 2, ] / df
    
    MSB <- h * (Mu[, 1] - Mu[, 2])^2
    MSW <- (S11 + S22) / 2
    
    good <- which(MSB > MSW)
    n_good <- length(good)
    
    if (n_good > 0L) {
      zero_streak <- 0L
      
      # Paired lnM* uses Delta/(n) inside the log
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
    
    total <- total + chunk
    
    # ---- adaptive chunking (same logic as independent) ----
    acc <- if (total > 0) kept / total else 0
    
    if (is_positive(acc)) {
      next_needed <- min(chunk_max, ceiling(max(0, min_kept - kept) / acc))
    } else {
      next_needed <- chunk * 2L
    }
    
    chunk <- as.integer(max(chunk_init, min(chunk_max, next_needed)))
  }
  
  list(
    point = if (kept >= 2) mean(lnM_star) else NA_real_,
    var   = if (kept >= 2) var(lnM_star)  else NA_real_,
    kept  = kept,
    total = total,
    attempts = attempts,
    status = status
  )
}


# ================================================================
#  Comparison wrapper: Delta-1 vs SAFE (paired or independent)
# ================================================================
#
# INPUTS
#   x1, x2, s1, s2 : summary stats
#   n1, n2         : if n2 is NULL, treat as paired (n1 = number of pairs)
#   r              : within-pair correlation required for paired
#   min_kept       : target number of usable SAFE draws
#   chunk_init, chunk_max, max_draws : SAFE controls
#
# OUTPUT
#   Prints a small table and SAFE diagnostics; invisibly returns the table.
#
compare_methods <- function(x1, x2, s1, s2, n1, n2 = NULL, r = NULL,
                            min_kept = 1e4, B = NULL,
                            chunk_init = 4e3, chunk_max = 2e6, max_draws = Inf) {
  
  # B is currently unused (kept for backward compatibility / future extension)
  # You can delete it if you prefer, but leaving it avoids breaking old code.
  
  if (is.null(n2)) {
    # Paired
    if (is.null(r)) stop("For paired data supply 'r' (within-pair correlation).")
    
    d1   <- lnM_delta1_dep(x1, x2, s1, s2, n1, r)
    safe <- safe_lnM_dep(x1, x2, s1, s2, n1, r,
                         min_kept = min_kept,
                         chunk_init = chunk_init,
                         chunk_max = chunk_max,
                         max_draws = max_draws)
    
  } else {
    # Independent
    d1   <- lnM_delta1_indep(x1, x2, s1, s2, n1, n2)
    safe <- safe_lnM_indep(x1, x2, s1, s2, n1, n2,
                           min_kept = min_kept,
                           chunk_init = chunk_init,
                           chunk_max = chunk_max,
                           max_draws = max_draws)
  }
  
  out <- dplyr::bind_rows(
    data.frame(method = "Delta-1",        lnM = d1["point"],    Var = d1["var"]),
    data.frame(method = "SAFE bootstrap", lnM = safe$point,     Var = safe$var)
  )
  
  print(out, row.names = FALSE, digits = 4)
  
  # SAFE diagnostics: “usable” = MSB > MSW
  pct <- 100 * safe$kept / max(1L, safe$total)
  cat(sprintf("\nSAFE kept %d of %d draws (%.2f%% usable) across %d chunk(s) [%s]\n",
              safe$kept, safe$total, pct, safe$attempts, safe$status))
  
  invisible(out)
}

# ================================================================
#  Quick demos (commented out)
# ================================================================
# set.seed(123)
#
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
