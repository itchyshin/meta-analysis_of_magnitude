# ================================================================
#  lnM: Delta-1 vs SAFE bootstrap (independent & paired)
#  - SAFE uses Chi-square / Wishart (no MVN for variances)
#  - Single delta version (delta-1) only
#  - Bounded while-loops via `max_chunks`
#  - Aligns with MS Equations (1)–(7) and SAFE §2.3
# ================================================================

library(MASS)    # for mvrnorm (paired means)
library(dplyr)

# ---------- tiny helpers ----------------------------------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(gap) ifelse(gap <= 0, NA_real_, gap)

# ---------- Delta-1 (independent)  -- MS Eqs (1)–(4), Var Eq (6) ----------
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2)
{
  h    <- n1 * n2 / (n1 + n2)                         # harmonic component
  MSB  <- h * (x1bar - x2bar)^2                       # Eq (1)
  MSW  <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)   # Eq (2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA, var = NA, se = NA))
  
  # lnM: Eq (4); here n0 = 2h for independent case
  lnM <- 0.5 * (log(Delta / (2 * h)) - log(MSW))
  
  # Delta-1 variance: MS Eq (6)
  sigmaD2 <- s1^2 / n1 + s2^2 / n2
  delta   <- x1bar - x2bar
  Var_B   <- h^2 * (2 * sigmaD2^2 + 4 * sigmaD2 * delta^2)
  Var_W   <- 2 * MSW^2 / (n1 + n2 - 2)
  
  gB   <- 0.5 / Delta
  gW   <- -0.5 * MSB / (Delta * MSW)
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ---------- Delta-1 (paired)  -- MS Eqs (3)–(4) specialized; Var Eq (7) ----
lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r)
{
  h     <- n / 2
  MSB   <- h * (x1bar - x2bar)^2                       # paired MSB
  MSW   <- (s1^2 + s2^2) / 2                           # paired MSW
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(point = NA, var = NA, se = NA))
  
  # lnM: Eq (4); here n0 = n for paired case
  lnM <- 0.5 * (log(Delta / n) - log(MSW))
  
  # Delta-1 variance: MS Eq (7)
  sigmaD2 <- s1^2 + s2^2 - 2 * r * s1 * s2
  delta   <- x1bar - x2bar
  Var_B   <- h^2 * (2 * sigmaD2^2 / n^2 + 4 * delta^2 * sigmaD2 / n)
  Var_W   <- (s1^4 + s2^4 + 2 * r^2 * s1^2 * s2^2) / (2 * (n - 1))
  
  gB   <- 0.5 / Delta
  gW   <- -0.5 * MSB / (Delta * MSW)
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ---------- SAFE bootstrap : independent (Chi-square) -----------------------
# Step 1–4 as in MS §2.3; means Normal, variances scaled Chi-square.
safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           B = 1e5, chunk = 5e3, max_chunks = 50)
{
  df1 <- n1 - 1L
  df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)   # same h as above
  
  lnM_star <- numeric(0)
  total    <- 0L      # attempted draws
  kept     <- 0L      # usable lnM*
  attempts <- 0L      # chunk iterations
  
  while (length(lnM_star) < B && attempts < max_chunks) {
    attempts <- attempts + 1L
    
    # Step 2: means ~ Normal; variances ~ Chi-square (always positive)
    m1 <- rnorm(chunk, mean = x1bar, sd = s1 / sqrt(n1))
    m2 <- rnorm(chunk, mean = x2bar, sd = s2 / sqrt(n2))
    v1 <- s1^2 * stats::rchisq(chunk, df = df1) / df1
    v2 <- s2^2 * stats::rchisq(chunk, df = df2) / df2
    
    total <- total + chunk
    
    # Step 3: transform
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
  
  list(point    = if (length(lnM_star)) mean(lnM_star) else NA_real_,
       var      = if (length(lnM_star)) var(lnM_star)  else NA_real_,
       kept     = kept,
       total    = total,
       attempts = attempts,
       status   = status)
}

# ---------- SAFE bootstrap : paired (Wishart) ------------------------------
# Step 1–4 as in MS §2.3; means bivariate Normal with Sigma/n; S ~ Wishart.
safe_lnM_dep <- function(x1bar, x2bar, s1, s2, n, r,
                         B = 1e5, chunk = 5e3, max_chunks = 50)
{
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
    
    # Step 2: means (μ*) and covariances (S*) draws
    Mu <- MASS::mvrnorm(n = chunk, mu = c(x1bar, x2bar), Sigma = Sig / n)
    W  <- stats::rWishart(n = chunk, df = df, Sigma = Sig)   # 2x2xchunk array
    S11 <- W[1,1,] / df
    S22 <- W[2,2,] / df
    
    total <- total + chunk
    
    # Step 3: transform
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
  
  list(point    = if (length(lnM_star)) mean(lnM_star) else NA_real_,
       var      = if (length(lnM_star)) var(lnM_star)  else NA_real_,
       kept     = kept,
       total    = total,
       attempts = attempts,
       status   = status)
}

# ---------- comparison wrapper (Delta-1 vs SAFE only) ----------------------
compare_methods <- function(x1, x2, s1, s2,
                            n1, n2 = NULL, r = NULL,
                            B = 1e5, chunk = 5e3, max_chunks = 50)
{
  if (is.null(n2)) {                       # paired design
    if (is.null(r)) stop("For paired data supply 'r' (within-pair correlation).")
    d1   <- lnM_delta1_dep  (x1, x2, s1, s2, n1, r)
    safe <- safe_lnM_dep    (x1, x2, s1, s2, n1, r, B, chunk, max_chunks)
  } else {                                 # independent design
    d1   <- lnM_delta1_indep (x1, x2, s1, s2, n1, n2)
    safe <- safe_lnM_indep   (x1, x2, s1, s2, n1, n2, B, chunk, max_chunks)
  }
  
  out <- dplyr::bind_rows(
    data.frame(method = "Delta-1",         lnM = d1["point"],  Var = d1["var"]),
    data.frame(method = "SAFE bootstrap",  lnM = safe$point,   Var = safe$var)
  )
  
  print(out, row.names = FALSE, digits = 4)
  pct <- 100 * safe$kept / max(1L, safe$total)
  cat(sprintf("\nSAFE kept %d of %d draws (%.2f %% usable) in %d chunk(s) [%s]\n",
              safe$kept, safe$total, pct, safe$attempts, safe$status))
  invisible(out)
}

# ---------- quick demo ------------------------------------------
set.seed(123)

# Independent example
compare_methods(x1 = 5, x2 = 5,
                s1 = 1,  s2 = 1,
                n1 = 40, n2 = 40,
                B  = 5e4, chunk = 4000, max_chunks = 50)


compare_methods(x1 = 7, x2 = 5,
                s1 = 1,  s2 = 1,
                n1 = 40, n2 = 40,
                B  = 5e4, chunk = 4000, max_chunks = 50)

compare_methods(x1 = 8, x2 = 5,
                s1 = 2,  s2 = 2,
                n1 = 40, n2 = 40,
                B  = 5e4, chunk = 4000, max_chunks = 50)

# Paired example (r = 0.8)
compare_methods(x1 = 5, x2 = 5,
                s1 = 1,  s2 = 1,
                n1 = 40, r  = 0.8,
                B  = 5e4, chunk = 4000, max_chunks = 50)

compare_methods(x1 = 7, x2 = 5,
                s1 = 1,  s2 = 1,
                n1 = 40, r  = 0.8,
                B  = 5e4, chunk = 4000, max_chunks = 50)

compare_methods(x1 = 8, x2 = 5,
                s1 = 2,  s2 = 2,
                n1 = 40, r  = 0.8,
                B  = 5e4, chunk = 4000, max_chunks = 50)
