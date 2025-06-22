#===============================================================
#  lnM: Δ-method vs SAFE bootstrap (independent groups only)
#===============================================================
#     Shinichi-style SAFE utilities, June 2025
#---------------------------------------------------------------
#  Requires no external packages except knitr for nicer tables
#  (remove knitr::kable if you prefer plain data.frame printing)
#===============================================================

# ---------- Analytic first-order Δ-method variance -------------
lnM_delta <- function(x1bar, x2bar, s1, s2, n1, n2) {
  
  # pooled within-group variance
  nu   <- n1 + n2 - 2
  sW2  <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / nu
  
  # between-group variance (Lynch & Walsh 1998, eq. 25-3)
  n0   <- 2 * n1 * n2 / (n1 + n2)
  MSB  <- (n1 * n2 / (n1 + n2)) * (x1bar - x2bar)^2
  sB2  <- (MSB - sW2) / n0
  if (sB2 <= 0)
    stop("sB^2 ≤ 0; lnM undefined for these data")
  
  # point estimate
  lnM  <- 0.5 * log(sB2) - 0.5 * log(sW2)
  
  # Δ-method variance components
  c_   <- 1 / n1 + 1 / n2
  var1 <- ((x1bar - x2bar)^2 * sW2) / (4 * sB2^2) * c_
  var2 <- (sW2^2 / (2 * nu)) * (1 / (n0 * sB2) + 1 / sW2)^2
  var  <- var1 + var2
  
  list(lnM = lnM, var = var)
}

# ---------- SAFE parametric bootstrap -------------------------
lnM_safe <- function(x1bar, x2bar, s1, s2, n1, n2, B = 20000) {
  
  nu   <- n1 + n2 - 2
  n0   <- 2 * n1 * n2 / (n1 + n2)
  lnM_vec <- numeric(B)
  
  for (b in seq_len(B)) {
    
    # (1) bootstrap means
    m1 <- rnorm(1, x1bar, s1 / sqrt(n1))
    m2 <- rnorm(1, x2bar, s2 / sqrt(n2))
    
    # (2) bootstrap variances
    s1b2 <- (n1 - 1) * s1^2 / rchisq(1, df = n1 - 1)
    s2b2 <- (n2 - 1) * s2^2 / rchisq(1, df = n2 - 1)
    
    # (3) derived quantities
    sW2b <- ((n1 - 1) * s1b2 + (n2 - 1) * s2b2) / nu
    MSBb <- (n1 * n2 / (n1 + n2)) * (m1 - m2)^2
    sB2b <- (MSBb - sW2b) / n0
    
    lnM_vec[b] <- if (sB2b > 0)
      0.5 * log(sB2b) - 0.5 * log(sW2b) else NA_real_
  }
  
  lnM_vec <- lnM_vec[!is.na(lnM_vec)]
  
  list(lnM = mean(lnM_vec),
       var = var(lnM_vec),
       draws_kept = length(lnM_vec))
}

# ---------- Side-by-side comparator ---------------------------
compare_methods <- function(x1bar, x2bar, s1, s2, n1, n2,
                            B = 20000, digits = 4) {
  
  delta <- lnM_delta(x1bar, x2bar, s1, s2, n1, n2)
  safe  <- lnM_safe (x1bar, x2bar, s1, s2, n1, n2, B)
  
  out <- data.frame(
    method = c("Delta", "SAFE bootstrap"),
    lnM    = c(delta$lnM, safe$lnM),
    Var    = c(delta$var,  safe$var)
  )
  
  if (requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(out, digits = digits))
  } else {
    print(round(out, digits))
  }
  
  cat(sprintf(
    "\nSAFE kept %d / %d draws (%.2f %% usable)\n\n",
    safe$draws_kept, B, 100 * safe$draws_kept / B))
  
  invisible(out)
}

#===============================================================
# Example usage  (delete / replace with your own study numbers)
#===============================================================
set.seed(1)  # reproducible SAFE draws

n1 <- 120; x1 <- 10.3; s1 <- 1.9
n2 <- 140; x2 <-  7.8; s2 <- 2.4

compare_methods(x1, x2, s1, s2, n1, n2, B = 20000)