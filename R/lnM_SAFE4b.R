## ==========================================================
##  Robust helpers for ln M   (guaranteed no NaN)
## ==========================================================

# ---------- 0.  House-keeping helpers ----------------------

stop_if_any_na <- function(...) {
  x <- unlist(list(...), use.names = FALSE)
  if (any(is.na(x)))
    stop("One or more inputs are NA; cannot proceed.")
}

# Positive gap check (used twice)
safe_gap <- function(gap, name = "gap") {
  if (gap <= 0) {
    warning(
      sprintf(
        "%s: MS_B ≤ MS_W → ln M undefined, returning NA", 
        name          # <<–– placeholder filled here
      ),
      call. = FALSE
    )
    return(NA_real_)
  }
  gap
}
# Fast, silent pmax for numerical noise
posify <- function(x, eps = 1e-12) pmax(x, eps)

# ---------- 1.  delta-1 estimator & variance ---------------

lnM_delta1 <- function(x1bar, x2bar, s1, s2, n1, n2) {
  
  ## -- basic sanity ----------------------------------------------------------
  stop_if_any_na(x1bar, x2bar, s1, s2, n1, n2)
  if (min(n1, n2) <= 1L)
    stop("Need n1 > 1 and n2 > 1 for variance calculations.")
  if (min(s1, s2) <= 0)
    stop("Both group SDs must be strictly positive.")
  
  ## -- ANOVA building blocks -------------------------------------------------
  h   <- n1 * n2 / (n1 + n2)
  MSB <- h * (x1bar - x2bar)^2
  MSW <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  
  Delta <- safe_gap(MSB - MSW, "delta-1")
  if (is.na(Delta))
    return(c(point = NA_real_, var = NA_real_, se = NA_real_))
  
  n0  <- 2 * h                         # same algebraically
  lnM <- 0.5 * (log(Delta / n0) - log(MSW))
  
  ## -- first-order delta variance -------------------------------------------
  #   Var(MSB)
  sigmaD2 <- s1^2 / n1 + s2^2 / n2          # Var(X̄₁ − X̄₂)
  delta   <- x1bar - x2bar
  Var_MSB <- h^2 * (2 * sigmaD2^2 + 4 * sigmaD2 * delta^2)
  
  #   Var(MSW)  [exact for χ²]
  Var_MSW <- 2 * MSW^2 / (n1 + n2 - 2)
  
  #   Cov(MSB, MSW)  ≈ 0 under independence of numerator & denominator
  #   (delta-1 assumes it, see Sokal & Rohlf 1995, Box 13.5)
  
  gB <-  0.5 / Delta
  gW <- -0.5 * MSB / (Delta * MSW)
  
  Var1 <- posify(gB^2 * Var_MSB + gW^2 * Var_MSW)
  
  c(point = lnM,
    var   = Var1,
    se    = sqrt(Var1))
}

# ---------- 2.  SAFE parametric bootstrap -----------------

safe_lnM <- function(x1bar, x2bar, s1, s2, n1, n2,
                     B       = 3e5,
                     chunk   = 5e4,
                     max_iter = 500) {
  
  stop_if_any_na(x1bar, x2bar, s1, s2, n1, n2)
  if (min(n1, n2) <= 1L)
    stop("Need n1 > 1 and n2 > 1 for bootstrap.")
  if (min(s1, s2) <= 0)
    stop("Both group SDs must be strictly positive.")
  
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package ‘MASS’ is required for mvrnorm().")
  
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- diag(c(s1^2 / n1,
                s2^2 / n2,
                2 * s1^4 / (n1 - 1),
                2 * s2^4 / (n2 - 1)))
  
  lnM_star <- numeric(0)
  iter <- 0L
  
  while (length(lnM_star) < B && iter < max_iter) {
    iter <- iter + 1L
    draws <- MASS::mvrnorm(chunk, mu, Sig)
    
    ## keep only physically admissible draws -------------------------------
    keep <- draws[, 3] > 0 & draws[, 4] > 0
    if (!any(keep)) next
    draws <- draws[keep, , drop = FALSE]
    
    m1 <- draws[, 1]; m2 <- draws[, 2]
    v1 <- draws[, 3]; v2 <- draws[, 4]
    
    MSB <- (n1 * n2) / (n1 + n2) * (m1 - m2)^2
    MSW <- ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    
    good <- MSB > MSW                      # ln M well-defined
    if (!any(good)) next
    
    n0 <- 2 * n1 * n2 / (n1 + n2)
    lnM_draws <- 0.5 * (log((MSB[good] - MSW[good]) / n0) -
                          log(MSW[good]))
    
    lnM_star <- c(lnM_star, lnM_draws)
  }
  
  if (length(lnM_star) < B) {
    warning("SAFE: reached iteration limit before collecting B draws; ",
            "returning NA.", call. = FALSE)
    return(c(point = NA_real_, var = NA_real_, se = NA_real_))
  }
  
  lnM_star <- lnM_star[seq_len(B)]         # truncate to exactly B
  se_SAFE  <- sd(lnM_star)
  
  ## ---------- bias anchoring ---------------------------------------------
  d1 <- lnM_delta1(x1bar, x2bar, s1, s2, n1, n2)
  lnM_hat <- if (is.na(d1["point"])) mean(lnM_star) else unname(d1["point"])
  
  lnM_BC <- 2 * lnM_hat - mean(lnM_star)   # bias-corrected
  
  c(point = lnM_BC,
    var   = se_SAFE^2,
    se    = se_SAFE)
}


set.seed(1)
## -------- 2.  Worked example ------------------------------------------
x1bar <- 20.3;  s1 <- 2.7;  n1 <- 200
x2bar <- 13.1;  s2 <- 2.1;  n2 <- 200

delta1 <- lnM_delta1(x1bar,x2bar,s1,s2,n1,n2)
safe   <- safe_lnM   (x1bar,x2bar,s1,s2,n1,n2, B = 1e6)  # half-million draws

rbind(delta1  = round(delta1,  4),
      SAFE_BC = round(safe[1:3],4))