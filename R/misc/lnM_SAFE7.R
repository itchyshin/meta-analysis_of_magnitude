# ================================================================
#  SAFE delta-1 vs. bootstrap for ln M   (independent & paired)
#  ––– now counts usable vs. attempted bootstrap draws
# ================================================================

library(MASS)
library(dplyr)

# ---------- tiny helpers ----------------------------------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(gap) ifelse(gap <= 0, NA_real_, gap)

# ---------- delta-1 : independent --------------------------------
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2)
{
  h   <- n1 * n2 / (n1 + n2)
  MSB <- h * (x1bar - x2bar)^2
  MSW <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  Δ   <- safe_gap(MSB - MSW)
  if (is.na(Δ)) return(c(point = NA, var = NA, se = NA))
  
  lnM <- 0.5 * (log(Δ / (2 * h)) - log(MSW))
  
  σD2   <- s1^2 / n1 + s2^2 / n2
  δ     <- x1bar - x2bar
  Var_B <- h^2 * (2 * σD2^2 + 4 * σD2 * δ^2)
  Var_W <- 2 * MSW^2 / (n1 + n2 - 2)
  
  gB <- 0.5 / Δ
  gW <- -0.5 * MSB / (Δ * MSW)
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ---------- delta-1 : paired ------------------------------------
lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r)
{
  h   <- n / 2
  MSB <- h * (x1bar - x2bar)^2
  MSW <- (s1^2 + s2^2) / 2
  Δ   <- safe_gap(MSB - MSW)
  if (is.na(Δ)) return(c(point = NA, var = NA, se = NA))
  
  lnM <- 0.5 * (log(Δ / n) - log(MSW))
  
  σD2   <- s1^2 + s2^2 - 2 * r * s1 * s2
  δ     <- x1bar - x2bar
  Var_B <- h^2 * (2 * σD2^2 / n^2 + 4 * δ^2 * σD2 / n)
  Var_W <- (s1^4 + s2^4 + 2 * r^2 * s1^2 * s2^2) / (2 * (n - 1))
  
  gB <- 0.5 / Δ
  gW <- -0.5 * MSB / (Δ * MSW)
  Var1 <- posify(gB^2 * Var_B + gW^2 * Var_W)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ---------- SAFE bootstrap : independent ------------------------
safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           B = 1e5, chunk = 5e3)
{
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Σ   <- diag(c(s1^2/n1, s2^2/n2,
                2 * s1^4 / (n1 - 1), 2 * s2^4 / (n2 - 1)))
  
  lnM_star <- numeric(0)
  total    <- 0L     # generated draws
  kept     <- 0L     # draws that gave usable ln M
  
  while (length(lnM_star) < B) {
    
    draws0     <- MASS::mvrnorm(chunk, mu, Σ)
    total      <- total + nrow(draws0)
    good_var   <- draws0[,3] > 0 & draws0[,4] > 0
    draws      <- draws0[good_var, , drop = FALSE]
    if (!nrow(draws)) next
    
    m1 <- draws[,1]; m2 <- draws[,2]
    v1 <- draws[,3]; v2 <- draws[,4]
    
    MSB  <- (n1 * n2) / (n1 + n2) * (m1 - m2)^2
    MSW  <- ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    good <- MSB > MSW
    kept <- kept + sum(good)
    
    lnM_star <- c(lnM_star,
                  0.5 * (log((MSB[good] - MSW[good]) /
                               (2 * n1 * n2 / (n1 + n2))) -
                           log(MSW[good])))
  }
  
  lnM_star <- lnM_star[seq_len(B)]
  list(point = mean(lnM_star),
       var   = var(lnM_star),
       kept  = kept,
       total = total)
}

# ---------- SAFE bootstrap : paired -----------------------------
safe_lnM_dep <- function(x1bar, x2bar, s1, s2, n, r,
                         B = 1e5, chunk = 5e3)
{
  mu <- c(x1bar, x2bar, s1^2, s2^2)
  Σ  <- matrix(0, 4, 4)
  Σ[1,1] <- s1^2 / n
  Σ[2,2] <- s2^2 / n
  Σ[1,2] <- Σ[2,1] <- r * s1 * s2 / n
  Σ[3,3] <- 2 * s1^4 / (n - 1)
  Σ[4,4] <- 2 * s2^4 / (n - 1)
  Σ[3,4] <- Σ[4,3] <- 2 * r^2 * s1^2 * s2^2 / (n - 1)
  
  lnM_star <- numeric(0)
  total    <- 0L
  kept     <- 0L
  
  while (length(lnM_star) < B) {
    
    draws0   <- MASS::mvrnorm(chunk, mu, Σ)
    total    <- total + nrow(draws0)
    good_var <- draws0[,3] > 0 & draws0[,4] > 0
    draws    <- draws0[good_var, , drop = FALSE]
    if (!nrow(draws)) next
    
    m1 <- draws[,1]; m2 <- draws[,2]
    v1 <- draws[,3]; v2 <- draws[,4]
    
    MSB  <- (n / 2) * (m1 - m2)^2
    MSW  <- (v1 + v2) / 2
    good <- MSB > MSW
    kept <- kept + sum(good)
    
    lnM_star <- c(lnM_star,
                  0.5 * (log((MSB[good] - MSW[good]) / n) -
                           log(MSW[good])))
  }
  
  lnM_star <- lnM_star[seq_len(B)]
  list(point = mean(lnM_star),
       var   = var(lnM_star),
       kept  = kept,
       total = total)
}

# ---------- delta-SD (shared-df, first order) --------------------
# *same* point-estimate as before, but uses the simple
# Var(lnM) = 1/df      (independent)
#             = (1-r^2)/df  (paired)

lnM_deltaSD_indep <- function(x1bar, x2bar, s1, s2, n1, n2)
{
  h   <- n1 * n2 / (n1 + n2)
  MSB <- h * (x1bar - x2bar)^2
  MSW <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  Δ   <- if (MSB > MSW) (MSB - MSW) else NA_real_
  if (is.na(Δ)) return(c(point = NA, var = NA, se = NA))
  
  lnM <- 0.5 * (log(Δ / (2 * h)) - log(MSW))
  df  <- n1 + n2 - 2
  Var <- 1 / df
  c(point = lnM, var = Var, se = sqrt(Var))
}

lnM_deltaSD_dep <- function(x1bar, x2bar, s1, s2, n, r)
{
  MSB <- (n / 2) * (x1bar - x2bar)^2
  MSW <- (s1^2 + s2^2) / 2
  Δ   <- if (MSB > MSW) (MSB - MSW) else NA_real_
  if (is.na(Δ)) return(c(point = NA, var = NA, se = NA))
  
  lnM <- 0.5 * (log(Δ / n) - log(MSW))
  df  <- n - 2
  Var <- (1 - r^2) / df           # shared-df, paired version
  c(point = lnM, var = Var, se = sqrt(Var))
}

# ---------- comparison wrapper  (now prints three methods) -------
compare_methods <- function(x1, x2, s1, s2,
                            n1, n2 = NULL, r = NULL,
                            B = 1e5, chunk = 5e3)
{
  if (is.null(n2)) {                       # paired design
    if (is.null(r))
      stop("For paired data supply 'r' (within-pair correlation).")
    dSD  <- lnM_deltaSD_dep (x1, x2, s1, s2, n1, r)
    d1   <- lnM_delta1_dep  (x1, x2, s1, s2, n1, r)
    safe <- safe_lnM_dep    (x1, x2, s1, s2, n1, r, B, chunk)
  } else {                                 # independent design
    dSD  <- lnM_deltaSD_indep(x1, x2, s1, s2, n1, n2)
    d1   <- lnM_delta1_indep (x1, x2, s1, s2, n1, n2)
    safe <- safe_lnM_indep   (x1, x2, s1, s2, n1, n2, B, chunk)
  }
  
  out <- dplyr::bind_rows(
    data.frame(method = "Delta-SD",        lnM = dSD["point"], Var = dSD["var"]),
    data.frame(method = "Delta-1",         lnM = d1["point"],  Var = d1["var"]),
    data.frame(method = "SAFE bootstrap",  lnM = safe$point,   Var = safe$var)
  )
  
  print(out, row.names = FALSE, digits = 4)
  pct <- 100 * safe$kept / safe$total
  cat(sprintf("\nSAFE kept %d of %d draws (%.2f %% usable)\n",
              safe$kept, safe$total, pct))
  invisible(out)
}

# ---------- quick demo ------------------------------------------
set.seed(123)

# Independent example
compare_methods(x1 = 5, x2 = 5,
                s1 = 1,  s2 = 1,
                n1 = 40,   n2 = 40,
                B  = 1e5, chunk = 4000)

# Paired example (r = 0.8)
compare_methods(x1 = 5, x2 = 5,
                s1 = 1,  s2 = 1,
                n1 = 40,   r  = 0.8,
                B  = 1e5, chunk = 4000)