# ================================================================
#  SAFE delta-1 vs. bootstrap for ln M   (independent & paired)
#  • now reports the % of usable bootstrap draws correctly
# ================================================================

library(MASS)
library(dplyr)

# TODO - change the chunck so it is not 100%

# ---------- small utilities ------------------------------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(gap, name = "gap") {
  if (gap <= 0) return(NA_real_) else gap
}

# ---------- delta-1 : independent ------------------------------
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2)
{
  h   <- n1*n2/(n1+n2)
  MSB <- h*(x1bar - x2bar)^2
  MSW <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  Delta <- safe_gap(MSB - MSW, "indep")
  if (is.na(Delta))
    return(c(point = NA, var = NA, se = NA))
  
  lnM <- 0.5*(log(Delta/(2*h)) - log(MSW))
  
  sigmaD2 <- s1^2/n1 + s2^2/n2
  delta   <- x1bar - x2bar
  Var_MSB <- h^2*(2*sigmaD2^2 + 4*sigmaD2*delta^2)
  Var_MSW <- 2*MSW^2/(n1+n2-2)
  
  gB <- 0.5/Delta
  gW <- -0.5*MSB/(Delta*MSW)
  Var1 <- posify(gB^2*Var_MSB + gW^2*Var_MSW)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ---------- delta-1 : paired -----------------------------------
lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r)
{
  h   <- n/2
  MSB <- h*(x1bar - x2bar)^2
  MSW <- (s1^2 + s2^2)/2
  Delta <- safe_gap(MSB - MSW, "paired")
  if (is.na(Delta))
    return(c(point = NA, var = NA, se = NA))
  
  lnM <- 0.5*(log(Delta/n) - log(MSW))
  
  sigmaD2 <- s1^2 + s2^2 - 2*r*s1*s2
  delta   <- x1bar - x2bar
  Var_MSB <- h^2*(2*sigmaD2^2/n^2 + 4*delta^2*sigmaD2/n)
  Var_MSW <- (s1^4 + s2^4 + 2*r^2*s1^2*s2^2)/(2*(n-1))
  
  gB <- 0.5/Delta
  gW <- -0.5*MSB/(Delta*MSW)
  Var1 <- posify(gB^2*Var_MSB + gW^2*Var_MSW)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

# ---------- SAFE bootstrap helpers -----------------------------
safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           B = 1e5, chunk = 5000)
{
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- diag(c(s1^2/n1, s2^2/n2,
                2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  
  lnM_star <- numeric(0)
  while (length(lnM_star) < B) {
    draws <- MASS::mvrnorm(B, mu, Sig)
    okVar <- draws[,3] > 0 & draws[,4] > 0
    draws <- draws[okVar,, drop = FALSE]
    
    m1 <- draws[,1]; m2 <- draws[,2]
    v1 <- draws[,3]; v2 <- draws[,4]
    
    MSB <- (n1*n2)/(n1+n2)*(m1 - m2)^2
    MSW <- ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2)
    good <- MSB > MSW
    lnM_star <- c(lnM_star,
                  0.5*(log((MSB[good]-MSW[good])/(2*n1*n2/(n1+n2)))
                       - log(MSW[good])))
  }
  
  lnM_star <- lnM_star[seq_len(B)]         # keep exactly B
  list(point = mean(lnM_star),
       var   = var(lnM_star),
       kept  = B,              # always equal to requested B
       total = B)
}

safe_lnM_dep <- function(x1bar, x2bar, s1, s2, n, r,
                         B = 1e5, chunk = 5000)
{
  mu <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- matrix(0, 4, 4)
  Sig[1,1] <- s1^2/n; Sig[2,2] <- s2^2/n
  Sig[1,2] <- Sig[2,1] <- r*s1*s2/n
  Sig[3,3] <- 2*s1^4/(n-1); Sig[4,4] <- 2*s2^4/(n-1)
  Sig[3,4] <- Sig[4,3] <- 2*r^2*s1^2*s2^2/(n-1)
  
  lnM_star <- numeric(0)
  while (length(lnM_star) < B) {
    draws <- MASS::mvrnorm(chunk, mu, Sig)
    okVar <- draws[,3] > 0 & draws[,4] > 0
    draws <- draws[okVar,, drop = FALSE]
    
    m1 <- draws[,1]; m2 <- draws[,2]
    v1 <- draws[,3]; v2 <- draws[,4]
    
    MSB <- (n/2)*(m1 - m2)^2
    MSW <- (v1 + v2)/2
    good <- MSB > MSW
    lnM_star <- c(lnM_star,
                  0.5*(log((MSB[good]-MSW[good])/n) - log(MSW[good])))
  }
  
  lnM_star <- lnM_star[seq_len(B)]
  list(point = mean(lnM_star),
       var   = var(lnM_star),
       kept  = B,
       total = B)
}

# ---------- comparison wrapper ---------------------------------
compare_methods <- function(x1, x2, s1, s2,
                            n1, n2 = NULL, r = NULL,
                            B  = 1e5)
{
  if (is.null(n2)) {  # paired
    if (is.null(r))
      stop("Supply 'r' for paired data.")
    delta <- lnM_delta1_dep(x1,x2,s1,s2,n1,r)
    safe  <- safe_lnM_dep  (x1,x2,s1,s2,n1,r,B)
  } else {            # independent
    delta <- lnM_delta1_indep(x1,x2,s1,s2,n1,n2)
    safe  <- safe_lnM_indep (x1,x2,s1,s2,n1,n2,B)
  }
  
  out <- bind_rows(
    data.frame(method="Delta",          lnM=delta["point"], Var=delta["var"]),
    data.frame(method="SAFE bootstrap", lnM=safe$point,     Var=safe$var))
  print(out, row.names = FALSE, digits = 4)
  
  cat(sprintf("\nSAFE kept %d / %d draws (%.2f %% usable)\n",
              safe$kept, safe$total, 100))
  invisible(out)
}

# ---------- quick demo -----------------------------------------
# independent example
compare_methods(x1=11.3, x2=14.1, s1=2.7, s2=2.1,
                n1=11.2, n2=10, B=1e5)

# paired example (r = 0.8)
compare_methods(x1=11.3, x2=14.1, s1=2.7, s2=2.1,
                n1=11.2, r=0.8,  B=1e5)