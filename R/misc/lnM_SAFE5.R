## ==========================================================
##  lnM  (log SD–ratio) : delta-1, delta-2, SAFE bootstrap
##  last update: 19-Jun-2025
## ==========================================================

## -------- 0.  Delta-1 helper (correct) --------------------
lnM_delta1 <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h   <- n1*n2/(n1+n2)
  MB  <- h*(x1bar - x2bar)^2
  MW  <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  n0  <- 2*n1*n2/(n1+n2)
  lnM <- log( sqrt((MB-MW)/n0) ) - log( sqrt(MW) )
  
  sigmaD2 <- s1^2/n1 + s2^2/n2
  delta   <- x1bar - x2bar
  Var_MB  <- h^2 * ( 2*sigmaD2^2 + 4*sigmaD2*delta^2 )
  Var_MW  <- 2*MW^2/(n1+n2-2)
  
  gB <-  1/(2*(MB-MW))
  gW <- -MB/(2*(MB-MW)*MW)
  Var1 <- gB^2*Var_MB + gW^2*Var_MW
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

## -------- 1.  Delta-2 helper (bias-corrected) -------------
lnM_delta2 <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h   <- n1*n2/(n1+n2)
  MB  <- h*(x1bar - x2bar)^2
  MW  <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  n0  <- 2*n1*n2/(n1+n2)
  lnM <- log( sqrt((MB-MW)/n0) ) - log( sqrt(MW) )
  
  sigmaD2 <- s1^2/n1 + s2^2/n2
  delta   <- x1bar - x2bar
  Var_MB  <- h^2 * ( 2*sigmaD2^2 + 4*sigmaD2*delta^2 )
  Var_MW  <- 2*MW^2/(n1+n2-2)
  
  Delta <- MB - MW
  gB  <-  1/(2*Delta)
  gW  <- -MB/(2*Delta*MW)
  gBB <- -0.5/Delta^2
  gWW <-  0.5*(1/Delta^2 + 1/MW^2)
  
  ## bias correction
  bias <- 0.25*((Var_MW - Var_MB)/Delta^2 + Var_MW/MW^2)
  lnM_BC <- lnM - bias
  
  ## variance (delta-2)
  Var1 <- gB^2*Var_MB + gW^2*Var_MW
  Var2 <- Var1 + 0.5*( gBB^2*Var_MB^2 + gWW^2*Var_MW^2 )
  
  c(point = lnM_BC, var = Var2, se = sqrt(Var2))
}

## -------- SAFE bootstrap helper (fixed) -------------------
safe_lnM <- function(x1bar, x2bar, s1, s2, n1, n2,
                     B = 2e5, chunk = 5e5) {
  
  out <- numeric(0)
  
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- diag(c(s1^2/n1, s2^2/n2,
                2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  library(MASS)
  
  while (length(out) < B) {
    draws <- mvrnorm(chunk, mu, Sig)
    pos   <- draws[,3] > 0 & draws[,4] > 0
    if (!any(pos)) next
    draws <- draws[pos,]
    
    m1 <- draws[,1]; m2 <- draws[,2]
    v1 <- draws[,3]; v2 <- draws[,4]
    
    MB <- (n1*n2)/(n1+n2) * (m1 - m2)^2
    MW <- ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2)
    good <- MB > MW + tiny
    if (!any(good)) next
    
    MB <- MB[good];  MW <- MW[good]
    n0 <- 2*n1*n2/(n1+n2)
    
    out <- c(out,
             log( sqrt((MB - MW)/n0) ) - log( sqrt(MW) ))
  }
  
  lnM_star <- out[1:B]
  se_SAFE  <- sd(lnM_star)
  
  ## -----------   bias correction (names stripped) ---------
  lnM_hat <- unname( lnM_delta1(
    x1bar,x2bar,s1,s2,n1,n2)["point"] )
  bias    <- mean(lnM_star) - lnM_hat
  lnM_BC  <- lnM_hat - bias           # numeric, no name
  
  c(point = lnM_BC,
    var   = se_SAFE^2,
    se    = se_SAFE)
}

## -------- 3.  Example data ---------------------------------
x1bar <- 15.3;  s1 <- 2.7;  n1 <- 10
x2bar <- 17.1;  s2 <- 2.1;  n2 <- 4

## -------- 4.  Run all three estimators ---------------------
delta1 <- lnM_delta1(x1bar,x2bar,s1,s2,n1,n2)
delta2 <- lnM_delta2(x1bar,x2bar,s1,s2,n1,n2)
safe   <- safe_lnM(x1bar,x2bar,s1,s2,n1,n2, B = 1e6)

## -------- 5.  Tidy comparison table ------------------------
library(dplyr)
bind_rows(delta1  = delta1,
          delta2  = delta2,
          SAFE_BC = safe,
          .id = "method") |>
  mutate(across(-method, round, 5))