## ==========================================================
##  lnM:  independent vs. dependent (paired) designs
##  delta-1 (first-order) variance   •    SAFE bootstrap
## ==========================================================

library(MASS)     # for mvrnorm()
library(dplyr)    # tidy output

## ----------------------------------------------------------
##  Helper to guarantee positive variances
posify <- function(x, eps = 1e-12) pmax(x, eps)

## ----------------------------------------------------------
##  1)  FIRST-ORDER (DELTA-1) FORMULAS
## ----------------------------------------------------------

## ---------- 1a. independent groups ------------------------
lnM_delta1_indep <- function(x1bar, x2bar, s1, s2, n1, n2)
{
  h   <- n1*n2/(n1+n2)
  MSB <- h*(x1bar - x2bar)^2
  MSW <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1 + n2 - 2)
  Delta <- MSB - MSW
  lnM   <- 0.5*(log(Delta/(2*h)) - log(MSW))
  
  ## delta-method variance
  sigmaD2 <- s1^2/n1 + s2^2/n2
  delta   <- x1bar - x2bar
  Var_MSB <- h^2*(2*sigmaD2^2 + 4*sigmaD2*delta^2)
  Var_MSW <- 2*MSW^2/(n1 + n2 - 2)
  
  gB <- 0.5/Delta
  gW <- -0.5*MSB/(Delta*MSW)
  Var1 <- posify(gB^2*Var_MSB + gW^2*Var_MSW)
  
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

## ---------- 1b.  CORRECT delta-1 for paired groups ----------
lnM_delta1_dep <- function(x1bar, x2bar, s1, s2, n, r)
{
  h   <- n/2
  MSB <- h * (x1bar - x2bar)^2
  MSW <- (s1^2 + s2^2) / 2
  Delta <- MSB - MSW
  lnM   <- 0.5 * (log(Delta / n) - log(MSW))
  
  ## correct Var(MSB)
  sigmaD2 <- s1^2 + s2^2 - 2 * r * s1 * s2
  delta   <- x1bar - x2bar
  Var_MSB <- h^2 * ( 2 * sigmaD2^2 / n^2 + 4 * delta^2 * sigmaD2 / n )
  
  ## Var(MSW) unchanged
  Var_MSW <- (s1^4 + s2^4 + 2 * r^2 * s1^2 * s2^2) / (2 * (n - 1))
  
  gB <- 0.5 / Delta
  gW <- -0.5 * MSB / (Delta * MSW)
  
  Var1 <- posify(gB^2 * Var_MSB + gW^2 * Var_MSW)
  c(point = lnM, var = Var1, se = sqrt(Var1))
}

## ----------------------------------------------------------
##  2)  SAFE PARAMETRIC BOOTSTRAPS
## ----------------------------------------------------------

## ---------- 2a. independent groups ------------------------
safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           B = 2e5, chunk = 5e4)
{
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- diag(c(s1^2/n1, s2^2/n2,
                2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  
  lnM_star <- numeric(0)
  while (length(lnM_star) < B) {
    draws <- mvrnorm(chunk, mu, Sig)
    goodV <- draws[,3] > 0 & draws[,4] > 0
    if (!any(goodV)) next
    d <- draws[goodV,, drop=FALSE]
    m1 <- d[,1]; m2 <- d[,2]; v1 <- d[,3]; v2 <- d[,4]
    MSB <- (n1*n2)/(n1+n2)*(m1 - m2)^2
    MSW <- ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2)
    good <- MSB > MSW
    if (!any(good)) next
    lnM_star <- c(lnM_star,
                  0.5*(log((MSB[good]-MSW[good])/(2*n1*n2/(n1+n2))) -
                         log(MSW[good])))
  }
  lnM_star <- lnM_star[seq_len(B)]
  c(point = mean(lnM_star), var = var(lnM_star), se = sd(lnM_star))
}

## ---------- 2b. dependent / paired groups -----------------
safe_lnM_dep <- function(x1bar, x2bar, s1, s2, n, r,
                         B = 2e5, chunk = 5e4)
{
  ## 4-vector: (mean1, mean2, var1, var2)
  mu <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- matrix(0, 4, 4)
  Sig[1,1] <- s1^2/n
  Sig[2,2] <- s2^2/n
  Sig[1,2] <- Sig[2,1] <- r*s1*s2/n
  Sig[3,3] <- 2*s1^4/(n-1)
  Sig[4,4] <- 2*s2^4/(n-1)
  Sig[3,4] <- Sig[4,3] <- 2*r^2*s1^2*s2^2/(n-1)
  
  lnM_star <- numeric(0)
  while (length(lnM_star) < B) {
    draws <- mvrnorm(chunk, mu, Sig)
    goodV <- draws[,3] > 0 & draws[,4] > 0
    if (!any(goodV)) next
    d <- draws[goodV,, drop=FALSE]
    m1 <- d[,1]; m2 <- d[,2]; v1 <- d[,3]; v2 <- d[,4]
    MSB <- (n/2)*(m1 - m2)^2
    MSW <- (v1 + v2)/2
    good <- MSB > MSW
    if (!any(good)) next
    lnM_star <- c(lnM_star,
                  0.5*(log((MSB[good]-MSW[good])/n) - log(MSW[good])))
  }
  lnM_star <- lnM_star[seq_len(B)]
  c(point = mean(lnM_star), var = var(lnM_star), se = sd(lnM_star))
}

## ----------------------------------------------------------
##  3)  WORKED EXAMPLE
## ----------------------------------------------------------

## independent design
x1_i <- 14.3;  s1_i <- 2.7;  n1 <- 10
x2_i <- 14.1;  s2_i <- 2.1;  n2 <- 10

## dependent design (same summary numbers, add correlation)
x1_d <- 14.3;  s1_d <- 2.7;  n  <- 10
x2_d <- 14.1;  s2_d <- 2.1;  r  <- 0.8   # choose any  |r|<1

## run all four estimators ----------------------------------

indep_delta <- lnM_delta1_indep(x1_i,x2_i,s1_i,s2_i,n1,n2)
indep_safe  <- safe_lnM_indep (x1_i,x2_i,s1_i,s2_i,n1,n2, B=1e6)

dep_delta   <- lnM_delta1_dep (x1_d,x2_d,s1_d,s2_d,n ,r )
dep_safe    <- safe_lnM_dep   (x1_d,x2_d,s1_d,s2_d,n ,r , B=1e6)

bind_rows(
  data.frame(method="delta-1 (indep)", t(indep_delta)),
  data.frame(method="SAFE    (indep)", t(indep_safe )),
  data.frame(method="delta-1 (dep  )", t(dep_delta  )),
  data.frame(method="SAFE    (dep  )", t(dep_safe   ))
) |>
  mutate(across(-method, round, 4))