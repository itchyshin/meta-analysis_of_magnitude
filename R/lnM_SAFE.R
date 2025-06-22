#==============================================================================  
# 1. Define safe_lnM() via 4‐variate sampling, with Delta1, Delta2, and SAFE
#==============================================================================  
safe_lnM <- function(x1bar, x2bar, s1, s2, n1, n2, B = 50000, min_surv = 100) {
  # 1a. degrees of freedom & scaling
  dfW <- n1 + n2 - 2
  n0  <- 2 * n1 * n2 / (n1 + n2)
  
  # 1b. plug‐in quantities
  W_hat <- ((n1-1)*s1^2 + (n2-1)*s2^2) / dfW
  SSB   <- (n1*n2)/(n1+n2) * (x1bar - x2bar)^2
  B_hat <- (SSB - W_hat) / n0
  if (W_hat <= 0 || B_hat <= 0) {
    warning("Non-positive W_hat or B_hat; returning all NA")
    return(data.frame(
      lnM = rep(NA_real_, 3),
      SE  = rep(NA_real_, 3),
      row.names = c("Delta1","Delta2","SAFE")
    ))
  }
  
  # first‐order plug‐in
  lnM   <- 0.5 * (log(B_hat) - log(W_hat))
  var1  <-  1/(n1 + n2 - 2)
  
  # second‐order analytic correction
  lnM2  <- lnM 
  #var2  <- var1 + 1/(2*n0) + 1/(2*dfW^2)
  var2 <- 1/(n1 + n2 - 2) + 1/(n1 + n2 - 2)^2
  
  # 2. SAFE bootstrap via 4‐variate normal
  requireNamespace("MASS", quietly=TRUE)
  mu    <- c(x1bar, x2bar, s1^2, s2^2)
  Sigma <- diag(c(
    s1^2/n1,
    s2^2/n2,
    2*s1^4/(n1-1),
    2*s2^4/(n2-1)
  ))
  sims <- MASS::mvrnorm(B, mu, Sigma)
  x1s  <- sims[,1]; x2s  <- sims[,2]
  v1s  <- sims[,3]; v2s  <- sims[,4]
  
  W_s   <- ((n1-1)*v1s + (n2-1)*v2s) / dfW
  SSB_s <- (n1*n2)/(n1+n2) * (x1s - x2s)^2
  B_s   <- (SSB_s - W_s) / n0
  
  ok <- which(W_s > 0 & B_s > 0)
  if (length(ok) < min_surv) {
    warning("Too few valid bootstraps: ", length(ok), " < min_surv=", min_surv)
    lnM_SAFE <- NA_real_
    SE_SAFE  <- NA_real_
  } else {
    lnM_star <- 0.5 * (log(B_s[ok]) - log(W_s[ok]))
    bias     <- mean(lnM_star) - lnM
    lnM_SAFE <- lnM - bias
    SE_SAFE  <- sd(lnM_star)
  }
  
  # 3. return all three estimates
  data.frame(
    lnM = c(lnM, lnM2, lnM_SAFE),
    SE  = c(sqrt(var1), sqrt(var2), SE_SAFE),
    row.names = c("Delta1","Delta2","SAFE")
  )
}

#==============================================================================  
# 2. Worked example
#==============================================================================  
set.seed(123)
res <- safe_lnM(
  x1bar = 15, x2bar = 30,
  s1    =  2,  s2    =  3,
  n1    = 2000,  n2    = 2000,
  B     = 50000
)

print(res)