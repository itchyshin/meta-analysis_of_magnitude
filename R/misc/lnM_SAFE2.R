safe_lnM <- function(x1bar, x2bar, s1, s2, n1, n2, r = 0, B = 20000, min_surv = 100) {
  # degrees of freedom
  dfW <- n1 + n2 - 2  # = n-2
  n0  <- 2 * n1 * n2 / (n1 + n2)
  
  # 1. plug‐in lnM
  MS_B <- (n1*n2)/(n1+n2) * (x1bar - x2bar)^2
  MS_W <- ((n1-1)*s1^2 + (n2-1)*s2^2) / dfW
  sB2  <- (MS_B - MS_W) / n0
  if (MS_W <= 0 || sB2 <= 0) {
    stop("Non‐positive MS_W or sB2; cannot compute lnM")
  }
  lnM <- log(sqrt(sB2) / sqrt(MS_W))
  
  # 2. analytic SEs
  # Delta1: independent case
  var2 <-  1/dfW + 1/(dfW^2)
  # Delta2: dependent case (full Eq. VlnM2d)
  var1 <- 1/dfW - r^2/dfW + 1/(dfW^2) +
    (r^4 * (s1^8 + s2^8)) /
    (2 * dfW^2 * s1^4 * s2^4)
  
  # 3. SAFE bootstrap (4‐variate)
  requireNamespace("MASS", quietly = TRUE)
  mu    <- c(x1bar, x2bar, s1^2, s2^2)
  Sigma <- matrix(0, 4, 4)
  Sigma[1,1] <- s1^2/n1
  Sigma[2,2] <- s2^2/n2
  Sigma[1,2] <- Sigma[2,1] <- r * s1 * s2 / sqrt(n1 * n2)
  Sigma[3,3] <- 2 * s1^4 / (n1 - 1)
  Sigma[4,4] <- 2 * s2^4 / (n2 - 1)
  Sigma[3,4] <- Sigma[4,3] <- 
    2 * r^2 * s1^2 * s2^2 / sqrt((n1 - 1)*(n2 - 1))
  
  sims <- MASS::mvrnorm(B, mu, Sigma)
  x1s  <- sims[,1]; x2s  <- sims[,2]
  v1s  <- sims[,3]; v2s  <- sims[,4]
  
  ok1 <- which(v1s > 0 & v2s > 0)
  if (length(ok1) < min_surv) {
    stop("Too few valid draws after var>0 filter: ", length(ok1))
  }
  MS_B_s <- (n1*n2)/(n1+n2) * (x1s[ok1] - x2s[ok1])^2
  MS_W_s <- ((n1-1)*v1s[ok1] + (n2-1)*v2s[ok1]) / dfW
  sB2_s  <- (MS_B_s - MS_W_s) / n0
  
  ok2 <- which(MS_W_s > 0 & sB2_s > 0)
  if (length(ok2) < min_surv) {
    stop("Too few valid draws after MS>0 filter: ", length(ok2))
  }
  lnM_s <- log(sqrt(sB2_s[ok2]) / sqrt(MS_W_s[ok2]))
  
  bias     <- mean(lnM_s) - lnM
  lnM_SAFE <- lnM - bias
  var_SAFE <- var(lnM_s)
  
  # 4. assemble results
  df <- data.frame(
    method = c("Delta1","Delta2","SAFE"),
    lnM    = c(lnM, lnM, lnM_SAFE),
    SE     = sqrt(c(var1, var2, var_SAFE)),
    stringsAsFactors = FALSE
  )
  df
}

# Worked example (paired, r = 0.5)
set.seed(42)
out <- safe_lnM(
  x1bar = 15, x2bar = 10,
  s1    =  2,  s2    = 2.2,
  n1    = 40,  n2    = 40,
  r     = 0,
  B     = 300000
)

# Round only numeric columns for display
out[ , c("lnM","SE")] <- lapply(out[ , c("lnM","SE")], round, 4)
print(out)
