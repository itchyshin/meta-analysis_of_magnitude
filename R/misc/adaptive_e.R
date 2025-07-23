########################################################################
##  quick-and-dirty ε calibration for SAFE-BC
##  – pick ε that minimises absolute bias at θ = 1
########################################################################

library(MASS)

## --- a helper that returns the SAFE-BC bias for *one* ε ---------------
bias_one_eps <- function(eps, n1, n2, paired = FALSE,   # design
                         Bboot = 5e3, K = 300,          # speed!
                         theta = 1, sigma = 1, rho = .8,
                         trim = .05)                    # same trim as main code
{
  lnMfun <- function(Δ, MSW, n0) 0.5*(log(Δ) - log(n0) - log(MSW))
  pos  <- function(x) pmax(x, eps)                      # ε *here*
  
  ## reference truth ----------------------------------------------------
  lnM_true <- if (!paired) {
    msb <- (n1*n2)/(n1+n2)*theta^2
    msw <- sigma^2
    lnMfun(msb-msw, msw, 2*n1*n2/(n1+n2))
  } else {
    msb <- (n1/2)*theta^2
    msw <- sigma^2
    lnMfun(msb-msw, msw, n1)
  }
  
  ## one replicate ------------------------------------------------------
  one_rep <- function(){
    if (!paired){
      x1 <- rnorm(n1, 0, sigma)
      x2 <- rnorm(n2, theta, sigma)
      rho_hat <- 0
    } else {
      S <- matrix(c(sigma^2, rho*sigma^2,
                    rho*sigma^2, sigma^2), 2, 2)
      xy <- mvrnorm(n1, c(0, theta), S)
      x1 <- xy[,1]; x2 <- xy[,2]
      rho_hat <- cor(x1,x2)
    }
    x1bar <- mean(x1); s1 <- sd(x1)
    x2bar <- mean(x2); s2 <- sd(x2)
    
    ## sample ln M (with ε)
    if (!paired){
      h   <- n1*n2/(n1+n2)
      MSB <- h*(x1bar-x2bar)^2
      MSW <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
      n0  <- 2*n1*n2/(n1+n2)
    } else {
      MSB <- (n1/2)*(x1bar-x2bar)^2
      MSW <- (s1^2+s2^2)/2
      n0  <- n1
    }
    z_raw <- lnMfun(pos(MSB-MSW), MSW, n0)
    
    ## bootstrap cloud --------------------------------------------------
    mu <- c(x1bar,x2bar,s1^2,s2^2)
    if (!paired){
      Sig <- diag(c(s1^2/n1, s2^2/n2,
                    2*s1^4/(n1-1), 2*s2^4/(n2-1)))
    } else {
      Sig <- matrix(0,4,4)
      Sig[1,1] <- s1^2/n1; Sig[2,2] <- s2^2/n1
      Sig[1,2] <- Sig[2,1] <- rho_hat*s1*s2/n1
      Sig[3,3] <- 2*s1^4/(n1-1); Sig[4,4] <- 2*s2^4/(n1-1)
      Sig[3,4] <- Sig[4,3] <- 2*rho_hat^2*s1^2*s2^2/(n1-1)
    }
    
    z_star <- numeric(Bboot)
    k <- 0
    while(k < Bboot){
      d <- mvrnorm(Bboot, mu, Sig)
      good <- d[,3] > 0 & d[,4] > 0
      d <- d[good,,drop = FALSE]
      need <- Bboot - k
      if (nrow(d) < need) next
      d <- d[seq_len(need),]
      
      if (!paired){
        MSB <- (n1*n2)/(n1+n2)*(d[,1]-d[,2])^2
        MSW <- ((n1-1)*d[,3] + (n2-1)*d[,4])/(n1+n2-2)
      } else {
        MSB <- (n1/2)*(d[,1]-d[,2])^2
        MSW <- (d[,3]+d[,4])/2
      }
      Δ <- pos(MSB-MSW)
      z_star[(k+1):(k+need)] <- lnMfun(Δ, MSW, n0)
      k <- k + need
    }
    
    ## bias-corrected ln M
    z_bc <- z_raw - (mean(z_star, trim = trim) - z_raw)
    z_bc
  }
  
  mean(replicate(K, one_rep())) - lnM_true
}

## ---- grid search for ONE design (n1 = n2 = 10 here) ------------------
grid_eps <- 10^seq(-12, -2, by = 1)          # 1e-12 … 1e-2
biases   <- sapply(grid_eps, bias_one_eps,
                   n1 = 10, n2 = 10, paired = FALSE,
                   K = 300, Bboot = 3000)     # quick!

data.frame(log10eps = log10(grid_eps),
           bias     = biases)