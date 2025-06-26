########################################################################
## lnM  – simulation, summary & graphics (Δ-method vs SAFE-BC)
## now includes Monte-Carlo standard errors (MCSE)
## 24 independent + 12 paired designs × 18 θ  = 648 rows
## last tested: 25-Jun-2025  (R ≥ 4.3, MASS 7.3-60, ggplot2 3.5-0)
########################################################################

library(MASS)        # mvrnorm()
library(ggplot2)     # plotting
library(dplyr)       # used later for facet ordering
theme_set(theme_bw(11))

## ---------- 0. globals & helpers -------------------------------------
maxVar   <- 20                        # clamp huge Δ‐variances
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)

## lnM kernel
.lnM <- function(Δ, MSW, n0) 0.5 * (log(Δ) - log(n0) - log(MSW))

## ---------- 0a. true lnM ---------------------------------------------
lnM_true_ind <- function(theta, n1, n2, sigma = 1) {
  msb <- (n1*n2)/(n1+n2)*theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  .lnM(msb-msw, msw, 2*n1*n2/(n1+n2))
}
lnM_true_dep <- function(theta, n, sigma = 1) {
  msb <- (n/2)*theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  .lnM(msb-msw, msw, n)
}

## ---------- 1. Δ-method plug-in --------------------------------------
#  (independent and paired versions – unchanged)
lnM_delta_ind <- function(x1, x2, s1, s2, n1, n2) {
  h   <- n1 * n2 / (n1 + n2)
  MSB <- h * (x1 - x2)^2
  MSW <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2)
  Δ   <- safe_gap(MSB - MSW)
  if (is.na(Δ)) return(c(pt = NA, var = NA))
  pt  <- .lnM(Δ, MSW, 2*h)
  sD2 <- s1^2/n1 + s2^2/n2
  dif <- x1 - x2
  vB  <- h^2 * (2*sD2^2 + 4*sD2*dif^2)
  vW  <- 2*MSW^2 / (n1 + n2 - 2)
  g1  <- 0.5/Δ
  g2  <- -0.5*MSB/(Δ*MSW)
  var <- posify(g1^2*vB + g2^2*vW)
  c(pt = pt, var = pmin(var, maxVar))
}

lnM_delta_dep <- function(x1, x2, s1, s2, n, rho) {
  MSB <- (n/2)*(x1 - x2)^2
  MSW <- (s1^2 + s2^2)/2
  Δ   <- safe_gap(MSB - MSW)
  if (is.na(Δ)) return(c(pt = NA, var = NA))
  pt  <- .lnM(Δ, MSW, n)
  sD2 <- s1^2 + s2^2 - 2*rho*s1*s2
  dif <- x1 - x2
  vB  <- (n/2)^2 * (2*sD2^2/n^2 + 4*dif^2*sD2/n)
  vW  <- (s1^4 + s2^4 + 2*rho^2*s1^2*s2^2) / (2*(n-1))
  g1  <- 0.5/Δ
  g2  <- -0.5*MSB/(Δ*MSW)
  var <- posify(g1^2*vB + g2^2*vW)
  c(pt = pt, var = pmin(var, maxVar))
}

## ---------- 2. SAFE-bootstrap (indep & paired) ------------------------
# ...  (safe_ind and safe_dep exactly as before – omitted for brevity) ...
## ---------- 2. SAFE-bootstrap (independent & paired) -----------------
# … unchanged (safe_ind, safe_dep) …

## ---------- 2. SAFE-BC bootstrap (pivot CI, double truncation) --------
safe_ind <- function(x1bar, x2bar, s1, s2, n1, n2,
                     B = 1e4, chunk = 5e3) {
  h    <- n1 * n2 / (n1 + n2)
  n0   <- 2 * h
  MSB0 <- h * (x1bar - x2bar)^2
  MSW0 <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2)
  if (MSB0 <= MSW0)
    return(list(pt=NA, var=NA, lo=NA, hi=NA, kept=0L, total=B))
  z_raw <- .lnM(MSB0-MSW0, MSW0, n0)
  mu  <- c(x1bar,x2bar,s1^2,s2^2)
  Sig <- diag(c(s1^2/n1, s2^2/n2,
                2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  cloud <- numeric(B); kept <- 0L; k <- 0L
  while(k < B) {
    d  <- mvrnorm(chunk, mu, Sig)
    ok <- d[,3]>0 & d[,4]>0
    if(!any(ok)) next
    m1  <- d[ok,1]; m2 <- d[ok,2]
    v1  <- d[ok,3]; v2 <- d[ok,4]
    MSB <- h*(m1-m2)^2
    MSW <- ((n1-1)*v1+(n2-1)*v2)/(n1+n2-2)
    use <- MSB>MSW
    if(!any(use)) next
    kept <- kept + sum(use)
    vals <- .lnM(MSB[use]-MSW[use], MSW[use], n0)
    take_len <- min(length(vals), B-k)
    cloud[(k+1):(k+take_len)] <- vals[1:take_len]
    k <- k + take_len
  }
  cloud <- cloud[1:B]
  m_boot <- mean(cloud); pt <- 2*z_raw - m_boot; v_est <- var(cloud)
  cen <- cloud - m_boot
  qs  <- quantile(cen, c(0.025,0.975))
  lo <- pt + qs[1]; hi <- pt + qs[2]
  list(pt=pt, var=v_est, lo=lo, hi=hi, kept=kept, total=B)
}

safe_dep <- function(x1bar, x2bar, s1, s2, n, rho,
                     B = 1e4, chunk = 5e3) {
  MSB0 <- (n/2)*(x1bar-x2bar)^2
  MSW0 <- (s1^2+s2^2)/2
  if (MSB0 <= MSW0)
    return(list(pt=NA, var=NA, lo=NA, hi=NA, kept=0L, total=B))
  z_raw <- .lnM(MSB0-MSW0, MSW0, n)
  mu  <- c(x1bar,x2bar,s1^2,s2^2)
  Sig <- matrix(0,4,4)
  Sig[1,1]<-s1^2/n; Sig[2,2]<-s2^2/n
  Sig[1,2]<-Sig[2,1]<-rho*s1*s2/n
  Sig[3,3]<-2*s1^4/(n-1); Sig[4,4]<-2*s2^4/(n-1)
  Sig[3,4]<-Sig[4,3]<-2*rho^2*s1^2*s2^2/(n-1)
  cloud<-numeric(B); kept<-0L; k<-0L
  while(k < B) {
    d<-mvrnorm(chunk, mu, Sig)
    ok<-d[,3]>0 & d[,4]>0
    if(!any(ok)) next
    m1<-d[ok,1]; m2<-d[ok,2]
    v1<-d[ok,3]; v2<-d[ok,4]
    MSB<-(n/2)*(m1-m2)^2
    MSW<-(v1+v2)/2
    use<-MSB>MSW
    if(!any(use)) next
    kept<-kept+sum(use)
    vals<-.lnM(MSB[use]-MSW[use], MSW[use], n)
    take_len <- min(length(vals), B-k)
    cloud[(k+1):(k+take_len)] <- vals[1:take_len]
    k<-k+take_len
  }
  cloud<-cloud[1:B]
  m_boot<-mean(cloud); pt<-2*z_raw-m_boot; v_est<-var(cloud)
  cen <- cloud - m_boot
  qs  <- quantile(cen, c(0.025,0.975))
  lo <- pt + qs[1]; hi <- pt + qs[2]
  list(pt=pt, var=v_est, lo=lo, hi=hi, kept=kept, total=B)
}

## ---------- 3. one replicate -----------------------------------------
one_rep <- function(mu1, mu2, sd1, sd2,
                    n1, n2 = NULL, rho = 0, B = 1e4) {
  
  if (is.null(n2)) {                    # paired / dependent
    Sigma <- matrix(c(sd1^2, rho*sd1*sd2,
                      rho*sd1*sd2, sd2^2), 2)
    xy  <- mvrnorm(n1, c(mu1,mu2), Sigma)
    x1  <- xy[,1]; x2 <- xy[,2]; rho_hat <- cor(x1,x2)
  } else {                              # independent
    x1 <- rnorm(n1, mu1, sd1)
    x2 <- rnorm(n2, mu2, sd2); rho_hat <- 0
  }
  x1bar <- mean(x1); s1 <- sd(x1)
  x2bar <- mean(x2); s2 <- sd(x2)
  
  d <- if (is.null(n2))
    lnM_delta_dep(x1bar,x2bar,s1,s2,n1,rho_hat)
  else
    lnM_delta_ind(x1bar,x2bar,s1,s2,n1,n2)
  
  b <- if (is.null(n2))
    safe_dep (x1bar,x2bar,s1,s2,n1,rho_hat,B)
  else
    safe_ind (x1bar,x2bar,s1,s2,n1,n2,B)
  
  c(delta_pt  = d["pt"],  delta_var = d["var"],
    safe_pt   = b$pt,     safe_var  = b$var,
    safe_lo   = b$lo,     safe_hi   = b$hi,
    safe_kept = b$kept,   safe_total= b$total)
}

## ---------- 4. parameter grid ----------------------------------------
theta_vals <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
                1,1.2,1.4,1.6,1.8,2,2.5,3,4,5)

pairs_ind <- data.frame(n1=c(5,10,20,100,3,6,12,60),
                        n2=c(5,10,20,100,7,14,28,140))

grid_ind <- expand.grid(theta=theta_vals, idx=seq_len(nrow(pairs_ind)))
grid_ind$n1 <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2 <- pairs_ind$n2[grid_ind$idx]
grid_ind$design <- "indep"; grid_ind$idx <- NULL

grid_dep <- expand.grid(theta=theta_vals, n=c(5,10,20,100))
grid_dep$n1 <- grid_dep$n2 <- grid_dep$n; grid_dep$n <- NULL
grid_dep$design <- "paired"

param_grid <- rbind(grid_ind, grid_dep)

## ---------- 5. simulation driver (MCSE) ------------------------------
set.seed(20250625)
K_repl <- 1000                       # per parameter point
B_boot <- 1e4                        # bootstrap size

mcse_mean <- function(x) sqrt(var(x) / length(x))
mcse_prop <- function(x) sqrt(mean(x)*(1-mean(x)) / length(x))

runner <- function(i) {
  p       <- param_grid[i, ]
  sd0     <- 1
  true_ln <- if (p$design=="indep")
    lnM_true_ind(p$theta,p$n1,p$n2,sd0)
  else
    lnM_true_dep(p$theta,p$n1,sd0)
  
  ## replicate() now returns a matrix with row-names intact
  M <- replicate(
    K_repl,
    if (p$design=="indep")
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,n2=p$n2,B=B_boot)
    else
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,rho=0.8,B=B_boot),
    simplify = "matrix")
  
  ok  <- !is.na(M["delta_pt", ])
  Mok <- M[, ok, drop = FALSE]
  
  tv_s <- var(Mok["safe_pt", ])        # MC “truth”
  
  cover_d <- abs(Mok["delta_pt",]-true_ln) <= 1.96*sqrt(Mok["delta_var",])
  cover_s <- abs(Mok["safe_pt", ]-true_ln) <= 1.96*sqrt(Mok["safe_var", ])
  
  ## MCSEs
  mcse_delta_bias   <- mcse_mean(Mok["delta_pt", ] - true_ln)
  mcse_safe_bias    <- mcse_mean(Mok["safe_pt",  ] - true_ln)
  mcse_delta_varbar <- mcse_mean(Mok["delta_var", ])
  mcse_safe_varbar  <- mcse_mean(Mok["safe_var",  ])
  mcse_delta_cover  <- mcse_prop(cover_d)
  mcse_safe_cover   <- mcse_prop(cover_s)
  
  data.frame(
    theta          = p$theta,
    design         = p$design,
    n1             = p$n1,
    n2             = ifelse(p$design=="indep", p$n2, p$n1),
    true_lnM       = true_ln,
    
    # summaries
    delta_mean     = mean(Mok["delta_pt", ]),
    safe_mean      = mean(Mok["safe_pt",  ]),
    delta_bias     = mean(Mok["delta_pt", ]) - true_ln,
    safe_bias      = mean(Mok["safe_pt",  ]) - true_ln,
    mean_var_delta = mean(Mok["delta_var", ]),
    mean_var_safe  = mean(Mok["safe_var",  ]),
    relbias_delta  = 100*(mean(Mok["delta_var", ]) / tv_s - 1),
    relbias_safe   = 100*(mean(Mok["safe_var",  ]) / tv_s - 1),
    rmse_delta     = sqrt(mean((Mok["delta_pt", ] - true_ln)^2)),
    rmse_safe      = sqrt(mean((Mok["safe_pt",  ] - true_ln)^2)),
    cover_delta    = mean(cover_d),
    cover_safe     = mean(cover_s),
    boot_keep      = sum(Mok["safe_kept", ]),
    boot_total     = sum(Mok["safe_total", ]),
    
    # MCSEs
    mcse_bias_delta   = mcse_delta_bias,
    mcse_bias_safe    = mcse_safe_bias,
    mcse_varbar_delta = mcse_delta_varbar,
    mcse_varbar_safe  = mcse_safe_varbar,
    mcse_cover_delta  = mcse_delta_cover,
    mcse_cover_safe   = mcse_safe_cover
  )
}

results <- do.call(rbind, lapply(seq_len(nrow(param_grid)), runner))



## ---------- 6. plots (unchanged; you can use `results` directly) ------
## ... plotting section stays exactly as you had it ...