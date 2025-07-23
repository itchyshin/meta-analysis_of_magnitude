########################################################################
##  lnM  – simulation, summary & graphics
##         Δ-method   vs   SAFE-1   vs   SAFE-2
##  24 independent + 12 paired designs  (36 total)
##  ASCII-only code – runs on Linux | macOS | Windows
##  tested: 27-Jun-2025   (R ≥ 4.3, MASS 7.3-60, ggplot2 3.5-0)
########################################################################

library(MASS)       # mvrnorm()
library(parallel)   # mclapply / parLapply
library(ggplot2)    # plots
theme_set(theme_bw(11))

## ───────────────────────────────────────────────────────────────
##  0.   HELPERS
## ───────────────────────────────────────────────────────────────
posify   <- function(x, eps = 1e-12) pmax(x, eps)          # force ≥ 0
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)        # NA if MSB≤MSW
lnM_from_F <- function(Fdraw, n0)
  0.5*(log(pmax(Fdraw - 1, 1e-12)) - log(n0))              # key identity

## ---------- “true” ln M (only defined for MSB>MSW) ----------
lnM_true_ind <- function(theta, n1, n2, sigma = 1) {
  msb <- (n1*n2)/(n1+n2)*theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  0.5*(log(msb-msw) - log(2*n1*n2/(n1+n2)) - log(msw))
}
lnM_true_dep <- function(theta, n, sigma = 1) {
  msb <- (n/2)*theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  0.5*(log(msb-msw) - log(n) - log(msw))
}

## ───────────────────────────────────────────────────────────────
##  1.   Δ-METHOD (plug-in) ESTIMATORS
## ───────────────────────────────────────────────────────────────
lnM_delta_ind <- function(x1bar, x2bar, s1, s2, n1, n2) {
  h   <- n1*n2/(n1+n2)
  MSB <- h*(x1bar-x2bar)^2
  MSW <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  n0  <- 2*n1*n2/(n1+n2)
  Δ   <- safe_gap(MSB-MSW);  if (is.na(Δ)) return(c(pt=NA, var=NA))
  
  ## point estimate –  equivalently  lnM_from_F(MSB/MSW, n0)
  zhat <- 0.5*(log(Δ) - log(n0) - log(MSW))
  
  ## delta-method variance
  sD2  <- s1^2/n1 + s2^2/n2
  difM <- x1bar - x2bar
  VarB <- h^2*(2*sD2^2 + 4*sD2*difM^2)
  VarW <- 2*MSW^2/(n1+n2-2)
  gB <- 0.5/Δ;   gW <- -0.5*MSB/(Δ*MSW)
  vhat <- posify(gB^2*VarB + gW^2*VarW)
  
  c(pt = zhat, var = vhat)
}

lnM_delta_dep <- function(x1bar, x2bar, s1, s2, n, rho) {
  MSB <- (n/2)*(x1bar-x2bar)^2
  MSW <- (s1^2+s2^2)/2
  Δ   <- safe_gap(MSB-MSW);  if (is.na(Δ)) return(c(pt=NA, var=NA))
  
  zhat <- 0.5*(log(Δ) - log(n) - log(MSW))
  
  sD2  <- s1^2 + s2^2 - 2*rho*s1*s2
  difM <- x1bar - x2bar
  VarB <- (n/2)^2*(2*sD2^2/n^2 + 4*difM^2*sD2/n)
  VarW <- (s1^4+s2^4+2*rho^2*s1^2*s2^2)/(2*(n-1))
  gB <- 0.5/Δ;   gW <- -0.5*MSB/(Δ*MSW)
  vhat <- posify(gB^2*VarB + gW^2*VarW)
  
  c(pt = zhat, var = vhat)
}

## ───────────────────────────────────────────────────────────────
##  2.   SAFE-1  (parametric bootstrap in **data** space)
## ───────────────────────────────────────────────────────────────
safe1_ind <- function(x1bar,x2bar,s1,s2,n1,n2,B=1e4,chunk=5e3){
  mu <- c(x1bar,x2bar,s1^2,s2^2)
  Σ  <- diag(c(s1^2/n1, s2^2/n2, 2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  cloud <- numeric(); kept <- 0L
  while (length(cloud) < B) {
    d <- mvrnorm(chunk, mu, Σ)
    okV <- d[,3]>0 & d[,4]>0
    d   <- d[okV,,drop=FALSE];  if(!nrow(d)) next
    m1<-d[,1]; m2<-d[,2]; v1<-d[,3]; v2<-d[,4]
    MSB <- (n1*n2)/(n1+n2)*(m1-m2)^2
    MSW <- ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2)
    ok  <- MSB > MSW
    kept <- kept + sum(ok);  if(!any(ok)) next
    Δ <- MSB[ok] - MSW[ok]
    cloud <- c(cloud,
               0.5*(log(Δ) - log(2*n1*n2/(n1+n2)) - log(MSW[ok])))
  }
  cloud <- cloud[seq_len(B)]
  list(pt = mean(cloud), var = var(cloud), kept = kept, total = B)
}

safe1_dep <- function(x1bar,x2bar,s1,s2,n,rho,B=1e4,chunk=5e3){
  mu <- c(x1bar,x2bar,s1^2,s2^2)
  Σ  <- matrix(0,4,4)
  Σ[1,1]<-s1^2/n; Σ[2,2]<-s2^2/n; Σ[1,2]<-Σ[2,1]<-rho*s1*s2/n
  Σ[3,3]<-2*s1^4/(n-1); Σ[4,4]<-2*s2^4/(n-1)
  Σ[3,4]<-Σ[4,3]<-2*rho^2*s1^2*s2^2/(n-1)
  cloud <- numeric(); kept <- 0L
  while (length(cloud) < B) {
    d <- mvrnorm(chunk, mu, Σ)
    okV <- d[,3]>0 & d[,4]>0
    d   <- d[okV,,drop=FALSE];  if(!nrow(d)) next
    m1<-d[,1]; m2<-d[,2]; v1<-d[,3]; v2<-d[,4]
    MSB <- (n/2)*(m1-m2)^2
    MSW <- (v1+v2)/2
    ok  <- MSB > MSW
    kept <- kept + sum(ok);  if(!any(ok)) next
    Δ <- MSB[ok] - MSW[ok]
    cloud <- c(cloud,
               0.5*(log(Δ) - log(n) - log(MSW[ok])))
  }
  cloud <- cloud[seq_len(B)]
  list(pt = mean(cloud), var = var(cloud), kept = kept, total = B)
}

## ───────────────────────────────────────────────────────────────
##  3.   SAFE-2  (bootstrap of the **F statistic** only)
##       ――――――――――――――――――――――――――――――――――――――――
##  • Point estimate: uses the *observed* F-ratio.
##  • Variance: draws B replicates from the *central* F(1,df2)
##    and plugs into lnM_from_F.
##  This stays valid because lnM depends only on (F-1).
## ───────────────────────────────────────────────────────────────
safe2_ind <- function(x1bar,x2bar,s1,s2,n1,n2,B=1e4){
  h   <- n1*n2/(n1+n2)
  MSB <- h*(x1bar-x2bar)^2
  MSW <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  if (MSB <= MSW) return(c(pt = NA, var = NA))
  
  Fobs <- MSB / MSW
  n0   <- 2*n1*n2/(n1+n2)
  df2  <- n1 + n2 - 2
  
  pt   <- lnM_from_F(Fobs, n0)
  var  <- var( lnM_from_F(rf(B, 1, df2), n0) )
  
  c(pt = pt, var = var)
}

safe2_dep <- function(x1bar,x2bar,s1,s2,n,B=1e4){
  MSB <- (n/2)*(x1bar-x2bar)^2
  MSW <- (s1^2+s2^2)/2
  if (MSB <= MSW) return(c(pt = NA, var = NA))
  
  Fobs <- MSB / MSW
  n0   <- n
  df2  <- n - 1                    # numerator df = 1
  
  pt   <- lnM_from_F(Fobs, n0)
  var  <- var( lnM_from_F(rf(B, 1, df2), n0) )
  
  c(pt = pt, var = var)
}

## ───────────────────────────────────────────────────────────────
##  4.   ONE REPLICATE
## ───────────────────────────────────────────────────────────────
one_rep <- function(mu1,mu2,sd1,sd2,n1,n2=NULL,rho=0,
                    B1=1e4,B2=1e4){
  if (is.null(n2)) {                       # paired
    Σ <- matrix(c(sd1^2, rho*sd1*sd2,
                  rho*sd1*sd2, sd2^2), 2, 2)
    xy <- mvrnorm(n1, c(mu1,mu2), Σ)
    x1 <- xy[,1];  x2 <- xy[,2];  rho_hat <- cor(x1,x2)
  } else {                                 # independent
    x1 <- rnorm(n1, mu1, sd1)
    x2 <- rnorm(n2, mu2, sd2)
    rho_hat <- 0
  }
  
  x1bar <- mean(x1);  x2bar <- mean(x2)
  s1 <- sd(x1);  s2 <- sd(x2)
  
  delta <- if (is.null(n2))
    lnM_delta_dep (x1bar,x2bar,s1,s2,n1,rho_hat)
  else lnM_delta_ind (x1bar,x2bar,s1,s2,n1,n2)
  
  safe1 <- if (is.null(n2))
    safe1_dep (x1bar,x2bar,s1,s2,n1,rho_hat,B1)
  else safe1_ind (x1bar,x2bar,s1,s2,n1,n2,      B1)
  
  safe2 <- if (is.null(n2))
    safe2_dep (x1bar,x2bar,s1,s2,n1,           B2)
  else safe2_ind (x1bar,x2bar,s1,s2,n1,n2,      B2)
  
  c(delta_pt = delta["pt"], delta_var = delta["var"],
    safe1_pt = safe1["pt"], safe1_var = safe1["var"],
    safe2_pt = safe2["pt"], safe2_var = safe2["var"])
}

## ───────────────────────────────────────────────────────────────
##  5.   PARAMETER GRID  (same 168 rows as before)
## ───────────────────────────────────────────────────────────────
theta_vals <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,4,5)
pairs_ind  <- data.frame(n1=c(5,10,20,100, 3,6,12,60),
                         n2=c(5,10,20,100, 7,14,28,140))

grid_ind <- expand.grid(theta=theta_vals, idx=seq_len(nrow(pairs_ind)),
                        KEEP.OUT.ATTRS=FALSE)
grid_ind$n1 <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2 <- pairs_ind$n2[grid_ind$idx]
grid_ind$rho <- 0; grid_ind$design <- "indep"; grid_ind$idx <- NULL

grid_dep <- expand.grid(theta=theta_vals, n=c(5,10,20,100),
                        KEEP.OUT.ATTRS=FALSE)
grid_dep$rho <- 0.8; grid_dep$design <- "paired"
grid_dep$n1 <- grid_dep$n; grid_dep$n2 <- grid_dep$n; grid_dep$n <- NULL
grid_dep <- grid_dep[,names(grid_ind)]

param_grid <- rbind(grid_ind, grid_dep)
stopifnot(nrow(param_grid) == 168)

## ───────────────────────────────────────────────────────────────
##  6.   SIMULATION DRIVER
## ───────────────────────────────────────────────────────────────
set.seed(20250627)
K_repl  <- 1000
B1 <- 1e4;   B2 <- 1e4
n_cores <- max(1, detectCores() - 1)
use_fork <- (.Platform$OS.type != "windows")

runner <- function(i){
  p   <- param_grid[i, ];  sd0 <- 1
  true_lnM <- if (p$design=="indep")
    lnM_true_ind(p$theta,p$n1,p$n2,sd0)
  else lnM_true_dep(p$theta,p$n1,sd0)
  
  reps <- replicate(
    K_repl,
    if (p$design=="indep")
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,n2=p$n2,
              B1=B1,B2=B2)
    else  one_rep(0,p$theta,sd0,sd0,n1=p$n1,rho=p$rho,
                  B1=B1,B2=B2),
    simplify = FALSE)
  M <- do.call(cbind, reps)
  rownames(M) <- c("d_pt","d_var","s1_pt","s1_var","s2_pt","s2_var")
  
  truVar_d  <- var(M["d_pt", ], na.rm=TRUE)
  truVar_s1 <- var(M["s1_pt",], na.rm=TRUE)
  truVar_s2 <- var(M["s2_pt",], na.rm=TRUE)
  
  ok <- !is.na(M["d_pt",])      # succeed ↔ MSB>MSW
  Mok <- M[, ok, drop=FALSE]
  
  data.frame(
    theta = p$theta, design=p$design, n1=p$n1, n2=p$n2, rho=p$rho,
    true_lnM = true_lnM,
    
    delta_bias = mean(Mok["d_pt", ]) - true_lnM,
    safe1_bias = mean(Mok["s1_pt",]) - true_lnM,
    safe2_bias = mean(Mok["s2_pt",]) - true_lnM,
    
    relbias_var_delta = 100*(mean(Mok["d_var", ])  / truVar_d  - 1),
    relbias_var_s1    = 100*(mean(Mok["s1_var", ]) / truVar_s1 - 1),
    relbias_var_s2    = 100*(mean(Mok["s2_var", ]) / truVar_s2 - 1),
    
    cover_delta = mean(abs(Mok["d_pt", ] - true_lnM) <=
                         1.96*sqrt(Mok["d_var", ])),
    cover_s1    = mean(abs(Mok["s1_pt",] - true_lnM) <=
                         1.96*sqrt(Mok["s1_var",])),
    cover_s2    = mean(abs(Mok["s2_pt",] - true_lnM) <=
                         1.96*sqrt(Mok["s2_var",]))
  )
}

if (use_fork) {
  results <- do.call(rbind,
                     mclapply(seq_len(nrow(param_grid)),
                              runner, mc.cores = n_cores))
} else {
  cl <- makeCluster(n_cores)
  results <- do.call(rbind,
                     parLapply(cl, seq_len(nrow(param_grid)), runner))
  stopCluster(cl)
}

## ───────────────────────────────────────────────────────────────
##  7.   GRAPHICS
## ───────────────────────────────────────────────────────────────
pal <- c(delta="firebrick", `SAFE-1`="steelblue", `SAFE-2`="darkgreen")

## A. Bias
bias_df <- rbind(
  data.frame(results[,1:4], estimator="delta",  bias=results$delta_bias),
  data.frame(results[,1:4], estimator="SAFE-1", bias=results$safe1_bias),
  data.frame(results[,1:4], estimator="SAFE-2", bias=results$safe2_bias))

ggplot(bias_df, aes(theta,bias,colour=estimator,group=estimator))+
  geom_hline(yintercept=0,lty=2,col="grey50")+
  geom_line(size=.8)+
  facet_wrap(~paste0(design,"  n1=",n1,"  n2=",n2),ncol=4)+
  scale_colour_manual(values=pal)+
  labs(x=expression(theta), y="Bias (estimate − true lnM)", colour=NULL)

## B. Variance relative bias
relvar_df <- rbind(
  data.frame(results[,1:4], estimator="delta",
             rb=results$relbias_var_delta),
  data.frame(results[,1:4], estimator="SAFE-1",
             rb=results$relbias_var_s1),
  data.frame(results[,1:4], estimator="SAFE-2",
             rb=results$relbias_var_s2))

ggplot(relvar_df, aes(theta,rb,colour=estimator,group=estimator))+
  geom_hline(yintercept=0,lty=2,col="grey50")+
  geom_line(size=.8)+
  facet_wrap(~paste0(design,"  n1=",n1,"  n2=",n2),ncol=4)+
  scale_colour_manual(values=pal)+
  labs(x=expression(theta), y="Relative bias of Var (%)", colour=NULL)

## C. Coverage
cov_df <- rbind(
  data.frame(results[,1:4], estimator="delta",  cover=results$cover_delta),
  data.frame(results[,1:4], estimator="SAFE-1", cover=results$cover_s1),
  data.frame(results[,1:4], estimator="SAFE-2", cover=results$cover_s2))

ggplot(cov_df, aes(theta,cover,colour=estimator,group=estimator))+
  geom_hline(yintercept=.95,lty=2,col="grey50")+
  geom_line(size=.8)+
  facet_wrap(~paste0(design,"  n1=",n1,"  n2=",n2),ncol=4)+
  scale_colour_manual(values=pal)+
  labs(x=expression(theta), y="Empirical coverage", colour=NULL)