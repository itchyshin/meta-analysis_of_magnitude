########################################################################
## lnM – Δ‐method vs SAFE‐BC (bootstrap, bias‐corrected pivot CI)
## 24 independent + 12 paired sample‐size sets × 14 θ = 168 rows
## tested 25‐Jun‐2025 (R ≥ 4.3, MASS 7.3‐60, parallel, ggplot2 3.5‐0)
########################################################################

library(MASS)       # for mvrnorm()
library(parallel)   # for mclapply / parLapply
library(ggplot2)    # for plotting
theme_set(theme_bw(11))

## ---------------------------------------------------------------------
## 0.  helpers
## ---------------------------------------------------------------------
posify   <- function(x, eps = 1e-12) pmax(x, eps)          # force > 0
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)       # drop non‐pos Δ
lnM      <- function(Δ, MSW, n0) 0.5*(log(Δ) - log(n0) - log(MSW))

## ---------------------------------------------------------------------
## 0a. “true” lnM
## ---------------------------------------------------------------------
lnM_true_ind <- function(theta, n1, n2, sigma = 1) {
  msb <- n1*n2/(n1+n2)*theta^2; msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  lnM(msb-msw, msw, 2*n1*n2/(n1+n2))
}
lnM_true_dep <- function(theta, n, sigma = 1) {
  msb <- (n/2)*theta^2; msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  lnM(msb-msw, msw, n)
}

## ---------------------------------------------------------------------
## 1. Δ‐method plug‐in
## ---------------------------------------------------------------------
lnM_delta_ind <- function(x1bar,x2bar,s1,s2,n1,n2){
  h   <- n1*n2/(n1+n2)
  MSB <- h*(x1bar-x2bar)^2
  MSW <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  Δ   <- safe_gap(MSB-MSW); if (is.na(Δ)) return(c(pt=NA,var=NA))
  n0  <- 2*h
  pt  <- lnM(Δ, MSW, n0)
  sD2 <- s1^2/n1 + s2^2/n2; dif<-x1bar-x2bar
  vB  <- h^2*(2*sD2^2 + 4*sD2*dif^2)
  vW  <- 2*MSW^2/(n1+n2-2)
  g1  <- 0.5/Δ; g2 <- -0.5*MSB/(Δ*MSW)
  c(pt=pt, var=posify(g1^2*vB + g2^2*vW))
}

lnM_delta_dep <- function(x1bar,x2bar,s1,s2,n,rho){
  MSB <- (n/2)*(x1bar-x2bar)^2
  MSW <- (s1^2+s2^2)/2
  Δ   <- safe_gap(MSB-MSW); if (is.na(Δ)) return(c(pt=NA,var=NA))
  pt  <- lnM(Δ, MSW, n)
  sD2 <- s1^2+s2^2-2*rho*s1*s2; dif<-x1bar-x2bar
  vB  <- (n/2)^2*(2*sD2^2/n^2 + 4*dif^2*sD2/n)
  vW  <- (s1^4+s2^4+2*rho^2*s1^2*s2^2)/(2*(n-1))
  g1  <- 0.5/Δ; g2 <- -0.5*MSB/(Δ*MSW)
  c(pt=pt, var=posify(g1^2*vB + g2^2*vW))
}

## ---------------------------------------------------------------------
## 2. SAFE‐BC with pivot CI (single truncation only)
## ---------------------------------------------------------------------
safeBC_ind <- function(x1bar,x2bar,s1,s2,n1,n2,
                       B=1e4, chunk=5e3){
  h    <- n1*n2/(n1+n2); n0<-2*h
  MSB0 <- h*(x1bar-x2bar)^2
  MSW0 <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  if(MSB0<=MSW0) return(c(pt=NA,var=NA,lo=NA,hi=NA))
  z_raw <- lnM(MSB0-MSW0, MSW0, n0)
  mu  <- c(x1bar,x2bar,s1^2,s2^2)
  Sig <- diag(c(s1^2/n1, s2^2/n2, 2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  z_star <- numeric(B); k<-0L
  while(k < B){
    d <- mvrnorm(chunk, mu, Sig)
    ok<- d[,3]>0 & d[,4]>0
    if(!any(ok)) next
    m1<-d[ok,1]; m2<-d[ok,2]; v1<-d[ok,3]; v2<-d[ok,4]
    MSB<-h*(m1-m2)^2
    MSW<-((n1-1)*v1+(n2-1)*v2)/(n1+n2-2)
    keep<-MSB>MSW
    if(!any(keep)) next
    vals<-lnM(MSB[keep]-MSW[keep], MSW[keep], n0)
    take<-seq_len(min(length(vals), B-k))
    z_star[(k+1):(k+length(take))]<-vals[take]
    k<-k+length(take)
  }
  z_star<-z_star[1:B]
  m_boot <- mean(z_star)
  pt     <- 2*z_raw - m_boot
  var    <- var(z_star)
  centered <- z_star - m_boot
  qs <- quantile(centered, c(0.025,0.975))
  lo <- pt + qs[1]; hi <- pt + qs[2]
  c(pt=pt, var=var, lo=lo, hi=hi)
}

safeBC_dep <- function(x1bar,x2bar,s1,s2,n,rho,
                       B=1e4, chunk=5e3){
  MSB0 <- (n/2)*(x1bar-x2bar)^2
  MSW0 <- (s1^2+s2^2)/2
  if(MSB0<=MSW0) return(c(pt=NA,var=NA,lo=NA,hi=NA))
  z_raw <- lnM(MSB0-MSW0, MSW0, n)
  mu  <- c(x1bar,x2bar,s1^2,s2^2)
  Sig <- matrix(0,4,4)
  Sig[1,1]<-s1^2/n; Sig[2,2]<-s2^2/n
  Sig[1,2]<-Sig[2,1]<-rho*s1*s2/n
  Sig[3,3]<-2*s1^4/(n-1); Sig[4,4]<-2*s2^4/(n-1)
  Sig[3,4]<-Sig[4,3]<-2*rho^2*s1^2*s2^2/(n-1)
  z_star<-numeric(B); k<-0L
  while(k < B){
    d <- mvrnorm(chunk, mu, Sig)
    ok<-d[,3]>0 & d[,4]>0
    if(!any(ok)) next
    m1<-d[ok,1]; m2<-d[ok,2]; v1<-d[ok,3]; v2<-d[ok,4]
    MSB<-(n/2)*(m1-m2)^2
    MSW<-(v1+v2)/2
    keep<-MSB>MSW
    if(!any(keep)) next
    vals<-lnM(MSB[keep]-MSW[keep], MSW[keep], n)
    take<-seq_len(min(length(vals),B-k))
    z_star[(k+1):(k+length(take))]<-vals[take]
    k<-k+length(take)
  }
  z_star<-z_star[1:B]
  m_boot <- mean(z_star)
  pt     <- 2*z_raw - m_boot
  var    <- var(z_star)
  centered <- z_star - m_boot
  qs <- quantile(centered, c(0.025,0.975))
  lo <- pt + qs[1]; hi <- pt + qs[2]
  c(pt=pt, var=var, lo=lo, hi=hi)
}

## ---------------------------------------------------------------------
## 3. one_rep() – generate data & return both estimators + CI
## ---------------------------------------------------------------------
one_rep <- function(mu1,mu2,sd1,sd2,n1,n2=NULL,rho=0,B=1e4){
  if(is.null(n2)){
    Σ <- matrix(c(sd1^2,rho*sd1*sd2,
                  rho*sd1*sd2,sd2^2),2,2)
    xy <- mvrnorm(n1,c(mu1,mu2),Σ)
    x1<-xy[,1]; x2<-xy[,2]; rho_hat<-cor(x1,x2)
  } else {
    x1<-rnorm(n1,mu1,sd1); x2<-rnorm(n2,mu2,sd2); rho_hat<-0
  }
  x1bar<-mean(x1); s1<-sd(x1)
  x2bar<-mean(x2); s2<-sd(x2)
  d <- if(is.null(n2))
    lnM_delta_dep(x1bar,x2bar,s1,s2,n1,rho_hat)
  else
    lnM_delta_ind(x1bar,x2bar,s1,s2,n1,n2)
  b <- if(is.null(n2))
    safeBC_dep(x1bar,x2bar,s1,s2,n1,rho_hat,B)
  else
    safeBC_ind(x1bar,x2bar,s1,s2,n1,n2,B)
  c(delta_pt=d["pt"], delta_var=d["var"],
    SAFE_pt  =b["pt"],  SAFE_var =b["var"],
    SAFE_lo  =b["lo"],  SAFE_hi  =b["hi"])
}

## ---------------------------------------------------------------------
## 4. parameter grid
## ---------------------------------------------------------------------
theta_vals <- c(0.2,0.4,0.6,0.8,1,1.5,2,3,4,5)
pairs_ind  <- data.frame(n1=c(5,10,20,100,3,6,12,60),
                         n2=c(5,10,20,100,7,14,28,140))

grid_ind <- expand.grid(theta=theta_vals,
                        idx=seq_len(nrow(pairs_ind)),
                        KEEP.OUT.ATTRS=FALSE)
grid_ind$n1    <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2    <- pairs_ind$n2[grid_ind$idx]
grid_ind$rho   <- 0
grid_ind$design<-"indep"
grid_ind$idx   <- NULL

grid_dep <- expand.grid(theta=theta_vals,
                        n=c(5,10,20,100),
                        KEEP.OUT.ATTRS=FALSE)
grid_dep$rho   <- 0.8
grid_dep$design<- "paired"
grid_dep$n1    <- grid_dep$n; grid_dep$n2<-grid_dep$n; grid_dep$n<-NULL
grid_dep       <- grid_dep[,names(grid_ind)]

param_grid <- rbind(grid_ind, grid_dep)

## ---------------------------------------------------------------------
## 5. simulation driver
## ---------------------------------------------------------------------
set.seed(20250625)
K_repl  <- 1000
B_boot  <- 1e5
n_cores <- max(1, detectCores()-1)
use_fork<- (.Platform$OS.type!="windows")

runner <- function(i){
  p <- param_grid[i,]; sd0<-1
  true_ln <- if(p$design=="indep")
    lnM_true_ind(p$theta,p$n1,p$n2,sd0)
  else
    lnM_true_dep(p$theta,p$n1,sd0)
  reps <- replicate(
    K_repl,
    if(p$design=="indep")
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,n2=p$n2,B=B_boot)
    else
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,rho=p$rho,B=B_boot),
    simplify=FALSE
  )
  M <- do.call(cbind, reps)
  rownames(M) <- c("delta_pt","delta_var",
                   "SAFE_pt","SAFE_var",
                   "SAFE_lo","SAFE_hi")
  ## summarise function
  summarise <- function(pref){
    pts <- M[paste0(pref,"_pt"), ]
    vs  <- M[paste0(pref,"_var"), ]
    lo  <- M[paste0(pref,"_lo"), ]
    hi  <- M[paste0(pref,"_hi"), ]
    ok  <- !is.na(pts)
    pts<-pts[ok]; vs<-vs[ok]; lo<-lo[ok]; hi<-hi[ok]
    c(bias  = mean(pts)-true_ln,
      rvar  = 100*(mean(vs)/var(pts)-1),
      cover = mean(lo<=true_ln & true_ln<=hi))
  }
  d <- summarise("delta")
  b <- summarise("SAFE")
  data.frame(theta         = p$theta,
             design        = p$design,
             n1            = p$n1,
             n2            = p$n2,
             delta_bias    = d["bias"],
             SAFE_bias     = b["bias"],
             rb_delta      = d["rvar"],
             rb_SAFE       = b["rvar"],
             cover_delta   = d["cover"],
             cover_SAFE    = b["cover"])
}

if(use_fork){
  results_list <- mclapply(seq_len(nrow(param_grid)),
                           runner, mc.cores=n_cores)
} else {
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, library(MASS))
  clusterExport(cl, varlist=c("param_grid","runner",
                              "lnM_true_ind","lnM_true_dep",
                              "lnM_delta_ind","lnM_delta_dep",
                              "safeBC_ind","safeBC_dep",
                              "one_rep"), envir=environment())
  results_list <- parLapply(cl, seq_len(nrow(param_grid)), runner)
  stopCluster(cl)
}
results <- do.call(rbind, results_list)

## ---------------------------------------------------------------------
## 6. graphics: bias, rel‐bias, coverage
## ---------------------------------------------------------------------
pal <- c(delta="firebrick", SAFE="steelblue")

# A) Bias
bias_df <- rbind(
  data.frame(results[,1:4], estimator="delta", bias=results$delta_bias),
  data.frame(results[,1:4], estimator="SAFE",  bias=results$SAFE_bias )
)
p_bias <- ggplot(bias_df,
                 aes(theta, bias, colour=estimator, group=estimator)) +
  geom_hline(yintercept=0,lty=2,col="grey50") +
  geom_line(size=0.8) +
  facet_wrap(~ paste0(design," n1=",n1," n2=",n2), ncol=4) +
  scale_colour_manual(values=pal) +
  labs(x=expression(theta),
       y="Bias (estimate−true lnM)") 
print(p_bias)

# B) Relative bias of variance
rb_df <- rbind(
  data.frame(results[,1:4], estimator="delta", relbias=results$rb_delta),
  data.frame(results[,1:4], estimator="SAFE",  relbias=results$rb_SAFE )
)
p_var <- ggplot(rb_df,
                aes(theta, relbias, colour=estimator, group=estimator)) +
  geom_hline(yintercept=0,lty=2,col="grey50") +
  geom_line(size=0.8) +
  facet_wrap(~ paste0(design," n1=",n1," n2=",n2), ncol=4) +
  scale_colour_manual(values=pal) +
  labs(x=expression(theta),
       y="Relative bias of Var (%)")
print(p_var)

# C) Coverage
cov_df <- rbind(
  data.frame(results[,1:4], estimator="delta", cover=results$cover_delta),
  data.frame(results[,1:4], estimator="SAFE",  cover=results$cover_SAFE )
)
p_cov <- ggplot(cov_df,
                aes(theta, cover, colour=estimator, group=estimator)) +
  geom_hline(yintercept=0.95,lty=2,col="grey50") +
  geom_line(size=0.8) +
  facet_wrap(~ paste0(design," n1=",n1," n2=",n2), ncol=4) +
  scale_colour_manual(values=pal) +
  labs(x=expression(theta),
       y="Empirical coverage")
print(p_cov)