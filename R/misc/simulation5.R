########################################################################
##  lnM  –  Δ-method      vs      SAFE-T      vs      SAFE-BC (adaptive ε)
##  8 independent + 4 paired sample-size designs  ×  14 θ’s  =  168 rows
##  tested: R ≥ 4.3   (MASS 7.3-60, parallel, ggplot2 3.5-0)
########################################################################

library(MASS)       # mvrnorm()
library(parallel)   # mclapply / parLapply
library(ggplot2)    # plots
theme_set(theme_bw(11))

## ---------------------------------------------------------------------
##  Data-adaptive ε  (single global κ*, local Δlo)  <<<<<<<<<<<<<<<<<<<
## ---------------------------------------------------------------------
KAPPA <- 0.15        # <-  one-time calibration constant (see note)

## helper that *always* keeps Δ positive with the local ε
mk_pos <- function(eps) function(x) pmax(x, eps)

lnM <- function(Δ, MSW, n0) 0.5*(log(Δ) - log(n0) - log(MSW))

## ---------------------------------------------------------------------
##  “True” ln M  for independent and paired designs
## ---------------------------------------------------------------------
lnM_true_ind <- function(theta, n1, n2, sigma = 1){
  msb <- (n1*n2)/(n1+n2) * theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  lnM(msb-msw, msw, 2*n1*n2/(n1+n2))
}
lnM_true_dep <- function(theta, n, sigma = 1){
  msb <- (n/2)*theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  lnM(msb-msw, msw, n)
}

## ---------------------------------------------------------------------
##  Δ-method plug-in variance – independent / paired   (unchanged)
## ---------------------------------------------------------------------
lnM_delta_ind <- function(x1bar,x2bar,s1,s2,n1,n2){
  h   <- n1*n2/(n1+n2)
  MSB <- h*(x1bar-x2bar)^2
  MSW <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  Δ   <- MSB-MSW;  if (Δ <= 0) return(c(pt=NA,var=NA))
  n0  <- 2*n1*n2/(n1+n2);  z <- lnM(Δ, MSW, n0)
  
  sD2 <- s1^2/n1 + s2^2/n2
  dif <- x1bar-x2bar
  varB <- h^2*(2*sD2^2 + 4*sD2*dif^2)
  varW <- 2*MSW^2/(n1+n2-2)
  gB <- 0.5/Δ;  gW <- -0.5*MSB/(Δ*MSW)
  v  <- pmax(0, gB^2*varB + gW^2*varW)
  c(pt=z, var=v)
}

lnM_delta_dep <- function(x1bar,x2bar,s1,s2,n,rho){
  MSB <- (n/2)*(x1bar-x2bar)^2
  MSW <- (s1^2+s2^2)/2
  Δ   <- MSB-MSW;  if (Δ <= 0) return(c(pt=NA,var=NA))
  z   <- lnM(Δ, MSW, n)
  
  sD2 <- s1^2 + s2^2 - 2*rho*s1*s2
  dif <- x1bar-x2bar
  varB <- (n/2)^2*(2*sD2^2/n^2 + 4*dif^2*sD2/n)
  varW <- (s1^4+s2^4+2*rho^2*s1^2*s2^2)/(2*(n-1))
  gB <- 0.5/Δ;  gW <- -0.5*MSB/(Δ*MSW)
  v  <- pmax(0, gB^2*varB + gW^2*varW)
  c(pt=z, var=v)
}

## ---------------------------------------------------------------------
##  SAFE-T  (original single-truncation) – *unchanged* code omitted
##  …  keep your previous safeT_ind() / safeT_dep()
## ---------------------------------------------------------------------

## -----------------------------------------------------------------
##  SAFE-BC-lite     (paired / dependent groups)
## -----------------------------------------------------------------
safeBC_dep <- function(x1bar, x2bar, s1, s2, n, rho,
                       B = 1e4, chunk = 5e3,
                       trim = 0.05)
{
  ## ---------- sample (raw) ----------
  MSB0 <- (n/2)*(x1bar - x2bar)^2
  MSW0 <- (s1^2 + s2^2)/2
  if (MSB0 <= MSW0)
    return(c(pt = NA, var = NA))
  z_raw <- lnM(MSB0 - MSW0, MSW0, n)
  
  ## ---------- bootstrap cloud -------
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- matrix(0,4,4)
  Sig[1,1] <- s1^2/n;           Sig[2,2] <- s2^2/n
  Sig[1,2] <- Sig[2,1] <- rho*s1*s2/n
  Sig[3,3] <- 2*s1^4/(n-1);     Sig[4,4] <- 2*s2^4/(n-1)
  Sig[3,4] <- Sig[4,3] <- 2*rho^2*s1^2*s2^2/(n-1)
  
  z <- numeric()
  while (length(z) < B) {
    d   <- mvrnorm(chunk, mu, Sig)
    okV <- d[,3] > 0 & d[,4] > 0
    d   <- d[okV,,drop = FALSE];  if(!nrow(d)) next
    
    m1  <- d[,1];  m2 <- d[,2]
    v1  <- d[,3];  v2 <- d[,4]
    MSB <- (n/2)*(m1 - m2)^2
    MSW <- (v1 + v2)/2
    ok  <- MSB > MSW
    if(!any(ok)) next
    
    z <- c(z, lnM(MSB[ok] - MSW[ok], MSW[ok], n))
  }
  z <- z[seq_len(B)]
  
  z_bc <- z_raw - (mean(z, trim = trim) - z_raw)
  c(pt = z_bc, var = var(z))
}
## ---------------------------------------------------------------------
##  one_rep()  – generate data & compute three estimators  (UNCHANGED)
## ---------------------------------------------------------------------
one_rep <- function(mu1,mu2,sd1,sd2,n1,n2=NULL,rho=0,B=1e4){
  if(is.null(n2)){                       # paired
    Σ  <- matrix(c(sd1^2, rho*sd1*sd2,
                   rho*sd1*sd2, sd2^2), 2,2)
    xy <- mvrnorm(n1, c(mu1,mu2), Σ)
    x1<-xy[,1]; x2<-xy[,2]; rho_hat<-cor(x1,x2)
  } else {                               # independent
    x1<-rnorm(n1,mu1,sd1);  x2<-rnorm(n2,mu2,sd2); rho_hat<-0
  }
  
  x1bar<-mean(x1); s1<-sd(x1)
  x2bar<-mean(x2); s2<-sd(x2)
  
  delta <- if(is.null(n2))
    lnM_delta_dep(x1bar,x2bar,s1,s2,n1,rho_hat)
  else lnM_delta_ind(x1bar,x2bar,s1,s2,n1,n2)
  
  sT <- if(is.null(n2))
    safeT_dep (x1bar,x2bar,s1,s2,n1,rho_hat,B)
  else safeT_ind(x1bar,x2bar,s1,s2,n1,n2,B)
  
  sB <- if(is.null(n2))
    safeBC_dep(x1bar,x2bar,s1,s2,n1,rho_hat,B)
  else safeBC_ind(x1bar,x2bar,s1,s2,n1,n2,B)
  
  c(delta_pt=delta["pt"], delta_var=delta["var"],
    sT_pt=sT["pt"],       sT_var=sT["var"],
    sB_pt=sB["pt"],       sB_var=sB["var"])
}

## ---------------------------------------------------------------------
##  Parameter grid: 8 indep + 4 paired  × 14 θ’s  = 168 rows
## ---------------------------------------------------------------------
theta_vals <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,4,5)
pairs_ind  <- data.frame(n1=c(5,10,20,100, 3,6,12,60),
                         n2=c(5,10,20,100, 7,14,28,140))

grid_ind <- expand.grid(theta=theta_vals,
                        idx   = seq_len(nrow(pairs_ind)),
                        KEEP.OUT.ATTRS=FALSE)
grid_ind$n1     <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2     <- pairs_ind$n2[grid_ind$idx]
grid_ind$rho    <- 0
grid_ind$design <- "indep"
grid_ind$idx    <- NULL

grid_dep <- expand.grid(theta=theta_vals,
                        n=c(5,10,20,100),
                        KEEP.OUT.ATTRS=FALSE)
grid_dep$rho    <- 0.8
grid_dep$design <- "paired"
grid_dep$n1     <- grid_dep$n
grid_dep$n2     <- grid_dep$n
grid_dep$n      <- NULL
grid_dep        <- grid_dep[,names(grid_ind)]

param_grid <- rbind(grid_ind, grid_dep)

## ---------------------------------------------------------------------
##  Simulation driver  (unchanged apart from function names)
## ---------------------------------------------------------------------
set.seed(20250627)
K_repl  <- 1000
B_boot  <- 1e4
n_cores <- max(1, detectCores()-1)
use_fork<- (.Platform$OS.type!="windows")

runner <- function(i){
  p <- param_grid[i,]; sd0 <- 1
  true_ln <- if(p$design=="indep")
    lnM_true_ind(p$theta,p$n1,p$n2,sd0) else
      lnM_true_dep(p$theta,p$n1,sd0)
  
  reps <- replicate(
    K_repl,
    if(p$design=="indep")
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,n2=p$n2,B=B_boot) else
        one_rep(0,p$theta,sd0,sd0,n1=p$n1,rho=p$rho,B=B_boot),
    simplify=FALSE)
  M <- do.call(cbind,reps)
  rownames(M)<-c("d_pt","d_var","t_pt","t_var","b_pt","b_var")
  ok <- !is.na(M["d_pt",]); Mok<-M[,ok,drop=FALSE]
  
  summarise <- function(pref){
    pts <- Mok[paste0(pref,"_pt"),]
    vars<- Mok[paste0(pref,"_var"),]
    bias <- mean(pts)-true_ln
    rvar <- 100*(mean(vars)/var(pts)-1)
    cover<- mean(abs(pts-true_ln)<=1.96*sqrt(vars))
    c(bias,rvar,cover)
  }
  d<-summarise("d"); t<-summarise("t"); b<-summarise("b")
  data.frame(theta=p$theta,design=p$design,n1=p$n1,n2=p$n2,
             delta_bias=d[1],SAFE_T_bias=t[1],SAFE_BC_bias=b[1],
             rb_delta=d[2],rb_SAFE_T=t[2],rb_SAFE_BC=b[2],
             cover_delta=d[3],cover_SAFE_T=t[3],cover_SAFE_BC=b[3])
}

if(use_fork){
  results <- do.call(rbind,
                     mclapply(seq_len(nrow(param_grid)),
                              runner, mc.cores=n_cores))
} else {
  cl <- makeCluster(n_cores)
  results <- do.call(rbind,
                     parLapply(cl, seq_len(nrow(param_grid)), runner))
  stopCluster(cl)
}

## ---------------------------------------------------------------------
##  quick visual check – bias panel only
## ---------------------------------------------------------------------
pal <- c(delta="firebrick", `SAFE-T`="darkgreen", `SAFE-BC`="steelblue")
bias_df <- rbind(
  data.frame(results[,1:4], est="delta",   bias=results$delta_bias),
  data.frame(results[,1:4], est="SAFE-T",  bias=results$SAFE_T_bias),
  data.frame(results[,1:4], est="SAFE-BC", bias=results$SAFE_BC_bias)
)

ggplot(bias_df,
       aes(theta, bias, colour=est, group=est))+
  geom_hline(yintercept=0,lty=2,col="grey50")+
  geom_line(size=.8)+
  facet_wrap(~paste0(design,"  n1=",n1,"  n2=",n2), ncol=4)+
  scale_colour_manual(values=pal)+
  labs(x=expression(theta),
       y="Bias (estimate – true lnM)",
       colour=NULL)

## -----------------------------------------------------------------
##  SAFE-BC-lite     (independent groups)
##  – identical to SAFE-BC except: draws with Δ ≤ 0 are discarded
## -----------------------------------------------------------------
safeBC_ind <- function(x1bar, x2bar, s1, s2,
                       n1, n2, B = 1e4, chunk = 5e3,
                       trim = 0.05)
{
  h   <- n1*n2/(n1+n2)
  n0  <- 2*n1*n2/(n1+n2)
  
  ## ---------- sample (raw) ----------
  MSB0 <- h*(x1bar-x2bar)^2
  MSW0 <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  if (MSB0 <= MSW0)          # data give Δ ≤ 0  ->  no estimate
    return(c(pt = NA, var = NA))
  z_raw <- lnM(MSB0 - MSW0, MSW0, n0)
  
  ## ---------- bootstrap cloud -------
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- diag(c(s1^2/n1, s2^2/n2,
                2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  z <- numeric()
  
  while (length(z) < B) {
    d   <- mvrnorm(chunk, mu, Sig)
    okV <- d[,3] > 0 & d[,4] > 0
    d   <- d[okV,,drop = FALSE];   if(!nrow(d)) next
    
    m1  <- d[,1];  m2 <- d[,2]
    v1  <- d[,3];  v2 <- d[,4]
    MSB <- h*(m1 - m2)^2
    MSW <- ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2)
    ok  <- MSB > MSW                 # <-  **single truncation**
    if(!any(ok)) next
    
    z <- c(z, lnM(MSB[ok] - MSW[ok], MSW[ok], n0))
  }
  z <- z[seq_len(B)]
  
  z_bc <- z_raw - (mean(z, trim = trim) - z_raw)
  c(pt = z_bc, var = var(z))
}

