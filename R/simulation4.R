########################################################################
##  lnM  –  Δ-method   vs   SAFE-T   vs   SAFE-BC  (adaptive ε)
##  24 independent + 12 paired designs  =  36 total
##  tested: 28-Jun-2025   R ≥ 4.3   MASS 7.3-60   ggplot2 3.5-0
########################################################################

library(MASS)          # mvrnorm()
library(parallel)
library(ggplot2)
theme_set(theme_bw(11))

## ------------------------------------------------------------------ ##
##  helpers                                                           ##
## ------------------------------------------------------------------ ##
eps0  <- 1e-6                          # only for Δ-method var-guard
pos0  <- function(x, e = eps0) pmax(x, e)
adapt_eps <- function(n_tot) 10^(round(-3 - log10(n_tot)))
lnM <- function(Δ, MSW, n0) 0.5*(log(Δ) - log(n0) - log(MSW))

## robust wrapper around MASS::mvrnorm --------------------------------
safe_mvrnorm <- function(n, mu, Sigma) {
  repeat {
    out <- tryCatch(
      mvrnorm(n, mu, Sigma + diag(1e-12, nrow(Sigma))),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }
}

## true ln M -----------------------------------------------------------
lnM_true_ind <- function(theta,n1,n2,sigma=1){
  MSB <- (n1*n2)/(n1+n2)*theta^2
  MSW <- sigma^2
  if (MSB <= MSW) return(NA_real_)
  lnM(MSB-MSW, MSW, 2*n1*n2/(n1+n2))
}
lnM_true_dep <- function(theta,n,sigma=1){
  MSB <- (n/2)*theta^2
  MSW <- sigma^2
  if (MSB <= MSW) return(NA_real_)
  lnM(MSB-MSW, MSW, n)
}

## Δ-method (independent & paired) ------------------------------------
lnM_delta_ind <- function(x1bar,x2bar,s1,s2,n1,n2){
  h   <- n1*n2/(n1+n2)
  MSB <- h*(x1bar-x2bar)^2
  MSW <- ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  Δ   <- MSB-MSW
  if (Δ <= 0) return(c(pt=NA,var=NA))
  n0  <- 2*n1*n2/(n1+n2); z <- lnM(Δ, MSW, n0)
  
  sD2 <- s1^2/n1 + s2^2/n2
  dif <- x1bar-x2bar
  VarB<- h^2*(2*sD2^2 + 4*sD2*dif^2)
  VarW<- 2*MSW^2/(n1+n2-2)
  gB  <- 0.5/Δ; gW <- -0.5*MSB/(Δ*MSW)
  v   <- pos0(gB^2*VarB + gW^2*VarW)
  c(pt=z, var=v)
}

lnM_delta_dep <- function(x1bar,x2bar,s1,s2,n,rho){
  MSB <- (n/2)*(x1bar-x2bar)^2
  MSW <- (s1^2+s2^2)/2
  Δ   <- MSB-MSW
  if (Δ <= 0) return(c(pt=NA,var=NA))
  z <- lnM(Δ, MSW, n)
  
  sD2 <- s1^2 + s2^2 - 2*rho*s1*s2
  dif <- x1bar-x2bar
  VarB<- (n/2)^2*(2*sD2^2/n^2 + 4*dif^2*sD2/n)
  VarW<- (s1^4+s2^4+2*rho^2*s1^2*s2^2)/(2*(n-1))
  gB  <- 0.5/Δ; gW <- -0.5*MSB/(Δ*MSW)
  v   <- pos0(gB^2*VarB + gW^2*VarW)
  c(pt=z, var=v)
}

## SAFE-T --------------------------------------------------------------
safeT_ind <- function(x1bar,x2bar,s1,s2,n1,n2,B=1e4,chunk=5e3){
  mu  <- c(x1bar,x2bar,s1^2,s2^2)
  Sig <- diag(c(s1^2/n1, s2^2/n2,
                2*s1^4/(n1-1), 2*s2^4/(n2-1)))
  n0  <- 2*n1*n2/(n1+n2);  z <- numeric()
  while (length(z) < B){
    d <- safe_mvrnorm(chunk, mu, Sig)
    d <- d[d[,3]>0 & d[,4]>0,,drop=FALSE]; if(!nrow(d)) next
    MSB <- (n1*n2)/(n1+n2)*(d[,1]-d[,2])^2
    MSW <- ((n1-1)*d[,3] + (n2-1)*d[,4])/(n1+n2-2)
    ok  <- MSB>MSW; if(!any(ok)) next
    z   <- c(z, lnM(MSB[ok]-MSW[ok], MSW[ok], n0))
  }
  z <- z[1:B]; c(pt=mean(z), var=var(z))
}

safeT_dep <- function(x1bar,x2bar,s1,s2,n,rho,B=1e4,chunk=5e3){
  mu  <- c(x1bar,x2bar,s1^2,s2^2)
  Sig <- matrix(0,4,4)
  Sig[1,1]<-s1^2/n; Sig[2,2]<-s2^2/n
  Sig[1,2]<-Sig[2,1]<-rho*s1*s2/n
  Sig[3,3]<-2*s1^4/(n-1); Sig[4,4]<-2*s2^4/(n-1)
  Sig[3,4]<-Sig[4,3]<-2*rho^2*s1^2*s2^2/(n-1)
  z <- numeric()
  while (length(z) < B){
    d <- safe_mvrnorm(chunk, mu, Sig)
    d <- d[d[,3]>0 & d[,4]>0,,drop=FALSE]; if(!nrow(d)) next
    MSB <- (n/2)*(d[,1]-d[,2])^2
    MSW <- (d[,3]+d[,4])/2
    ok  <- MSB>MSW; if(!any(ok)) next
    z   <- c(z, lnM(MSB[ok]-MSW[ok], MSW[ok], n))
  }
  z <- z[1:B]; c(pt=mean(z), var=var(z))
}

## SAFE-BC -------------------------------------------------------------
bias_correct <- function(z_raw,z_star,trim=.05)
  z_raw - (mean(z_star,trim=trim) - z_raw)

safeBC_ind <- function(x1bar,x2bar,s1,s2,n1,n2,B=1e4,chunk=5e3,trim=.05){
  n_tot<-n1+n2; eps<-adapt_eps(n_tot); pos<-function(x)pmax(x,eps)
  h <- n1*n2/(n1+n2)
  MSB0<-h*(x1bar-x2bar)^2
  MSW0<-((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2)
  if (MSB0<=MSW0) return(c(pt=NA,var=NA))
  n0<-2*n1*n2/(n1+n2); z_raw<-lnM(MSB0-MSW0,MSW0,n0)
  
  mu<-c(x1bar,x2bar,s1^2,s2^2)
  Sig<-diag(c(s1^2/n1,s2^2/n2,2*s1^4/(n1-1),2*s2^4/(n2-1)))
  z<-numeric()
  while(length(z)<B){
    d<-safe_mvrnorm(chunk,mu,Sig)
    d<-d[d[,3]>0 & d[,4]>0,,drop=FALSE]; if(!nrow(d)) next
    MSB<-h*(d[,1]-d[,2])^2
    MSW<-((n1-1)*d[,3] + (n2-1)*d[,4])/(n1+n2-2)
    ok <- MSB>MSW; if(!any(ok)) next
    z<-c(z, lnM(pos(MSB[ok]-MSW[ok]), MSW[ok], n0))
  }
  z<-z[1:B]; c(pt=bias_correct(z_raw,z,trim), var=var(z))
}

safeBC_dep <- function(x1bar,x2bar,s1,s2,n,rho,B=1e4,chunk=5e3,trim=.05){
  eps<-adapt_eps(n); pos<-function(x)pmax(x,eps)
  MSB0<-(n/2)*(x1bar-x2bar)^2
  MSW0<-(s1^2+s2^2)/2
  if (MSB0<=MSW0) return(c(pt=NA,var=NA))
  z_raw<-lnM(MSB0-MSW0,MSW0,n)
  
  mu<-c(x1bar,x2bar,s1^2,s2^2)
  Sig<-matrix(0,4,4)
  Sig[1,1]<-s1^2/n; Sig[2,2]<-s2^2/n
  Sig[1,2]<-Sig[2,1]<-rho*s1*s2/n
  Sig[3,3]<-2*s1^4/(n-1); Sig[4,4]<-2*s2^4/(n-1)
  Sig[3,4]<-Sig[4,3]<-2*rho^2*s1^2*s2^2/(n-1)
  
  z<-numeric()
  while(length(z)<B){
    d<-safe_mvrnorm(chunk,mu,Sig)
    d<-d[d[,3]>0 & d[,4]>0,,drop=FALSE]; if(!nrow(d)) next
    MSB<-(n/2)*(d[,1]-d[,2])^2
    MSW<-(d[,3]+d[,4])/2
    ok <- MSB>MSW; if(!any(ok)) next
    z<-c(z, lnM(pos(MSB[ok]-MSW[ok]), MSW[ok], n))
  }
  z<-z[1:B]; c(pt=bias_correct(z_raw,z,trim), var=var(z))
}

## one replicate -------------------------------------------------------
one_rep <- function(mu1,mu2,sd1,sd2,n1,n2=NULL,rho=0,B=1e4){
  if(is.null(n2)){                 # paired
    Sig<-matrix(c(sd1^2,rho*sd1*sd2,rho*sd1*sd2,sd2^2),2,2)
    xy <- safe_mvrnorm(n1,c(mu1,mu2),Sig)
    x1<-xy[,1]; x2<-xy[,2]; rho_hat<-cor(x1,x2)
  } else {                         # independent
    x1<-rnorm(n1,mu1,sd1); x2<-rnorm(n2,mu2,sd2); rho_hat<-0
  }
  x1bar<-mean(x1); s1<-sd(x1)
  x2bar<-mean(x2); s2<-sd(x2)
  
  d <- if(is.null(n2))
    lnM_delta_dep(x1bar,x2bar,s1,s2,n1,rho_hat)
  else lnM_delta_ind(x1bar,x2bar,s1,s2,n1,n2)
  t <- if(is.null(n2))
    safeT_dep(x1bar,x2bar,s1,s2,n1,rho_hat,B)
  else safeT_ind(x1bar,x2bar,s1,s2,n1,n2,B)
  b <- if(is.null(n2))
    safeBC_dep(x1bar,x2bar,s1,s2,n1,rho_hat,B)
  else safeBC_ind(x1bar,x2bar,s1,s2,n1,n2,B)
  
  c(delta_pt=d["pt"],delta_var=d["var"],
    sT_pt=t["pt"],   sT_var=t["var"],
    sB_pt=b["pt"],   sB_var=b["var"])
}

## parameter grid ------------------------------------------------------
theta_vals <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,4,5)
pairs_ind  <- data.frame(n1=c(5,10,20,100, 3,6,12,60),
                         n2=c(5,10,20,100, 7,14,28,140))
grid_ind <- expand.grid(theta=theta_vals,
                        idx=seq_len(nrow(pairs_ind)),
                        KEEP.OUT.ATTRS=FALSE)
grid_ind$n1<-pairs_ind$n1[grid_ind$idx]
grid_ind$n2<-pairs_ind$n2[grid_ind$idx]
grid_ind$rho<-0; grid_ind$design<-"indep"; grid_ind$idx<-NULL
grid_dep <- expand.grid(theta=theta_vals,
                        n=c(5,10,20,100),
                        KEEP.OUT.ATTRS=FALSE)
grid_dep$rho<-0.8; grid_dep$design<-"paired"
grid_dep$n1<-grid_dep$n; grid_dep$n2<-grid_dep$n; grid_dep$n<-NULL
grid_dep<-grid_dep[,names(grid_ind)]
param_grid <- rbind(grid_ind,grid_dep)

## simulation driver ---------------------------------------------------
set.seed(20250628)
K_repl <- 1000; B_boot <- 1e4
n_cores <- max(1, detectCores()-1)
use_fork <- (.Platform$OS.type!="windows")

runner <- function(i){
  p<-param_grid[i,]; sd0<-1
  true_ln <- if(p$design=="indep")
    lnM_true_ind(p$theta,p$n1,p$n2,sd0)
  else lnM_true_dep(p$theta,p$n1,sd0)
  
  reps <- replicate(
    K_repl,
    if(p$design=="indep")
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,n2=p$n2,B=B_boot)
    else
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,rho=p$rho,B=B_boot),
    simplify=FALSE)
  
  M <- do.call(cbind,reps)
  if (length(M)==0)                           # nothing worked
    return(data.frame(theta=p$theta,design=p$design,n1=p$n1,n2=p$n2,
                      delta_bias=NA,SAFE_T_bias=NA,SAFE_BC_bias=NA,
                      rb_delta=NA,rb_SAFE_T=NA,rb_SAFE_BC=NA,
                      cover_delta=NA,cover_SAFE_T=NA,cover_SAFE_BC=NA))
  
  ok <- !is.na(M["delta_pt",]); Mok <- M[,ok,drop=FALSE]
  if(ncol(Mok)==0)                            # no non-NA replicate
    return(data.frame(theta=p$theta,design=p$design,n1=p$n1,n2=p$n2,
                      delta_bias=NA,SAFE_T_bias=NA,SAFE_BC_bias=NA,
                      rb_delta=NA,rb_SAFE_T=NA,rb_SAFE_BC=NA,
                      cover_delta=NA,cover_SAFE_T=NA,cover_SAFE_BC=NA))
  
  summarise <- function(pref){
    pts <- Mok[paste0(pref,"_pt"),]
    vars<- Mok[paste0(pref,"_var"),]
    c(bias = mean(pts) - true_ln,
      rvar = 100*(mean(vars)/var(pts)-1),
      cover= mean(abs(pts-true_ln) <= 1.96*sqrt(vars)))
  }
  d<-summarise("delta"); t<-summarise("sT"); b<-summarise("sB")
  
  data.frame(theta=p$theta,design=p$design,n1=p$n1,n2=p$n2,
             delta_bias=d["bias"],SAFE_T_bias=t["bias"],
             SAFE_BC_bias=b["bias"],rb_delta=d["rvar"],
             rb_SAFE_T=t["rvar"],rb_SAFE_BC=b["rvar"],
             cover_delta=d["cover"],cover_SAFE_T=t["cover"],
             cover_SAFE_BC=b["cover"])
}

if(use_fork){
  results <- do.call(rbind,
                     mclapply(seq_len(nrow(param_grid)),
                              runner, mc.cores=n_cores))
} else {
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, library(MASS))
  clusterExport(cl, ls(), envir=environment())
  results <- do.call(rbind,
                     parLapply(cl, seq_len(nrow(param_grid)), runner))
  stopCluster(cl)
}

## plots ---------------------------------------------------------------
pal <- c(delta="firebrick",`SAFE-T`="darkgreen",`SAFE-BC`="steelblue")

bias_df <- rbind(
  data.frame(results[,1:4], est="delta",  bias=results$delta_bias),
  data.frame(results[,1:4], est="SAFE-T", bias=results$SAFE_T_bias),
  data.frame(results[,1:4], est="SAFE-BC",bias=results$SAFE_BC_bias))

ggplot(bias_df,aes(theta,bias,colour=est,group=est))+
  geom_hline(yintercept=0,lty=2,col="grey50")+
  geom_line(size=.8)+
  facet_wrap(~paste0(design,"  n1=",n1,"  n2=",n2), ncol=4)+
  scale_colour_manual(values=pal)+
  labs(x=expression(theta),
       y="Bias (estimate − true lnM)", colour=NULL)

rb_df <- rbind(
  data.frame(results[,1:4], est="delta",  rb=results$rb_delta),
  data.frame(results[,1:4], est="SAFE-T", rb=results$rb_SAFE_T),
  data.frame(results[,1:4], est="SAFE-BC",rb=results$rb_SAFE_BC))

ggplot(rb_df,aes(theta,rb,colour=est,group=est))+
  geom_hline(yintercept=0,lty=2,col="grey50")+
  geom_line(size=.8)+
  facet_wrap(~paste0(design,"  n1=",n1,"  n2=",n2), ncol=4)+
  scale_colour_manual(values=pal)+
  labs(x=expression(theta),
       y="Relative bias of Var (%)", colour=NULL)

cov_df <- rbind(
  data.frame(results[,1:4], est="delta",  cover=results$cover_delta),
  data.frame(results[,1:4], est="SAFE-T", cover=results$cover_SAFE_T),
  data.frame(results[,1:4], est="SAFE-BC",cover=results$cover_SAFE_BC))

ggplot(cov_df,aes(theta,cover,colour=est,group=est))+
  geom_hline(yintercept=.95,lty=2,col="grey50")+
  geom_line(size=.8)+
  facet_wrap(~paste0(design,"  n1=",n1,"  n2=",n2), ncol=4)+
  scale_colour_manual(values=pal)+
  labs(x=expression(theta),
       y="Empirical coverage", colour=NULL)