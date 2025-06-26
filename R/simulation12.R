########################################################################
## lnM  – simulation, summary & graphics (Δ-method vs SAFE-BC)
## 24 independent + 12 paired designs × 14 θ  = 168 rows
## last tested: 25-Jun-2025  (R ≥ 4.3, MASS 7.3-60, ggplot2 3.5-0)
########################################################################

library(MASS)       # for mvrnorm()
library(ggplot2)    # for plotting
theme_set(theme_bw(11))

## ---------- 0. globals & helpers -------------------------------------
maxVar   <- 20       # cap for Δ-method variance
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)

## small internal lnM for bootstrap
.lnM <- function(Δ, MSW, n0) {
  0.5 * ( log(Δ) - log(n0) - log(MSW) )
}

## ---------- 0a. true-lnM ------------------------------------------------
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

## ---------- 1. Δ-method plug-in ----------------------------------------
lnM_delta_ind <- function(x1, x2, s1, s2, n1, n2) {
  h   <- n1 * n2 / (n1 + n2)
  MSB <- h * (x1 - x2)^2
  MSW <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2)
  Δ   <- safe_gap(MSB - MSW)
  if (is.na(Δ)) return(c(pt = NA, var = NA))
  n0  <- 2 * h
  pt  <- 0.5 * (log(Δ) - log(n0) - log(MSW))
  ## delta-method variance
  sD2 <- s1^2/n1 + s2^2/n2
  dif <- x1 - x2
  vB  <- h^2 * (2*sD2^2 + 4*sD2*dif^2)
  vW  <- 2*MSW^2 / (n1 + n2 - 2)
  g1  <- 0.5/Δ
  g2  <- -0.5*MSB/(Δ*MSW)
  var <- posify(g1^2 * vB + g2^2 * vW)
  ## cap any extreme variance:
  var <- pmin(var, maxVar)
  c(pt = pt, var = var)
}

lnM_delta_dep <- function(x1, x2, s1, s2, n, rho) {
  MSB <- (n/2) * (x1 - x2)^2
  MSW <- (s1^2 + s2^2)/2
  Δ   <- safe_gap(MSB - MSW)
  if (is.na(Δ)) return(c(pt = NA, var = NA))
  pt  <- 0.5 * (log(Δ) - log(n) - log(MSW))
  ## delta-method variance
  sD2 <- s1^2 + s2^2 - 2*rho*s1*s2
  dif <- x1 - x2
  vB  <- (n/2)^2 * (2*sD2^2/n^2 + 4*dif^2*sD2/n)
  vW  <- (s1^4 + s2^4 + 2*rho^2*s1^2*s2^2) / (2*(n-1))
  g1  <- 0.5/Δ
  g2  <- -0.5*MSB/(Δ*MSW)
  var <- posify(g1^2 * vB + g2^2 * vW)
  ## cap any extreme variance:
  var <- pmin(var, maxVar)
  c(pt = pt, var = var)
}

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

## ---------- 3.–6. rest of your simulation & plotting -----------------
## (unchanged from your existing code)
## ---------- 3. one replicate (collect both estimators) -------------
one_rep <- function(mu1, mu2, sd1, sd2,
                    n1, n2 = NULL, rho = 0, B = 1e4) {
  
  if (is.null(n2)) {
    Σ  <- matrix(c(sd1^2, rho*sd1*sd2,
                   rho*sd1*sd2, sd2^2), 2, 2)
    xy <- mvrnorm(n1, c(mu1, mu2), Σ)
    x1 <- xy[,1]; x2 <- xy[,2]
    rho_hat <- cor(x1, x2)
  } else {
    x1 <- rnorm(n1, mu1, sd1)
    x2 <- rnorm(n2, mu2, sd2)
    rho_hat <- 0
  }
  
  x1bar <- mean(x1); s1 <- sd(x1)
  x2bar <- mean(x2); s2 <- sd(x2)
  
  d <- if (is.null(n2))
    lnM_delta_dep(x1bar, x2bar, s1, s2, n1, rho_hat)
  else
    lnM_delta_ind(x1bar, x2bar, s1, s2, n1, n2)
  
  b <- if (is.null(n2))
    safe_dep(x1bar, x2bar, s1, s2, n1, rho_hat, B)
  else
    safe_ind(x1bar, x2bar, s1, s2, n1, n2, B)
  
  c(delta_pt    = d["pt"],
    delta_var   = d["var"],
    safe_pt     = b$pt,
    safe_var    = b$var,
    safe_lo     = b$lo,
    safe_hi     = b$hi,
    safe_kept   = b$kept,
    safe_total  = b$total)
}

## ---------- 4. build parameter grid -------------------------------
theta_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4, 5)
pairs_ind  <- data.frame(
  n1 = c(5,10,20,100,3,6,12,60),
  n2 = c(5,10,20,100,7,14,28,140)
)

grid_ind <- expand.grid(theta = theta_vals,
                        idx   = seq_len(nrow(pairs_ind)),
                        KEEP.OUT.ATTRS = FALSE)
grid_ind$n1     <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2     <- pairs_ind$n2[grid_ind$idx]
grid_ind$rho    <- 0
grid_ind$design <- "indep"
grid_ind$idx    <- NULL

grid_dep <- expand.grid(theta = theta_vals,
                        n     = c(5,10,20,100),
                        KEEP.OUT.ATTRS = FALSE)
grid_dep$rho    <- 0.8
grid_dep$design <- "paired"
grid_dep$n1     <- grid_dep$n
grid_dep$n2     <- grid_dep$n
grid_dep$n      <- NULL
grid_dep        <- grid_dep[, names(grid_ind)]

param_grid <- rbind(grid_ind, grid_dep)

## ---------- 5. simulation driver ------------------------------------
set.seed(20250625)
K_repl <- 1000   # replicates per parameter point (true one 1e5)
B_boot <- 1e4    # bootstrap size (true one 1e6)

runner <- function(i) {
  p       <- param_grid[i, ]
  sd0     <- 1
  true_ln <- if (p$design == "indep")
    lnM_true_ind(p$theta, p$n1, p$n2, sd0)
  else
    lnM_true_dep(p$theta, p$n1, sd0)
  
  reps <- replicate(
    K_repl,
    if (p$design=="indep")
      one_rep(0, p$theta, sd0, sd0, n1=p$n1, n2=p$n2, B=B_boot)
    else
      one_rep(0, p$theta, sd0, sd0, n1=p$n1, rho=p$rho, B=B_boot),
    simplify = FALSE
  )
  
  M <- do.call(cbind, reps)
  rownames(M) <- c("delta_pt","delta_var",
                   "safe_pt","safe_var",
                   "safe_lo","safe_hi",
                   "safe_kept","safe_total")
  
  ok  <- !is.na(M["delta_pt", ])
  Mok <- M[, ok, drop=FALSE]
  
  # true MC variances
  tv_d <- var(Mok["delta_pt", ])
  tv_s <- var(Mok["safe_pt",  ])
  
  # coverage
  cov_d <- mean(abs(Mok["delta_pt",]-true_ln) <=
                  1.96*sqrt(Mok["delta_var",]), na.rm=TRUE)
  #cov_s <- mean(Mok["safe_lo",] <= true_ln &
  #                true_ln <= Mok["safe_hi",], na.rm=TRUE)
  cov_s <- mean(
    abs(Mok["safe_pt",] - true_ln)
    <= 1.96*sqrt(Mok["safe_var",]),
    na.rm=TRUE
  )
  
  
  # relative bias of variance estimators
  # the best estimator's MC (simulated variance)
  rb_d <- 100 * (mean(Mok["delta_var", ]) / tv_s - 1) # changed tv_d
  rb_s <- 100 * (mean(Mok["safe_var",  ]) / tv_s - 1)
  
  data.frame(
    theta          = p$theta,
    design         = p$design,
    n1             = p$n1,
    n2             = p$n2,
    true_lnM       = true_ln,
    true_var_delta = tv_d,
    true_var_safe  = tv_s,
    delta_mean     = mean(Mok["delta_pt", ]),
    safe_mean      = mean(Mok["safe_pt",  ]),
    delta_bias     = mean(Mok["delta_pt", ]) - true_ln,
    safe_bias      = mean(Mok["safe_pt",  ]) - true_ln,
    mean_var_delta = mean(Mok["delta_var", ]),
    mean_var_safe  = mean(Mok["safe_var",  ]),
    relbias_delta  = rb_d,
    relbias_safe   = rb_s,
    rmse_delta     = sqrt(mean((Mok["delta_pt", ] - true_ln)^2)),
    rmse_safe      = sqrt(mean((Mok["safe_pt",  ] - true_ln)^2)),
    cover_delta    = cov_d,
    cover_safe     = cov_s,
    boot_keep      = sum(Mok["safe_kept", ]),
    boot_total     = sum(Mok["safe_total", ])
  )
}

res_list <- lapply(seq_len(nrow(param_grid)), runner)
results  <- do.call(rbind, res_list)

## ---------- 6. graphics ------------------------------------------------
# A. Generate facet label
results$facet_label <- with(results, paste0(design, " n1=", n1, " n2=", n2))

# A.1 Derive desired order from param_grid
param_grid$facet_label <- with(param_grid, paste0(design, " n1=", n1, " n2=", n2))

# A.2 Create facet order: indep balanced -> indep unbalanced -> paired
facet_order <- param_grid |>
  dplyr::distinct(facet_label, design, n1, n2) |>
  dplyr::mutate(group = dplyr::case_when(
    design == "paired"            ~ 3,
    design == "indep" & n1 == n2  ~ 1,
    TRUE                          ~ 2
  )) |>
  dplyr::arrange(group, n1 + n2) |>
  dplyr::pull(facet_label)

# A.3 Convert to factor with desired order
results$facet_label <- factor(results$facet_label, levels = facet_order)

# A.4 Create bias_df using facet_label
bias_df <- rbind(
  data.frame(results[, c("theta", "facet_label")],
             estimator = "delta", bias = results$delta_bias),
  data.frame(results[, c("theta", "facet_label")],
             estimator = "SAFE",  bias = results$safe_bias)
)

# A.5 Plot Bias
p_bias <- ggplot(bias_df,
                 aes(theta, bias, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4) +
  labs(x = expression(theta), y = "Bias (estimate − true lnM)") +
  #scale_colour_manual(values = c(delta = "firebrick", SAFE = "steelblue")) +
  theme_bw(11)

#print(p_bias)



# A.6 Optional: Log-scaled theta axis
p_bias_log10 <- p_bias +
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5),
                labels = c("0.2", "0.5", "1", "2", "5"))

#p_bias_log10

p_bias_inset_topright <-
  p_bias_log10 +
  
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(
      fill   = "white",
      colour = "grey80"
    ),
    legend.key.size      = unit(0.8, "lines")
  )

p_bias_inset_topright + scale_colour_manual(values = c(delta = "firebrick",
                                                              SAFE  = "steelblue"),
                                                   labels = c(delta = "PI",  # displayed text
                                                              SAFE  = "SAFE")) 


# B. relative bias of variance

# 1. build a single label column, in the exact order of param_grid
param_grid$lab <- with(param_grid,
                       paste0(design, "   n1=", n1, "   n2=", n2)
)

# capture the unique labels *in the order they appear* in param_grid
lab_levels <- unique(param_grid$lab)

# make it a factor with those levels
param_grid$lab <- factor(param_grid$lab, levels = lab_levels)

# 2. when you build your rb_df, carry that lab column along
rb_df <- rbind(
  data.frame(
    results[, c("theta","n1","n2","design")],
    estimator = "delta",
    relbias   = results$relbias_delta
  ),
  data.frame(
    results[, c("theta","n1","n2","design")],
    estimator = "SAFE",
    relbias   = results$relbias_safe
  )
)
# join the lab factor back on
rb_df$lab <- with(rb_df,
                  factor(
                    paste0(design, "   n1=", n1, "   n2=", n2),
                    levels = lab_levels
                  )
)

# 3. facet by lab
p_rb_log10 <- ggplot(rb_df, aes(theta, relbias,
                                colour=estimator, group=estimator)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_line() +
  facet_wrap(~ lab, ncol = 4) +
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5),
                labels = c("0.2","0.5","1","2","5")) +
  labs(x = expression(theta),
       y = "Relative bias of Var (%)") +
  scale_colour_manual(values = c(delta="firebrick", SAFE="steelblue")) +
  theme_bw(11)

#print(p_rb_log10)

p_rb_log10_inset_topright <-
  p_rb_log10 +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(
      fill   = "white",
      colour = "grey80"
    ),
    legend.key.size      = unit(0.8, "lines")
  )

print(p_rb_log10_inset_topright)


# C. coverage
cov_df <- rbind(
  data.frame(results[,c("theta","design","n1","n2")],
             estimator="delta", cover=results$cover_delta),
  data.frame(results[,c("theta","design","n1","n2")],
             estimator="SAFE",  cover=results$cover_safe)
)
p_cov <- ggplot(cov_df,
                aes(theta, cover, colour=estimator, group=estimator)) +
  geom_hline(yintercept=0.95, linetype=2, colour="grey50") +
  geom_line() +
  facet_wrap(~ paste0(design," n1=",n1," n2=",n2), ncol=4) +
  labs(x=expression(theta), y="Empirical coverage") +
  scale_colour_manual(values=c(delta="firebrick", SAFE="steelblue"))
#print(p_cov)

# D. RMSE
rmse_df <- rbind(
  data.frame(results[,c("theta","design","n1","n2")],
             estimator="delta", rmse=results$rmse_delta),
  data.frame(results[,c("theta","design","n1","n2")],
             estimator="SAFE",  rmse=results$rmse_safe)
)
p_rmse <- ggplot(rmse_df,
                 aes(theta, rmse, colour=estimator, group=estimator)) +
  geom_line(size=0.8) +
  facet_wrap(~ paste0(design," n1=",n1," n2=",n2), ncol=4) +
  scale_colour_manual(values=c(delta="firebrick", SAFE="steelblue")) +
  labs(x=expression(theta),
       y="RMSE ( ln M̂ − ln M )",
       colour=NULL) +
  theme_bw(11)
#print(p_rmse)

# average abs(bias)
mean(abs(results$delta_bias), na.rm = TRUE)  # average Δ-method bias
mean(abs(results$safe_bias), na.rm = TRUE)  # average SAFE-BC bias

mean(abs(results$relbias_delta), na.rm = TRUE)  # average Δ-method bias
mean(abs(results$relbias_safe) , na.rm = TRUE)  # average SAFE-BC bias
