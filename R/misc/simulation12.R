########################################################################
## lnM  – simulation, summary & graphics (Δ-method vs SAFE-BC)
## now includes Monte-Carlo standard errors (MCSE) and
## PROGRESS MESSAGES  (outer + inner loops)
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

## ---------- 1. Δ-method plug-in  (independent & paired) --------------
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
##  (safe_ind and safe_dep unchanged – omitted for brevity)

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
pairs_ind  <- data.frame(n1=c(5,10,20,100,3,6,12,60),
                         n2=c(5,10,20,100,7,14,28,140))

grid_ind <- expand.grid(theta=theta_vals, idx=seq_len(nrow(pairs_ind)))
grid_ind$n1 <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2 <- pairs_ind$n2[grid_ind$idx]
grid_ind$design <- "indep"; grid_ind$idx <- NULL

grid_dep <- expand.grid(theta=theta_vals, n=c(5,10,20,100))
grid_dep$n1 <- grid_dep$n2 <- grid_dep$n; grid_dep$n <- NULL
grid_dep$design <- "paired"

param_grid <- rbind(grid_ind, grid_dep)

## ---------- 5. simulation driver (with PROGRESS) ---------------------
set.seed(20250625)
K_repl <- 1e5                       # per parameter point
B_boot <- 1e5                        # bootstrap size

outer_verbose <- TRUE      # print which parameter row finished
inner_every   <- 100       # print every 100th replicate

mcse_mean <- function(x) sqrt(var(x) / length(x))
mcse_prop <- function(x) sqrt(mean(x)*(1-mean(x)) / length(x))

runner <- function(i) {
  p       <- param_grid[i, ]
  sd0     <- 1
  true_ln <- if (p$design=="indep")
    lnM_true_ind(p$theta,p$n1,p$n2,sd0)
  else
    lnM_true_dep(p$theta,p$n1,sd0)
  
  ## --- inner loop with progress -----
  M <- matrix(NA_real_, 8, K_repl,
              dimnames = list(
                c("delta_pt","delta_var","safe_pt","safe_var",
                  "safe_lo","safe_hi","safe_kept","safe_total"),
                NULL))
  for (k in seq_len(K_repl)) {
    M[,k] <- if (p$design=="indep")
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,n2=p$n2,B=B_boot)
    else
      one_rep(0,p$theta,sd0,sd0,n1=p$n1,rho=0.8,B=B_boot)
    if (k %% inner_every == 0)
      message(sprintf("[row %3d] replicate %4d / %d done",
                      i, k, K_repl))
  }
  
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
  
  if (outer_verbose)
    message(sprintf("Finished row %3d  (%s  θ=%s  n1=%d n2=%d)",
                    i, p$design, p$theta, p$n1,
                    ifelse(p$design=="indep", p$n2, p$n1)))
  
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

## optional — nice outer progress bar
pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)

results <- vector("list", nrow(param_grid))
for (i in seq_len(nrow(param_grid))) {
  results[[i]] <- runner(i)
  setTxtProgressBar(pb, i)
}
close(pb)

results <- do.call(rbind, results)

## ---------- 6. plotting section (unchanged) --------------------------
## … use `results` as in your previous script …

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

mean(abs(results$relbias_delta), na.rm = TRUE) 
mean(abs(results$relbias_safe) , na.rm = TRUE) 

mean(abs(results$mcse_bias_delta), na.rm = TRUE)  
mean(abs(results$mcse_bias_safe) , na.rm = TRUE)  

mean(abs(results$mcse_varbar_delta), na.rm = TRUE) 
mean(abs(results$mcse_varbar_safe) , na.rm = TRUE) 

