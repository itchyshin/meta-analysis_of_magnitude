########################################################################
## lnM – simulation, summary & graphics  (Delta-method vs SAFE-BC)
## * ASCII-only, multicore outer loop, verbose progress, disk persist *
## 24 independent + 12 paired designs × 18 theta  = 648 rows
## tested: 26-Jun-2025   (R ≥ 4.3, MASS 7.3-60, ggplot2 3.5-0)
########################################################################

library(MASS)        # mvrnorm()
library(ggplot2)     # plotting
library(dplyr)       # facet ordering
library(parallel)    # multicore helpers
library(pbapply)     # progress bars for *apply
library(here)       # file paths relative to this script
theme_set(theme_bw(11))

## -------- 0. globals & helpers ---------------------------------------
maxVar   <- 20
posify   <- function(x, eps = 1e-12) pmax(x, eps)
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)

## lnM kernel (ASCII)
lnM_core <- function(Delta, MSW, n0)
  0.5 * (log(Delta) - log(n0) - log(MSW))

## -------- 0a. true lnM ------------------------------------------------
lnM_true_ind <- function(theta, n1, n2, sigma = 1) {
  msb <- (n1 * n2) / (n1 + n2) * theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  lnM_core(msb - msw, msw, 2 * n1 * n2 / (n1 + n2))
}
lnM_true_dep <- function(theta, n, sigma = 1) {
  msb <- (n / 2) * theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  lnM_core(msb - msw, msw, n)
}

## -------- 1. Delta-method plug-in (indep & paired) --------------------
lnM_delta_ind <- function(x1, x2, s1, s2, n1, n2) {
  h   <- n1 * n2 / (n1 + n2)
  MSB <- h * (x1 - x2)^2
  MSW <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA, var = NA))
  pt  <- lnM_core(Delta, MSW, 2 * h)
  sD2 <- s1^2 / n1 + s2^2 / n2
  dif <- x1 - x2
  vB  <- h^2 * (2 * sD2^2 + 4 * sD2 * dif^2)
  vW  <- 2 * MSW^2 / (n1 + n2 - 2)
  g1  <- 0.5 / Delta
  g2  <- -0.5 * MSB / (Delta * MSW)
  var <- posify(g1^2 * vB + g2^2 * vW)
  c(pt = pt, var = pmin(var, maxVar))
}

lnM_delta_dep <- function(x1, x2, s1, s2, n, rho) {
  MSB <- (n / 2) * (x1 - x2)^2
  MSW <- (s1^2 + s2^2) / 2
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA, var = NA))
  pt  <- lnM_core(Delta, MSW, n)
  sD2 <- s1^2 + s2^2 - 2 * rho * s1 * s2
  dif <- x1 - x2
  vB  <- (n / 2)^2 * (2 * sD2^2 / n^2 + 4 * dif^2 * sD2 / n)
  vW  <- (s1^4 + s2^4 + 2 * rho^2 * s1^2 * s2^2) / (2 * (n - 1))
  g1  <- 0.5 / Delta
  g2  <- -0.5 * MSB / (Delta * MSW)
  var <- posify(g1^2 * vB + g2^2 * vW)
  c(pt = pt, var = pmin(var, maxVar))
}

## -------- 2. SAFE-BC bootstrap (indep & paired) -----------------------
## (only the kernel name changed → lnM_core)
safe_ind <- function(x1bar, x2bar, s1, s2, n1, n2,
                     B = 1e4, chunk = 5e3) {
  h  <- n1 * n2 / (n1 + n2); n0 <- 2 * h
  MSB0 <- h * (x1bar - x2bar)^2
  MSW0 <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  if (MSB0 <= MSW0)
    return(list(pt = NA, var = NA, lo = NA, hi = NA,
                kept = 0L, total = B))
  z_raw <- lnM_core(MSB0 - MSW0, MSW0, n0)
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- diag(c(s1^2 / n1, s2^2 / n2,
                2 * s1^4 / (n1 - 1), 2 * s2^4 / (n2 - 1)))
  cloud <- numeric(B); kept <- 0L; k <- 0L
  while (k < B) {
    d  <- mvrnorm(chunk, mu, Sig)
    ok <- d[, 3] > 0 & d[, 4] > 0
    if (!any(ok)) next
    m1 <- d[ok, 1]; m2 <- d[ok, 2]
    v1 <- d[ok, 3]; v2 <- d[ok, 4]
    MSB <- h * (m1 - m2)^2
    MSW <- ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    use <- MSB > MSW
    if (!any(use)) next
    kept <- kept + sum(use)
    vals <- lnM_core(MSB[use] - MSW[use], MSW[use], n0)
    take <- min(length(vals), B - k)
    cloud[(k + 1):(k + take)] <- vals[1:take]; k <- k + take
  }
  cloud <- cloud[1:B]
  m_boot <- mean(cloud); pt <- 2 * z_raw - m_boot; v_est <- var(cloud)
  cen <- cloud - m_boot; qs <- quantile(cen, c(0.025, 0.975))
  list(pt = pt, var = v_est,
       lo = pt + qs[1], hi = pt + qs[2],
       kept = kept, total = B)
}
safe_dep <- function(x1bar, x2bar, s1, s2, n, rho,
                     B = 1e4, chunk = 5e3) {
  MSB0 <- (n / 2) * (x1bar - x2bar)^2
  MSW0 <- (s1^2 + s2^2) / 2
  if (MSB0 <= MSW0)
    return(list(pt = NA, var = NA, lo = NA, hi = NA,
                kept = 0L, total = B))
  z_raw <- lnM_core(MSB0 - MSW0, MSW0, n)
  mu  <- c(x1bar, x2bar, s1^2, s2^2)
  Sig <- matrix(0, 4, 4)
  Sig[1, 1] <- s1^2 / n;  Sig[2, 2] <- s2^2 / n
  Sig[1, 2] <- Sig[2, 1] <- rho * s1 * s2 / n
  Sig[3, 3] <- 2 * s1^4 / (n - 1); Sig[4, 4] <- 2 * s2^4 / (n - 1)
  Sig[3, 4] <- Sig[4, 3] <- 2 * rho^2 * s1^2 * s2^2 / (n - 1)
  cloud <- numeric(B); kept <- 0L; k <- 0L
  while (k < B) {
    d <- mvrnorm(chunk, mu, Sig)
    ok <- d[, 3] > 0 & d[, 4] > 0
    if (!any(ok)) next
    m1 <- d[ok, 1]; m2 <- d[ok, 2]
    v1 <- d[ok, 3]; v2 <- d[ok, 4]
    MSB <- (n / 2) * (m1 - m2)^2
    MSW <- (v1 + v2) / 2
    use <- MSB > MSW
    if (!any(use)) next
    kept <- kept + sum(use)
    vals <- lnM_core(MSB[use] - MSW[use], MSW[use], n)
    take <- min(length(vals), B - k)
    cloud[(k + 1):(k + take)] <- vals[1:take]; k <- k + take
  }
  cloud <- cloud[1:B]
  m_boot <- mean(cloud); pt <- 2 * z_raw - m_boot; v_est <- var(cloud)
  cen <- cloud - m_boot; qs <- quantile(cen, c(0.025, 0.975))
  list(pt = pt, var = v_est,
       lo = pt + qs[1], hi = pt + qs[2],
       kept = kept, total = B)
}

## -------- 3. one replicate -------------------------------------------
one_rep <- function(mu1, mu2, sd1, sd2,
                    n1, n2 = NULL, rho = 0, B = 1e4) {
  if (is.null(n2)) {                      # paired
    Sigma <- matrix(c(sd1^2, rho * sd1 * sd2,
                      rho * sd1 * sd2, sd2^2), 2)
    xy <- mvrnorm(n1, c(mu1, mu2), Sigma)
    x1 <- xy[, 1]; x2 <- xy[, 2]; rho_hat <- cor(x1, x2)
  } else {                                # independent
    x1 <- rnorm(n1, mu1, sd1)
    x2 <- rnorm(n2, mu2, sd2); rho_hat <- 0
  }
  x1bar <- mean(x1); s1 <- sd(x1)
  x2bar <- mean(x2); s2 <- sd(x2)
  d <- if (is.null(n2))
    lnM_delta_dep(x1bar, x2bar, s1, s2, n1, rho_hat)
  else
    lnM_delta_ind(x1bar, x2bar, s1, s2, n1, n2)
  b <- if (is.null(n2))
    safe_dep (x1bar, x2bar, s1, s2, n1, rho_hat, B)
  else
    safe_ind (x1bar, x2bar, s1, s2, n1, n2, B)
  c(delta_pt  = d["pt"],  delta_var = d["var"],
    safe_pt   = b$pt,     safe_var  = b$var,
    safe_lo   = b$lo,     safe_hi   = b$hi,
    safe_kept = b$kept,   safe_total = b$total)
}

## -------- 4. parameter grid ------------------------------------------
theta_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4, 5)

pairs_ind <- data.frame(n1 = c(5, 10, 20, 100, 3, 6, 12, 60),
                        n2 = c(5, 10, 20, 100, 7, 14, 28, 140))

grid_ind <- expand.grid(theta = theta_vals,
                        idx = seq_len(nrow(pairs_ind)))
grid_ind$n1     <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2     <- pairs_ind$n2[grid_ind$idx]
grid_ind$design <- "indep"
grid_ind$idx    <- NULL

grid_dep <- expand.grid(theta = theta_vals,
                        n = c(5, 10, 20, 100))
grid_dep$n1     <- grid_dep$n
grid_dep$n2     <- grid_dep$n
grid_dep$n      <- NULL
grid_dep$design <- "paired"

param_grid <- rbind(grid_ind, grid_dep)

## -------- 5. simulation driver ---------------------------------------
set.seed(20250625)
K_repl <- 1e4      # <-- dial back from 1e5 for demo / speed
B_boot <- 1e5

outer_verbose <- TRUE
inner_every   <- 100

mcse_mean <- function(x) sqrt(var(x) / length(x))
mcse_prop <- function(x) sqrt(mean(x) * (1 - mean(x)) / length(x))

runner <- function(i) {
  p   <- param_grid[i, ]
  sd0 <- 1
  true_ln <- if (p$design == "indep")
    lnM_true_ind(p$theta, p$n1, p$n2, sd0)
  else
    lnM_true_dep(p$theta, p$n1, sd0)
  
  M <- matrix(NA_real_, 8, K_repl,
              dimnames = list(
                c("delta_pt", "delta_var", "safe_pt", "safe_var",
                  "safe_lo", "safe_hi", "safe_kept", "safe_total"),
                NULL))
  for (k in seq_len(K_repl)) {
    M[, k] <- if (p$design == "indep")
      one_rep(0, p$theta, sd0, sd0, n1 = p$n1, n2 = p$n2, B = B_boot)
    else
      one_rep(0, p$theta, sd0, sd0, n1 = p$n1, rho = 0.8, B = B_boot)
    if (k %% inner_every == 0)
      message(sprintf("[row %3d] replicate %4d / %d done",
                      i, k, K_repl))
  }
  
  ok  <- !is.na(M["delta_pt", ])
  Mok <- M[, ok, drop = FALSE]
  
  tv_s <- var(Mok["safe_pt", ])
  
  cover_d <- abs(Mok["delta_pt", ] - true_ln) <=
    1.96 * sqrt(Mok["delta_var", ])
  cover_s <- abs(Mok["safe_pt", ]  - true_ln) <=
    1.96 * sqrt(Mok["safe_var", ])
  
  out <- data.frame(
    theta          = p$theta,
    design         = p$design,
    n1             = p$n1,
    n2             = ifelse(p$design == "indep", p$n2, p$n1),
    true_lnM       = true_ln,
    delta_mean     = mean(Mok["delta_pt", ]),
    safe_mean      = mean(Mok["safe_pt",  ]),
    delta_bias     = mean(Mok["delta_pt", ]) - true_ln,
    safe_bias      = mean(Mok["safe_pt",  ]) - true_ln,
    mean_var_delta = mean(Mok["delta_var", ]),
    mean_var_safe  = mean(Mok["safe_var",  ]),
    relbias_delta  = 100 * (mean(Mok["delta_var", ]) / tv_s - 1),
    relbias_safe   = 100 * (mean(Mok["safe_var",  ]) / tv_s - 1),
    rmse_delta     = sqrt(mean((Mok["delta_pt", ] - true_ln)^2)),
    rmse_safe      = sqrt(mean((Mok["safe_pt",  ] - true_ln)^2)),
    cover_delta    = mean(cover_d),
    cover_safe     = mean(cover_s),
    boot_keep      = sum(Mok["safe_kept", ]),
    boot_total     = sum(Mok["safe_total", ]),
    mcse_bias_delta   = mcse_mean(Mok["delta_pt", ] - true_ln),
    mcse_bias_safe    = mcse_mean(Mok["safe_pt",  ] - true_ln),
    mcse_varbar_delta = mcse_mean(Mok["delta_var", ]),
    mcse_varbar_safe  = mcse_mean(Mok["safe_var", ]),
    mcse_cover_delta  = mcse_prop(cover_d),
    mcse_cover_safe   = mcse_prop(cover_s)
  )
  
  ## attach raw replicate matrix so caller can persist it
  attr(out, "raw_M") <- M
  
  if (outer_verbose)
    message(sprintf("Finished row %3d  (%s  theta=%s  n1=%d n2=%d)",
                    i, p$design, p$theta, p$n1,
                    ifelse(p$design == "indep", p$n2, p$n1)))
  out
}

## -------- 5a. PARALLEL outer loop ------------------------------------
## use all cores except two (never less than one)
n_cores_use <- max(1L, detectCores() - 92)

pbop <- pbapply::pboptions(type = "txt")   # progress bar style

if (.Platform$OS.type == "windows") {
  ## ---- Windows: PSOCK cluster ---------------------------------------
  cl <- makeCluster(n_cores_use)
  clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)
  results_list <- pbapply::pblapply(
    seq_len(nrow(param_grid)), runner, cl = cl)
  stopCluster(cl)
} else {
  ## ---- Unix/macOS: forked workers -----------------------------------
  results_list <- pbapply::pblapply(
    seq_len(nrow(param_grid)), runner, cl = n_cores_use)
}

pbapply::pboptions(pbop)  # restore previous pbapply settings

## -------- 5c. COLLAPSE + SERIALISE -----------------------------------
results <- do.call(rbind, results_list)

## summary RDS
summary_file <- sprintf("lnM_summary_%s.rds", Sys.Date())
saveRDS(results, summary_file)
message("Saved overall summary to: ", summary_file)

## (optional) save each raw replicate matrix
save_raw <- TRUE     # flip to FALSE to skip heavy files
if (save_raw) {
  dir.create("raw_runs", showWarnings = FALSE)
  for (i in seq_along(results_list)) {
    raw_i <- attr(results_list[[i]], "raw_M")
    if (is.null(raw_i)) next
    saveRDS(raw_i,
            sprintf("raw_runs/row_%03d.rds", i),
            compress = "xz")
  }
  message("Saved ", length(results_list),
          " raw replicate files into raw_runs/")
}

## also CSV for spreadsheets
write.csv(results,
          file = sprintf("lnM_summary_%s.csv", Sys.Date()),
          row.names = FALSE)

## -------- 6. quick Bias plot example ---------------------------------
results <- readRDS(here("Rdata", "lnM_summary_2025-06-26.rds"))

# 2. Make a combined facet label ----------------------------------------

results <- results %>%
  mutate(facet_label = paste0(design, " n1=", n1, " n2=", n2))

# 3. Define the new ordering: 3 rows × 4 cols --------------------------

facet_info <- results %>%
  distinct(facet_label, design, n1, n2)

balanced_ind   <- facet_info %>%
  filter(design == "indep", n1 == n2) %>%
  arrange(n1) %>%
  pull(facet_label)

unbalanced_ind <- facet_info %>%
  filter(design == "indep", n1 != n2) %>%
  arrange(n1) %>%
  pull(facet_label)

paired_balanced <- facet_info %>%
  filter(design == "paired") %>%
  arrange(n1) %>%
  pull(facet_label)

new_levels <- c(balanced_ind, unbalanced_ind, paired_balanced)

results <- results %>%
  mutate(facet_label = factor(facet_label, levels = new_levels))

# 4. Absolute‐bias plot -------------------------------------------------

bias_df <- bind_rows(
  results %>% select(theta, facet_label, delta_bias) %>%
    rename(bias = delta_bias) %>% mutate(estimator = "delta"),
  results %>% select(theta, facet_label, safe_bias) %>%
    rename(bias = safe_bias)  %>% mutate(estimator = "SAFE")
)

p_bias <- ggplot(bias_df,
                 aes(theta, bias, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
  labs(x = "theta",
       y = "Bias (estimate \u2212 true lnM)",
       colour = "Estimator") +
  scale_colour_manual(
    values = c(delta = "firebrick", SAFE = "steelblue"),
    labels = c(delta = "PI", SAFE = "SAFE")
  ) +
  theme_bw(11)

print(p_bias)

# 5. Relative‐bias of variance plot ------------------------------------

rb_df <- bind_rows(
  results %>% select(theta, facet_label, relbias_delta) %>%
    rename(relbias = relbias_delta) %>% mutate(estimator = "delta"),
  results %>% select(theta, facet_label, relbias_safe) %>%
    rename(relbias = relbias_safe)  %>% mutate(estimator = "SAFE")
)

p_relbias <- ggplot(rb_df,
                    aes(theta, relbias, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
  labs(x = "theta",
       y = "Relative bias of Var (%)",
       colour = "Estimator") +
  scale_colour_manual(
    values = c(delta = "firebrick", SAFE = "steelblue"),
    labels = c(delta = "PI", SAFE = "SAFE")
  ) +
  theme_bw(11)

print(p_relbias)