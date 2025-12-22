########################################################################
## lnM – simulation, summary & graphics  (Delta-method vs SAFE-BC)
## * ASCII-only, multicore outer loop, verbose progress, disk persist *
##
## This version uses SAFE_fun.R (safe_lnM_indep / safe_lnM_dep) and
## targets a minimum number of ACCEPTED bootstrap draws (min_kept).
##
## Relative bias of variance: baseline = MC Var of SAFE-BC point estimator.
##   (i.e., MC variance of M["safe_pt", ] after applying SAFE-BC correction)
########################################################################

library(MASS)        # mvrnorm()
library(ggplot2)     # plotting
library(dplyr)       # facet ordering
library(parallel)    # multicore helpers
library(pbapply)     # progress bars for *apply
library(here)        # file paths relative to this script
theme_set(theme_bw(11))

## -------- 0. globals & helpers ---------------------------------------
maxVar   <- 20
posify   <- function(x, eps = 1e-12) pmax(x, eps) # for flowing point issue
safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)

## lnM kernel (ASCII)
lnM_core <- function(Delta, MSW, n0)
  0.5 * (log(Delta) - log(n0) - log(MSW))

## -------- 0b. load SAFE_fun.R ----------------------------------------
SAFE_FILE <- here("R", "SAFE_fun.R")
if (!file.exists(SAFE_FILE)) stop("SAFE_fun.R not found at: ", SAFE_FILE)
source(SAFE_FILE)
if (!exists("safe_lnM_indep") || !exists("safe_lnM_dep")) {
  stop("SAFE_fun.R must define safe_lnM_indep() and safe_lnM_dep().")
}

## A small helper so we do not care about exact list field names
safe_get <- function(S, name, default = NA) {
  if (!is.list(S)) return(default)
  if (!is.null(S[[name]])) return(S[[name]])
  alt <- switch(name,
                pt     = c("point"),
                var    = c("v"),
                kept   = c("keep","accepted"),
                tried  = c("total","draws"),
                status = c("state"),
                character(0))
  for (a in alt) if (!is.null(S[[a]])) return(S[[a]])
  default
}

## -------- 0c. true lnM ------------------------------------------------
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
## returns (pt, var, capped)
lnM_delta_ind <- function(x1, x2, s1, s2, n1, n2) {
  h   <- n1 * n2 / (n1 + n2)
  MSB <- h * (x1 - x2)^2
  MSW <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA_real_, var = NA_real_, capped = NA_real_))
  
  pt  <- lnM_core(Delta, MSW, 2 * h)
  
  sD2 <- s1^2 / n1 + s2^2 / n2
  dif <- x1 - x2
  vB  <- h^2 * (2 * sD2^2 + 4 * sD2 * dif^2)
  vW  <- 2 * MSW^2 / (n1 + n2 - 2)
  g1  <- 0.5 / Delta
  g2  <- -0.5 * MSB / (Delta * MSW)
  var_raw <- posify(g1^2 * vB + g2^2 * vW)
  
  capped <- as.integer(is.finite(var_raw) && var_raw > maxVar)
  c(pt = pt, var = pmin(var_raw, maxVar), capped = capped)
}

lnM_delta_dep <- function(x1, x2, s1, s2, n, rho) {
  MSB <- (n / 2) * (x1 - x2)^2
  MSW <- (s1^2 + s2^2) / 2
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA_real_, var = NA_real_, capped = NA_real_))
  
  pt  <- lnM_core(Delta, MSW, n)
  
  sD2 <- s1^2 + s2^2 - 2 * rho * s1 * s2
  dif <- x1 - x2
  vB  <- (n / 2)^2 * (2 * sD2^2 / n^2 + 4 * dif^2 * sD2 / n)
  vW  <- (s1^4 + s2^4 + 2 * rho^2 * s1^2 * s2^2) / (2 * (n - 1))
  g1  <- 0.5 / Delta
  g2  <- -0.5 * MSB / (Delta * MSW)
  var_raw <- posify(g1^2 * vB + g2^2 * vW)
  
  capped <- as.integer(is.finite(var_raw) && var_raw > maxVar)
  c(pt = pt, var = pmin(var_raw, maxVar), capped = capped)
}

## -------- 2. SAFE-BC via SAFE_fun.R ----------------------------------
safe_ind_fun <- function(x1bar, x2bar, s1, s2, n1, n2,
                         min_kept = 2000, chunk_init = 4000,
                         chunk_max = 2e6, max_draws = Inf,
                         patience_noaccept = 5) {
  
  S <- safe_lnM_indep(x1bar, x2bar, s1, s2, n1, n2,
                      min_kept = min_kept,
                      chunk_init = chunk_init,
                      chunk_max = chunk_max,
                      max_draws = max_draws,
                      patience_noaccept = patience_noaccept)
  
  list(pt     = safe_get(S, "pt", NA_real_),
       var    = safe_get(S, "var", NA_real_),
       kept   = safe_get(S, "kept", NA_integer_),
       tried  = safe_get(S, "tried", NA_integer_),
       status = safe_get(S, "status", NA_character_))
}

safe_dep_fun <- function(x1bar, x2bar, s1, s2, n, rho,
                         min_kept = 2000, chunk_init = 4000,
                         chunk_max = 2e6, max_draws = Inf,
                         patience_noaccept = 5) {
  
  S <- safe_lnM_dep(x1bar, x2bar, s1, s2, n, rho,
                    min_kept = min_kept,
                    chunk_init = chunk_init,
                    chunk_max = chunk_max,
                    max_draws = max_draws,
                    patience_noaccept = patience_noaccept)
  
  list(pt     = safe_get(S, "pt", NA_real_),
       var    = safe_get(S, "var", NA_real_),
       kept   = safe_get(S, "kept", NA_integer_),
       tried  = safe_get(S, "tried", NA_integer_),
       status = safe_get(S, "status", NA_character_))
}

## -------- 3. one replicate -------------------------------------------
one_rep <- function(mu1, mu2, sd1, sd2,
                    n1, n2 = NULL, rho = 0,
                    min_kept = 2000, chunk_init = 4000,
                    chunk_max = 2e6, max_draws = Inf,
                    patience_noaccept = 5) {
  
  if (is.null(n2)) {                      # paired
    Sigma <- matrix(c(sd1^2, rho * sd1 * sd2,
                      rho * sd1 * sd2, sd2^2), 2)
    xy <- mvrnorm(n1, c(mu1, mu2), Sigma)
    x1 <- xy[, 1]; x2 <- xy[, 2]
    rho_hat <- suppressWarnings(cor(x1, x2))
    if (!is.finite(rho_hat)) rho_hat <- rho
  } else {                                # independent
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
  
  s <- if (is.null(n2))
    safe_dep_fun(x1bar, x2bar, s1, s2, n1, rho_hat,
                 min_kept, chunk_init, chunk_max, max_draws, patience_noaccept)
  else
    safe_ind_fun(x1bar, x2bar, s1, s2, n1, n2,
                 min_kept, chunk_init, chunk_max, max_draws, patience_noaccept)
  
  delta_pt  <- unname(d["pt"])
  safe_raw  <- s$pt
  
  ## SAFE-BC: 2*PI - SAFE_mean (when PI exists)
  safe_bc <- if (is.na(delta_pt)) safe_raw else 2 * delta_pt - safe_raw
  
  c(delta_pt    = delta_pt,
    delta_var   = unname(d["var"]),
    delta_cap   = unname(d["capped"]),
    safe_pt     = safe_bc,              # SAFE-BC point estimator
    safe_var    = s$var,                # SAFE variance estimator
    safe_kept   = s$kept,
    safe_tried  = s$tried,
    safe_ok     = as.integer(isTRUE(s$status == "ok")),
    delta_fail  = as.integer(is.na(delta_pt)),
    safe_fail   = as.integer(is.na(safe_bc)))
}

## -------- 4. parameter grid ------------------------------------------
theta_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4, 5)

pairs_ind <- data.frame(n1 = c(5, 10, 20, 100, 3, 6, 12, 40),
                        n2 = c(5, 10, 20, 100, 7, 14, 28, 160))

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

K_repl     <- as.integer(Sys.getenv("K_REPL", "200"))     # demo default
MIN_KEPT   <- as.integer(Sys.getenv("MIN_KEPT", "2000"))   # accepted SAFE draws target
CHUNK_INIT <- as.integer(Sys.getenv("CHUNK_INIT", "5000"))
CHUNK_MAX  <- as.numeric(Sys.getenv("CHUNK_MAX", "2000000"))
MAX_DRAWS  <- as.numeric(Sys.getenv("MAX_DRAWS", "Inf"))
PATIENCE   <- as.integer(Sys.getenv("PATIENCE", "5"))

outer_verbose <- TRUE
inner_every   <- 100

mcse_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sqrt(var(x) / length(x))
}
mcse_prop <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  p <- mean(x)
  sqrt(p * (1 - p) / length(x))
}

runner <- function(i) {
  p   <- param_grid[i, ]
  sd0 <- 1
  
  true_ln <- if (p$design == "indep")
    lnM_true_ind(p$theta, p$n1, p$n2, sd0)
  else
    lnM_true_dep(p$theta, p$n1, sd0)
  
  M <- matrix(NA_real_, 10, K_repl,
              dimnames = list(
                c("delta_pt","delta_var","delta_cap",
                  "safe_pt","safe_var",
                  "safe_kept","safe_tried","safe_ok",
                  "delta_fail","safe_fail"),
                NULL))
  
  for (k in seq_len(K_repl)) {
    M[, k] <- if (p$design == "indep")
      one_rep(0, p$theta, sd0, sd0, n1 = p$n1, n2 = p$n2, rho = 0,
              min_kept = MIN_KEPT, chunk_init = CHUNK_INIT,
              chunk_max = CHUNK_MAX, max_draws = MAX_DRAWS,
              patience_noaccept = PATIENCE)
    else
      one_rep(0, p$theta, sd0, sd0, n1 = p$n1, n2 = NULL, rho = 0.8,
              min_kept = MIN_KEPT, chunk_init = CHUNK_INIT,
              chunk_max = CHUNK_MAX, max_draws = MAX_DRAWS,
              patience_noaccept = PATIENCE)
    
    if (k %% inner_every == 0)
      message(sprintf("[row %3d] replicate %7d / %d done", i, k, K_repl))
  }
  
  ## keep a conservative 'ok' set for delta-variance validity
  ok_d <- is.finite(M["delta_pt", ]) & is.finite(M["delta_var", ]) & (M["delta_var", ] > 0)
  Mok  <- M[, ok_d, drop = FALSE]
  
  ## BASELINE: MC Var of SAFE-BC point estimator (computed on same Mok subset)
  ok_safe_pt <- is.finite(Mok["safe_pt", ])
  Var_MC_SAFEpt <- if (sum(ok_safe_pt) >= 2) var(Mok["safe_pt", ok_safe_pt]) else NA_real_
  
  ## coverage only meaningful if true_ln is finite
  if (is.finite(true_ln) && sum(ok_d) > 0) {
    cover_d <- abs(Mok["delta_pt", ] - true_ln) <= 1.96 * sqrt(Mok["delta_var", ])
    
    ok_s <- is.finite(Mok["safe_pt", ]) & is.finite(Mok["safe_var", ]) & (Mok["safe_var", ] > 0)
    cover_s <- abs(Mok["safe_pt", ok_s] - true_ln) <= 1.96 * sqrt(Mok["safe_var", ok_s])
    
    cover_delta <- mean(cover_d)
    cover_safe  <- if (length(cover_s)) mean(cover_s) else NA_real_
  } else {
    cover_delta <- NA_real_
    cover_safe  <- NA_real_
  }
  
  delta_fail_prop <- mean(M["delta_fail", ], na.rm = TRUE)
  safe_fail_prop  <- mean(M["safe_fail",  ], na.rm = TRUE)
  
  boot_keep  <- sum(M["safe_kept", ],  na.rm = TRUE)
  boot_tried <- sum(M["safe_tried", ], na.rm = TRUE)
  boot_accept_prop <- if (boot_tried > 0) boot_keep / boot_tried else NA_real_
  
  delta_cap_rate <- mean(M["delta_cap", ][ok_d], na.rm = TRUE)
  delta_cap_n    <- sum(M["delta_cap", ][ok_d] == 1, na.rm = TRUE)
  
  mean_var_delta <- mean(M["delta_var", ], na.rm = TRUE)
  mean_var_safe  <- mean(M["safe_var",  ], na.rm = TRUE)
  
  ## Relative bias of variance: vs MC(SAFE-BC point)
  relbias_delta <- if (is.finite(Var_MC_SAFEpt) && Var_MC_SAFEpt > 0)
    100 * (mean_var_delta / Var_MC_SAFEpt - 1) else NA_real_
  relbias_safe  <- if (is.finite(Var_MC_SAFEpt) && Var_MC_SAFEpt > 0)
    100 * (mean_var_safe  / Var_MC_SAFEpt - 1) else NA_real_
  
  out <- data.frame(
    theta      = p$theta,
    design     = p$design,
    n1         = p$n1,
    n2         = ifelse(p$design == "indep", p$n2, p$n1),
    true_lnM   = true_ln,
    
    delta_mean = mean(M["delta_pt", ], na.rm = TRUE),
    safe_mean  = mean(M["safe_pt",  ], na.rm = TRUE),
    
    delta_bias = if (is.finite(true_ln)) mean(M["delta_pt", ], na.rm = TRUE) - true_ln else NA_real_,
    safe_bias  = if (is.finite(true_ln)) mean(M["safe_pt",  ], na.rm = TRUE) - true_ln else NA_real_,
    
    mean_var_delta = mean_var_delta,
    mean_var_safe  = mean_var_safe,
    
    Var_MC_SAFEpt  = Var_MC_SAFEpt,
    relbias_delta  = relbias_delta,
    relbias_safe   = relbias_safe,
    
    rmse_delta     = if (is.finite(true_ln)) sqrt(mean((M["delta_pt", ] - true_ln)^2, na.rm = TRUE)) else NA_real_,
    rmse_safe      = if (is.finite(true_ln)) sqrt(mean((M["safe_pt",  ] - true_ln)^2, na.rm = TRUE)) else NA_real_,
    
    cover_delta    = cover_delta,
    cover_safe     = cover_safe,
    
    delta_fail_prop   = delta_fail_prop,
    safe_fail_prop    = safe_fail_prop,
    
    delta_cap_rate = delta_cap_rate,
    delta_cap_n    = delta_cap_n,
    maxVar         = maxVar,
    
    boot_keep        = boot_keep,
    boot_tried       = boot_tried,
    boot_accept_prop = boot_accept_prop,
    SAFE_ok_rate     = mean(M["safe_ok", ], na.rm = TRUE),
    
    mcse_bias_delta   = if (is.finite(true_ln)) mcse_mean(M["delta_pt", ] - true_ln) else NA_real_,
    mcse_bias_safe    = if (is.finite(true_ln)) mcse_mean(M["safe_pt",  ] - true_ln) else NA_real_,
    mcse_varbar_delta = mcse_mean(M["delta_var", ]),
    mcse_varbar_safe  = mcse_mean(M["safe_var", ]),
    mcse_cover_delta  = if (is.finite(true_ln) && sum(ok_d) > 0) mcse_prop(cover_d) else NA_real_,
    mcse_cover_safe   = if (is.finite(true_ln) && exists("cover_s") && length(cover_s)) mcse_prop(cover_s) else NA_real_
  )
  
  attr(out, "raw_M") <- M
  
  if (outer_verbose)
    message(sprintf("Finished row %3d  (%s  theta=%s  n1=%d n2=%d)",
                    i, p$design, p$theta, p$n1,
                    ifelse(p$design == "indep", p$n2, p$n1)))
  out
}

## -------- 5a. PARALLEL outer loop ------------------------------------
n_cores_use <- suppressWarnings(as.integer(Sys.getenv("N_CORES", "")))
if (is.na(n_cores_use) || n_cores_use <= 0) n_cores_use <- max(1L, detectCores() - 184)

pbop <- pbapply::pboptions(type = "txt")

if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_cores_use)
  clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)
  results_list <- pbapply::pblapply(seq_len(nrow(param_grid)), runner, cl = cl)
  stopCluster(cl)
} else {
  results_list <- pbapply::pblapply(seq_len(nrow(param_grid)), runner, cl = n_cores_use)
}

pbapply::pboptions(pbop)

## -------- 5c. COLLAPSE (The SAFE way) ------------------------

# 1. Identify which results are actual data frames (successful)
# and which are character strings (errors)
is_df <- sapply(results_list, is.data.frame)

# 2. Only combine the successful ones
results <- do.call(rbind, results_list[is_df])

# 3. Force numeric types just in case (Safety first)
results$theta <- as.numeric(as.character(results$theta))
results$n1    <- as.numeric(as.character(results$n1))

# 4. Report any failures so you know which designs to check
if(any(!is_df)) {
  message(sprintf("Warning: %d parameter sets failed and were excluded.", sum(!is_df)))
  message("Failed indices: ", paste(which(!is_df), collapse = ", "))
}

## -------- SAVE OUTPUTS -----------------------------------------------

# 1) Save overall summary (RDS + CSV)
summary_rds <- sprintf("lnM_summary_SAFEfun_%s.rds", Sys.Date())
saveRDS(results, summary_rds)
message("Saved overall summary (RDS) to: ", summary_rds)

summary_csv <- sprintf("lnM_summary_SAFEfun_%s.csv", Sys.Date())
write.csv(results, file = summary_csv, row.names = FALSE)
message("Saved overall summary (CSV) to: ", summary_csv)

# 2) Optionally save raw replicate results (one file per param-grid row)
#    Saved format: data.frame with runs as rows, stats as columns
save_raw <- TRUE

if (isTRUE(save_raw)) {
  raw_dir <- "raw_runs"
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  
  n_saved <- 0L
  for (i in seq_along(results_list)) {
    
    raw_M <- attr(results_list[[i]], "raw_M")
    if (is.null(raw_M)) next
    
    # raw_M is stats x K_REPL -> convert to K_REPL x stats
    raw_df <- as.data.frame(t(raw_M))
    
    # keep replicate id
    raw_df$rep <- seq_len(nrow(raw_df))
    
    # put rep first
    raw_df <- raw_df[, c("rep", setdiff(names(raw_df), "rep"))]
    
    raw_path <- file.path(raw_dir, sprintf("row_%03d.rds", i))
    saveRDS(raw_df, raw_path, compress = "xz")
    n_saved <- n_saved + 1L
  }
  
  message("Saved ", n_saved, " raw replicate file(s) into: ", raw_dir)
}


## -------- 6. plots (same style as your code) --------------------------
results <- readRDS(here("lnM_summary_SAFEfun_2025-12-21.rds"))

results <- results %>%
  mutate(facet_label = paste0(design, " n1=", n1, " n2=", n2))

facet_info <- results %>% distinct(facet_label, design, n1, n2)

balanced_ind <- facet_info %>%
  filter(design == "indep", n1 == n2) %>% arrange(n1) %>% pull(facet_label)
unbalanced_ind <- facet_info %>%
  filter(design == "indep", n1 != n2) %>% arrange(n1) %>% pull(facet_label)
paired_balanced <- facet_info %>%
  filter(design == "paired") %>% arrange(n1) %>% pull(facet_label)

new_levels <- c(balanced_ind, unbalanced_ind, paired_balanced)
results <- results %>% mutate(facet_label = factor(facet_label, levels = new_levels))

## Bias plot
bias_df <- bind_rows(
  results %>% select(theta, facet_label, delta_bias) %>%
    rename(bias = delta_bias) %>% mutate(estimator = "delta"),
  results %>% select(theta, facet_label, safe_bias) %>%
    rename(bias = safe_bias)  %>% mutate(estimator = "SAFE")
)

p_bias <- ggplot(bias_df, aes(theta, bias, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
  labs(x = "theta", y = "Bias (estimate - true lnM)", colour = "Estimator") +
  scale_colour_manual(values = c(delta = "firebrick", SAFE = "steelblue"),
                      labels = c(delta = "PI", SAFE = "SAFE-BC"))
print(p_bias)

## Relative-bias of variance plot (baseline = Var_MC_SAFEpt)
rb_df <- bind_rows(
  results %>% select(theta, facet_label, relbias_delta) %>%
    rename(relbias = relbias_delta) %>% mutate(estimator = "delta"),
  results %>% select(theta, facet_label, relbias_safe) %>%
    rename(relbias = relbias_safe)  %>% mutate(estimator = "SAFE")
)

p_relbias <- ggplot(rb_df, aes(theta, relbias, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3) +
  labs(x = "theta", y = "Relative bias of Var (%) vs MC(SAFE-BC point)", colour = "Estimator") +
  scale_colour_manual(values = c(delta = "firebrick", SAFE = "steelblue"),
                      labels = c(delta = "PI", SAFE = "SAFE-var"))
print(p_relbias)

## Coverage
cov_df <- rbind(
  data.frame(results[,c("theta","design","n1","n2")], estimator="delta", cover=results$cover_delta),
  data.frame(results[,c("theta","design","n1","n2")], estimator="SAFE",  cover=results$cover_safe)
)

p_cov <- ggplot(cov_df, aes(theta, cover, colour=estimator, group=estimator)) +
  geom_hline(yintercept=0.95, linetype=2, colour="grey50") +
  geom_line() +
  facet_wrap(~ paste0(design," n1=",n1," n2=",n2), ncol=4) +
  labs(x=expression(theta), y="Empirical coverage") +
  scale_colour_manual(values=c(delta="firebrick", SAFE="steelblue"),
                      labels=c(delta="PI", SAFE="SAFE-BC"))
print(p_cov)

## RMSE
rmse_df <- rbind(
  data.frame(results[,c("theta","design","n1","n2")], estimator="delta", rmse=results$rmse_delta),
  data.frame(results[,c("theta","design","n1","n2")], estimator="SAFE",  rmse=results$rmse_safe)
)

p_rmse <- ggplot(rmse_df, aes(theta, rmse, colour=estimator, group=estimator)) +
  geom_line() +
  facet_wrap(~ paste0(design," n1=",n1," n2=",n2), ncol=4) +
  scale_colour_manual(values=c(delta="firebrick", SAFE="steelblue"),
                      labels=c(delta="PI", SAFE="SAFE-BC")) +
  labs(x=expression(theta), y="RMSE ( ln Mhat - ln M )", colour=NULL) +
  theme_bw(11)
print(p_rmse)

## Quick diagnostics
message("Mean abs(delta_bias): ", mean(abs(results$delta_bias), na.rm = TRUE))
message("Mean abs(safe_bias):  ", mean(abs(results$safe_bias),  na.rm = TRUE))
message("Mean abs(relbias_delta): ", mean(abs(results$relbias_delta), na.rm = TRUE))
message("Mean abs(relbias_safe):  ", mean(abs(results$relbias_safe),  na.rm = TRUE))
message("Mean delta_cap_rate: ", mean(results$delta_cap_rate, na.rm = TRUE))
message("Mean SAFE_ok_rate:   ", mean(results$SAFE_ok_rate,   na.rm = TRUE))
message("Mean boot_accept_prop (kept/tried): ", mean(results$boot_accept_prop, na.rm = TRUE))

########################################################################
## How to run (examples)
## Linux/macOS:
##   N_CORES=48 K_REPL=2000 MIN_KEPT=2000 CHUNK_INIT=4000 Rscript simulation16_SAFEfun_style.R
## Windows (cmd):
##   set N_CORES=48 & set K_REPL=2000 & set MIN_KEPT=2000 & set CHUNK_INIT=4000 & Rscript simulation16_SAFEfun_style.R
########################################################################
## ================================================================
## 7) Monte-Carlo error diagnostics from disk (raw_runs/*.rds)
##    (NO dependency on results_list)
## ================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(pbapply)
library(here)

## ---- 7a. Load summary results + reconstruct facet labels -------------
## (Use your own filename here)
summary_file <- here("lnM_summary_SAFEfun_2025-12-21.rds")
results <- readRDS(summary_file)

results <- results %>%
  mutate(facet_label = paste0(design, " n1=", n1, " n2=", n2))

## Keep your facet ordering logic (optional)
facet_info <- results %>% distinct(facet_label, design, n1, n2)

balanced_ind <- facet_info %>%
  filter(design == "indep", n1 == n2) %>% arrange(n1) %>% pull(facet_label)
unbalanced_ind <- facet_info %>%
  filter(design == "indep", n1 != n2) %>% arrange(n1) %>% pull(facet_label)
paired_balanced <- facet_info %>%
  filter(design == "paired") %>% arrange(n1) %>% pull(facet_label)

new_levels <- c(balanced_ind, unbalanced_ind, paired_balanced)
results <- results %>% mutate(facet_label = factor(facet_label, levels = new_levels))

## ---- 7b. Locate raw files -------------------------------------------
raw_dir <- here("raw_runs")
raw_files <- list.files(raw_dir, pattern = "^row_[0-9]{3}\\.rds$", full.names = TRUE)
if (!length(raw_files)) stop("No raw files found in: ", raw_dir)

## row index extracted from filename row_XXX.rds
row_index_from_path <- function(p) as.integer(sub("^row_([0-9]{3})\\.rds$", "\\1", basename(p)))

## ---- 7c. Helper: compute diagnostics from ONE raw_df -----------------
## raw_df format: runs as rows, stats as columns, includes rep column
## columns expected: delta_pt, delta_var, safe_pt, safe_var (at least)
one_raw_diag <- function(raw_df) {
  
  # defensive: allow rep or not
  if ("rep" %in% names(raw_df)) raw_df <- raw_df %>% select(-rep)
  
  # ok_d filter as in the main script
  ok_d <- is.finite(raw_df$delta_pt) & is.finite(raw_df$delta_var) & (raw_df$delta_var > 0)
  
  if (!any(ok_d)) {
    return(data.frame(
      n_ok_d = 0L,
      n_ok_safe = 0L,
      RSE_VarMC_SAFEpt = NA_real_,
      # MCSE of bias cannot be computed without true_lnM; filled later
      sd_PI = NA_real_,
      sd_SAFE = NA_real_
    ))
  }
  
  df_ok <- raw_df[ok_d, , drop = FALSE]
  
  # usable SAFE points within ok_d subset (matches your Var_MC_SAFEpt definition)
  ok_safe_pt <- is.finite(df_ok$safe_pt)
  n_ok_safe <- sum(ok_safe_pt)
  
  # RSE of sample variance (approx)
  RSE_varmc <- if (n_ok_safe >= 3) sqrt(2 / (n_ok_safe - 1)) else NA_real_
  
  # sd needed later to compute MCSE(bias) = sd(est - true)/sqrt(n_eff)
  # (true_lnM is constant, so sd(est-true)=sd(est))
  sd_PI   <- sd(df_ok$delta_pt, na.rm = TRUE)
  sd_SAFE <- sd(df_ok$safe_pt,  na.rm = TRUE)
  
  data.frame(
    n_ok_d = nrow(df_ok),
    n_ok_safe = n_ok_safe,
    RSE_VarMC_SAFEpt = RSE_varmc,
    sd_PI = sd_PI,
    sd_SAFE = sd_SAFE
  )
}

## ---- 7d. Compute diagnostics for ALL rows from disk ------------------
## We join by "row id" = file index, which should match param_grid row number
diag_list <- pbapply::pblapply(raw_files, function(f) {
  rid <- row_index_from_path(f)
  raw_df <- readRDS(f)
  d <- one_raw_diag(raw_df)
  d$row_id <- rid
  d
})

diag_df <- bind_rows(diag_list)

## Join diagnostics onto summary results
## IMPORTANT: this assumes row_id == row number in param_grid used when saving files.
## If you saved with the same i used in runner loop, this is correct.
results2 <- results %>%
  mutate(row_id = seq_len(nrow(results))) %>%
  left_join(diag_df, by = "row_id")

## ---- 7e. Recompute MCSE(bias) from disk (recommended) ----------------
## MCSE(bias) = sd(estimator)/sqrt(n_eff)
## n_eff should match the subset used for bias curves in your summaries.
## Here we use:
##  - PI: n_ok_d
##  - SAFE: n_ok_safe (within ok_d)
results2 <- results2 %>%
  mutate(
    mcse_bias_PI_disk   = ifelse(is.finite(sd_PI)   & n_ok_d >= 2,   sd_PI   / sqrt(n_ok_d),   NA_real_),
    mcse_bias_SAFE_disk = ifelse(is.finite(sd_SAFE) & n_ok_safe >= 2, sd_SAFE / sqrt(n_ok_safe), NA_real_)
  )

## If you prefer to keep the MCSE computed during the run, keep:
##   mcse_bias_delta and mcse_bias_safe
## But for “from disk only”, use *_disk below.

## ---- 7f. Plot diagnostics (like your figure) -------------------------
mc_long <- results2 %>%
  transmute(
    theta, facet_label,
    mcse_bias_PI = mcse_bias_PI_disk,
    mcse_bias_SAFE = mcse_bias_SAFE_disk,
    RSE_VarMC_SAFEpt = RSE_VarMC_SAFEpt
  ) %>%
  pivot_longer(cols = c(mcse_bias_PI, mcse_bias_SAFE, RSE_VarMC_SAFEpt),
               names_to = "metric", values_to = "value")

p_mc <- ggplot(mc_long, aes(theta, value, colour = metric, group = metric)) +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4) +
  labs(x = expression(theta),
       y = "MC error metric",
       title = "Monte-Carlo error diagnostics (per grid cell; computed from raw_runs)") +
  theme(legend.position = "bottom")
print(p_mc)

## ---- 7g. “Worst offenders” tables -----------------------------------
worst_mcse <- results2 %>%
  mutate(worst_mcse = pmax(abs(mcse_bias_PI_disk), abs(mcse_bias_SAFE_disk), na.rm = TRUE)) %>%
  arrange(desc(worst_mcse)) %>%
  select(row_id, design, n1, n2, theta, facet_label,
         n_ok_d, n_ok_safe,
         mcse_bias_PI_disk, mcse_bias_SAFE_disk, RSE_VarMC_SAFEpt) %>%
  head(20)

print(worst_mcse, row.names = FALSE)

worst_rse <- results2 %>%
  arrange(desc(RSE_VarMC_SAFEpt)) %>%
  select(row_id, design, n1, n2, theta, facet_label,
         n_ok_safe, RSE_VarMC_SAFEpt) %>%
  head(20)

print(worst_rse, row.names = FALSE)

## ================================================================
## Plot MCSE of the *average* variance estimates:
##   - mcse_varbar_delta : MCSE of mean(delta_var)
##   - mcse_varbar_safe  : MCSE of mean(safe_var)
## Assumes you already have: results (or results2) + facet_label
## ================================================================

var_mcse_df <- results %>%
  transmute(theta, facet_label,
            PI_var_MCSE   = mcse_varbar_delta,
            SAFE_var_MCSE = mcse_varbar_safe) %>%
  pivot_longer(cols = c(PI_var_MCSE, SAFE_var_MCSE),
               names_to = "estimator", values_to = "mcse") %>%
  mutate(estimator = recode(estimator,
                            PI_var_MCSE   = "PI (delta-var)",
                            SAFE_var_MCSE = "SAFE (safe-var)"))

p_var_mcse <- ggplot(var_mcse_df, aes(theta, mcse, colour = estimator, group = estimator)) +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4) +
  labs(x = expression(theta),
       y = "MCSE of mean(variance estimate)",
       colour = NULL,
       title = "Monte-Carlo error of variance estimators (MCSE of the mean)") +
  theme(legend.position = "bottom")

print(p_var_mcse)

## Optional: same plot on log scale (useful if MCSE spans orders of magnitude)
p_var_mcse_log <- p_var_mcse + scale_y_log10()
print(p_var_mcse_log)

## ---- 7h. Optional: convergence trace from ONE raw file (robust) ------
TRACE_ROW <- as.integer(Sys.getenv("TRACE_ROW", "1"))
trace_path <- file.path(raw_dir, sprintf("row_%03d.rds", TRACE_ROW))

if (!file.exists(trace_path)) {
  message("Trace file not found: ", trace_path)
} else {
  
  raw_df <- readRDS(trace_path)
  
  # baseline ok filter (as in main script) for PI-related validity
  ok_d <- is.finite(raw_df$delta_pt) & is.finite(raw_df$delta_var) & (raw_df$delta_var > 0)
  
  # separate usable streams
  df_PI   <- raw_df[ok_d & is.finite(raw_df$delta_pt), , drop = FALSE]
  df_SAFE <- raw_df[ok_d & is.finite(raw_df$safe_pt),  , drop = FALSE]
  
  message(sprintf("TRACE_ROW=%d: ok_d=%d, PI usable=%d, SAFE usable=%d",
                  TRACE_ROW, sum(ok_d), nrow(df_PI), nrow(df_SAFE)))
  
  run_mean <- function(x) cumsum(x) / seq_along(x)
  
  true_ln <- results2$true_lnM[results2$row_id == TRACE_ROW][1]
  
  # If true lnM is undefined for this grid cell, plot running means (not bias)
  do_bias <- is.finite(true_ln)
  
  tr_df <- bind_rows(
    data.frame(
      k = seq_len(nrow(df_PI)),
      estimator = "PI",
      value = run_mean(df_PI$delta_pt) - if (do_bias) true_ln else 0
    ),
    data.frame(
      k = seq_len(nrow(df_SAFE)),
      estimator = "SAFE-BC",
      value = run_mean(df_SAFE$safe_pt) - if (do_bias) true_ln else 0
    )
  )
  
  p_trace <- ggplot(tr_df, aes(k, value, colour = estimator)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_line() +
    labs(
      x = "Usable replicates (after filters)",
      y = if (do_bias) "Running bias" else "Running mean (centered at 0)",
      title = if (do_bias)
        sprintf("Convergence trace (row %03d): running bias", TRACE_ROW)
      else
        sprintf("Convergence trace (row %03d): true lnM is NA, showing running means", TRACE_ROW),
      colour = NULL
    )
  
  print(p_trace)
}