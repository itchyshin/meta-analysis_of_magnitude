########################################################################
## |d| simulation, summary  (naive absolute value vs folded-normal plug-in)
## -----------------------------------------------------------------------
## What this script does
##   This is a stand-alone simulation driver parallel in spirit to the
##   attached `sim_study.R` file used for lnM, but here the focus is on the
##   absolute standardized mean difference, |d|.
##
##   In other words, this script asks:
##     “If one wants to meta-analyse the magnitude of a standardized mean
##      difference by taking its absolute value, how well does that behave?”
##
##   The script compares two approaches:
##
##   (1) Naive |d| approach
##       - point estimate:  |d_hat|
##       - variance used:   Var(d_hat)
##         (that is, the usual sampling-variance formula for the signed SMD,
##          inherited unchanged after taking the absolute value)
##
##   (2) Folded-normal (FN) plug-in approach
##       - treat d_hat as approximately Normal(mu, sigma^2)
##       - then |d_hat| is approximately Folded Normal
##       - point estimate:  E(|X|) under the FN with mu = d_hat,
##                          sigma^2 = Var(d_hat)
##       - variance used:   Var(|X|) under that same FN plug-in
##
## Design mirrors sim_study.R
##   As in the attached lnM script, the design includes both:
##     - independent two-group comparisons (unpaired)
##     - paired (dependent) two-condition comparisons with correlation rho
##
## Main outputs per parameter set
##   For each combination of theta, sample sizes, and design, we record:
##     - bias and RMSE of point estimators relative to true |d|
##     - mean reported variances
##     - relative bias of variance estimators against Monte Carlo variance
##     - coverage of nominal 95% Wald intervals
##     - failure rates
##
## Important conceptual notes
##   - Folded-normal calculations rely on the approximation
##         d_hat ~ Normal(mu, sigma^2)
##     which is only approximate, especially in small samples.
##
##   - The "naive" approach is included because this is the common practice:
##       compute |d_hat| and then analyse it as though it were just another
##       Gaussian outcome with the usual SMD variance.
##
##   - For paired data, we use a Cohen-type repeated-measures d based on the
##     average / pooled within-condition SD. This choice keeps the population
##     target aligned with the theta grid used in the lnM simulation when
##     sigma1 = sigma2 = 1.
##
## Relationship to the attached sim_study.R
##   The overall structure, sectioning, helper philosophy, and summary logic
##   are deliberately kept parallel to the attached lnM simulation file, so
##   that the two scripts are easy to compare side by side.
########################################################################

library(MASS)        # mvrnorm() for paired-data simulation
library(parallel)    # multicore / cluster helpers for outer-loop parallelism
library(pbapply)     # progress bars for parallel *apply

## -------- 0. globals & helpers ---------------------------------------
##
## This section contains small utility functions used throughout the script.
## These mirror the kind of “infrastructure helpers” used in sim_study.R.

posify <- function(x, eps = 1e-12) pmax(x, eps)
# posify() is a numerical safety helper.
# In theory, quantities such as variances should be non-negative, but in
# finite-precision arithmetic tiny negative values can occasionally arise
# from round-off. We therefore truncate to a very small positive number.

mcse_mean <- function(x) {
  # Monte Carlo standard error (MCSE) of a sample mean.
  # This quantifies simulation uncertainty in estimated means/biases.
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sqrt(stats::var(x) / length(x))
}

mcse_prop <- function(x) {
  # Monte Carlo standard error (MCSE) of a sample proportion.
  # This is used for quantities such as coverage.
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  p <- mean(x)
  sqrt(p * (1 - p) / length(x))
}

## -------- 1. d and variance formulas ---------------------------------
##
## This section defines the signed standardized mean difference estimators
## and the corresponding large-sample sampling-variance formulas used in
## the simulation.
##
## Independent groups:
##   d = (x1bar - x2bar) / s_pooled
##   Var(d) ≈ (n1+n2)/(n1*n2) + d^2 / {2(n1+n2-2)}
##
## Paired groups:
##   We use a repeated-measures / paired d with the average within-condition
##   SD in the denominator:
##     d_av = (x1bar - x2bar) / sqrt((s1^2 + s2^2)/2)
##
##   Approximate variance:
##     Var(d_av) ≈ 2(1-rho)/n + d_av^2 / {2(n-1)}
##
## Why this paired denominator?
##   The goal here is not to explore every possible paired-SMD definition,
##   but to keep the target aligned with the lnM simulation when sigma1=sigma2=1.
##   Under that setup, the population denominator is 1, so true signed d = theta.

d_ind <- function(x1bar, x2bar, s1, s2, n1, n2) {
  # Independent-groups Cohen's d:
  #   d = (mean difference) / pooled within-group SD
  
  sp2 <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  # sp2 = pooled variance estimate
  
  sp  <- sqrt(posify(sp2))
  # pooled SD, guarded by posify() for numerical stability
  
  d   <- (x1bar - x2bar) / sp
  # signed standardized mean difference
  
  vd  <- (n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2 - 2))
  # common large-sample approximation to Var(d)
  
  c(d = d, vd = posify(vd))
}

d_dep <- function(x1bar, x2bar, s1, s2, n, rho) {
  # Paired / repeated-measures d using the average within-condition SD:
  #   d = (mean difference) / sqrt((s1^2 + s2^2)/2)
  
  sav2 <- (s1^2 + s2^2) / 2
  # average within-condition variance
  
  sav  <- sqrt(posify(sav2))
  # average within-condition SD
  
  d    <- (x1bar - x2bar) / sav
  # signed standardized mean difference for paired data
  
  vd   <- 2 * (1 - rho) / n + d^2 / (2 * (n - 1))
  # approximate variance formula that decreases as within-pair correlation
  # increases, since highly correlated pairs provide more information
  
  c(d = d, vd = posify(vd))
}

## -------- 2. Folded-normal functions ---------------------------------
##
## If X ~ Normal(mu, sigma^2), then Y = |X| follows a Folded Normal
## distribution. This section implements the mean and variance formulas of
## that distribution.
##
## Specifically:
##
##   E(Y) = sigma * sqrt(2/pi) * exp(-mu^2 / (2 sigma^2))
##          + |mu| * {2 Phi(|mu| / sigma) - 1}
##
##   Var(Y) = mu^2 + sigma^2 - [E(Y)]^2
##
## We use these in two ways:
##   (i)  plug-in estimator based on observed d_hat and v_d_hat
##   (ii) theoretical benchmark based on true theta and true v_d
##
## The folded-normal idea is important here because taking an absolute value
## changes the sampling distribution even if the signed estimator is roughly
## Normal.

fn_mean <- function(mu, var) {
  # Mean of |X| when X ~ Normal(mu, var)
  
  var <- posify(var)
  sig <- sqrt(var)
  a   <- abs(mu) / sig
  
  sig * sqrt(2 / pi) * exp(-0.5 * a^2) +
    abs(mu) * (2 * pnorm(a) - 1)
}

fn_var <- function(mu, var) {
  # Variance of |X| when X ~ Normal(mu, var)
  
  var <- posify(var)
  m   <- fn_mean(mu, var)
  posify(mu^2 + var - m^2)
}

## -------- 3. True target and theoretical FN benchmark ----------------
##
## In this simulation design, sigma1 = sigma2 = 1, so true signed d = theta.
## Therefore the target magnitude is:
##   true_absd = |theta|
##
## We also compute the theoretical folded-normal mean/variance implied by the
## large-sample Normal approximation to signed d_hat.
##
## These “true” variance formulas are not used to construct the estimators.
## Instead, they are used as external benchmarks for comparison.

vd_true_ind <- function(theta, n1, n2) {
  # Theoretical large-sample variance for signed d in the independent design
  (n1 + n2) / (n1 * n2) + theta^2 / (2 * (n1 + n2 - 2))
}

vd_true_dep <- function(theta, n, rho) {
  # Theoretical large-sample variance for signed d in the paired design
  2 * (1 - rho) / n + theta^2 / (2 * (n - 1))
}

## -------- 4. One replicate -------------------------------------------
##
## This function simulates ONE dataset (one replicate) for either the
## independent or paired design, then computes:
##
##   - signed d and its variance
##   - naive |d| and inherited variance
##   - folded-normal plug-in point estimate and variance
##
## It returns a named vector so the outer simulation loop can collect all
## replicates into a matrix.

one_rep_absd <- function(mu1, mu2, sd1, sd2,
                         n1, n2 = NULL, rho = 0) {
  
  if (is.null(n2)) {   # paired
    # If n2 is NULL, we interpret the design as paired / dependent.
    # We generate bivariate Normal observations with correlation rho.
    
    Sigma <- matrix(c(sd1^2, rho * sd1 * sd2,
                      rho * sd1 * sd2, sd2^2), 2)
    # 2x2 covariance matrix for the paired outcomes
    
    xy <- mvrnorm(n1, c(mu1, mu2), Sigma)
    x1 <- xy[, 1]
    x2 <- xy[, 2]
    
    rho_hat <- suppressWarnings(cor(x1, x2))
    # use the sample within-pair correlation in the variance formula,
    # matching the idea that the analyst would usually only observe the sample
    
    if (!is.finite(rho_hat)) rho_hat <- rho
    # guard against degenerate rare cases
    
  } else {             # independent
    # If n2 is provided, we simulate two independent groups.
    x1 <- rnorm(n1, mu1, sd1)
    x2 <- rnorm(n2, mu2, sd2)
    rho_hat <- 0
  }
  
  # Sample summary statistics used by the estimators
  x1bar <- mean(x1)
  x2bar <- mean(x2)
  s1    <- sd(x1)
  s2    <- sd(x2)
  
  # Compute signed d and its approximate variance
  est <- if (is.null(n2)) {
    d_dep(x1bar, x2bar, s1, s2, n1, rho_hat)
  } else {
    d_ind(x1bar, x2bar, s1, s2, n1, n2)
  }
  
  d_hat  <- unname(est["d"])
  vd_hat <- unname(est["vd"])
  
  ## Naive absolute-value approach
  absd_pt  <- abs(d_hat)
  absd_var <- vd_hat
  # The point estimate is |d_hat|, but the variance is simply inherited from
  # the signed d formula. This is exactly the “naive” practice under study.
  
  ## Folded-normal plug-in approach
  fn_pt  <- fn_mean(d_hat, vd_hat)
  fn_v   <- fn_var(d_hat, vd_hat)
  # Here we treat the signed d_hat as approximately Normal(d_hat, vd_hat)
  # and then transform to the corresponding folded-normal mean/variance.
  
  c(
    d_pt      = d_hat,
    d_var     = vd_hat,
    absd_pt   = absd_pt,
    absd_var  = absd_var,
    fn_pt     = fn_pt,
    fn_var    = fn_v,
    fail      = 0
  )
}

## -------- 5. parameter grid ------------------------------------------
##
## This section defines the simulation design.
## It deliberately mirrors the design used in the attached sim_study.R file,
## so that the |d| results and lnM results are directly comparable.

theta_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4, 5)
# theta is the true signed standardized mean difference because sigma1=sigma2=1

pairs_ind <- data.frame(
  n1 = c(5, 10, 20, 100, 3, 6, 12, 40),
  n2 = c(5, 10, 20, 100, 7, 14, 28, 160)
)
# independent-group sample-size combinations:
# balanced (5,5), (10,10), ...
# and unbalanced (3,7), (6,14), ...

grid_ind <- expand.grid(
  theta = theta_vals,
  idx   = seq_len(nrow(pairs_ind))
)
grid_ind$n1     <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2     <- pairs_ind$n2[grid_ind$idx]
grid_ind$design <- "indep"
grid_ind$idx    <- NULL

grid_dep <- expand.grid(
  theta = theta_vals,
  n     = c(5, 10, 20, 100)
)
grid_dep$n1     <- grid_dep$n
grid_dep$n2     <- grid_dep$n
grid_dep$n      <- NULL
grid_dep$design <- "paired"

param_grid <- rbind(grid_ind, grid_dep)
# final parameter grid: one row per combination of theta, design, and sample size

## -------- 6. simulation driver ---------------------------------------
##
## These are the run-time settings and the main per-row simulation function.
##
## As in sim_study.R:
##   - K_repl controls the number of Monte Carlo replicates per parameter set
##   - RHO_PAIRED is the population within-pair correlation for paired data
##   - inner_every controls progress messages within each row

set.seed(20260321)

K_repl        <- as.integer(Sys.getenv("K_REPL", "100000"))
RHO_PAIRED    <- as.numeric(Sys.getenv("RHO_PAIRED", "0.8"))
outer_verbose <- TRUE
inner_every   <- 100
# inner_every = print progress every 100 replicates within a parameter set

runner_absd <- function(i) {
  # Run the full Monte Carlo simulation for one row of param_grid
  
  p   <- param_grid[i, ]
  sd0 <- 1
  
  true_absd <- abs(p$theta)
  # target magnitude, since sigma1=sigma2=1
  
  true_vd <- if (p$design == "indep") {
    vd_true_ind(p$theta, p$n1, p$n2)
  } else {
    vd_true_dep(p$theta, p$n1, RHO_PAIRED)
  }
  
  ## Theoretical folded-normal moments under Normal approx to signed d_hat
  true_fn_mean <- fn_mean(p$theta, true_vd)
  true_fn_var  <- fn_var(p$theta, true_vd)
  # these are “benchmark” FN quantities under the population truth
  
  M <- matrix(
    NA_real_, 7, K_repl,
    dimnames = list(
      c("d_pt", "d_var", "absd_pt", "absd_var", "fn_pt", "fn_var", "fail"),
      NULL
    )
  )
  # replicate-level storage for this parameter set
  
  for (k in seq_len(K_repl)) {
    M[, k] <- if (p$design == "indep") {
      one_rep_absd(0, p$theta, sd0, sd0, n1 = p$n1, n2 = p$n2, rho = 0)
    } else {
      one_rep_absd(0, p$theta, sd0, sd0, n1 = p$n1, n2 = NULL, rho = RHO_PAIRED)
    }
    
    if (k %% inner_every == 0) {
      message(sprintf("[row %3d] replicate %7d / %d done", i, k, K_repl))
    }
  }
  
  ## Monte Carlo variances of the point estimators
  Var_MC_absd <- if (sum(is.finite(M["absd_pt", ])) >= 2) {
    var(M["absd_pt", ], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  Var_MC_fn <- if (sum(is.finite(M["fn_pt", ])) >= 2) {
    var(M["fn_pt", ], na.rm = TRUE)
  } else {
    NA_real_
  }
  # These Monte Carlo variances are treated as the empirical “truth”
  # for evaluating the variance estimators.
  
  ## 95% Wald coverage against true |d| target
  cover_absd <- abs(M["absd_pt", ] - true_absd) <= 1.96 * sqrt(M["absd_var", ])
  cover_fn   <- abs(M["fn_pt",   ] - true_absd) <= 1.96 * sqrt(M["fn_var",   ])
  # Wald intervals are constructed on the transformed scale for each estimator
  
  out <- data.frame(
    theta   = p$theta,
    design  = p$design,
    n1      = p$n1,
    n2      = ifelse(p$design == "indep", p$n2, p$n1),
    rho     = ifelse(p$design == "paired", RHO_PAIRED, 0),
    
    ## Targets / theoretical benchmarks
    true_absd    = true_absd,
    true_vd      = true_vd,
    true_fn_mean = true_fn_mean,
    true_fn_var  = true_fn_var,
    
    ## Mean point estimates
    d_mean      = mean(M["d_pt", ], na.rm = TRUE),
    absd_mean   = mean(M["absd_pt", ], na.rm = TRUE),
    fn_mean_est = mean(M["fn_pt", ], na.rm = TRUE),
    
    ## Bias relative to the target |d_true| = |theta|
    absd_bias = mean(M["absd_pt", ], na.rm = TRUE) - true_absd,
    fn_bias   = mean(M["fn_pt",   ], na.rm = TRUE) - true_absd,
    
    ## Also compare MC means to the theoretical folded-normal benchmark
    absd_vs_fn_mean = mean(M["absd_pt", ], na.rm = TRUE) - true_fn_mean,
    fn_vs_fn_mean   = mean(M["fn_pt",   ], na.rm = TRUE) - true_fn_mean,
    
    ## Mean reported variances
    mean_var_absd = mean(M["absd_var", ], na.rm = TRUE),
    mean_var_fn   = mean(M["fn_var",   ], na.rm = TRUE),
    
    ## MC variances of the point estimators
    Var_MC_absd = Var_MC_absd,
    Var_MC_fn   = Var_MC_fn,
    
    ## Relative bias of the variance estimators vs their own MC variance
    relbias_var_absd = if (is.finite(Var_MC_absd) && Var_MC_absd > 0) {
      100 * (mean(M["absd_var", ], na.rm = TRUE) / Var_MC_absd - 1)
    } else {
      NA_real_
    },
    
    relbias_var_fn = if (is.finite(Var_MC_fn) && Var_MC_fn > 0) {
      100 * (mean(M["fn_var", ], na.rm = TRUE) / Var_MC_fn - 1)
    } else {
      NA_real_
    },
    
    ## Extra diagnostic: compare reported variances with theoretical FN variance
    relbias_var_absd_vsFN = if (is.finite(true_fn_var) && true_fn_var > 0) {
      100 * (mean(M["absd_var", ], na.rm = TRUE) / true_fn_var - 1)
    } else {
      NA_real_
    },
    
    relbias_var_fn_vsFN = if (is.finite(true_fn_var) && true_fn_var > 0) {
      100 * (mean(M["fn_var", ], na.rm = TRUE) / true_fn_var - 1)
    } else {
      NA_real_
    },
    
    ## RMSE relative to true |d|
    rmse_absd = sqrt(mean((M["absd_pt", ] - true_absd)^2, na.rm = TRUE)),
    rmse_fn   = sqrt(mean((M["fn_pt",   ] - true_absd)^2, na.rm = TRUE)),
    
    ## Coverage
    cover_absd = mean(cover_absd, na.rm = TRUE),
    cover_fn   = mean(cover_fn,   na.rm = TRUE),
    
    ## Failures
    fail_prop = mean(M["fail", ], na.rm = TRUE),
    
    ## MCSEs
    mcse_bias_absd   = mcse_mean(M["absd_pt", ] - true_absd),
    mcse_bias_fn     = mcse_mean(M["fn_pt",   ] - true_absd),
    mcse_varbar_absd = mcse_mean(M["absd_var", ]),
    mcse_varbar_fn   = mcse_mean(M["fn_var",   ]),
    mcse_cover_absd  = mcse_prop(cover_absd),
    mcse_cover_fn    = mcse_prop(cover_fn)
  )
  
  attr(out, "raw_M") <- M
  # store raw replicate-level matrix as an attribute so it can optionally be
  # saved later, matching the philosophy of sim_study.R
  
  if (outer_verbose) {
    message(sprintf(
      "Finished row %3d  (%s  theta=%s  n1=%d n2=%d)",
      i, p$design, p$theta, p$n1,
      ifelse(p$design == "indep", p$n2, p$n1)
    ))
  }
  
  out
}

## -------- 7. PARALLEL outer loop -------------------------------------
##
## As in sim_study.R, we parallelise over parameter sets (outer loop),
## not over replicates within a parameter set. This keeps each row self-
## contained and makes progress reporting easier to interpret.

n_cores_use <- suppressWarnings(as.integer(Sys.getenv("N_CORES", "")))
if (is.na(n_cores_use) || n_cores_use <= 0) {
  n_cores_use <- max(1L, detectCores() / 2)
}
# default to roughly half the available cores if not specified

pbop <- pbapply::pboptions(type = "txt")

if (.Platform$OS.type == "windows") {
  # On Windows, pbapply uses an explicit PSOCK cluster
  cl <- makeCluster(n_cores_use)
  clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)
  results_list <- pbapply::pblapply(seq_len(nrow(param_grid)), runner_absd, cl = cl)
  stopCluster(cl)
} else {
  # On Unix-like systems, pbapply can use multicore directly
  results_list <- pbapply::pblapply(seq_len(nrow(param_grid)), runner_absd, cl = n_cores_use)
}

pbapply::pboptions(pbop)

## -------- 8. COLLAPSE -------------------------------------------------
##
## Convert the list of per-row data frames into one combined summary table.

is_df <- sapply(results_list, is.data.frame)
results <- do.call(rbind, results_list[is_df])

results$theta <- as.numeric(as.character(results$theta))
results$n1    <- as.numeric(as.character(results$n1))
results$n2    <- as.numeric(as.character(results$n2))

if (any(!is_df)) {
  message(sprintf("Warning: %d parameter sets failed and were excluded.", sum(!is_df)))
  message("Failed indices: ", paste(which(!is_df), collapse = ", "))
}

## -------- 9. SAVE OUTPUTS --------------------------------------------
##
## Save the combined summary table in both RDS and CSV form.
## Optionally also save raw replicate-level outputs per parameter row.

summary_rds <- sprintf("absd_summary_%s.rds", Sys.Date())
saveRDS(results, summary_rds)
message("Saved overall summary (RDS) to: ", summary_rds)

summary_csv <- sprintf("absd_summary_%s.csv", Sys.Date())
write.csv(results, file = summary_csv, row.names = FALSE)
message("Saved overall summary (CSV) to: ", summary_csv)

# saveRDS(
#   results_list,
#   file = sprintf("absd_results_list_%s.rds", Sys.Date()),
#   compress = "xz"
# )
# Optional: save the full list object directly if desired.

save_raw <- TRUE

if (isTRUE(save_raw)) {
  raw_dir <- "raw_runs_absd"
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  
  n_saved <- 0L
  for (i in seq_along(results_list)) {
    raw_M <- attr(results_list[[i]], "raw_M")
    if (is.null(raw_M)) next
    
    raw_df <- as.data.frame(t(raw_M))
    raw_df$rep <- seq_len(nrow(raw_df))
    raw_df <- raw_df[, c("rep", setdiff(names(raw_df), "rep"))]
    
    raw_path <- file.path(raw_dir, sprintf("row_%03d.rds", i))
    saveRDS(raw_df, raw_path, compress = "xz")
    n_saved <- n_saved + 1L
  }
  
  message("Saved ", n_saved, " raw replicate file(s) into: ", raw_dir)
}