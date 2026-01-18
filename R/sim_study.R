########################################################################
## lnM – simulation, summary  (Delta-method vs SAFE-BC)
## --------------------------------------------------------------------
## What this script does
##   This is the main simulation driver used in the lnM manuscript to
##   compare:
##     (1) Delta-method plug-in estimator (“Delta-1”) for lnM + its variance
##     (2) SAFE parametric bootstrap estimator for lnM, with bias-correction
##         (“SAFE-BC”) and its variance estimated from accepted bootstrap draws
##
##   The design includes both:
##     - independent two-group comparisons (unpaired)
##     - paired (dependent) two-condition comparisons with correlation rho
##
## Key idea: lnM is only defined when a “gap” is positive:
##   Delta = MSB - MSW  must be > 0
## where
##   MSB = between-group mean square  (signal)
##   MSW = within-group mean square   (noise)
## If MSB <= MSW then Delta <= 0 and lnM is undefined (log of non-positive).
## We treat those cases as NA (for both Delta-method and SAFE).
##
## SAFE parametric bootstrap in this project
##   SAFE draws replicate (m1*, m2*, v1*, v2*) from the finite-sample sampling
##   distributions implied by observed summary statistics. It then computes
##   lnM* and KEEPS ONLY accepted draws with MSB* > MSW*.
##
##   Importantly, this version targets a minimum number of ACCEPTED draws
##   (min_kept), using adaptive chunking (see SAFE_fun.R). The total number of
##   attempted draws varies depending on acceptance probability.
##
## Bias correction used (SAFE-BC)
##   SAFE returns the mean of accepted lnM* draws (call it mean(lnM*)).
##   When Delta-method point estimate exists (PI), we apply:
##       lnM_SAFE-BC = 2 * lnM_PI - mean(lnM*)
##   If PI is undefined (NA), we fall back to the SAFE mean (no BC).
##
## Performance metrics recorded per parameter set include:
##   - bias and RMSE of point estimators
##   - mean reported variances (Delta-method variance, SAFE variance)
##   - coverage of nominal 95% intervals (when true lnM is defined)
##   - failure rates (when estimators are undefined)
##   - SAFE acceptance rates (accepted / attempted)
##   - rate at which Delta-method variance exceeds a cap (maxVar)
##
## Baseline for variance “relative bias”
##   We treat the Monte Carlo variance of the SAFE-BC point estimator as the
##   reference (“truth”) for variance evaluation:
##     Var_MC_SAFEpt = var( lnM_SAFE-BC across replicates )
##   and compute relative bias (%) of each variance estimator against it.
##
## Reproducibility and runtime notes
##   - outer loop is parallelised (pbapply + multicore/cluster)
##   - verbose progress inside each parameter set (every inner_every reps)
##   - results saved as RDS+CSV; optionally raw replicate-level outputs saved
##
## This version uses SAFE_fun.R (safe_lnM_indep / safe_lnM_dep).
########################################################################

library(MASS)        # mvrnorm()
library(ggplot2)     # plotting
library(dplyr)       # facet ordering
library(parallel)    # multicore helpers
library(pbapply)     # progress bars for *apply
library(here)        # file paths relative to this script

## -------- 0. globals & helpers ---------------------------------------
maxVar   <- 20
# maxVar is a pragmatic cap on Delta-method variance to avoid a tiny number of
# extreme values dominating summaries/plots when Delta is near 0 and gradients
# explode. We track how often capping occurs (delta_cap_rate / delta_cap_n).

posify   <- function(x, eps = 1e-12) pmax(x, eps)
# posify() is a “floating point safety” helper: occasionally expressions that
# should be non-negative can become slightly negative (e.g., -1e-16) due to
# numerical round-off. posify() truncates to a small positive epsilon.

safe_gap <- function(g) ifelse(g <= 0, NA_real_, g)
# lnM is only defined when Delta = MSB - MSW > 0. We use safe_gap() to map
# non-positive values to NA, consistently across methods.

## lnM core transformation
lnM_core <- function(Delta, MSW, n0)
  0.5 * (log(Delta) - log(n0) - log(MSW))
# This is the common lnM formula once we have Delta (>0), MSW (>0), and the
# design-specific scaling n0:
#   lnM = 0.5 * [ log(Delta) - log(n0) - log(MSW) ]
# n0 differs between independent and paired designs (see lnM_true_* and delta
# estimators below). For independent data, n0 = 2h with h = n1*n2/(n1+n2).

## -------- 0b. load SAFE_fun.R ----------------------------------------
SAFE_FILE <- here("R", "SAFE_fun.R")
# We keep SAFE bootstrap machinery in a separate file to:
#   (i)   keep this script focused on simulation design + summarisation
#   (ii)  allow SAFE to be unit-tested independently
#   (iii) reduce code duplication across projects / manuscripts

if (!file.exists(SAFE_FILE)) stop("SAFE_fun.R not found at: ", SAFE_FILE)
source(SAFE_FILE)

# Sanity check: confirm the expected entry points exist after sourcing
if (!exists("safe_lnM_indep") || !exists("safe_lnM_dep")) {
  stop("SAFE_fun.R must define safe_lnM_indep() and safe_lnM_dep().")
}

## A small helper so we do not care about exact list field names
safe_get <- function(S, name, default = NA) {
  # SAFE_fun.R functions may return slightly different field names across
  # versions (e.g., 'point' vs 'pt', 'total' vs 'tried'). safe_get() provides
  # a robust interface so this simulation script does not break when field
  # names are updated.
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
# “True” lnM here means the deterministic value under the population model used
# in the simulation (Normal with sigma = 1, and mean difference theta).
# We use these functions for bias, RMSE, and coverage calculations.
#
# Important: true lnM itself may be undefined (NA) when MSB <= MSW at the
# population level. In those regions, it is not meaningful to compute coverage
# of confidence intervals around an undefined target. We therefore set coverage
# to NA when true_lnM is NA.

lnM_true_ind <- function(theta, n1, n2, sigma = 1) {
  # Independent design “population” ANOVA quantities:
  # MSB = h * theta^2, where h = n1*n2/(n1+n2)
  # MSW = sigma^2
  msb <- (n1 * n2) / (n1 + n2) * theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  # For independent design, n0 = 2h
  lnM_core(msb - msw, msw, 2 * n1 * n2 / (n1 + n2))
}

lnM_true_dep <- function(theta, n, sigma = 1) {
  # Paired design “population” quantities (simplified two-condition repeated
  # measures setup):
  # MSB = (n/2) * theta^2
  # MSW = sigma^2
  msb <- (n / 2) * theta^2
  msw <- sigma^2
  if (msb <= msw) return(NA_real_)
  # For paired design, n0 = n
  lnM_core(msb - msw, msw, n)
}


## -------- 1. Delta-method plug-in (indep & paired) --------------------
## Returns a named vector: (pt, var, capped)
##
## Delta-method estimators use observed sample summaries:
##   x1bar, x2bar, s1, s2
## and plug them into the lnM formula plus a first-order delta-method variance
## approximation.
##
## lnM is undefined when Delta <= 0; we return NA in that case.
## Variance capping:
##   We compute var_raw, then:
##     capped = 1 if var_raw > maxVar, else 0
##     var    = min(var_raw, maxVar)
##
## Note: capping is ONLY applied to the Delta-method variance estimate here.
## SAFE variance comes from the empirical variance of accepted lnM* draws and is
## not capped in this script.

lnM_delta_ind <- function(x1, x2, s1, s2, n1, n2) {
  # h enters MSB (between-group mean square)
  h   <- n1 * n2 / (n1 + n2)
  
  # Between- and within- mean squares based on sample summaries
  MSB <- h * (x1 - x2)^2
  MSW <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  
  # Gap must be positive
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA_real_, var = NA_real_, capped = NA_real_))
  
  # Point estimate: lnM_core with n0 = 2h for independent groups
  pt  <- lnM_core(Delta, MSW, 2 * h)
  
  # delta-method variance components:
  # sigmaD2 is Var(x1bar - x2bar) under independence
  sD2 <- s1^2 / n1 + s2^2 / n2
  dif <- x1 - x2
  
  # Approx Var(MSB) and Var(MSW) under finite-sample normal theory
  vB  <- h^2 * (2 * sD2^2 + 4 * sD2 * dif^2)
  vW  <- 2 * MSW^2 / (n1 + n2 - 2)
  
  # Gradients of lnM wrt MSB and MSW (chain rule through Delta and logs)
  g1  <- 0.5 / Delta
  g2  <- -0.5 * MSB / (Delta * MSW)
  
  # Combine; enforce positivity for numeric stability
  var_raw <- posify(g1^2 * vB + g2^2 * vW)
  
  # Apply variance cap
  capped <- as.integer(is.finite(var_raw) && var_raw > maxVar)
  c(pt = pt, var = pmin(var_raw, maxVar), capped = capped)
}

lnM_delta_dep <- function(x1, x2, s1, s2, n, rho) {
  # Paired design mean squares:
  # MSB uses h = n/2; MSW is average of the two within-condition variances.
  MSB <- (n / 2) * (x1 - x2)^2
  MSW <- (s1^2 + s2^2) / 2
  
  Delta <- safe_gap(MSB - MSW)
  if (is.na(Delta)) return(c(pt = NA_real_, var = NA_real_, capped = NA_real_))
  
  # Point estimate: paired uses n0 = n
  pt  <- lnM_core(Delta, MSW, n)
  
  # Variance of paired mean difference accounting for correlation rho
  sD2 <- s1^2 + s2^2 - 2 * rho * s1 * s2
  dif <- x1 - x2
  
  # Approx Var(MSB) and Var(MSW) for paired setting
  vB  <- (n / 2)^2 * (2 * sD2^2 / n^2 + 4 * dif^2 * sD2 / n)
  vW  <- (s1^4 + s2^4 + 2 * rho^2 * s1^2 * s2^2) / (2 * (n - 1))
  
  # Gradients
  g1  <- 0.5 / Delta
  g2  <- -0.5 * MSB / (Delta * MSW)
  
  var_raw <- posify(g1^2 * vB + g2^2 * vW)
  
  capped <- as.integer(is.finite(var_raw) && var_raw > maxVar)
  c(pt = pt, var = pmin(var_raw, maxVar), capped = capped)
}


## -------- 2. SAFE-BC via SAFE_fun.R ----------------------------------
## These wrappers standardise the SAFE output fields and keep this script
## independent of small changes in SAFE_fun.R list element names.
##
## SAFE_lnM_* functions implement:
##   - parametric draws of means and variances (or MVN + Wishart for paired)
##   - acceptance rule: keep draws with MSB* > MSW*
##   - adaptive chunking until min_kept accepted draws are reached (or fail)
##
## We return:
##   pt   : SAFE estimate of lnM (typically mean of accepted lnM*)
##   var  : variance of accepted lnM* draws (empirical)
##   kept : number of accepted draws kept
##   tried: total number of draws attempted (accepted + rejected)
##   status: "ok" if completed, or an informative failure status

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
## one_rep() generates one simulated dataset under a given parameter set, then
## computes:
##   - Delta-method point + variance (with variance capping)
##   - SAFE point + variance (then applies SAFE-BC correction)
## and returns a fixed vector of summary quantities (length 10).
##
## Inputs:
##   mu1, mu2      : population means (mu2 - mu1 = theta in the grid)
##   sd1, sd2      : population SDs (fixed at 1 in this study)
##   n1, n2        : sample sizes (n2 = NULL indicates paired design)
##   rho           : population correlation in paired design (0.8 in grid)
##
## SAFE tuning parameters:
##   min_kept      : number of accepted bootstrap draws we target
##   chunk_*       : adaptive chunking controls
##   max_draws     : hard cap on total attempted draws (safety)
##   patience_*    : stop if too many consecutive chunks have zero acceptance
##
## Important detail for paired simulations:
##   We *generate* data with population rho, but then estimate rho_hat from the
##   simulated sample and feed rho_hat into both delta and SAFE estimators.
##   This mirrors real applications, where we often only have an estimated within-
##   study correlation rather than the “true” population rho.

one_rep <- function(mu1, mu2, sd1, sd2,
                    n1, n2 = NULL, rho = 0,
                    min_kept = 2000, chunk_init = 4000,
                    chunk_max = 2e6, max_draws = Inf,
                    patience_noaccept = 5) {
  
  if (is.null(n2)) {                      # paired (dependent)
    # Construct the population covariance matrix
    Sigma <- matrix(c(sd1^2, rho * sd1 * sd2,
                      rho * sd1 * sd2, sd2^2), 2)
    
    # Generate paired observations (x1, x2) ~ MVN((mu1, mu2), Sigma)
    xy <- mvrnorm(n1, c(mu1, mu2), Sigma)
    x1 <- xy[, 1]; x2 <- xy[, 2]
    
    # Estimate the within-pair correlation from the sample
    rho_hat <- suppressWarnings(cor(x1, x2))
    # In edge cases (e.g., zero variance), cor() can be NA/NaN; fallback to rho
    if (!is.finite(rho_hat)) rho_hat <- rho
    
  } else {                                # independent groups
    # Generate independent samples from N(mu, sd)
    x1 <- rnorm(n1, mu1, sd1)
    x2 <- rnorm(n2, mu2, sd2)
    rho_hat <- 0
  }
  
  # Sample summaries (the estimators only use summaries)
  x1bar <- mean(x1); s1 <- sd(x1)
  x2bar <- mean(x2); s2 <- sd(x2)
  
  # ---- Delta-method plug-in estimator ----
  d <- if (is.null(n2))
    lnM_delta_dep(x1bar, x2bar, s1, s2, n1, rho_hat)
  else
    lnM_delta_ind(x1bar, x2bar, s1, s2, n1, n2)
  
  # ---- SAFE bootstrap estimator (mean of accepted lnM*) ----
  s <- if (is.null(n2))
    safe_dep_fun(x1bar, x2bar, s1, s2, n1, rho_hat,
                 min_kept, chunk_init, chunk_max, max_draws, patience_noaccept)
  else
    safe_ind_fun(x1bar, x2bar, s1, s2, n1, n2,
                 min_kept, chunk_init, chunk_max, max_draws, patience_noaccept)
  
  delta_pt  <- unname(d["pt"])
  safe_raw  <- s$pt
  
  ## SAFE-BC: bias-correct SAFE using the Delta-method plug-in when available
  ##   lnM_BC = 2 * lnM_PI - mean(lnM*)
  ## If Delta-method point estimate is NA (undefined), use SAFE raw mean.
  safe_bc <- if (is.na(delta_pt)) safe_raw else 2 * delta_pt - safe_raw
  
  # Return a consistent vector; these are later stored in matrix M
  c(delta_pt    = delta_pt,
    delta_var   = unname(d["var"]),
    delta_cap   = unname(d["capped"]),
    
    safe_pt     = safe_bc,              # SAFE-BC point estimator
    safe_var    = s$var,                # empirical Var(lnM*) over accepted draws
    
    safe_kept   = s$kept,               # accepted bootstrap draws
    safe_tried  = s$tried,              # total attempted draws
    safe_ok     = as.integer(isTRUE(s$status == "ok")),
    
    delta_fail  = as.integer(is.na(delta_pt)),   # PI undefined due to Delta <= 0
    safe_fail   = as.integer(is.na(safe_bc)))    # SAFE may fail if no usable draws
}


## -------- 4. parameter grid ------------------------------------------
## We evaluate a range of mean differences (theta) and sample sizes for both
## independent and paired designs.
##
## theta_vals spans small effects (near boundary where lnM may be undefined) to
## large effects (where acceptance is high and both methods should behave well).
##
## For independent designs we include both equal-n and unequal-n configurations.
## For paired designs we vary n (number of pairs) and fix rho=0.8 in DGP.

theta_vals <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 4, 5)

# Independent sample-size pairs:
#   (5,5), (10,10), (20,20), (100,100) plus some imbalanced designs
pairs_ind <- data.frame(n1 = c(5, 10, 20, 100, 3, 6, 12, 40),
                        n2 = c(5, 10, 20, 100, 7, 14, 28, 160))

# Build independent grid: theta x each (n1,n2) pair
grid_ind <- expand.grid(theta = theta_vals,
                        idx = seq_len(nrow(pairs_ind)))
grid_ind$n1     <- pairs_ind$n1[grid_ind$idx]
grid_ind$n2     <- pairs_ind$n2[grid_ind$idx]
grid_ind$design <- "indep"
grid_ind$idx    <- NULL

# Build paired grid: theta x n (number of pairs)
grid_dep <- expand.grid(theta = theta_vals,
                        n = c(5, 10, 20, 100))
grid_dep$n1     <- grid_dep$n
grid_dep$n2     <- grid_dep$n
grid_dep$n      <- NULL
grid_dep$design <- "paired"

# Combined parameter grid across both designs
param_grid <- rbind(grid_ind, grid_dep)


## -------- 5. simulation driver ---------------------------------------
## Settings are controlled by environment variables, so the same code can be run
## on laptops (small demo) and HPC/cluster runs (large K, large MIN_KEPT) without
## editing the script.
##
## Defaults shown here are large (K=1e5, MIN_KEPT=1e5) and are intended for
## production runs. For quick tests, set K_REPL and MIN_KEPT smaller.

set.seed(20250625)  # fixed seed for reproducibility

K_repl     <- as.integer(Sys.getenv("K_REPL", "100000"))     # Monte Carlo replicates per parameter set
MIN_KEPT   <- as.integer(Sys.getenv("MIN_KEPT", "100000"))   # accepted SAFE draws per replicate
CHUNK_INIT <- as.integer(Sys.getenv("CHUNK_INIT", "5000"))   # initial SAFE chunk size
CHUNK_MAX  <- as.numeric(Sys.getenv("CHUNK_MAX", "2000000")) # maximum SAFE chunk size
MAX_DRAWS  <- as.numeric(Sys.getenv("MAX_DRAWS", "Inf"))     # cap on total attempted SAFE draws
PATIENCE   <- as.integer(Sys.getenv("PATIENCE", "5"))        # stop SAFE after this many zero-accept chunks

outer_verbose <- TRUE  # print per-parameter-set completion messages
inner_every   <- 100   # print progress every N inner replicates

# Monte Carlo standard error (MCSE) helpers for reporting uncertainty of Monte
# Carlo summaries (bias, mean variance, coverage, etc.).
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

## runner() executes the full Monte Carlo simulation for ONE parameter set (one
## row of param_grid), returning a one-row data.frame of summary metrics.
runner <- function(i) {
  p   <- param_grid[i, ]
  sd0 <- 1  # fixed population SD in this study
  
  # True lnM depends on design; may be NA in boundary regions (msb <= msw)
  true_ln <- if (p$design == "indep")
    lnM_true_ind(p$theta, p$n1, p$n2, sd0)
  else
    lnM_true_dep(p$theta, p$n1, sd0)
  
  # Storage for replicate-level outputs:
  # M is (statistic x replicate). This is fast and compact.
  M <- matrix(NA_real_, 10, K_repl,
              dimnames = list(
                c("delta_pt","delta_var","delta_cap",
                  "safe_pt","safe_var",
                  "safe_kept","safe_tried","safe_ok",
                  "delta_fail","safe_fail"),
                NULL))
  
  # Inner Monte Carlo loop: generate dataset -> compute Delta + SAFE
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
    
    # Optional progress output for long runs
    if (k %% inner_every == 0)
      message(sprintf("[row %3d] replicate %7d / %d done", i, k, K_repl))
  }
  
  ## keep a conservative 'ok' set for delta-variance validity
  # We define ok_d based on Delta-method having a defined point estimate and a
  # positive variance estimate. We use this mask to create Mok, a restricted
  # matrix used for some downstream evaluations.
  ok_d <- is.finite(M["delta_pt", ]) & is.finite(M["delta_var", ]) & (M["delta_var", ] > 0)
  Mok  <- M[, ok_d, drop = FALSE]
  
  ## BASELINE: MC Var of SAFE-BC point estimator (computed on same Mok subset)
  # For variance evaluation we define a baseline “truth” as the Monte Carlo
  # variance of the SAFE-BC point estimator, using the SAME subset Mok so that
  # comparisons are made on a common set of replicates.
  ok_safe_pt <- is.finite(Mok["safe_pt", ])
  Var_MC_SAFEpt <- if (sum(ok_safe_pt) >= 2) var(Mok["safe_pt", ok_safe_pt]) else NA_real_
  
  ## coverage only meaningful if true_ln is finite
  # Coverage checks whether the nominal 95% interval contains the true target.
  # We only compute coverage if true lnM is defined (finite) and there are
  # enough valid replicates to evaluate.
  if (is.finite(true_ln) && sum(ok_d) > 0) {
    # Delta-method coverage using its own variance estimate
    cover_d <- abs(Mok["delta_pt", ] - true_ln) <= 1.96 * sqrt(Mok["delta_var", ])
    
    # SAFE coverage: needs SAFE point and SAFE variance defined and >0
    ok_s <- is.finite(Mok["safe_pt", ]) & is.finite(Mok["safe_var", ]) & (Mok["safe_var", ] > 0)
    cover_s <- abs(Mok["safe_pt", ok_s] - true_ln) <= 1.96 * sqrt(Mok["safe_var", ok_s])
    
    cover_delta <- mean(cover_d)
    cover_safe  <- if (length(cover_s)) mean(cover_s) else NA_real_
  } else {
    cover_delta <- NA_real_
    cover_safe  <- NA_real_
  }
  
  # Failure rates: proportion of replicates where estimator is NA
  delta_fail_prop <- mean(M["delta_fail", ], na.rm = TRUE)
  safe_fail_prop  <- mean(M["safe_fail",  ], na.rm = TRUE)
  
  # SAFE acceptance rate aggregated across replicates:
  # total accepted draws / total attempted draws (both sums across K replicates)
  boot_keep  <- sum(M["safe_kept", ],  na.rm = TRUE)
  boot_tried <- sum(M["safe_tried", ], na.rm = TRUE)
  boot_accept_prop <- if (boot_tried > 0) boot_keep / boot_tried else NA_real_
  
  # Delta variance capping frequency among ok_d replicates
  delta_cap_rate <- mean(M["delta_cap", ][ok_d], na.rm = TRUE)
  delta_cap_n    <- sum(M["delta_cap", ][ok_d] == 1, na.rm = TRUE)
  
  # Average reported variances
  mean_var_delta <- mean(M["delta_var", ], na.rm = TRUE)
  mean_var_safe  <- mean(M["safe_var",  ], na.rm = TRUE)
  
  ## Relative bias of variance: vs MC(SAFE-BC point)
  # Interpreted as:
  #   +20% means variance estimator is 20% larger than MC baseline variance
  #   -20% means variance estimator is 20% smaller than MC baseline variance
  relbias_delta <- if (is.finite(Var_MC_SAFEpt) && Var_MC_SAFEpt > 0)
    100 * (mean_var_delta / Var_MC_SAFEpt - 1) else NA_real_
  relbias_safe  <- if (is.finite(Var_MC_SAFEpt) && Var_MC_SAFEpt > 0)
    100 * (mean_var_safe  / Var_MC_SAFEpt - 1) else NA_real_
  
  # Summarise into a one-row data.frame
  out <- data.frame(
    theta      = p$theta,
    design     = p$design,
    n1         = p$n1,
    n2         = ifelse(p$design == "indep", p$n2, p$n1),
    true_lnM   = true_ln,
    
    # Mean point estimates
    delta_mean = mean(M["delta_pt", ], na.rm = TRUE),
    safe_mean  = mean(M["safe_pt",  ], na.rm = TRUE),
    
    # Bias relative to true lnM (only meaningful if true lnM is defined)
    delta_bias = if (is.finite(true_ln)) mean(M["delta_pt", ], na.rm = TRUE) - true_ln else NA_real_,
    safe_bias  = if (is.finite(true_ln)) mean(M["safe_pt",  ], na.rm = TRUE) - true_ln else NA_real_,
    
    # Average variance estimates
    mean_var_delta = mean_var_delta,
    mean_var_safe  = mean_var_safe,
    
    # MC baseline variance for SAFE-BC point estimator
    Var_MC_SAFEpt  = Var_MC_SAFEpt,
    relbias_delta  = relbias_delta,
    relbias_safe   = relbias_safe,
    
    # RMSE (only if true is defined)
    rmse_delta     = if (is.finite(true_ln)) sqrt(mean((M["delta_pt", ] - true_ln)^2, na.rm = TRUE)) else NA_real_,
    rmse_safe      = if (is.finite(true_ln)) sqrt(mean((M["safe_pt",  ] - true_ln)^2, na.rm = TRUE)) else NA_real_,
    
    # Coverage
    cover_delta    = cover_delta,
    cover_safe     = cover_safe,
    
    # Failure rates
    delta_fail_prop   = delta_fail_prop,
    safe_fail_prop    = safe_fail_prop,
    
    # Variance capping diagnostics
    delta_cap_rate = delta_cap_rate,
    delta_cap_n    = delta_cap_n,
    maxVar         = maxVar,
    
    # SAFE acceptance diagnostics
    boot_keep        = boot_keep,
    boot_tried       = boot_tried,
    boot_accept_prop = boot_accept_prop,
    SAFE_ok_rate     = mean(M["safe_ok", ], na.rm = TRUE),
    
    # MC standard errors for key summaries
    mcse_bias_delta   = if (is.finite(true_ln)) mcse_mean(M["delta_pt", ] - true_ln) else NA_real_,
    mcse_bias_safe    = if (is.finite(true_ln)) mcse_mean(M["safe_pt",  ] - true_ln) else NA_real_,
    mcse_varbar_delta = mcse_mean(M["delta_var", ]),
    mcse_varbar_safe  = mcse_mean(M["safe_var", ]),
    mcse_cover_delta  = if (is.finite(true_ln) && sum(ok_d) > 0) mcse_prop(cover_d) else NA_real_,
    mcse_cover_safe   = if (is.finite(true_ln) && exists("cover_s") && length(cover_s)) mcse_prop(cover_s) else NA_real_
  )
  
  # Attach replicate-level matrix so we can save raw runs later (optional)
  attr(out, "raw_M") <- M
  
  # Outer-loop progress message
  if (outer_verbose)
    message(sprintf("Finished row %3d  (%s  theta=%s  n1=%d n2=%d)",
                    i, p$design, p$theta, p$n1,
                    ifelse(p$design == "indep", p$n2, p$n1)))
  out
}


## -------- 5a. PARALLEL outer loop ------------------------------------
## We parallelise across parameter sets (rows of param_grid), because each row
## is independent and expensive (K_repl simulations per row).
##
## Windows uses PSOCK clusters; Unix/macOS can use forked multicore via pbapply.
## We default to half available cores unless N_CORES is set.

n_cores_use <- suppressWarnings(as.integer(Sys.getenv("N_CORES", "")))
if (is.na(n_cores_use) || n_cores_use <= 0) n_cores_use <- max(1L, detectCores()/2)

pbop <- pbapply::pboptions(type = "txt")

if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_cores_use)
  # Export everything in the global environment so worker processes see all
  # functions/objects defined above (runner, one_rep, grids, etc.)
  clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)
  results_list <- pbapply::pblapply(seq_len(nrow(param_grid)), runner, cl = cl)
  stopCluster(cl)
} else {
  # On Unix/macOS, pbapply can use multicore by passing an integer core count.
  results_list <- pbapply::pblapply(seq_len(nrow(param_grid)), runner, cl = n_cores_use)
}

pbapply::pboptions(pbop)

## -------- 5c. COLLAPSE (The SAFE way) ------------------------
## The outer loop returns a list. In rare cases, a parameter set can fail and
## return a non-data.frame object (e.g., an error message). We:
##   1) keep only successful data.frames
##   2) rbind them into one summary table
##   3) report which indices failed, so they can be investigated

is_df <- sapply(results_list, is.data.frame)
results <- do.call(rbind, results_list[is_df])

# Type-safety: ensure key fields are numeric (useful if expand.grid coerced)
results$theta <- as.numeric(as.character(results$theta))
results$n1    <- as.numeric(as.character(results$n1))

if(any(!is_df)) {
  message(sprintf("Warning: %d parameter sets failed and were excluded.", sum(!is_df)))
  message("Failed indices: ", paste(which(!is_df), collapse = ", "))
}

## -------- SAVE OUTPUTS -----------------------------------------------
## We save:
##   - summary table (RDS + CSV) for analysis/plots/manuscript tables
##   - full results_list (RDS) so we retain raw matrices as attributes
##   - optionally, raw replicate-level tables per parameter set in /raw_runs
##     (one file per parameter row), to enable detailed diagnostics later.

summary_rds <- sprintf("lnM_summary_SAFEfun_%s.rds", Sys.Date())
saveRDS(results, summary_rds)
message("Saved overall summary (RDS) to: ", summary_rds)

summary_csv <- sprintf("lnM_summary_SAFEfun_%s.csv", Sys.Date())
write.csv(results, file = summary_csv, row.names = FALSE)
message("Saved overall summary (CSV) to: ", summary_csv)

# Save full list (including failures + raw matrices in attributes)
saveRDS(results_list,
        file = sprintf("lnM_results_list_SAFEfun_%s.rds", Sys.Date()),
        compress = "xz")

## Optional: save raw replicate results (one file per param-grid row)
## Saved format: data.frame with K_repl rows (replicates) x stats columns
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
    
    # keep replicate id for traceability
    raw_df$rep <- seq_len(nrow(raw_df))
    
    # put rep first
    raw_df <- raw_df[, c("rep", setdiff(names(raw_df), "rep"))]
    
    raw_path <- file.path(raw_dir, sprintf("row_%03d.rds", i))
    saveRDS(raw_df, raw_path, compress = "xz")
    n_saved <- n_saved + 1L
  }
  
  message("Saved ", n_saved, " raw replicate file(s) into: ", raw_dir)
}