# ================================================================
# Visualising lnM simulation results (Delta-method vs SAFE-BC)
# ================================================================
#
# Files produced by the simulation script (simulation16_SAFEfun_style.R):
#   - results/lnM_summary_SAFEfun_YYYY-MM-DD.rds  (one row per grid cell)
#   - raw_runs/row_XXX.rds                        (optional; replicate-level outputs)
#
# What each figure is intended to show
#   Fig 1: point estimate bias
#          Bias = E[ estimator ] - true lnM
#
#   Fig 2: relative bias of variance estimate
#          RelBias(Var) = 100 * ( mean(Var_hat) / Var_MC_SAFEpt - 1 )
#          where Var_MC_SAFEpt = Monte-Carlo variance of SAFE-BC point estimates
#          (computed within each grid cell).
#
#   Fig S1: coverage of nominal 95% intervals
#          Coverage = proportion of intervals containing true lnM
#
#   Fig S2: RMSE of point estimates
#          RMSE = sqrt( E[ (estimate - true)^2 ] )
#
#   Fig S3: MCSE of bias (from raw_runs)
#          MCSE(bias) approx = sd(estimator) / sqrt(n_eff)
#          (computed per grid cell using replicate-level files)
#
#   Fig S4: MCSE of mean variance estimates (from run summary)
#          MCSE(mean Var_hat) = sqrt( var(Var_hat) / n_eff )
#
# Notes on terminology
#   - "PI" here refers to the Delta-method plug-in estimator (point estimate)
#   - "Delta" in Fig 2 refers to the Delta-method *variance* estimator
#   - "SAFE" refers to SAFE bootstrap outputs; "SAFE-BC" is the bias-corrected
#     SAFE point estimate used throughout the main manuscript comparisons.
# ================================================================

#library(MASS)      # often used elsewhere (be aware it masks select(); we use dplyr::select explicitly)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(here)
library(kableExtra)
library(grid)      # needed for manual legend placement
library(pbapply)   # for pblapply in the raw-run diagnostics section

results <- readRDS(here("results","lnM_summary_SAFEfun_2025-12-21.rds"))

# ------------------------- load summary results -----------------------------
# This RDS is the collapsed output from the simulation driver:
# one row per grid cell (theta x design x (n1,n2)).
results <- readRDS(here("results","lnM_summary_SAFEfun_2025-12-21.rds"))

# ------------------------- facet ordering -------------------------------------
# We build a single facet label so each panel corresponds to a design + sample sizes.
# This is purely for plotting convenience.
results <- results %>%
  mutate(facet_label = paste0(design, " n1=", n1, " n2=", n2))

# Create a helper data frame listing unique facets and their attributes
facet_info <- results %>% distinct(facet_label, design, n1, n2)

# We want the facet order to be:
#   1) independent & balanced (n1==n2), increasing n
#   2) independent & unbalanced, increasing n1
#   3) paired (always balanced), increasing n
balanced_ind <- facet_info %>%
  filter(design == "indep", n1 == n2) %>% arrange(n1) %>% pull(facet_label)

unbalanced_ind <- facet_info %>%
  filter(design == "indep", n1 != n2) %>% arrange(n1) %>% pull(facet_label)

paired_balanced <- facet_info %>%
  filter(design == "paired") %>% arrange(n1) %>% pull(facet_label)

new_levels <- c(balanced_ind, unbalanced_ind, paired_balanced)
results <- results %>% mutate(facet_label = factor(facet_label, levels = new_levels))

# A nicer strip label for the facets: "indep   n1=5   n2=5"
strip_labeller <- function(x) gsub(" n1=", "   n1=", gsub(" n2=", "   n2=", x))

# Base theme shared across plots so the style is consistent
base_theme <- theme_bw(11) +
  theme(
    legend.title = element_text(),
    legend.position = "right",
    strip.background = element_rect(fill = "grey85", colour = "grey40"),
    strip.text = element_text(face = "plain"),
    panel.grid.minor = element_blank()
  )

# ------------------------- legend-inside helpers ------------------------------
# We draw the legend inside the faceted plot (useful when there are many panels).
# NOTE: This is a “manual” grid approach—good for final figures, but a bit fragile:
# you may need to tweak x/y depending on device size and facet layout.
extract_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[idx]]
}

# Draw a faceted ggplot and overlay its legend inside the panels.
# x,y,w,h are in npc (0..1) relative coordinates on the page.
draw_with_legend_inside <- function(p, x = 0.80, y = 0.73, w = 0.18, h = 0.26) {
  leg <- extract_legend(p + theme(legend.position = "right"))
  p0  <- p + theme(legend.position = "none")
  
  g <- ggplotGrob(p0)
  grid.newpage()
  grid.draw(g)
  
  vp <- viewport(x = x, y = y, width = w, height = h, just = c("left", "bottom"))
  pushViewport(vp)
  grid.draw(leg)
  upViewport()
  
  invisible(NULL)
}

# ================================================================
# Fig 1) Bias of point estimates
# ================================================================
# bias_df stacks the bias columns into a “long” data frame:
#   estimator = PI (Delta plug-in) or SAFE (SAFE-BC point estimate)
#   bias      = estimate - true lnM
#
# Interpretation:
#   - Bias = 0 is ideal (horizontal dashed line).
#   - Large positive bias means systematic overestimation of lnM.
#   - Large negative bias means systematic underestimation of lnM.
# ================================================================
bias_df <- dplyr::bind_rows(
  results %>% dplyr::select(theta, facet_label, delta_bias) %>%
    dplyr::rename(bias = delta_bias) %>% dplyr::mutate(estimator = "PI"),
  results %>% dplyr::select(theta, facet_label, safe_bias) %>%
    dplyr::rename(bias = safe_bias)  %>% dplyr::mutate(estimator = "SAFE")
)

p_bias <- ggplot(bias_df, aes(theta, bias, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3,
             labeller = labeller(facet_label = strip_labeller)) +
  labs(x = expression(theta),
       y = "Bias (estimate \u2212 true lnM)",
       colour = "estimator") +
  scale_colour_manual(values = c("PI" = "firebrick", "SAFE" = "steelblue")) +
  base_theme

draw_with_legend_inside(p_bias, x = 0.80, y = 0.73, w = 0.18, h = 0.26)

# ================================================================
# Fig 2) Relative bias of variance estimates (percent)
# ================================================================
# Here the “truth” baseline is:
#   Var_MC_SAFEpt = Monte-Carlo variance of the SAFE-BC point estimates
# within each grid cell (computed in the simulation driver).
#
# relbias_delta and relbias_safe are:
#   100 * ( mean(Var_hat) / Var_MC_SAFEpt - 1 )
#
# Interpretation:
#   - 0% is ideal.
#   - +50% means the variance estimator is 50% too large on average
#     relative to the MC baseline.
#   - -50% means it is 50% too small (intervals may undercover).
# ================================================================
rb_df <- bind_rows(
  results %>% dplyr::select(theta, facet_label, relbias_delta) %>%
    dplyr::rename(relbias = relbias_delta) %>% dplyr::mutate(estimator = "Delta"),
  results %>% dplyr::select(theta, facet_label, relbias_safe) %>%
    dplyr::rename(relbias = relbias_safe)  %>% dplyr::mutate(estimator = "SAFE")
)

p_relbias <- ggplot(rb_df, aes(theta, relbias, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3,
             labeller = labeller(facet_label = strip_labeller)) +
  labs(x = expression(theta),
       y = "Relative bias of Var (%)",
       colour = "estimator") +
  scale_colour_manual(values = c("Delta" = "firebrick", "SAFE" = "steelblue")) +
  base_theme

draw_with_legend_inside(p_relbias, x = 0.80, y = 0.73, w = 0.18, h = 0.26)

# ================================================================
# Fig S1) Coverage
# ================================================================
# Coverage is only meaningful where true lnM is defined (not NA).
# In grid cells where true lnM is NA, coverage values may be NA as well.
#
# Interpretation:
#   - dashed line at 0.95 is nominal 95% target.
#   - below 0.95 indicates under-coverage (intervals too narrow)
#   - above 0.95 indicates over-coverage (intervals conservative)
# ================================================================
cov_df <- bind_rows(
  results %>% transmute(theta, facet_label, estimator = "PI",   cover = cover_delta),
  results %>% transmute(theta, facet_label, estimator = "SAFE", cover = cover_safe)
)

p_cov <- ggplot(cov_df, aes(theta, cover, colour = estimator, group = estimator)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3,
             labeller = labeller(facet_label = strip_labeller)) +
  labs(x = expression(theta), y = "Empirical coverage", colour = "estimator") +
  scale_colour_manual(values = c("PI" = "firebrick", "SAFE" = "steelblue")) +
  base_theme

draw_with_legend_inside(p_cov, x = 0.80, y = 0.73, w = 0.18, h = 0.26)

# ================================================================
# Fig S2) RMSE of point estimates
# ================================================================
# RMSE combines bias and variance:
#   RMSE = sqrt( mean( (estimate - true)^2 ) )
# Interpretation:
#   - lower is better
#   - RMSE often decreases as theta increases or n increases
# ================================================================
rmse_df <- bind_rows(
  results %>% transmute(theta, facet_label, estimator = "PI",   rmse = rmse_delta),
  results %>% transmute(theta, facet_label, estimator = "SAFE", rmse = rmse_safe)
)

p_rmse <- ggplot(rmse_df, aes(theta, rmse, colour = estimator, group = estimator)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ facet_label, ncol = 4, nrow = 3,
             labeller = labeller(facet_label = strip_labeller)) +
  labs(x = expression(theta),
       y = expression(RMSE~"("~hat(lnM)~"\u2212"~lnM~")"),
       colour = "estimator") +
  scale_colour_manual(values = c("PI" = "firebrick", "SAFE" = "steelblue")) +
  base_theme

draw_with_legend_inside(p_rmse, x = 0.80, y = 0.73, w = 0.18, h = 0.26)

# ------------------------------------------------------------------------------
# Quick diagnostics printed to console:
#   - average absolute bias across all grid cells
#   - average absolute relative variance bias across all grid cells
#   - average delta variance capping rate
#   - average SAFE "ok" rate and acceptance probability
# ------------------------------------------------------------------------------
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
## Fig S3) Monte-Carlo error diagnostics from disk (raw_runs/*.rds)
##        (NO dependency on results_list)
## ================================================================
#
# Why this exists:
#   The simulation driver computes MCSE quantities during the run, but if you
#   only have the saved summary results and raw replicate files, you can
#   reconstruct MCSE(bias) “from disk”.
#
# raw_runs/row_XXX.rds contains replicate-level rows:
#   delta_pt, delta_var, safe_pt, safe_var, ... plus rep id
#
# IMPORTANT:
#   This section assumes:
#     row_id extracted from row_XXX matches the ordering of param_grid used in
#     the simulation driver.
# ================================================================

raw_dir <- here("raw_runs")
raw_files <- list.files(raw_dir, pattern = "^row_[0-9]{3}\\.rds$", full.names = TRUE)
if (!length(raw_files)) stop("No raw files found in: ", raw_dir)

row_index_from_path <- function(p) as.integer(sub("^row_([0-9]{3})\\.rds$", "\\1", basename(p)))

# Compute diagnostics for a single grid cell (single raw_df)
one_raw_diag <- function(raw_df) {
  
  # Defensive: remove rep id if present (rep is not used in computations)
  if ("rep" %in% names(raw_df)) raw_df <- raw_df %>% dplyr::select(-rep)
  
  # Match the “ok_d” logic used in the simulation summary:
  # keep only replicates where Delta-method pt & var are defined and var>0
  ok_d <- is.finite(raw_df$delta_pt) & is.finite(raw_df$delta_var) & (raw_df$delta_var > 0)
  
  # If nothing survives the ok_d filter, return NAs
  if (!any(ok_d)) {
    return(data.frame(
      n_ok_d = 0L,
      n_ok_safe = 0L,
      RSE_VarMC_SAFEpt = NA_real_,
      sd_PI = NA_real_,
      sd_SAFE = NA_real_
    ))
  }
  
  df_ok <- raw_df[ok_d, , drop = FALSE]
  
  # Within the ok_d subset, check how many finite SAFE-BC point estimates exist
  ok_safe_pt <- is.finite(df_ok$safe_pt)
  n_ok_safe <- sum(ok_safe_pt)
  
  # Relative standard error (RSE) of the MC variance estimate of SAFE-BC point:
  #   RSE(Var_hat) approx sqrt(2/(n-1)) for sample variance under normality
  RSE_varmc <- if (n_ok_safe >= 3) sqrt(2 / (n_ok_safe - 1)) else NA_real_
  
  # SD of point estimates (needed for MCSE(bias) below)
  # Since true_lnM is constant within a grid cell, sd(est - true) = sd(est).
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

diag_list <- pbapply::pblapply(raw_files, function(f) {
  rid <- row_index_from_path(f)
  raw_df <- readRDS(f)
  d <- one_raw_diag(raw_df)
  d$row_id <- rid
  d
})

diag_df <- bind_rows(diag_list)

results2 <- results %>%
  mutate(row_id = seq_len(nrow(results))) %>%
  left_join(diag_df, by = "row_id")

# Recompute MCSE(bias) “from disk”:
#   MCSE(bias) ≈ sd(estimator)/sqrt(n_eff)
# where n_eff matches the number of usable replicates in each grid cell.
results2 <- results2 %>%
  mutate(
    mcse_bias_PI_disk   = ifelse(is.finite(sd_PI)   & n_ok_d >= 2,    sd_PI   / sqrt(n_ok_d),     NA_real_),
    mcse_bias_SAFE_disk = ifelse(is.finite(sd_SAFE) & n_ok_safe >= 2, sd_SAFE / sqrt(n_ok_safe),  NA_real_)
  )

mc_long <- results2 %>%
  transmute(
    theta, facet_label,
    mcse_bias_PI   = mcse_bias_PI_disk,
    mcse_bias_SAFE = mcse_bias_SAFE_disk
  ) %>%
  pivot_longer(
    cols = c(mcse_bias_PI, mcse_bias_SAFE),
    names_to = "metric", values_to = "value"
  )

p_mc <- ggplot(mc_long, aes(theta, value, colour = metric, group = metric)) +
  geom_line() +
  facet_wrap(~ facet_label, ncol = 4) +
  labs(
    x = expression(theta),
    y = "MCSE of bias",
    title = "Monte-Carlo error diagnostics: MCSE of bias (per grid cell; from raw_runs)"
  ) +
  base_theme + 
  theme(legend.position = "bottom")

print(p_mc)

## ================================================================
## Fig S4) MCSE of the mean variance estimates (from summary file)
## ================================================================
# These MCSE quantities were computed during the simulation run:
#   - mcse_varbar_delta : MCSE of mean(delta_var)
#   - mcse_varbar_safe  : MCSE of mean(safe_var)
#
# Interpretation:
#   - smaller MCSE means the average variance estimates are stable across MC reps
#   - larger MCSE indicates more MC replicates (K_repl) might be needed
# ================================================================

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
  base_theme + 
  theme(legend.position = "bottom")

print(p_var_mcse)

# Optional: log scale if MCSE spans orders of magnitude
# p_var_mcse_log <- p_var_mcse + scale_y_log10()
# print(p_var_mcse_log)

# ================================================================
# Table 1: Method feasibility / stability by facet (summarised over theta)
# Columns:
#   - delta_fail_prop: PI undefined rate
#   - safe_fail_prop:  SAFE-BC undefined rate
#   - SAFE_ok_rate:    proportion of SAFE runs with status == "ok"
#   - delta_cap_rate:  proportion of PI variance estimates capped at maxVar
# ================================================================

tab1_stability <- results %>%
  group_by(design, n1, n2, facet_label) %>%
  summarise(
    theta_min = min(theta, na.rm = TRUE),
    theta_max = max(theta, na.rm = TRUE),
    
    # "how often do we fail?"
    delta_fail_median = median(delta_fail_prop, na.rm = TRUE),
    safe_fail_median  = median(safe_fail_prop,  na.rm = TRUE),
    
    # "how often does SAFE return ok status?"
    SAFE_ok_median    = median(SAFE_ok_rate, na.rm = TRUE),
    
    # "how often does Delta variance get capped?"
    delta_cap_median  = median(delta_cap_rate, na.rm = TRUE),
    
    # also show worst-case over theta (useful for identifying boundary regions)
    delta_fail_max = max(delta_fail_prop, na.rm = TRUE),
    safe_fail_max  = max(safe_fail_prop,  na.rm = TRUE),
    delta_cap_max  = max(delta_cap_rate,  na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(factor(facet_label, levels = levels(results$facet_label))) %>%
  mutate(
    design = as.character(design),
    facet_label = as.character(facet_label)
  )

kable(tab1_stability, digits = 3,
      caption = "Table 1. Feasibility / stability diagnostics summarised across θ within each design × (n1,n2) facet.") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))

# ================================================================
# Table 2: SAFE bootstrap acceptance (kept/tried) by facet (summarised over theta)
# boot_accept_prop is already kept/tried aggregated across MC reps in each grid cell.
# We summarise across theta to get a compact facet-level picture.
# ================================================================

tab2_accept <- results %>%
  group_by(design, n1, n2, facet_label) %>%
  summarise(
    q0  = quantile(boot_accept_prop, 0.00, na.rm = TRUE),
    q25 = quantile(boot_accept_prop, 0.25, na.rm = TRUE),
    q50 = quantile(boot_accept_prop, 0.50, na.rm = TRUE),
    q75 = quantile(boot_accept_prop, 0.75, na.rm = TRUE),
    q100= quantile(boot_accept_prop, 1.00, na.rm = TRUE),
    mean = mean(boot_accept_prop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(facet_label, levels = levels(results$facet_label)))

kable(tab2_accept, digits = 4,
      caption = "Table 2. SAFE acceptance rate (kept/tried) summarised across θ for each facet. Quantiles highlight worst-case acceptance.") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))

# ================================================================
# Table 3: Worst offenders (grid-cell level)
#   A) lowest acceptance
#   B) highest delta cap rate
#   C) highest failure rates
# ================================================================

tab3_low_accept <- results %>%
  dplyr::select(design, n1, n2, theta, facet_label, boot_accept_prop, SAFE_ok_rate) %>%
  arrange(boot_accept_prop) %>%
  head(20)

tab3_high_cap <- results %>%
  dplyr::select(design, n1, n2, theta, facet_label, delta_cap_rate, maxVar) %>%
  arrange(desc(delta_cap_rate)) %>%
  head(20)

tab3_fail <- results %>%
  dplyr::select(design, n1, n2, theta, facet_label, delta_fail_prop, safe_fail_prop, SAFE_ok_rate) %>%
  mutate(worst_fail = pmax(delta_fail_prop, safe_fail_prop, na.rm = TRUE)) %>%
  arrange(desc(worst_fail)) %>%
  dplyr::select(-worst_fail) %>%
  head(20)

kable(tab3_low_accept, digits = 4,
      caption = "Table 3A. Lowest SAFE acceptance grid cells (smallest kept/tried).") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover","condensed"))

kable(tab3_high_cap, digits = 3,
      caption = "Table 3B. Highest Delta variance capping rates (fraction of replicates capped).") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover","condensed"))

kable(tab3_fail, digits = 3,
      caption = "Table 3C. Highest failure rates (PI undefined and/or SAFE-BC undefined).") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover","condensed"))

# ================================================================
# Table 4: MCSE diagnostics (run-computed) by facet
# Uses:
#   - mcse_bias_delta, mcse_bias_safe (if present in your results file)
#   - mcse_varbar_delta, mcse_varbar_safe (present)
# If mcse_bias_* are missing in your stored summary, comment those lines out.
# ================================================================

has_mcse_bias <- all(c("mcse_bias_delta","mcse_bias_safe") %in% names(results))

tab4_mcse <- results %>%
  group_by(design, n1, n2, facet_label) %>%
  summarise(
    # bias MCSE (optional)
    mcse_bias_PI_max   = if (has_mcse_bias) max(mcse_bias_delta, na.rm = TRUE) else NA_real_,
    mcse_bias_SAFE_max = if (has_mcse_bias) max(mcse_bias_safe,  na.rm = TRUE) else NA_real_,
    
    # variance MCSE (these are in your results)
    mcse_var_PI_max    = max(mcse_varbar_delta, na.rm = TRUE),
    mcse_var_SAFE_max  = max(mcse_varbar_safe,  na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(factor(facet_label, levels = levels(results$facet_label)))

kable(tab4_mcse, digits = 5,
      caption = "Table 4. Monte-Carlo precision (max MCSE over θ within each facet).") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover","condensed"))




