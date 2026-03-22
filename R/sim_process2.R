# ================================================================
# Process and visualise simulation results for |d| and lnM
# ================================================================
#
# Purpose
#   This script is the post-processing / plotting companion to the
#   simulation scripts for:
#     1) lnM  (Delta plug-in vs SAFE)
#     2) |d|  (Naive absolute value vs Folded-normal plug-in)
#
#   It plays the same role as sim_process.R for the main lnM study,
#   but extends that workflow so that we can:
#     - generate stand-alone figures for the supplementary |d| simulation
#     - regenerate the lnM figures in the same plotting style
#     - produce combined |d| + lnM comparison figures
#     - reconstruct Monte Carlo error diagnostics for |d| from raw replicate
#       files saved on disk
#     - save compact summary tables for later insertion into index.qmd
#
# Expected inputs
#
#   Collapsed summary files (in results/)
#     - absd_summary_2026-03-21.csv
#     - lnM_summary_SAFEfun_2025-12-21.csv
#
#   Raw replicate files for |d| MC diagnostics (in raw_runs_absd/)
#     - row_001.rds, row_002.rds, ..., row_XXX.rds
#
# Main figure outputs (saved to figs2/)
#   1)  bias_absd.png
#   2)  varbias_absd.png
#   3)  coverage_absd.png
#   4)  rmse_absd.png
#   5)  bias_lnm.png
#   6)  varbias_lnm.png
#   7)  coverage_lnm.png              [if present in lnM summary]
#   8)  rmse_lnm.png                  [if present in lnM summary]
#   9)  bias_combined_absd_lnm.png
#  10)  varbias_combined_absd_lnm.png
#  11)  mcse_bias_absd.png
#  12)  mcse_var_absd.png
#
# Table outputs (saved to data/)
#   - sim_overview_absd.csv
#   - sim_overview_lnm.csv
#   - sim_overview_combined.csv
#   - sim_table_absd_bias_var.csv
#   - sim_table_lnm_bias_var.csv
#   - sim_table_combined_bias_var.csv
#   - tab_mcse_absd.csv
#
# Notes
#   - The |d| results correspond to the supplementary simulation using:
#       * naive |d| = |d_hat|
#       * folded-normal plug-in point estimator
#       * folded-normal plug-in variance estimator
#
#   - The lnM results correspond to the main simulation using:
#       * Delta plug-in point estimator
#       * SAFE bias-corrected point estimator
#       * Delta variance estimator
#       * SAFE variance estimator
#
#   - Most plots and summary tables are built from the collapsed summary
#     files in results/.
#
#   - The two Monte Carlo error diagnostics for |d| use both sources:
#       * MCSE of bias is reconstructed from raw replicate files in
#         raw_runs_absd/
#       * MCSE of the mean variance estimator is read directly from the
#         collapsed |d| summary file
#
#   - This script assumes that the numbering of raw |d| files
#       row_001.rds, row_002.rds, ...
#     matches the row order of the parameter grid used in the original
#     |d| simulation driver.
# ================================================================

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(readr)
library(stringr)

# ------------------------------------------------------------
# 0) file paths and output folders
# ------------------------------------------------------------
#
# We assume the working project layout shown in your screenshots:
#   meta-analysis_of_magnitude/
#     results/
#     figs2/
#     data/
#
# If figs2/ or data/ do not exist, create them.

dir.create(here("figs2"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("data"),  showWarnings = FALSE, recursive = TRUE)

absd_file <- here("results", "absd_summary_2026-03-21.csv")
lnm_file  <- here("results", "lnM_summary_SAFEfun_2025-12-21.csv")

# ------------------------------------------------------------
# 1) read collapsed summary files
# ------------------------------------------------------------
#
# absd = supplementary |d| simulation summary
# lnm  = main lnM simulation summary

absd <- read_csv(absd_file, show_col_types = FALSE)
lnm  <- read_csv(lnm_file,  show_col_types = FALSE)

# ------------------------------------------------------------
# 2) helper: construct design groups and facet labels
# ------------------------------------------------------------
#
# We want the same facet order as in your earlier figures:
#   row 1: independent balanced
#   row 2: independent unbalanced
#   row 3: paired
#
# We therefore create:
#   design_group = indep_equal / indep_unequal / paired
#   facet_lab    = nice label for facet strips

make_design_labels <- function(dat) {
  dat %>%
    mutate(
      design_group = case_when(
        design == "paired" ~ "paired",
        design == "indep" & n1 == n2 ~ "indep_equal",
        design == "indep" & n1 != n2 ~ "indep_unequal",
        TRUE ~ NA_character_
      ),
      facet_lab = case_when(
        design_group == "indep_equal"   ~ paste0("indep: n1=n2=", n1),
        design_group == "indep_unequal" ~ paste0("indep: n1=", n1, ", n2=", n2),
        design_group == "paired"        ~ paste0("paired: n=", n1),
        TRUE ~ NA_character_
      )
    )
}

absd <- make_design_labels(absd)
lnm  <- make_design_labels(lnm)

facet_levels <- c(
  paste0("indep: n1=n2=", c(5, 10, 20, 100)),
  paste0("indep: n1=", c(3, 6, 12, 40), ", n2=", c(7, 14, 28, 160)),
  paste0("paired: n=", c(5, 10, 20, 100))
)

absd <- absd %>%
  mutate(facet_lab = factor(facet_lab, levels = facet_levels))

lnm <- lnm %>%
  mutate(facet_lab = factor(facet_lab, levels = facet_levels))

# ------------------------------------------------------------
# 3) helper: plotting theme
# ------------------------------------------------------------
#
# Keep style consistent with your existing simulation plots.

sim_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey50"),
    strip.text = element_text(face = "plain"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 4) reshape |d| results to long format
# ------------------------------------------------------------
#
# Point estimators:
#   - absd_bias  = bias of naive |d|
#   - fn_bias    = bias of folded-normal |d|
#
# Variance estimators:
#   - relbias_var_absd = relative bias of naive inherited variance
#   - relbias_var_fn   = relative bias of FN plug-in variance
#
# We also prepare coverage and RMSE long tables for supplementary use.

absd_bias_dat <- absd %>%
  select(theta, design_group, facet_lab, absd_bias, fn_bias) %>%
  pivot_longer(
    cols = c(absd_bias, fn_bias),
    names_to = "estimator",
    values_to = "bias"
  ) %>%
  mutate(
    estimator = recode(
      estimator,
      absd_bias = "Naive |d|",
      fn_bias   = "Folded-normal |d|"
    )
  )

absd_var_dat <- absd %>%
  select(theta, design_group, facet_lab, relbias_var_absd, relbias_var_fn) %>%
  pivot_longer(
    cols = c(relbias_var_absd, relbias_var_fn),
    names_to = "estimator",
    values_to = "rel_bias"
  ) %>%
  mutate(
    estimator = recode(
      estimator,
      relbias_var_absd = "Naive |d| variance",
      relbias_var_fn   = "Folded-normal |d| variance"
    )
  )

absd_cover_dat <- absd %>%
  select(theta, design_group, facet_lab, cover_absd, cover_fn) %>%
  pivot_longer(
    cols = c(cover_absd, cover_fn),
    names_to = "estimator",
    values_to = "coverage"
  ) %>%
  mutate(
    estimator = recode(
      estimator,
      cover_absd = "Naive |d|",
      cover_fn   = "Folded-normal |d|"
    )
  )

absd_rmse_dat <- absd %>%
  select(theta, design_group, facet_lab, rmse_absd, rmse_fn) %>%
  pivot_longer(
    cols = c(rmse_absd, rmse_fn),
    names_to = "estimator",
    values_to = "rmse"
  ) %>%
  mutate(
    estimator = recode(
      estimator,
      rmse_absd = "Naive |d|",
      rmse_fn   = "Folded-normal |d|"
    )
  )

# ------------------------------------------------------------
# 5) reshape lnM results to long format
# ------------------------------------------------------------
#
# Point estimators:
#   - delta_bias = bias of Delta plug-in lnM point estimate
#   - safe_bias  = bias of SAFE bias-corrected lnM point estimate
#
# Variance estimators:
#   - relbias_delta = relative bias of Delta variance estimator
#   - relbias_safe  = relative bias of SAFE variance estimator
#
# We also prepare coverage and RMSE long tables if they exist in the
# lnM summary file.

lnm_bias_dat <- lnm %>%
  select(theta, design_group, facet_lab, delta_bias, safe_bias) %>%
  pivot_longer(
    cols = c(delta_bias, safe_bias),
    names_to = "estimator",
    values_to = "bias"
  ) %>%
  mutate(
    estimator = recode(
      estimator,
      delta_bias = "Delta lnM",
      safe_bias  = "SAFE lnM"
    )
  )

lnm_var_dat <- lnm %>%
  select(theta, design_group, facet_lab, relbias_delta, relbias_safe) %>%
  pivot_longer(
    cols = c(relbias_delta, relbias_safe),
    names_to = "estimator",
    values_to = "rel_bias"
  ) %>%
  mutate(
    estimator = recode(
      estimator,
      relbias_delta = "Delta lnM variance",
      relbias_safe  = "SAFE lnM variance"
    )
  )

# coverage / RMSE are expected from the lnM summary generated by sim_study.R
lnm_cover_dat <- NULL
lnm_rmse_dat  <- NULL

if (all(c("cover_delta", "cover_safe") %in% names(lnm))) {
  lnm_cover_dat <- lnm %>%
    select(theta, design_group, facet_lab, cover_delta, cover_safe) %>%
    pivot_longer(
      cols = c(cover_delta, cover_safe),
      names_to = "estimator",
      values_to = "coverage"
    ) %>%
    mutate(
      estimator = recode(
        estimator,
        cover_delta = "Delta lnM",
        cover_safe  = "SAFE lnM"
      )
    )
}

if (all(c("rmse_delta", "rmse_safe") %in% names(lnm))) {
  lnm_rmse_dat <- lnm %>%
    select(theta, design_group, facet_lab, rmse_delta, rmse_safe) %>%
    pivot_longer(
      cols = c(rmse_delta, rmse_safe),
      names_to = "estimator",
      values_to = "rmse"
    ) %>%
    mutate(
      estimator = recode(
        estimator,
        rmse_delta = "Delta lnM",
        rmse_safe  = "SAFE lnM"
      )
    )
}

# ------------------------------------------------------------
# 6) helper: save both png and pdf if desired later
# ------------------------------------------------------------
#
# For now we save png only, because that is what you already used.

save_plot <- function(plot_obj, filename, width, height, dpi = 300) {
  ggsave(
    filename = here("figs2", filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

# ------------------------------------------------------------
# 7) |d| plots
# ------------------------------------------------------------

p_bias_absd <- ggplot(absd_bias_dat, aes(x = theta, y = bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Naive |d|" = "firebrick2",
      "Folded-normal |d|" = "dodgerblue3"
    )
  ) +
  labs(
    x = expression(theta ~ "(= true " * abs(d) * " here)"),
    y = "Bias of point estimator",
    title = expression("Bias of " * abs(d) * " estimators")
  ) +
  sim_theme

p_var_absd <- ggplot(absd_var_dat, aes(x = theta, y = rel_bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Naive |d| variance" = "firebrick2",
      "Folded-normal |d| variance" = "dodgerblue3"
    )
  ) +
  labs(
    x = expression(theta ~ "(= true " * abs(d) * " here)"),
    y = "Relative bias of variance estimator (%)",
    title = expression("Relative bias of variance estimators for " * abs(d))
  ) +
  sim_theme

p_cov_absd <- ggplot(absd_cover_dat, aes(x = theta, y = coverage, colour = estimator)) +
  geom_hline(yintercept = 0.95, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Naive |d|" = "firebrick2",
      "Folded-normal |d|" = "dodgerblue3"
    )
  ) +
  labs(
    x = expression(theta ~ "(= true " * abs(d) * " here)"),
    y = "Coverage of 95% Wald interval",
    title = expression("Coverage of " * abs(d) * " estimators")
  ) +
  sim_theme

p_rmse_absd <- ggplot(absd_rmse_dat, aes(x = theta, y = rmse, colour = estimator)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Naive |d|" = "firebrick2",
      "Folded-normal |d|" = "dodgerblue3"
    )
  ) +
  labs(
    x = expression(theta ~ "(= true " * abs(d) * " here)"),
    y = "RMSE",
    title = expression("RMSE of " * abs(d) * " estimators")
  ) +
  sim_theme

print(p_bias_absd)
print(p_var_absd)
print(p_cov_absd)
print(p_rmse_absd)

save_plot(p_bias_absd, "bias_absd.png",     width = 12, height = 8.5)
save_plot(p_var_absd,  "varbias_absd.png",  width = 12, height = 8.5)
save_plot(p_cov_absd,  "coverage_absd.png", width = 12, height = 8.5)
save_plot(p_rmse_absd, "rmse_absd.png",     width = 12, height = 8.5)

# ------------------------------------------------------------
# 8) lnM plots
# ------------------------------------------------------------

p_bias_lnm <- ggplot(lnm_bias_dat, aes(x = theta, y = bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Delta lnM" = "darkorange2",
      "SAFE lnM"  = "darkgreen"
    )
  ) +
  labs(
    x = expression(theta),
    y = "Bias of point estimator",
    title = "Bias of lnM estimators"
  ) +
  sim_theme

p_var_lnm <- ggplot(lnm_var_dat, aes(x = theta, y = rel_bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Delta lnM variance" = "darkorange2",
      "SAFE lnM variance"  = "darkgreen"
    )
  ) +
  labs(
    x = expression(theta),
    y = "Relative bias of variance estimator (%)",
    title = "Relative bias of variance estimators for lnM"
  ) +
  sim_theme

print(p_bias_lnm)
print(p_var_lnm)

save_plot(p_bias_lnm, "bias_lnm.png",    width = 12, height = 8.5)
save_plot(p_var_lnm,  "varbias_lnm.png", width = 12, height = 8.5)

if (!is.null(lnm_cover_dat)) {
  p_cov_lnm <- ggplot(lnm_cover_dat, aes(x = theta, y = coverage, colour = estimator)) +
    geom_hline(yintercept = 0.95, linetype = 2, linewidth = 0.4, colour = "grey40") +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ facet_lab, ncol = 4) +
    scale_colour_manual(
      values = c(
        "Delta lnM" = "darkorange2",
        "SAFE lnM"  = "darkgreen"
      )
    ) +
    labs(
      x = expression(theta),
      y = "Coverage of 95% interval",
      title = "Coverage of lnM estimators"
    ) +
    sim_theme
  
  print(p_cov_lnm)
  save_plot(p_cov_lnm, "coverage_lnm.png", width = 12, height = 8.5)
}

if (!is.null(lnm_rmse_dat)) {
  p_rmse_lnm <- ggplot(lnm_rmse_dat, aes(x = theta, y = rmse, colour = estimator)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ facet_lab, ncol = 4) +
    scale_colour_manual(
      values = c(
        "Delta lnM" = "darkorange2",
        "SAFE lnM"  = "darkgreen"
      )
    ) +
    labs(
      x = expression(theta),
      y = "RMSE",
      title = "RMSE of lnM estimators"
    ) +
    sim_theme
  
  print(p_rmse_lnm)
  save_plot(p_rmse_lnm, "rmse_lnm.png", width = 12, height = 8.5)
}

# ------------------------------------------------------------
# 9) combined |d| + lnM plots
# ------------------------------------------------------------
#
# These figures are useful for the supplement because they show, panel by
# panel, how the |d| estimators behave relative to the lnM estimators.

combined_bias <- bind_rows(
  absd_bias_dat %>% mutate(method_family = "|d|"),
  lnm_bias_dat  %>% mutate(method_family = "lnM")
)

combined_var <- bind_rows(
  absd_var_dat %>% mutate(method_family = "|d|"),
  lnm_var_dat  %>% mutate(method_family = "lnM")
)

p_bias_combined <- ggplot(combined_bias, aes(x = theta, y = bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Naive |d|"         = "firebrick2",
      "Folded-normal |d|" = "dodgerblue3",
      "Delta lnM"         = "darkorange2",
      "SAFE lnM"          = "darkgreen"
    )
  ) +
  labs(
    x = expression(theta),
    y = "Bias of point estimator",
    title = "Bias of |d| and lnM estimators"
  ) +
  sim_theme

p_var_combined <- ggplot(combined_var, aes(x = theta, y = rel_bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Naive |d| variance"         = "firebrick2",
      "Folded-normal |d| variance" = "dodgerblue3",
      "Delta lnM variance"         = "darkorange2",
      "SAFE lnM variance"          = "darkgreen"
    )
  ) +
  labs(
    x = expression(theta),
    y = "Relative bias of variance estimator (%)",
    title = "Relative bias of variance estimators for |d| and lnM"
  ) +
  sim_theme

print(p_bias_combined)
print(p_var_combined)

save_plot(p_bias_combined, "bias_combined_absd_lnm.png",    width = 13, height = 9)
save_plot(p_var_combined,  "varbias_combined_absd_lnm.png", width = 13, height = 9)

# ------------------------------------------------------------
# 10) summary tables for qmd
# ------------------------------------------------------------
#
# These tables are intentionally compact. They are not the full simulation
# outputs; instead they are “qmd-ready” summaries that can be read directly
# inside index.qmd later.
#
# A) Overview tables:
#    average absolute bias / variance-bias / coverage / RMSE by method
#
# B) Design-level summaries:
#    average values within design group (indep_equal / indep_unequal / paired)

sim_overview_absd <- tibble(
  method = c("Naive |d|", "Folded-normal |d|"),
  mean_abs_bias = c(
    mean(abs(absd$absd_bias), na.rm = TRUE),
    mean(abs(absd$fn_bias),   na.rm = TRUE)
  ),
  mean_abs_relbias_var = c(
    mean(abs(absd$relbias_var_absd), na.rm = TRUE),
    mean(abs(absd$relbias_var_fn),   na.rm = TRUE)
  ),
  mean_coverage = c(
    mean(absd$cover_absd, na.rm = TRUE),
    mean(absd$cover_fn,   na.rm = TRUE)
  ),
  mean_rmse = c(
    mean(absd$rmse_absd, na.rm = TRUE),
    mean(absd$rmse_fn,   na.rm = TRUE)
  )
)

sim_overview_lnm <- tibble(
  method = c("Delta lnM", "SAFE lnM"),
  mean_abs_bias = c(
    mean(abs(lnm$delta_bias), na.rm = TRUE),
    mean(abs(lnm$safe_bias),  na.rm = TRUE)
  ),
  mean_abs_relbias_var = c(
    mean(abs(lnm$relbias_delta), na.rm = TRUE),
    mean(abs(lnm$relbias_safe),  na.rm = TRUE)
  ),
  mean_coverage = c(
    if ("cover_delta" %in% names(lnm)) mean(lnm$cover_delta, na.rm = TRUE) else NA_real_,
    if ("cover_safe"  %in% names(lnm)) mean(lnm$cover_safe,  na.rm = TRUE) else NA_real_
  ),
  mean_rmse = c(
    if ("rmse_delta" %in% names(lnm)) mean(lnm$rmse_delta, na.rm = TRUE) else NA_real_,
    if ("rmse_safe"  %in% names(lnm)) mean(lnm$rmse_safe,  na.rm = TRUE) else NA_real_
  )
)

sim_overview_combined <- bind_rows(
  sim_overview_absd %>% mutate(family = "|d|"),
  sim_overview_lnm  %>% mutate(family = "lnM")
) %>%
  select(family, everything())

sim_table_absd_bias_var <- absd %>%
  group_by(design_group) %>%
  summarise(
    mean_abs_bias_naive = mean(abs(absd_bias), na.rm = TRUE),
    mean_abs_bias_fn    = mean(abs(fn_bias), na.rm = TRUE),
    mean_abs_var_naive  = mean(abs(relbias_var_absd), na.rm = TRUE),
    mean_abs_var_fn     = mean(abs(relbias_var_fn), na.rm = TRUE),
    mean_cover_naive    = mean(cover_absd, na.rm = TRUE),
    mean_cover_fn       = mean(cover_fn, na.rm = TRUE),
    mean_rmse_naive     = mean(rmse_absd, na.rm = TRUE),
    mean_rmse_fn        = mean(rmse_fn, na.rm = TRUE),
    .groups = "drop"
  )

sim_table_lnm_bias_var <- lnm %>%
  group_by(design_group) %>%
  summarise(
    mean_abs_bias_delta = mean(abs(delta_bias), na.rm = TRUE),
    mean_abs_bias_safe  = mean(abs(safe_bias), na.rm = TRUE),
    mean_abs_var_delta  = mean(abs(relbias_delta), na.rm = TRUE),
    mean_abs_var_safe   = mean(abs(relbias_safe), na.rm = TRUE),
    mean_cover_delta    = if ("cover_delta" %in% names(cur_data())) mean(cover_delta, na.rm = TRUE) else NA_real_,
    mean_cover_safe     = if ("cover_safe"  %in% names(cur_data())) mean(cover_safe,  na.rm = TRUE) else NA_real_,
    mean_rmse_delta     = if ("rmse_delta" %in% names(cur_data())) mean(rmse_delta, na.rm = TRUE) else NA_real_,
    mean_rmse_safe      = if ("rmse_safe"  %in% names(cur_data())) mean(rmse_safe,  na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

sim_table_combined_bias_var <- full_join(
  sim_table_absd_bias_var,
  sim_table_lnm_bias_var,
  by = "design_group"
)

# save qmd-ready tables
write_csv(sim_overview_absd,           here("data", "sim_overview_absd.csv"))
write_csv(sim_overview_lnm,            here("data", "sim_overview_lnm.csv"))
write_csv(sim_overview_combined,       here("data", "sim_overview_combined.csv"))
write_csv(sim_table_absd_bias_var,     here("data", "sim_table_absd_bias_var.csv"))
write_csv(sim_table_lnm_bias_var,      here("data", "sim_table_lnm_bias_var.csv"))
write_csv(sim_table_combined_bias_var, here("data", "sim_table_combined_bias_var.csv"))

# ------------------------------------------------------------
# 11) console diagnostics
# ------------------------------------------------------------
#
# These are useful quick checks so that, after running the script, you can
# immediately see if the processed results look sensible.

message("----- quick diagnostics: |d| -----")
message("Mean abs naive |d| bias: ", round(mean(abs(absd$absd_bias), na.rm = TRUE), 4))
message("Mean abs FN |d| bias:    ", round(mean(abs(absd$fn_bias),   na.rm = TRUE), 4))
message("Mean abs naive var RB:   ", round(mean(abs(absd$relbias_var_absd), na.rm = TRUE), 2))
message("Mean abs FN var RB:      ", round(mean(abs(absd$relbias_var_fn),   na.rm = TRUE), 2))

message("----- quick diagnostics: lnM -----")
message("Mean abs Delta lnM bias: ", round(mean(abs(lnm$delta_bias), na.rm = TRUE), 4))
message("Mean abs SAFE lnM bias:  ", round(mean(abs(lnm$safe_bias),  na.rm = TRUE), 4))
message("Mean abs Delta var RB:   ", round(mean(abs(lnm$relbias_delta), na.rm = TRUE), 2))
message("Mean abs SAFE var RB:    ", round(mean(abs(lnm$relbias_safe),  na.rm = TRUE), 2))

message("All figures saved to: ", here("figs2"))
message("All qmd-ready summary tables saved to: ", here("data"))


## ================================================================
## Fig Sx) Monte-Carlo error diagnostics for |d| from disk
##        (reconstructed from raw_runs_absd/*.rds)
## ================================================================
#
# Why this exists:
#   The |d| simulation driver computes MCSE quantities during the run, but if
#   we only have the saved summary file and raw replicate files, we can
#   reconstruct MCSE(bias) from disk.
#
# raw_runs_absd/row_XXX.rds contains replicate-level rows:
#   d_pt, d_var, absd_pt, absd_var, fn_pt, fn_var, ... plus rep id
#
# IMPORTANT:
#   This section assumes:
#     row_id extracted from row_XXX matches the ordering of param_grid used in
#     the |d| simulation driver.
## ================================================================

raw_dir_absd <- here("raw_runs_absd")
raw_files_absd <- list.files(raw_dir_absd, pattern = "^row_[0-9]{3}\\.rds$", full.names = TRUE)

if (!length(raw_files_absd)) {
  warning("No |d| raw files found in: ", raw_dir_absd)
}

row_index_from_path <- function(p) {
  as.integer(sub("^row_([0-9]{3})\\.rds$", "\\1", basename(p)))
}

# ------------------------------------------------------------
# Compute diagnostics for one |d| raw replicate file
# ------------------------------------------------------------
#
# For |d| there is no analogue of lnM's "ok_d" filter based on undefined PI.
# Here we only require finite point estimates.
#
# Since true |d| is constant within a grid cell:
#   sd(estimator - true_absd) = sd(estimator)
# so MCSE(bias) = sd(estimator) / sqrt(n_eff)

one_raw_diag_absd <- function(raw_df) {
  
  # drop rep column if present
  if ("rep" %in% names(raw_df)) raw_df <- raw_df %>% dplyr::select(-rep)
  
  ok_absd <- is.finite(raw_df$absd_pt)
  ok_fn   <- is.finite(raw_df$fn_pt)
  
  n_ok_absd <- sum(ok_absd)
  n_ok_fn   <- sum(ok_fn)
  
  sd_absd <- if (n_ok_absd >= 2) sd(raw_df$absd_pt[ok_absd], na.rm = TRUE) else NA_real_
  sd_fn   <- if (n_ok_fn   >= 2) sd(raw_df$fn_pt[ok_fn],     na.rm = TRUE) else NA_real_
  
  # Relative standard error (RSE) of the Monte Carlo variance estimate
  # under the usual normal-theory approximation:
  #   RSE(Var_hat) approx sqrt(2 / (n - 1))
  RSE_VarMC_absd <- if (n_ok_absd >= 3) sqrt(2 / (n_ok_absd - 1)) else NA_real_
  RSE_VarMC_fn   <- if (n_ok_fn   >= 3) sqrt(2 / (n_ok_fn   - 1)) else NA_real_
  
  data.frame(
    n_ok_absd = n_ok_absd,
    n_ok_fn   = n_ok_fn,
    sd_absd   = sd_absd,
    sd_fn     = sd_fn,
    RSE_VarMC_absdpt = RSE_VarMC_absd,
    RSE_VarMC_fnpt   = RSE_VarMC_fn
  )
}

diag_list_absd <- pbapply::pblapply(raw_files_absd, function(f) {
  rid <- row_index_from_path(f)
  raw_df <- readRDS(f)
  d <- one_raw_diag_absd(raw_df)
  d$row_id <- rid
  d
})

diag_df_absd <- bind_rows(diag_list_absd)

absd_results2 <- absd %>%
  mutate(row_id = seq_len(nrow(absd))) %>%
  left_join(diag_df_absd, by = "row_id") %>%
  mutate(
    mcse_bias_absd_disk = ifelse(is.finite(sd_absd) & n_ok_absd >= 2,
                                 sd_absd / sqrt(n_ok_absd), NA_real_),
    mcse_bias_fn_disk   = ifelse(is.finite(sd_fn)   & n_ok_fn >= 2,
                                 sd_fn / sqrt(n_ok_fn), NA_real_)
  )

mc_long_absd <- absd_results2 %>%
  transmute(
    theta, facet_lab,
    mcse_bias_absd = mcse_bias_absd_disk,
    mcse_bias_fn   = mcse_bias_fn_disk
  ) %>%
  pivot_longer(
    cols = c(mcse_bias_absd, mcse_bias_fn),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      mcse_bias_absd = "Naive |d|",
      mcse_bias_fn   = "Folded-normal |d|"
    )
  )

p_mc_absd <- ggplot(mc_long_absd,
                    aes(theta, value, colour = metric, group = metric)) +
  geom_line() +
  facet_wrap(~ facet_lab, ncol = 4) +
  labs(
    x = expression(theta),
    y = "MCSE of bias",
    colour = NULL,
    title = expression("Monte-Carlo error diagnostics: MCSE of bias for " * abs(d))
  ) +
  sim_theme +
  theme(legend.position = "bottom")

print(p_mc_absd)

ggsave(
  filename = here("figs2", "mcse_bias_absd.png"),
  plot = p_mc_absd,
  width = 12,
  height = 8.5,
  dpi = 300
)

## ================================================================
## Fig Sy) MCSE of the mean variance estimates for |d|
##        (read directly from summary file)
## ================================================================
#
# These MCSE quantities were computed during the |d| simulation run:
#   - mcse_varbar_absd : MCSE of mean(absd_var)
#   - mcse_varbar_fn   : MCSE of mean(fn_var)
#
# Interpretation:
#   - smaller MCSE means the average variance estimates are stable across MC reps
#   - larger MCSE indicates more MC replicates (K_repl) might be needed
## ================================================================

var_mcse_absd_df <- absd %>%
  transmute(
    theta, facet_lab,
    absd_var_MCSE = mcse_varbar_absd,
    fn_var_MCSE   = mcse_varbar_fn
  ) %>%
  pivot_longer(
    cols = c(absd_var_MCSE, fn_var_MCSE),
    names_to = "estimator", values_to = "mcse"
  ) %>%
  mutate(
    estimator = recode(
      estimator,
      absd_var_MCSE = "Naive |d| variance",
      fn_var_MCSE   = "Folded-normal |d| variance"
    )
  )

p_var_mcse_absd <- ggplot(var_mcse_absd_df,
                          aes(theta, mcse, colour = estimator, group = estimator)) +
  geom_line() +
  facet_wrap(~ facet_lab, ncol = 4) +
  labs(
    x = expression(theta),
    y = "MCSE of mean(variance estimate)",
    colour = NULL,
    title = expression("Monte-Carlo error of variance estimators for " * abs(d))
  ) +
  sim_theme +
  theme(legend.position = "bottom")

print(p_var_mcse_absd)

ggsave(
  filename = here("figs2", "mcse_var_absd.png"),
  plot = p_var_mcse_absd,
  width = 12,
  height = 8.5,
  dpi = 300
)

## ================================================================
## Table Sx: MCSE diagnostics for |d| by facet
## ================================================================

tab_mcse_absd <- absd_results2 %>%
  group_by(design, n1, n2, facet_lab) %>%
  summarise(
    mcse_bias_naive_max = max(mcse_bias_absd_disk, na.rm = TRUE),
    mcse_bias_fn_max    = max(mcse_bias_fn_disk,   na.rm = TRUE),
    mcse_var_naive_max  = max(mcse_varbar_absd,    na.rm = TRUE),
    mcse_var_fn_max     = max(mcse_varbar_fn,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(facet_lab, levels = levels(absd$facet_lab)))

write_csv(tab_mcse_absd, here("data", "tab_mcse_absd.csv"))