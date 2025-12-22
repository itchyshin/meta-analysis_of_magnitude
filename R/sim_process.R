# visualizing results

library(dplyr)
library(ggplot2)
library(tidyr)
library(here)
library(grid)

results <- readRDS(here("lnM_summary_SAFEfun_2025-12-21.rds"))

# ------------------------- facet ordering -------------------------------------
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

# a nicer strip label: "indep   n1=5   n2=5"
strip_labeller <- function(x) gsub(" n1=", "   n1=", gsub(" n2=", "   n2=", x))

base_theme <- theme_bw(11) +
  theme(
    legend.title = element_text(),
    legend.position = "right",
    strip.background = element_rect(fill = "grey85", colour = "grey40"),
    strip.text = element_text(face = "plain"),
    panel.grid.minor = element_blank()
  )

# ------------------------- legend-inside helpers ------------------------------
extract_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[idx]]
}

# Draw a faceted ggplot and overlay its legend inside the panels.
# x,y,w,h are in npc (0..1). Tweak these to land in your desired facet.
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

# ------------------------- 1) Bias: (estimate - true lnM) ---------------------
bias_df <- bind_rows(
  results %>% select(theta, facet_label, delta_bias) %>%
    rename(bias = delta_bias) %>% mutate(estimator = "PI"),
  results %>% select(theta, facet_label, safe_bias) %>%
    rename(bias = safe_bias)  %>% mutate(estimator = "SAFE")
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

# Print with legend inside (top-right-ish; tweak if needed)
draw_with_legend_inside(p_bias, x = 0.80, y = 0.73, w = 0.18, h = 0.26)

# If you also want the normal version:
# print(p_bias)

# ------------------------- 2) Relative bias of Var (%) ------------------------
rb_df <- bind_rows(
  results %>% select(theta, facet_label, relbias_delta) %>%
    rename(relbias = relbias_delta) %>% mutate(estimator = "Delta"),
  results %>% select(theta, facet_label, relbias_safe) %>%
    rename(relbias = relbias_safe)  %>% mutate(estimator = "SAFE")
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

# ------------------------- 3) Coverage (if you still want it) -----------------
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

# ------------------------- 4) RMSE (if you still want it) ---------------------
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
# Notes:
# - If the legend lands between panels, tweak x/y slightly.
# - Typical ranges for top-right panel in a 4x3 facet:
#   x ~ 0.78–0.86, y ~ 0.70–0.80, w ~ 0.14–0.22, h ~ 0.20–0.30
# ------------------------------------------------------------------------------

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

## ---- 7f. Plot diagnostics ----------------------------------------
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



################

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