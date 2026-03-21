## ============================================================
## Plot |d|, lnM, and combined simulation results
## ------------------------------------------------------------
## Files expected in working directory:
##   - absd_summary_2026-03-21.csv
##   - lnM_summary_SAFEfun_2025-12-21.csv
##
## Output:
##   1) bias_absd.png
##   2) varbias_absd.png
##   3) bias_lnm.png
##   4) varbias_lnm.png
##   5) bias_combined_absd_lnm.png
##   6) varbias_combined_absd_lnm.png
## ============================================================

library(tidyverse)
library(here)

## ------------------------------------------------------------
## 1) read data
## ------------------------------------------------------------
absd <- read_csv(here("results","absd_summary_2026-03-21.csv"), show_col_types = FALSE)
lnm  <- read_csv(here("results","lnM_summary_SAFEfun_2025-12-21.csv"), show_col_types = FALSE)

## ------------------------------------------------------------
## 2) helper: design labels / facet labels
## ------------------------------------------------------------
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

## ------------------------------------------------------------
## 3) long data: |d|
## ------------------------------------------------------------
absd_bias_dat <- absd %>%
  dplyr::select(theta, design_group, facet_lab, absd_bias, fn_bias) %>%
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
  dplyr::select(theta, design_group, facet_lab, relbias_var_absd, relbias_var_fn) %>%
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

## ------------------------------------------------------------
## 4) long data: lnM
## ------------------------------------------------------------
lnm_bias_dat <- lnm %>%
  dplyr::select(theta, design_group, facet_lab, delta_bias, safe_bias) %>%
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
  dplyr::select(theta, design_group, facet_lab, relbias_delta, relbias_safe) %>%
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

## ------------------------------------------------------------
## 5) plot theme
## ------------------------------------------------------------
sim_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", colour = "grey50"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

## ------------------------------------------------------------
## 6) |d| bias plot
## ------------------------------------------------------------
p_bias_absd <- ggplot(absd_bias_dat, aes(x = theta, y = bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Naive |d|" = "red3",
      "Folded-normal |d|" = "blue3"
    )
  ) +
  labs(
    x = expression(theta ~ "(= true " * abs(d) * " here)"),
    y = "Bias of point estimator",
    colour = NULL,
    title = expression("Bias of " * abs(d) * " estimators")
  ) +
  sim_theme

print(p_bias_absd)

ggsave("bias_absd.png", p_bias_absd, width = 12, height = 8.5, dpi = 300)

## ------------------------------------------------------------
## 7) |d| variance-bias plot
## ------------------------------------------------------------
p_var_absd <- ggplot(absd_var_dat, aes(x = theta, y = rel_bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Naive |d| variance" = "red3",
      "Folded-normal |d| variance" = "blue3"
    )
  ) +
  labs(
    x = expression(theta ~ "(= true " * abs(d) * " here)"),
    y = "Relative bias of variance estimator (%)",
    colour = NULL,
    title = expression("Relative bias of variance estimators for " * abs(d))
  ) +
  sim_theme

print(p_var_absd)

ggsave("varbias_absd.png", p_var_absd, width = 12, height = 8.5, dpi = 300)

## ------------------------------------------------------------
## 8) lnM bias plot
## ------------------------------------------------------------
p_bias_lnm <- ggplot(lnm_bias_dat, aes(x = theta, y = bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Delta lnM" = "red3",
      "SAFE lnM"  = "blue3"
    )
  ) +
  labs(
    x = expression(theta),
    y = "Bias of point estimator",
    colour = NULL,
    title = "Bias of lnM estimators"
  ) +
  sim_theme

print(p_bias_lnm)

ggsave("bias_lnm.png", p_bias_lnm, width = 12, height = 8.5, dpi = 300)

## ------------------------------------------------------------
## 9) lnM variance-bias plot
## ------------------------------------------------------------
p_var_lnm <- ggplot(lnm_var_dat, aes(x = theta, y = rel_bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Delta lnM variance" = "red3",
      "SAFE lnM variance"  = "blue3"
    )
  ) +
  labs(
    x = expression(theta),
    y = "Relative bias of variance estimator (%)",
    colour = NULL,
    title = "Relative bias of variance estimators for lnM"
  ) +
  sim_theme

print(p_var_lnm)

ggsave("varbias_lnm.png", p_var_lnm, width = 12, height = 8.5, dpi = 300)

## ------------------------------------------------------------
## 10) combined bias plot: |d| + lnM together
## ------------------------------------------------------------
combined_bias <- bind_rows(
  absd_bias_dat %>% mutate(method_family = "|d|"),
  lnm_bias_dat  %>% mutate(method_family = "lnM")
)

p_bias_combined <- ggplot(combined_bias, aes(x = theta, y = bias, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, colour = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Naive |d|"          = "firebrick2",
      "Folded-normal |d|"  = "dodgerblue3",
      "Delta lnM"          = "darkorange2",
      "SAFE lnM"           = "darkgreen"
    )
  ) +
  labs(
    x = expression(theta),
    y = "Bias of point estimator",
    colour = NULL,
    title = "Bias of |d| and lnM estimators"
  ) +
  sim_theme

print(p_bias_combined)

ggsave("bias_combined_absd_lnm.png", p_bias_combined, width = 13, height = 9, dpi = 300)

## ------------------------------------------------------------
## 11) combined variance-bias plot: |d| + lnM together
## ------------------------------------------------------------
combined_var <- bind_rows(
  absd_var_dat %>% mutate(method_family = "|d|"),
  lnm_var_dat  %>% mutate(method_family = "lnM")
)

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
    colour = NULL,
    title = "Relative bias of variance estimators for |d| and lnM"
  ) +
  sim_theme

print(p_var_combined)

ggsave("varbias_combined_absd_lnm.png", p_var_combined, width = 13, height = 9, dpi = 300)