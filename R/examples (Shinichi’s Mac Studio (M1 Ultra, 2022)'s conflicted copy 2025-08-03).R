# examples
library(tidyverse)
library(metafor)        # escalc()
library(here)
library(orchaRd)

source(here("R", "lnM_SAFE7.R"))   # defines lnM_delta1_indep()   &   safe_lnM_indep()

# function for folded normal


folded_norm <- function(mu, var) {
  sigma  <- sqrt(var)
  ## avoid division-by-zero warnings
  if (sigma == 0) return(c(point = abs(mu), var = 0))
  
  z      <- mu / sigma
  meanFN <- sigma * sqrt(2 / pi) * exp(-0.5 * z^2) +
    mu    * (2 * pnorm(z) - 1)
  varFN  <- mu^2 + var - meanFN^2
  
  c(point = meanFN, var = varFN)
}

# lnM - delta

## ── 4.  ln M  (non-SAFE delta-method)  ─────────────────────────────────────
get_lnM_raw <- function(m1, m2, s1, s2, n1, n2) {
  out <- lnM_delta1_indep(m1, m2, s1, s2, n1, n2)   # <- returns c(point, var)
  tibble(
    yi_lnM_raw = out["point"],
    vi_lnM_raw = out["var"]          # ← sampling variance here
  )
}

# lnM - SAFE
## ── 5.  SAFE ln M  (parametric bootstrap)  ─────────────────────────────────
##        B = 10 000 resamples is a good default.
get_lnM_safe <- function(m1, m2, s1, s2, n1, n2, B = 1e3) {
  out <- safe_lnM_indep(m1, m2, s1, s2, n1, n2, B = B)
  tibble(
    yi_lnM_safe = out$point,
    vi_lnM_safe = out$var,           # ← bootstrap sampling variance
    draws_kept  = out$kept,
    draws_total = out$total
  )
}



###########
# Example 1
###########
#############################################

# data - read in from CSV files

dat <- read_csv(here("data", "data_almeida.csv"), show_col_types = FALSE)

# phylogenetic correlation matrix
corMat.env <- as.matrix(read_csv(here("data","Almeida_et_al_Phylo_correlations_update_fulldataset_env.csv")))


# effect size 

dat <-  escalc(measure = "SMD", 
               m1i = m1i, 
               m2i = m2i, 
               sd1i = sd1i, 
               sd2i = sd2i, 
               n1i = n1i, 
               n2i = n2i,
               data = dat, append = TRUE,
               var.names = c("yi_g",   "vi_g"))

dat <- dat %>%                                   # <- your existing data frame
  mutate(
    abs_g   = map2_dfr(yi_g,   vi_g,   ~ {
      out <- folded_norm(.x, .y)
      tibble(yi_g_abs   = out["point"], vi_g_abs   = out["var"])
    })
  ) %>%                                           # flatten the tibbles
  unnest(c(abs_g))

# log(0.5 + yi)
dat <- dat %>% 
  mutate(
    log_half_d  = log(0.5 + yi_g),
    vi_log      = vi_g / (0.5 + yi_g)^2       #  Delta-method var
  )


# lnM

dat <- dat %>% mutate(
  lnM_raw = pmap_dfr(
    list(m1i, m2i, sd1i, sd2i, n1i, n2i),
    get_lnM_raw)
) %>% unnest(lnM_raw)

dat <- dat %>% mutate(
  lnM_safe = pmap_dfr(
    list(m1i, m2i, sd1i, sd2i, n1i, n2i),
    get_lnM_safe)
) %>% unnest(lnM_safe)

## ── 6.  Quick diagnostics ────────────────────────────────────────────────
summary_stats <- list(
  n_total               = nrow(dat),
  
  ## delta-method lnM
  n_lnM_raw_point_NA    = sum(is.na(dat$yi_lnM_raw)),
  n_lnM_raw_var_NA      = sum(is.na(dat$vi_lnM_raw)),
  
  ## SAFE lnM
  n_lnM_safe_point_NA   = sum(is.na(dat$yi_lnM_safe)),
  n_lnM_safe_var_NA     = sum(is.na(dat$vi_lnM_safe)),
  
  ## rescued comparisons (point estimate only)
  n_rescued_by_SAFE     = sum(is.na(dat$yi_lnM_raw) & !is.na(dat$yi_lnM_safe)),
  
  ## rescued variances
  n_var_rescued_by_SAFE = sum(is.na(dat$vi_lnM_raw) & !is.na(dat$vi_lnM_safe))
)

print(summary_stats)


# random efffects
rand_list <- list(
  ~ 1 | Study,
  ~ 1 | Species,
  ~ 1 | Species.no.phylo,
  ~ 1 | ES.ID
)
  
## 3a. |d|
mod_abs <- rma.mv(
  yi   = yi_g_abs, 
  V    = vi_abs,
  mods = ~ magnitude + duration + Recovery,
  random = rand_list,
  R = list(Species = corMat.env),
  Rscale = 0,
  data   = env_data_full,
  method = "REML")

# bubble plot
bubble_plot(mod_abs, mod = "magnitude", ylab = "abs(SMD)", group = "Study")
bubble_plot(mod_abs, mod = "duration", ylab = "abs(SMD)", group = "Study")
bubble_plot(mod_abs, mod = "Recovery", ylab = "abs(SMD)", group = "Study")

## 3b. log(0.5 + d)
mod_log <- rma.mv(
  yi   = log_half_d,
  V    = vi_g_abs,
  mods = ~ magnitude + duration + Recovery,
  random = rand_list,
  R = list(Species = corMat.env),
  Rscale = 0,
  data   = env_data_full,
  method = "REML")

## 3c. lnM
mod_lnM <- rma.mv(
  yi   = yi_lnM_safe,
  V    = vi_lnM_safe,
  mods = ~ magnitude + duration + Recovery,
  random = rand_list,
  R = list(Species = corMat.env),
  Rscale = 0,
  data   = env_data_full,
  method = "REML")



###########
# Example 2
###########




###########
# Example 3
###########