# examples
library(tidyverse)
library(metafor)        # escalc()
library(here)
library(orchaRd)

source(here("R", "lnM_SAFE8.R"))   # defines lnM_delta1_indep()   &   safe_lnM_indep()

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
corMat.env <- read.csv(here("data","Almeida_et_al_Phylo_correlations_update_fulldataset.csv"),sep=",",header=TRUE, row.names = 1)
#corMat.env <- as.matrix(read_csv(here("data","Almeida_et_al_Phylo_correlations_update_fulldataset_env.csv"),header=TRUE, row.names = 1))


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

  
# # meta-analysis
# 
# ma_abs <- rma.mvrma.mv(
#   yi   = yi_g=d_abs, 
#   V    = vi_g_abs,
#   mods = ~ magnitude + duration + Recovery,
#   random =  list(
#     ~ 1 | Study,
#     ~ 1 | Species,
#     ~ 1 | Species.no.phylo,
#     ~ 1 | ES.ID
#   ),
#   #R = list(Species = corMat.env),
#   #Rscale = 0,
#   data   = dat_cond,
#   method = "REML")
# 
# 
# 
# ## 3a. |d|

# remove . and replace with " " space 

dat$Species <- gsub("\\.", " ", dat$Species)

# filtering for just condition data

dat_cond <- dat %>% filter(Category == "Condition")
dim(dat_cond)


mod_abs <- rma.mv(
  yi   = yi_g_abs, 
  V    = vi_g_abs,
  mods = ~ magnitude + duration + Recovery,
  random =  list(
    ~ 1 | Study,
    ~ 1 | Species,
    ~ 1 | Species.no.phylo,
    ~ 1 | ES.ID
  ),
  #R = list(Species = corMat.env),
  #Rscale = 0,
  data   = dat_cond,
  method = "REML")
summary(mod_abs)
# bubble plot
p1 <- bubble_plot(mod_abs, mod = "magnitude", ylab = "abs(SMD)", group = "Study")
p2 <- bubble_plot(mod_abs, mod = "duration", ylab = "abs(SMD)", group = "Study")
p3 <- bubble_plot(mod_abs, mod = "Recovery", ylab = "abs(SMD)", group = "Study")

# ## 3b. log(0.5 + d)
# mod_log <- rma.mv(
#   yi   = log_half_d,
#   V    = vi_g_abs,
#   mods = ~ magnitude + duration + Recovery,
#   random = rand_list,
#  #R = list(Species = corMat.env),
#   #Rscale = 0,
#   data   = dat_cond,
#   method = "REML")
# 
# summary(mod_log)

## 3c. lnM
mod_lnM <- rma.mv(
  yi   = yi_lnM_safe,
  V    = vi_lnM_safe,
  mods = ~ magnitude + duration + Recovery,
  random =  list(
    ~ 1 | Study,
    ~ 1 | Species,
    ~ 1 | Species.no.phylo,
    ~ 1 | ES.ID
  ),
 # R = list(Species = corMat.env),
 # Rscale = 0,
  data   = dat_cond,
  method = "REML")

summary(mod_lnM)

# adjustment


## 3c. lnM

dat_cond$yi_lnM_safe2 <- dat_cond$yi_lnM_safe +log(sqrt(2))

mod_lnM <- rma.mv(
  yi   = yi_lnM_safe2,
  V    = vi_lnM_safe,
  mods = ~ magnitude + duration + Recovery,
  random =  list(
    ~ 1 | Study,
    ~ 1 | Species,
    ~ 1 | Species.no.phylo,
    ~ 1 | ES.ID
  ),
  # R = list(Species = corMat.env),
  # Rscale = 0,
  data   = dat_cond,
  method = "REML")

summary(mod_lnM)

# bubble plot
p4 <- bubble_plot(mod_lnM, mod = "magnitude", ylab = "lnM", group = "Study")
p5 <- bubble_plot(mod_lnM, mod = "duration", ylab = "lnM", group = "Study")
p6 <- bubble_plot(mod_lnM, mod = "Recovery", ylab = "lnM", group = "Study")

# want to plot p1 - p6 in two rows

library(patchwork)

(p1 + p2 + p3) / (p4 + p5 + p6) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16, face = "bold"))

# egger regressoin - dat$n0 <- (dat$sample.size.control * dat$sample.size.noise.1) / (dat$sample.size.control + dat$sample.size.noise.1)
#dat$nSE <- 1 / sqrt(dat$n0)
#dat$nV <- 1 / dat$n0

dat_cond$n0 <- (dat_cond$n1i * dat_cond$n2i) / (dat_cond$n1i + dat_cond$n2i)
dat_cond$nSE <- 1 / sqrt(dat_cond$n0)
dat_cond$nV <- 1 / dat_cond$n0

# the model

mod_egger <- rma.mv(
  yi   = yi_lnM_safe,
  V    = vi_lnM_safe,
  mods = ~ nSE,
  random =  list(
    ~ 1 | Study,
    ~ 1 | Species,
    ~ 1 | Species.no.phylo,
    ~ 1 | ES.ID
  ),
 # R = list(Species = corMat.env),
 # Rscale = 0,
  data   = dat_cond,
  method = "REML")

summary(mod_egger)

bubble_plot(mod_egger, mod = "nSE", ylab = "lnM", group = "Study")

# use yi_g_abs

mod_egger_g <- rma.mv(
  yi   = yi_g_abs,
  V    = vi_g_abs,
  mods = ~ nSE,
  random =  list(
    ~ 1 | Study,
    ~ 1 | Species,
    ~ 1 | Species.no.phylo,
    ~ 1 | ES.ID
  ),
 # R = list(Species = corMat.env),
 # Rscale = 0,
  data   = dat_cond,
  method = "REML")

summary(mod_egger_g)
bubble_plot(mod_egger_g, mod = "nSE", ylab = "abs(SMD)", group = "Study")

