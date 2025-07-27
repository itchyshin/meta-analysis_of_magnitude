## ── 1.  Set-up ────────────────────────────────────────────────────────────────
library(tidyverse)
library(metafor)        # escalc()
library(here)
library(orchaRd)

source(here("R", "lnM_SAFE7.R"))   # defines lnM_delta1_indep()   &   safe_lnM_indep()

## ── 2.  Load the PROCEED summary table and keep rows with full information ───
dat0 <- read_csv(here("data","PROCEED_Rates_DB_v5_for_upload_2021-11-25.csv"),
                 show_col_types = FALSE)

dat  <- dat0 %>%                                         # keep only usable rows
  filter(across(c(n1_r, n2_r, mean1_r, mean2_r, sd1_r, sd2_r),
                \(x) is.finite(x) & x > 0))

# write data to file for later use
#write.csv(dat, here("data", "proceed_complete.csv"), row.names = FALSE)

## ── 3.  Log response ratio (ROM) and Hedges’ g ───────────────────────────────
dat <- escalc(measure = "ROM",
              m1i = mean1_r, sd1i = sd1_r, n1i = n1_r,
              m2i = mean2_r, sd2i = sd2_r, n2i = n2_r,
              data = dat, append = TRUE,
              var.names = c("yi_rom", "vi_rom"))

dat <- escalc(measure = "SMDH",                     # bias-corrected Hedges g
              m1i = mean1_r, sd1i = sd1_r, n1i = n1_r,
              m2i = mean2_r, sd2i = sd2_r, n2i = n2_r,
              data = dat, append = TRUE,
              var.names = c("yi_g",   "vi_g"))

dat <- escalc(measure = "CVR",
              m1i = mean1_r, sd1i = sd1_r, n1i = n1_r,
              m2i = mean2_r, sd2i = sd2_r, n2i = n2_r,
              data = dat, append = TRUE,
              var.names = c("yi_cvr", "vi_cvr"))

dat <- escalc(measure = "VR",
              m1i = mean1_r, sd1i = sd1_r, n1i = n1_r,
              m2i = mean2_r, sd2i = sd2_r, n2i = n2_r,
              data = dat, append = TRUE,
              var.names = c("yi_vr",  "vi_vr"))

## ── folded-normal helper ──────────────────────────────────────────────────
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

## ── apply to ROM and g ────────────────────────────────────────────────────
dat <- dat %>%                                   # <- your existing data frame
  mutate(
    abs_rom = map2_dfr(yi_rom, vi_rom, ~ {
      out <- folded_norm(.x, .y)
      tibble(yi_rom_abs = out["point"], vi_rom_abs = out["var"])
    }),
    abs_g   = map2_dfr(yi_g,   vi_g,   ~ {
      out <- folded_norm(.x, .y)
      tibble(yi_g_abs   = out["point"], vi_g_abs   = out["var"])
    })
  ) %>%                                           # flatten the tibbles
  unnest(c(abs_rom, abs_g))
# 
# ## ── sanity check ──────────────────────────────────────────────────────────
# dat %>%
#   select(es_ID,
#          yi_rom,       vi_rom,
#          yi_rom_abs,   vi_rom_abs,
#          yi_g,         vi_g,
#          yi_g_abs,     vi_g_abs) %>%
#   slice_head(n = 5)

## ── 4.  ln M  (non-SAFE delta-method)  ─────────────────────────────────────
get_lnM_raw <- function(m1, m2, s1, s2, n1, n2) {
  out <- lnM_delta1_indep(m1, m2, s1, s2, n1, n2)   # <- returns c(point, var)
  tibble(
    yi_lnM_raw = out["point"],
    vi_lnM_raw = out["var"]          # ← sampling variance here
  )
}

dat <- dat %>% mutate(
  lnM_raw = pmap_dfr(
    list(mean1_r, mean2_r, sd1_r, sd2_r, n1_r, n2_r),
    get_lnM_raw)
) %>% unnest(lnM_raw)


## ── 5.  SAFE ln M  (parametric bootstrap)  ─────────────────────────────────
##        B = 10 000 resamples is a good default.
get_lnM_safe <- function(m1, m2, s1, s2, n1, n2, B = 1e4) {
  out <- safe_lnM_indep(m1, m2, s1, s2, n1, n2, B = B)
  tibble(
    yi_lnM_safe = out$point,
    vi_lnM_safe = out$var,           # ← bootstrap sampling variance
    draws_kept  = out$kept,
    draws_total = out$total
  )
}

dat <- dat %>% mutate(
  lnM_safe = pmap_dfr(
    list(mean1_r, mean2_r, sd1_r, sd2_r, n1_r, n2_r),
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


# draw fig to compare lnM_raw and lnM_safe

ggplot(dat, aes(x = yi_lnM_raw, y = yi_lnM_safe)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Comparison of ln M (raw vs SAFE)",
       x = "ln M (raw)",
       y = "ln M (SAFE bootstrap)") +
  theme_minimal() +
  coord_fixed(ratio = 1)

# draw fig to compare lnRR and SMD

ggplot(dat, aes(x = yi_rom, y = yi_g)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  labs(title = "Comparison of lnRR and Hedges' g",
       x = "lnRR",
       y = "Hedges' g") +
  theme_minimal() +
  coord_fixed(ratio = 1)

# draw fig to compare lnRR and lM (safe)

ggplot(dat, aes(x = abs(yi_rom), y = yi_lnM_safe)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "green", linetype = "dashed") +
  labs(title = "Comparison of lnRR and ln M (SAFE bootstrap)",
       x = "lnRR",
       y = "ln M (SAFE bootstrap)") +
  theme_minimal() #+
  #coord_fixed(ratio = 1)

# draw fig to compare SMD and lM (safe)

ggplot(dat, aes(x = abs(yi_g), y = yi_lnM_safe)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "purple", linetype = "dashed") +
  labs(title = "Comparison of Hedges' g and ln M (SAFE bootstrap)",
       x = "Hedges' g",
       y = "ln M (SAFE bootstrap)") +
  theme_minimal() #+
  #coord_fixed(ratio = 1)

################
# meta-analysis
###############3

# centre geneeration

dat$c_gen <- as.numeric(scale(dat$generations, center = TRUE, scale = FALSE))
dat$c_year <- as.numeric(scale(dat$years, center = TRUE, scale = FALSE))
dat$ln_year <- log(dat$years)
dat$urban_disturbance[is.na(dat$urban_disturbance) == TRUE] <- "Natural"

# ordering factor Natural as reference level
dat$disturbance <- factor(dat$urban_disturbance, 
                                levels = c("Natural", 
                                           "Biotic", 
                                           "Habitat Mod",
                                           "Heterogeneity", 
                                           "Novel",
                                           "Social"),
                                ordered = TRUE)
# ordering In situe natural variation as reference and merge 2 self introduction type stuff
dat$driver2 <- factor(dat$driver, 
                                levels = c("In situ natural variation", 
                                           "Climate change",
                                           "Harvest",
                                           "Introduction",
                                           "Introduction of predator/prey/host/competitor",
                                           "Pollution",
                                           "Range expansion after introduction",
                                           "Self induced range expansion",
                                           "Self induced range or host expansion" 
                                           ),
                                labels = c("Natural", 
                                           "Climate",
                                           "Harvest",
                                           "Introduction",
                                           "Competition",
                                           "Pollution",
                                           "Expantion",
                                           "Expantion",
                                           "Expantion" 
                                ),
                                ordered = TRUE)
# exclude p084 (ref_ID)

dat <- dat %>%
  filter(ref_ID != "p084") %>% # exclude p084
  mutate(sys_ID = factor(sys_ID), 
         ref_ID = factor(ref_ID),
         es_ID = factor(es_ID))

# Meta-analysis of ln M (SAFE bootstrap) - use rma.mv - random = ~ 1 | es_ID and sys_ID

ma_lnM <- rma.mv(yi = yi_lnM_safe, 
                      V = vi_lnM_safe, 
                      random = list(~ 1 | es_ID,
                               ~ 1 | ref_ID,
                               ~ 1 | sys_ID), 
                      sparse = TRUE,
                      data = dat)

summary(ma_lnM)

# meta-regression using c_gen

mr_lnM1 <- rma.mv(yi = yi_lnM_safe, 
                          V = vi_lnM_safe, 
                          mods = ~ ln_year,
                          random = list(~ 1 | es_ID,
                                        ~ 1 | ref_ID,
                                        ~ 1 | sys_ID), 
                          sparse = TRUE,
                          test = "t",
                          method = "REML", 
                          data = dat)

summary(mr_lnM1)

bubble_plot(mr_lnM1, mod = "c_year", xlab = "centered_year", group = "sys_ID")


mr_lnM2 <- rma.mv(yi = yi_lnM_safe, 
                  V = vi_lnM_safe, 
                  mods = ~ disturbance*ln_year,
                  random = list(~ 1 | es_ID,
                                ~ 1 | ref_ID,
                                ~ 1 | sys_ID), 
                  sparse = TRUE,
                  test = "t",
                  method = "REML", 
                  data = dat)

summary(mr_lnM2)

orchard_plot(mr_lnM2,  mod = "disturbance", xlab = "lnM",  group = "sys_ID")
bubble_plot(mr_lnM2, mod = "ln_year", xlab = "centered_year", group = "sys_ID")



mr_lnM3 <- rma.mv(yi = yi_lnM_safe, 
                  V = vi_lnM_safe, 
                  mods = ~ driver2*ln_year,
                  random = list(~ 1 | es_ID,
                                ~ 1 | ref_ID,
                                ~ 1 | sys_ID), 
                  
                  sparse = TRUE,
                  test = "t",
                  method = "REML", 
                  data = dat)

summary(mr_lnM3)

orchard_plot(mr_lnM3,  mod = "driver2", xlab = "lnM",  group = "sys_ID")
bubble_plot(mr_lnM3, mod = "ln_year", xlab = "centered_year", group = "sys_ID", by = "driver2")



mr_lnM4<- rma.mv(yi = yi_lnM_safe, 
                  V = vi_lnM_safe, 
                  mods = ~ driver2 + ln_year,
                  random = list(~driver2 | es_ID,
                                ~ 1 | ref_ID,
                                ~ 1 | sys_ID), 
                  struct = "DIAG",
                  sparse = TRUE,
                  test = "t",
                  method = "REML", 
                  data = dat[is.na(dat$driver2) == FALSE, ])

summary(mr_lnM4)
orchard_plot(mr_lnM4,  mod = "driver2", xlab = "lnM",  group = "sys_ID")
########
# lnCVR
#########
# 
# mr_lnCVR3 <- rma.mv(yi = yi_cvr, 
#                      V = vi_cvr, 
#                      mods = ~ driver2 + ln_year,
#                      random = list(~ 1 | es_ID,
#                                    ~ 1 | ref_ID,
#                                    ~ 1 | sys_ID), 
#                      sparse = TRUE,
#                      data = dat)
# summary(mr_lnCVR3)
# 
# orchard_plot(mr_lnCVR3,  mod = "driver2", xlab = "lnCVR",  group = "sys_ID")
# bubble_plot(mr_lnCVR3, mod = "ln_year", xlab = "centered_year", group = "sys_ID", by = "driver2")
# 
# #####
# # lnVR
# ######
# 
# mr_lnVR3 <- rma.mv(yi = yi_vr, 
#                      V = vi_vr, 
#                      mods = ~ driver2 + ln_year,
#                      random = list(~ 1 | es_ID,
#                                    ~ 1 | ref_ID,
#                                    ~ 1 | sys_ID), 
#                      sparse = TRUE,
#                      data = dat)
# 
# summary(mr_lnVR3)
# 
# orchard_plot(mr_lnVR3,  mod = "driver2", xlab = "lnVR",  group = "sys_ID")
# bubble_plot(mr_lnVR3, mod = "ln_year", xlab = "centered_year", group = "sys_ID", by = "driver2")
# 
# 

# Meta-analysis of yi_rom_abs

ma_lnRR <- rma.mv(yi = yi_rom_abs, 
                      V = vi_rom_abs, 
                  random = list(~ 1 | es_ID,
                                ~ 1 | ref_ID,
                                ~ 1 | sys_ID), 
                      sparse = TRUE,
                      data = dat)

summary(ma_lnRR)

# meta-regression using c_gen

mr_lnRR1 <- rma.mv(yi = yi_rom_abs, 
                     V = vi_rom_abs, 
                     mods = ~ c_year,
                     random = list(~ 1 | es_ID,
                                   ~ 1 | ref_ID,
                                   ~ 1 | sys_ID), 
                     sparse = TRUE,
                     data = dat)

summary(ma_lnRR_gen)

mr_lnRR2 <- rma.mv(yi = yi_rom_abs, 
                          V = vi_rom_abs, 
                          mods = ~ c_gen,
                          random = list(~ 1 | es_ID,
                                        ~ 1 | ref_ID,
                                        ~ 1 | sys_ID), 
                          sparse = TRUE,
                          data = dat)


# Meta-analysis of yi_g_abs

ma_SMD <- rma.mv(yi = yi_g_abs, 
                      V = vi_g_abs, 
                      random = list(~ 1 | es_ID,
                               ~ 1 | ref_ID,
                               ~ 1 | sys_ID), 
                      sparse = TRUE,
                      data = dat)

summary(ma_SMD)

# meta-regression using c_gen

ma_SMD_gen <- rma.mv(yi = yi_g_abs, 
                      V = vi_g_abs, 
                      mods = ~ c_year,
                      random = list(~ 1 | es_ID,
                                    ~ 1 | ref_ID,
                                    ~ 1 | sys_ID), 
                      sparse = TRUE,
                      data = dat)

summary(ma_SMD_gen)


#####
# Bayesian analysis

library(brms)
library(cmdstanr)  # makes sure CmdStan is available

# tell brms to use cmdstanr (required for reduce_sum / threading)
options(brms.backend = "cmdstanr")

# how many chains will run **in parallel**
options(mc.cores = 2)   

vcv <- diag(dat$vi_lnM_safe)
rownames(vcv) <- dat$es_ID
colnames(vcv) <- dat$es_ID

bf_lnM <- bf(
  yi_lnM_safe ~ 1 + driver2 + ln_year +
                (1|gr(es_ID, cov = vcv)) + (1 | ref_ID) + (1 | sys_ID),
                sigma ~ 1 + driver2 
)


# prior

prior_lnM <- default_prior(bf_lnM, 
                           data=dat, 
                           data2=list(vcv=vcv),
                           family=gaussian())

prior_lnM$prior[11] = "constant(1)"


mod_lnM <- brm(
  formula  = bf_lnM,
  data     = dat,
  data2.   =list(vcv=vcv),
  family   = gaussian(),
  backend  = "cmdstanr",
  prior    = prior_lnM,
  ## --- parallel settings -----------------------------------------------
  chains   = 2,                  # 2 chains run *between*-chain-parallel
  core     = 2,
  threads  = threading(8),       # 4 OpenMP threads *within* each chain
  iter     = 10000, 
  warmup   = 5000,
  control  = list(adapt_delta = 0.95,
                 max_treedepth = 15)
)


summary(mod_lnM)

# save the model as rds file

saveRDS(mod_lnM, here("data", "mod_lnM.rds"))
