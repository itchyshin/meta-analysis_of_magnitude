# example 3


################################################################################
## META-ANALYSES (multilevel: rma.mv) of
##   (A) |g|  – absolute Hedges’ g   (Dobler et al. 2014)
##   (B) lnM – ln geometric mean of control & treatment means (Inverts dataset)
## ---------------------------------------------------------------------------
## OUTPUT:
##   – list(|g| overall, |g|-by-kingdom, lnM overall, lnM-by-Stressor)  →  *.rds
##   – optional forest plots (commented ▹▹ if not needed)
################################################################################

library(tidyverse)   # data wrangling
library(data.table)  # fast fread()
library(metafor)     # rma.mv()
theme_set(theme_bw())

#### helper: folded-normal moments  -------------------------------------------
folded_norm <- function(mu, var){
  sigma <- sqrt(var)
  if (sigma == 0) return(c(mean = abs(mu), var = 0))
  z  <- mu / sigma
  m  <- sigma * sqrt(2/pi) * exp(-0.5*z^2) + mu*(2*pnorm(z) - 1)
  v  <- mu^2 + var - m^2
  c(mean = m, var = v)
}

################################################################################
## (A)  ABSOLUTE HEDGES’ g  – final_raw_data.csv
################################################################################
raw_g <- fread("final_raw_data.csv")

## calculate |g| on folded-normal scale
g_dat <- raw_g %>%
  mutate(across(c(g, gvar), as.numeric)) %>%
  mutate(
    tmp = map2(g, gvar, folded_norm),
    yi  = map_dbl(tmp, "mean"),
    vi  = map_dbl(tmp, "var")
  ) %>% select(-tmp) %>%
  mutate(ID = row_number())            # unique-effect identifier

## ── rma.mv models ────────────────────────────────────────────────────────────
## overall (random: publication + residual)
ma_g_overall <- rma.mv(yi, vi,
                       random = ~ 1 | study_nr/ID,
                       data   = g_dat,
                       method = "REML")

## by kingdom (animals vs plants) – same random structure inside each subset
ma_g_kingdom <- g_dat %>%
  group_split(kingdom) %>%
  set_names(unique(g_dat$kingdom)) %>%
  map(~ rma.mv(yi, vi,
               random = ~ 1 | study_nr/ID,
               data   = .x,
               method = "REML"))

################################################################################
## (B)  lnM  – DataInverts.csv  (metabolism example)
################################################################################
inv0 <- fread("DataInverts.csv")

inv <- inv0 %>%
  filter(across(c(ControlMeanTransformed, TreatmentMeanTransformed,
                  ControlSE, TreatmentSE), ~ is.finite(.x) & .x > 0)) %>%
  mutate(
    lnM     = 0.5*(log(ControlMeanTransformed) + 
                     log(TreatmentMeanTransformed)),
    var_lnM = 0.25*( (ControlSE/ControlMeanTransformed)^2 +
                       (TreatmentSE/TreatmentMeanTransformed)^2 ),
    PaperID = interaction(PaperAuthor, PaperYear, drop = TRUE),
    ID      = row_number()
  )

## overall multilevel model
ma_lnM_overall <- rma.mv(lnM, var_lnM,
                         random = ~ 1 | PaperID/ID,
                         data   = inv,
                         method = "REML")

## by Stressor (T, pH, TpH)  – same random structure per subset
ma_lnM_stressor <- inv %>%
  group_split(Stressor) %>%
  set_names(unique(inv$Stressor)) %>%
  map(~ rma.mv(lnM, var_lnM,
               random = ~ 1 | PaperID/ID,
               data   = .x,
               method = "REML"))

################################################################################
## SAVE ALL MODELS  -----------------------------------------------------------
saveRDS(list(
  abs_g_overall   = ma_g_overall,
  abs_g_kingdom   = ma_g_kingdom,
  lnM_overall     = ma_lnM_overall,
  lnM_stressor    = ma_lnM_stressor),
  file = "meta_multilevel_results.rds")

################################################################################
## OPTIONAL QUICK PLOTS  (uncomment if desired) -------------------------------
# pdf("forest_abs_g_overall.pdf", 7, 9)
# forest(ma_g_overall, slab = paste(raw_g$author, raw_g$year), xlab = "|g|")
# dev.off()
#
# pdf("forest_lnM_overall.pdf", 7, 9)
# forest(ma_lnM_overall, slab = paste(inv$PaperAuthor, inv$PaperYear), xlab = "ln M")
# dev.off()
################################################################################

## Session-info summary (good practice for reproducibility)
# sessionInfo()