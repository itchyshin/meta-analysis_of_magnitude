# example 2

################################################################################
## Meta-analysis of invertebrate data
##   – absolute lnRR  (|ln(Treatment/Control)|, folded-normal expectation)
##   – lnM            (½·[ln Control + ln Treatment], i.e. ln geometric mean)
## Author: (adapt as needed)
################################################################################

#### 1.  Set-up ----------------------------------------------------------------
library(tidyverse)   # dplyr, purrr, ggplot2, …
library(data.table)  # fast fread()
library(metafor)     # rma()
theme_set(theme_bw())

#### 2.  Load the data ---------------------------------------------------------
dat0 <- fread("DataInverts.csv", na.strings = c("", "NA"))

## keep rows with complete information for means & SEs (>0 and finite)
dat  <- dat0 %>%
  filter(
    is.finite(ControlMeanTransformed)  & ControlMeanTransformed  > 0,
    is.finite(TreatmentMeanTransformed) & TreatmentMeanTransformed > 0,
    is.finite(ControlSE),  is.finite(TreatmentSE)
  )

#### 3.  lnRR (re-computed) ----------------------------------------------------
dat <- dat %>%
  mutate(
    lnRR      = if_else(`Metric Type` == "negative",
                        log(ControlMeanTransformed / TreatmentMeanTransformed),
                        log(TreatmentMeanTransformed / ControlMeanTransformed)),
    var_lnRR  = (ControlSE^2   / ControlMeanTransformed^2) +
      (TreatmentSE^2 / TreatmentMeanTransformed^2)
  )

#### 3a.  Absolute lnRR via folded-normal moments -----------------------------
folded_norm <- function(mu, var){
  sigma <- sqrt(var)
  if (sigma == 0) return(c(point = abs(mu), var = 0))
  z      <- mu / sigma
  meanFN <- sigma * sqrt(2/pi) * exp(-0.5 * z^2) + mu * (2*pnorm(z) - 1)
  varFN  <- mu^2 + var - meanFN^2
  c(point = meanFN, var = varFN)
}

dat <- dat %>%
  mutate(
    abs_lnRR = map2_dfr(lnRR, var_lnRR,
                        ~{ out <- folded_norm(.x, .y);
                        tibble(yi_abs = out["point"],
                               vi_abs = out["var"])})
  ) %>% unnest(abs_lnRR)

#### 4.  lnM (geometric-mean scale) -------------------------------------------
## lnM  := ½·[ln(m_C) + ln(m_T)]               -- Eq. (1)
## Var   ≈ ¼·[ (SE_C / m_C)^2 + (SE_T / m_T)^2 ] -- delta method
dat <- dat %>%
  mutate(
    lnM     = 0.5 * (log(ControlMeanTransformed) +
                       log(TreatmentMeanTransformed)),
    var_lnM = 0.25 * ( (ControlSE   / ControlMeanTransformed )^2 +
                         (TreatmentSE / TreatmentMeanTransformed)^2 )
  )

#### 5.  Random-effects meta-analyses (REML) -----------------------------------
## overall
ma_abs_lnRR_total <- rma(yi_abs, vi_abs, data = dat, method = "REML")
ma_lnM_total      <- rma(lnM,    var_lnM, data = dat, method = "REML")

## by stressor
stressors <- sort(unique(dat$Stressor))
ma_abs_lnRR_st <- map(stressors,
                      ~ rma(yi_abs, vi_abs,
                            data = filter(dat, Stressor == .x),
                            method = "REML")) |> set_names(stressors)

ma_lnM_st <- map(stressors,
                 ~ rma(lnM, var_lnM,
                       data = filter(dat, Stressor == .x),
                       method = "REML")) |> set_names(stressors)

#### 6.  Inspect results -------------------------------------------------------
print(summary(ma_abs_lnRR_total))
print(summary(ma_lnM_total))

#### 7.  Quick forest plot of |lnRR| ------------------------------------------
pdf("forest_abs_lnRR.pdf", width = 8, height = 10)
forest(ma_abs_lnRR_total,
       slab = paste(dat$PaperAuthor, dat$PaperYear, sep = ", "),
       xlab = "|ln response ratio|")
dev.off()

#### 8.  Save the fitted objects ----------------------------------------------
saveRDS(list(
  total_abs_lnRR   = ma_abs_lnRR_total,
  total_lnM        = ma_lnM_total,
  stressor_abs_lnRR = ma_abs_lnRR_st,
  stressor_lnM      = ma_lnM_st),
  file = "inverts_meta_results.rds")

## ──────────────────────────────────────────────────────────────────────────────
## The *.rds* object can be loaded with readRDS("inverts_meta_results.rds")
## and the PDF forest plot will appear in your working directory.
################################################################################