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
library(here)        # here()
theme_set(theme_bw())

source(here("R", "lnM_SAFE8.R"))   # defines lnM_delta1_indep()   &   safe_lnM_indep()


#### functions

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



#### 2.  Load the data ---------------------------------------------------------
dat0 <- read.csv(here("data", "DataInverts.csv"), na.strings = c("", "NA"))




