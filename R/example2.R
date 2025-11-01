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
library(orchaRd)    # bubble_plot()
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
dat <- read.csv(here("data", "DataInverts.csv"), na.strings = c("", "NA"))

# need to fix data mean and SE to  fix escalc - use ControlMeanTransformed ControlSE ControlN  TreatmentMeanTransformed TreatmentSE TreatmentN

dat <- dat %>%
  mutate(
    m1i  = TreatmentMeanTransformed,
    sd1i = TreatmentSE * sqrt(TreatmentN),   # convert SE to SD
    n1i  = TreatmentN,
    
    m2i  = ControlMeanTransformed,
    sd2i = ControlSE * sqrt(ControlN),       # convert SE to SD
    n2i  = ControlN
  ) %>%
  filter(!is.na(m1i) & !is.na(m2i) & n1i > 2 & n2i > 2)  %>%  # need means and n > 2 for both groups
 # also get rid of zeros in mean and sd - none in this dataset
 filter(m1i > 0 & m2i > 0 & sd1i > 0 & sd2i > 0)   # none in this dataset


dat <-  escalc(measure = "ROM", 
               m1i = m1i,
               m2i = m2i,
               sd1i = sd1i,
               sd2i = sd2i,
               n1i = n1i,
               n2i = n2i,
               data = dat, append = TRUE,
               var.names = c("yi_rom",   "vi_rom"))

dat <- dat %>%                                   # <- your existing data frame
  mutate(
    abs_rom   = map2_dfr(yi_rom,   vi_rom,   ~ {
      out <- folded_norm(.x, .y)
      tibble(yi_rom_abs   = out["point"], vi_rom_abs   = out["var"])
    })
  ) %>%                                           # flatten the tibbles
  unnest(c(abs_rom))

# lnM - delta
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

# meta-analysis using rma.mv using both RoM and lnM

res_rom_abs     <- rma.mv(yi = yi_rom_abs, V = vi_rom_abs, 
                      random = ~ 1 | PaperAuthor, 
                      data = dat,
                      spares = TRUE,
                      control = list(optimizer = "Nelder-Mead"))
                      

summary(res_rom_abs)


res_rom_abs1     <- rma.mv(yi = yi_rom_abs, V = vi_rom_abs, 
                       mod = ~ Habitat,
                       random = ~ 1 | PaperAuthor, 
                       data = dat,
                       spares = TRUE,
                       control = list(optimizer = "Nelder-Mead"))

summary(res_rom_abs1)


# meta-regression - PaperYear

res_rom_abs2     <- rma.mv(yi = yi_rom_abs, V = vi_rom_abs, 
                       mod = ~ PaperYear,
                       random = ~ 1 | PaperAuthor, 
                       data = dat,
                       spares = TRUE,
                       control = list(optimizer = "Nelder-Mead"))

summary(res_rom_abs2)

# meta-regression - CO2Treatment

res_rom_abs3     <- rma.mv(yi = yi_rom_abs, V = vi_rom_abs, 
                       mod = ~ CO2Treatment - 1,
                       random = ~ 1 | PaperAuthor, 
                       data = dat,
                       spares = TRUE,
                       control = list(optimizer = "Nelder-Mead"))


summary(res_rom_abs3)


# use lmM - the same models but use safe

dat$yi_lnM_safe2 <- as.numeric(dat$yi_lnM_safe + log(sqrt(2))) 

res_lnM_safe <- rma.mv(yi = yi_lnM_safe2, 
                           V = vi_lnM_safe, 
                           random = ~ 1 | PaperAuthor, 
                           data = dat,
                           sparse= TRUE)

summary(res_lnM_safe)

orchard_plot(res_lnM_safe
             ,  xlab = "lnM (SAFE)", group = "PaperAuthor")

res_lnM_safe1     <- rma.mv(yi = yi_lnM_safe2, 
                            V = vi_lnM_safe, 
                            mod = ~  Habitat - 1,
                            random = ~ 1 | PaperAuthor, data = dat,
                            sparse= TRUE)

summary(res_lnM_safe1)

orchard_plot(res_lnM_safe1, mod = "Habitat", xlab = "lnM (SAFE)", group = "PaperAuthor")

# meta-regression - PaperYear

res_lnM_safe2     <- rma.mv(yi = yi_lnM_safe, 
                            V = vi_lnM_safe, 
                            mod = ~  PaperYear,
                            random = ~ 1 | PaperAuthor, data = dat,
                            sparse= TRUE)

summary(res_lnM_safe2)


bubble_plot(res_lnres_lnM_safe2M_safe5, mod = "PaperYear", ylab = "lnM (SAFE)", group = "PaperAuthor")

# meta-regression - CO2Treatment

res_lnM_safe3     <- rma.mv(yi = yi_lnM_safe, 
                            V = vi_lnM_safe, 
                            mod = ~  CO2Treatment - 1,
                            random = ~ 1 | PaperAuthor, data = dat,
                            sparse= TRUE)

summary(res_lnM_safe3) # something is wrong with this one

orchard_plot(res_lnM_safe3, mod = "CO2Treatment", xlab = "lnM (SAFE)", group = "PaperAuthor")


# meta-regression - sqrt(vi_lnM_safe) - to check for publication bias

res_lnM_safe4     <- rma.mv(yi = yi_lnM_safe, 
                            V = vi_lnM_safe, 
                            mod = ~  sqrt(vi_lnM_safe),
                            random = ~ 1 | PaperAuthor, data = dat,
                            sparse= TRUE)

summary(res_lnM_safe4)


# meta-regression - effective n

dat$effective_n <- 1 / (1/dat$n1i + 1/dat$n2i)
dat$sqrt_inverse_n <- sqrt(1/dat$effective_n )


res_lnM_safe5     <- rma.mv(yi = yi_lnM_safe, 
                            V = vi_lnM_safe, 
                            mod = ~ sqrt_inverse_n,
                            random = ~ 1 | PaperAuthor, data = dat,
                            sparse= TRUE)

summary(res_lnM_safe5)

bubble_plot(res_lnM_safe5, mod = "sqrt_inverse_n", ylab = "lnM (SAFE)", group = "PaperAuthor")


# meta-regression - log(n1 + n2)
dat$log_n <- log(dat$n1i + dat$n2i)

res_lnM_safe6     <- rma.mv(yi = yi_lnM_safe, 
                            V = vi_lnM_safe, 
                            mod = ~ log_n,
                            random = ~ 1 | PaperAuthor, data = dat,
                            sparse= TRUE)

summary(res_lnM_safe6)

bubble_plot(res_lnM_safe6, mod = "log_n", ylab = "lnM (SAFE)", group = "PaperAuthor")



# visualization

plot(dat$yi_rom, dat$vi_rom)

# filger vi_rom mover 100

dat_filt <- dat %>% filter(vi_rom < 20)

plot(dat_filt$yi_rom, dat_filt$vi_rom)


# cor between yi_rom_abs and yi_lnM_safe

cor.test(order(dat$yi_rom_abs), order(dat$yi_lnM_safe))

plot(dat$yi_rom_abs, dat$yi_lnM_safe)

# filtered

plot(dat_filt$yi_rom_abs, dat_filt$yi_lnM_safe)

cor.test(dat_filt$yi_rom_abs, dat_filt$yi_lnM_safe)

