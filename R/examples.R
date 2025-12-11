# examples
# TODO 
# figs can be improved etc

rm(list = ls())

# packages
library(tidyverse)
library(metafor)        # escalc()
library(here)
library(orchaRd)
library(multcomp)
library(patchwork)

source(here("R", "lnM_SAFE8.R")) 

# TODO 
# we need to be able to deal with NA for get_ln_safe
# 
# get_lnM_safe_or_na <- function(m1, m2, s1, s2, n1, n2, B = 1e3) {
#   if (any(is.na(c(m1, m2, s1, s2, n1, n2)))) {
#     # Return a 1-row tibble of NAs
#     return(tibble(
#       yi_lnM_safe = NA_real_,
#       vi_lnM_safe = NA_real_,
#       draws_kept  = NA_integer_,
#       draws_total = NA_integer_
#     ))
#   }
#   get_lnM_safe(m1, m2, s1, s2, n1, n2, B = B)
# }

# TODO
# d_eq conversion to SAFE


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

# data - read in from CSV files

dat <- read_csv(here("data", "data_almeida.csv"), show_col_types = FALSE)

# effect size calculations 

dat <-  escalc(measure = "SMD", 
               m1i = m1i, 
               m2i = m2i, 
               sd1i = sd1i, 
               sd2i = sd2i, 
               n1i = n1i, 
               n2i = n2i,
               data = dat, append = TRUE,
               var.names = c("yi_d",   "vi_d"))

dat <- dat %>%                                   # <- your existing data frame
  mutate(
    abs_d   = map2_dfr(yi_d,   vi_d,   ~ {
      out <- folded_norm(.x, .y)
      tibble(yi_d_abs   = out["point"], vi_d_abs   = out["var"])
    })
  ) %>%                                           # flatten the tibbles
  unnest(c(abs_d))


# lnM
dat <- dat %>% mutate(
  lnM_raw = pmap_dfr(
    list(m1i, m2i, sd1i, sd2i, n1i, n2i),
    get_lnM_raw)
) %>% unnest(lnM_raw)

# SAFE is stochastic so set a seed
set.seed(123)

dat <- dat %>% mutate(
  lnM_safe = pmap_dfr(
    list(m1i, m2i, sd1i, sd2i, n1i, n2i),
    get_lnM_safe)
) %>% unnest(lnM_safe)

## we need to report this online - 1/3 of effect sizes are saved
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

print(summary_stats) # 270 effect sizes saved


# filtering for just condition data
# just using a part of data

dat_cond <- dat %>% filter(Category == "Condition")
dim(dat_cond)

# meta-analysis using abs(d)

ma_abs <- rma.mv(
  yi   = yi_d_abs, 
  V    = vi_d_abs,
  random =  list(
    ~ 1 | Study,
    ~ 1 | ES.ID
  ),
  #R = list(Species = corMat.env),
  #Rscale = 0,
  data   = dat_cond,
  method = "REML")

summary(ma_abs)
i2_ml(ma_abs)

p1 <- orchard_plot(ma_abs,  xlab = "abs(SMD)", group = "Study")

# meta-analysis using safe

ma_safe <- rma.mv(
  yi   = yi_lnM_safe, 
  V    = vi_lnM_safe,
  random =  list(
    ~ 1 | Study,
    ~ 1 | ES.ID
  ),
  #R = list(Species = corMat.env),
  #Rscale = 0,
  data   = dat_cond,
  method = "REML",
  test = "t")

summary(ma_safe)
i2_ml(ma_safe)
ml_m1(ma_safe)

p2 <- orchard_plot(ma_safe,  xlab = "lnM", group = "Study")


# N = 265 as some of missing values in moderators

# meta-regression abs(d)
mod_abs <- rma.mv(
  yi   = yi_d_abs, 
  V    = vi_d_abs,
  mods = ~ magnitude + duration + Recovery,
  random =  list(
    ~ 1 | Study,
    ~ 1 | ES.ID
  ),
  data   = dat_cond,
  method = "REML",
  test = "t")
summary(mod_abs)
# bubble plot
p3 <- bubble_plot(mod_abs, mod = "Recovery", ylab = "abs(SMD)", xlab = "Recovery days", group = "Study")

# meta-regression lnM
mod_safe <- rma.mv(
  yi   = yi_lnM_safe, 
  V    = vi_lnM_safe,
  mods = ~ magnitude + duration + Recovery,
  random =  list(
    ~ 1 | Study,
    ~ 1 | ES.ID
  ),
  data   = dat_cond,
  method = "REML",
  test = "t")

summary(mod_safe)
# bubble plot
p4 <- bubble_plot(mod_safe, mod = "Recovery", ylab = "lnM", xlab = "Recovery days", group = "Study")

# Egger regression 
# Creaing a variable using n0
dat_cond$n0 <- (dat_cond$n1i * dat_cond$n2i) / (dat_cond$n1i + dat_cond$n2i)
dat_cond$nSE <- 1 / sqrt(dat_cond$n0)
dat_cond$nV <- 1 / dat_cond$n0

# using abs(d)
egger_abs <- rma.mv(
  yi   = yi_d_abs, 
  V    = vi_d_abs,
  mods = ~ nSE,
  random =  list(
    ~ 1 | Study,
    ~ 1 | ES.ID
  ),
  data   = dat_cond,
  method = "REML",
  test = "t")
summary(egger_abs)
# bubble plot
p5 <- bubble_plot(egger_abs, mod = "nSE", ylab = "abs(SMD)", xlab = "sqrt(inverse sample size)", group = "Study")

# using safe
egger_safe <- rma.mv(
  yi   = yi_lnM_safe, 
  V    = vi_lnM_safe,
  mods = ~ nSE,
  random =  list(
    ~ 1 | Study,
    ~ 1 | ES.ID
  ),
  data   = dat_cond,
  method = "REML",
  test = "t")
summary(egger_safe)
# bubble plot
p6 <-bubble_plot(egger_safe, mod = "nSE", ylab = "lnM", xlab = "sqrt(inverse sample size)", group = "Study")


egger_safe2 <- rma.mv(
  yi   = yi_lnM_safe, 
  V    = vi_lnM_safe,
  mods = ~ nV,
  random =  list(
    ~ 1 | Study,
    ~ 1 | ES.ID
  ),
  data   = dat_cond,
  method = "REML",
  test = "t")
summary(egger_safe2)
# bubble plot
bubble_plot(egger_safe, mod = "nSE", ylab = "lnM", xlab = "sqrt(inverse sample size)", group = "Study")

# making 2 x 3 panel fig
(p1 + p2) / (p3 + p4) / (p5 + p6) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16, face = "bold"))


# alternative model

vtilde <- 1/dat_cond$n0

Vf   <- diag(as.numeric(vtilde))
levs <- levels(factor(dat_cond$ES.ID))
rownames(Vf) <- levs
colnames(Vf) <- levs

ma_alt <- rma.mv(
  yi   = yi_lnM_safe, 
  V    = vi_lnM_safe,
  random =  list(
    ~ 1 | Study,
    ~ 1 | ES.ID
  ),
  data   = dat_cond,
  method = "REML",
  test = "t",
    R      = list(ES.ID = Vf),
    Rscale = FALSE)

summary(ma_alt)

###########
# Example 2
###########

# load data
data<-read.csv(here("data", "kunc_2019.csv"), na.strings=c("","NA"))
data<-as.data.frame(data)
data.2 <-data[which(data$species.latin!="NA"),]                ### exclude cases where no species were provided
str(data.2)

#  make sd.control and sd.noise numeric

data.2$sd.control<-as.numeric(as.character(data.2$sd.control))
data.2$sd.noise<-as.numeric(as.character(data.2$sd.noise))

# exclude NA in dat$sd.control

data.2 <- data.2[which(!is.na(data.2$sd.control)),]

# hetro-verion of SMD
dat <- escalc(measure="SMDH", 
              m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
              m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
              data=data.2,
              append = TRUE,
              var.names = c("yi_d",   "vi_d"))

# fitering out NA for yi_d and vi_d
dat <- dat[which(!is.na(dat$yi_d)),]

dat <- dat %>%                                   # <- your existing data frame
  mutate(
    abs_d   = map2_dfr(yi_d,   vi_d,   ~ {
      out <- folded_norm(.x, .y)
      tibble(yi_d_abs   = out["point"], vi_d_abs   = out["var"])
    })
  ) %>%                                           # flatten the tibbles
  unnest(c(abs_d))

# lnM
dat <- dat %>% mutate(
  lnM_raw = pmap_dfr(
    list(mean.control, mean.noise, sd.control, sd.noise, sample.size.control, sample.size.noise.1),
    get_lnM_raw)
) %>% unnest(lnM_raw)

# SAFE is stochastic so set a seed
set.seed(123)

dat <- dat %>% mutate(
  lnM_safe = pmap_dfr(
    list(mean.control, mean.noise, sd.control, sd.noise, sample.size.control, sample.size.noise.1),
    get_lnM_safe)
) %>% unnest(lnM_safe)


## we need to report this online - 1/3 of effect sizes are saved
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

print(summary_stats) # recovered 177


# restrict data just for birds, mammals, amphibians and fish
dat<-dat[which(dat$taxon.for.plot!="reptilia" & 
                 dat$taxon.for.plot!="mollusca" & 
                 dat$taxon.for.plot!="arthropoda"),]


# meta-analysis using abs(d)
ma_abs <- rma.mv(
  yi   = yi_d_abs, 
  V    = vi_d_abs,
  random =  list(
    ~ 1 | study,
    ~ 1 |case.nr
  ),
  data   = dat,
  method = "REML",
  test = "t")
summary(ma_abs)
i2_ml(ma_abs)
p7 <- orchard_plot(ma_abs,  xlab = "abs(SMD)", group = "study")

# meta-analysis using safe
ma_safe <- rma.mv(
  yi   = yi_lnM_safe,
  V    = vi_lnM_safe,
  random =  list(
    ~ 1 | study,
    ~ 1 |case.nr
  ),
  data   = dat,
  method = "REML",
  test = "t")
summary(ma_safe)
i2_ml(ma_safe)
p8 <- orchard_plot(ma_safe,  xlab = "lnM", group = "study")

# meta-regression abs(d) using taxon.for.plot
mod_abs <- rma.mv(
  yi   = yi_d_abs, 
  V    = vi_d_abs,
  mods = ~ taxon.for.plot - 1,
  random =  list(
    ~ 1 | study,
    ~ 1 |case.nr
  ),
  data   = dat,
  method = "REML",
  test = "t")
summary(mod_abs)
# orchard plot
p9 <- orchard_plot(mod_abs, mod = "taxon.for.plot", xlab = "abs(SMD)", group = "study")
# multiple comparison
summary(glht(mod_abs, linfct=cbind(contrMat(rep(1,4), type="Tukey"))), test=adjusted("none"))

# meta-regression safe using taxon.for.plot
mod_safe <- rma.mv(
  yi   = yi_lnM_safe,
  V    = vi_lnM_safe,
  mods = ~ taxon.for.plot -1,
  random =  list(
    ~ 1 | study,
    ~ 1 |case.nr
  ),
  data   = dat,
  method = "REML",
  test = "t")
summary(mod_safe)
# orchard plot
p10 <- orchard_plot(mod_safe, mod = "taxon.for.plot", xlab = "lnM", group = "study")
# multiple comparison
summary(glht(mod_safe, linfct=cbind(contrMat(rep(1,4), type="Tukey"))), test=adjusted("none"))

# Egger regression
# Creaing a variable using n0
dat$n0 <- (dat$sample.size.control * dat$sample.size.noise.1
                 ) / (dat$sample.size.control + dat$sample.size.noise.1)
dat$nSE <- 1 / sqrt(dat$n0)
dat$nV <- 1 / dat$n0

# using abs(d)
egger_abs <- rma.mv(
  yi   = yi_d_abs, 
  V    = vi_d_abs,
  mods = ~ nSE,
  random =  list(
    ~ 1 | study,
    ~ 1 |case.nr
  ),
  data   = dat,
  method = "REML",
  test = "t")
summary(egger_abs)
# bubble plot
p11 <- bubble_plot(egger_abs, mod = "nSE", ylab = "abs(SMD)", xlab = "sqrt(inverse sample size)", group = "study")

# using safe
egger_safe <- rma.mv(
  yi   = yi_lnM_safe,
  V    = vi_lnM_safe,
  mods = ~ nSE,
  random =  list(
    ~ 1 | study,
    ~ 1 |case.nr
  ),
  data   = dat,
  method = "REML",
  test = "t")
summary(egger_safe)
# bubble plot
p12 <- bubble_plot(egger_safe, mod = "nSE", ylab = "lnM", xlab = "sqrt(inverse sample size)", group = "study")


egger_safe2 <- rma.mv(
  yi   = yi_lnM_safe,
  V    = vi_lnM_safe,
  mods = ~ nV,
  random =  list(
    ~ 1 | study,
    ~ 1 |case.nr
  ),
  data   = dat,
  method = "REML",
  test = "t")
summary(egger_safe2)
# bubble plot
bubble_plot(egger_safe2, mod = "nV", ylab = "lnM", xlab = "sqrt(inverse sample size)", group = "study")

# making 2 x 3 panel fig
(p7 + p8) / (p9 + p10) / (p11 + p12) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16, face = "bold"))


# alternative model

vtilde <- 1/dat$n0 

Vf   <- diag(as.numeric(vtilde))
levs <- levels(factor(dat$case.nr ))
rownames(Vf) <- levs
colnames(Vf) <- levs

ma_alt2 <- rma.mv(
  yi   = yi_lnM_safe, 
  V    = vi_lnM_safe,
  random =  list(
    ~ 1 | study,
    ~ 1 |case.nr
  ),
  data   = dat,
  method = "REML",
  test = "t",
  R      = list(case.nr = Vf),
  Rscale = FALSE)

summary(ma_alt2)


