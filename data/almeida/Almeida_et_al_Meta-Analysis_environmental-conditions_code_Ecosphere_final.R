#######################################################################
######## Meta-Analysis of Legacies of Nutritional Restriction ##########
#######################################################################

### Effects of environmental conditions (magnitude, duration, recovery) ###

#R version 4.0.3




#29 March 2021
#The following is the R file used to run meta-regressions related to the effect of environmental conditions (magnitude, duration, recovery) 
## on the absolute value of the effect size for the full dataset as well as for the Condition-only, Reproduction-only, and Insect-only analyses.
#This includes running the meta-regression as well as plotting the results.


#Data files needed: ####
#  Data extracted from literature:   Almeida_et_al_Data_Ecosphere_final.csv
#  Correlations between species for full dataset where magnitude, duration, and recovery were reported:   Almeida_et_al_Phylo_correlations_update_fulldataset_env.csv
#  Covariance matrix for shared controls for full dataset where magnitude, duration, and recovery were reported:   Almeida_et_al_covarianceCC_fulldataset_env.csv

#  Correlations between species for condition-only dataset where magnitude, duration, and recovery were reported:   Almeida_et_al_Phylo_correlations_update_conddata_env.csv
#  Covariance matrix for shared controls for condition-only dataset where magnitude, duration, and recovery were reported:   Almeida_et_al_covarianceCC_conddata_env.csv

#  Correlations between species for reproduction-only dataset where magnitude, duration, and recovery were reported:   Almeida_et_al_Phylo_correlations_update_reprodata_env.csv
#  Covariance matrix for shared controls for reproduction-only dataset where magnitude, duration, and recovery were reported:   Almeida_et_al_covarianceCC_reprodata_env.csv

#  Correlations between species for insect-only dataset where magnitude, duration, and recovery were reported:   Almeida_et_al_Phylo_correlations_update_insectdata_env.csv
#  Covariance matrix for shared controls for insect-only dataset where magnitude, duration, and recovery were reported:   Almeida_et_al_covarianceCC_insectdata_env.csv







#Libraries needed
library(tidyverse) #used to arranged and rearrange data
library(rstudioapi) #used to set working directory
library(metafor) #used for meta-regressions
library(ggpubr) #used to plot figures together in different pannels
library(cowplot) #used to plot figures together in different pannels




setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets working directory to the same location that this R script is saved in













# Read in data
meta = read.csv("Almeida_et_al_Data_Ecosphere_final.csv",sep=",",header=TRUE)
head(meta)

# name variables
nTrt = meta$nTrt
nCtrl = meta$nCtrl
CtrlResponse = meta$CtrlResponse
TrtResponse = meta$TrtResponse
TrtSD = meta$TrtSD
CtrlSD = meta$CtrlSD

Category = meta$Category
Sub.cat = meta$Sub_category
Sp.group = meta$Sp_group
age = meta$age_treated
sex = meta$sex
Trt.bad = meta$Trt.bad
magnitude = 1 - meta$Trt_magnitude #Subtracting the value in the dataframe by 1 allows for greater differences between treatment and control to be represented as a larger magnitude
duration = meta$Trt_duration
Direction = meta$LoH
Species = meta$Species
Species.no.phylo = meta$Species
Recovery = meta$Recovery

Study = (meta$Study)
same.trt = (meta$same_trt)
same.ctrl = meta$Shared.Control

ES.ID = as.character(meta$ES.ID)


#assigning names that fit with escalc - how Hedge's g will be calculated below
m1i <- vector(length = length(CtrlResponse))
m2i <- vector(length = length(TrtResponse))
sd1i <- vector(length = length(CtrlResponse))
sd2i <- vector(length = length(TrtResponse))
n1i <- vector(length = length(CtrlResponse))
n2i <- vector(length = length(TrtResponse))


# Adjust for studies in which the lower value is more advantageous
for (i in 1:length(CtrlResponse)) {
  if(Direction[i] == "higher") {
    m1i[i] = TrtResponse[i]
  } else {
    m1i[i] = CtrlResponse[i]
    }
}

for (i in 1:length(CtrlResponse)) {
  if(Direction[i] == "higher") {
    m2i[i] = CtrlResponse[i]
  } else {
    m2i[i] = TrtResponse[i]
  }
}


for (i in 1:length(CtrlResponse)) {
  if(Direction[i] == "higher") {
    sd1i[i] = TrtSD[i]
  } else {
    sd1i[i] = CtrlSD[i]
  }
}

for (i in 1:length(CtrlResponse)) {
  if(Direction[i] == "higher") {
    sd2i[i] = CtrlSD[i]
  } else {
    sd2i[i] = TrtSD[i]
  }
}


for (i in 1:length(CtrlResponse)) {
  if(Direction[i] == "higher") {
    n1i[i] = nTrt[i]
  } else {
    n1i[i] = nCtrl[i]
  }
}

for (i in 1:length(CtrlResponse)) {
  if(Direction[i] == "higher") {
    n2i[i] = nCtrl[i]
  } else {
    n2i[i] = nTrt[i]
  }
}





#####################################################################################################################################################
################################################# Calculating Hedge's g & creating relevant matrices ################################################
#####################################################################################################################################################


All.data <- data.frame(Category, age, sex, magnitude, duration, Recovery, Sub.cat, same.trt, same.ctrl, Study, Species, Species.no.phylo, ES.ID, Sp.group, m1i, m2i, sd1i, sd2i, n1i, n2i, nCtrl)


Cond.only <- filter(All.data, Category == "Condition")
Rep.only <- filter(All.data, Category == "Reproduction")

Insect.only <- filter(All.data, Sp.group == "Insect")





#*************************************************************************************************************
# Calculate Hedge's g effect size for full dataset ####
#*************************************************************************************************************

hg <- escalc(measure = "SMD", m1i = m1i, m2i = m2i, sd1i = sd1i, sd2i = sd2i, n1i = n1i, n2i = n2i, slab = Study)

All.data.ES <- All.data %>%
  mutate(yi = hg$yi) %>%
  mutate(vi = hg$vi)

#write.csv(All.data.ES, "Cov_Common_Control.csv", row.names = FALSE)

#*************************************************************************************************************


#*************************************************************************************************************
# Calculate Hedge's g effect size for Condition dataset ####
#*************************************************************************************************************

hg.cond <- escalc(measure = "SMD", m1i = m1i, m2i = m2i, sd1i = sd1i, sd2i = sd2i, n1i = n1i, n2i = n2i, slab = Study, data = Cond.only)

Cond.data.ES <- Cond.only %>%
  mutate(yi = hg.cond$yi) %>%
  mutate(vi = hg.cond$vi)

#write.csv(All.data.ES, "Cov_Common_Control.csv", row.names = FALSE)

#*************************************************************************************************************


#*************************************************************************************************************
# Calculate Hedge's g effect size for Reproduction dataset ####
#*************************************************************************************************************

hg.repro <- escalc(measure = "SMD", m1i = m1i, m2i = m2i, sd1i = sd1i, sd2i = sd2i, n1i = n1i, n2i = n2i, slab = Study, data = Rep.only)

Repro.data.ES <- Rep.only %>%
  mutate(yi = hg.repro$yi) %>%
  mutate(vi = hg.repro$vi)

#write.csv(All.data.ES, "Cov_Common_Control.csv", row.names = FALSE)

#*************************************************************************************************************


#*************************************************************************************************************
# Calculate Hedge's g effect size for Insect dataset ####
#*************************************************************************************************************

hg.insect <- escalc(measure = "SMD", m1i = m1i, m2i = m2i, sd1i = sd1i, sd2i = sd2i, n1i = n1i, n2i = n2i, slab = Study, data = Insect.only)

Insect.data.ES <- Insect.only %>%
  mutate(yi = hg.insect$yi) %>%
  mutate(vi = hg.insect$vi)

#write.csv(All.data.ES, "Cov_Common_Control.csv", row.names = FALSE)

#*************************************************************************************************************











#*************************************************************************************************************
# Add matricies with phylogenetic correlations, common control variance-covariance, and multiple responses correlations ####
#*************************************************************************************************************

#----------------- Full dataset ----------------

#dataframe with only values that have magnitude, duration, and recovery reported
env_data_full <- All.data.ES %>%
  filter(!is.na(magnitude)) %>%
  filter(!is.na(duration)) %>%
  filter(!is.na(Recovery))

#correlation matrix specific to species included in the magnitude dataframe
corMat.env <- read.csv("Almeida_et_al_Phylo_correlations_update_fulldataset_env.csv",sep=",",header=TRUE, row.names = 1)

#Input variance-covariance matrix describing covariance between treatments that share a common control
cov.cc.env = read.csv("Almeida_et_al_covarianceCC_fulldataset_env.csv",sep=",",header=FALSE)
diag(cov.cc.env) = env_data_full$vi

#Create correlation matrix for assumption that all responses that come from the same treatment (multiple responses) are correlated at 0.5
mult.resp.env <- matrix(data = 0, nrow = length(env_data_full$ES.ID), ncol = length(env_data_full$ES.ID))

for (i in 1:length(env_data_full$ES.ID)) {
  for (j in 1: length(env_data_full$ES.ID)) {
    if (env_data_full$same.trt[i] == env_data_full$same.trt[j]) {
      mult.resp.env[i, j] = 0.5
    } 
  }
}

colnames(mult.resp.env) = env_data_full$ES.ID
row.names(mult.resp.env) = env_data_full$ES.ID

diag(mult.resp.env) <- c(rep(1, length(env_data_full$yi)))



#----------------- Condition dataset ----------------
#dataframe with only values that have magnitude, duration, and recovery reported
env_data_cond <- Cond.data.ES %>%
  filter(!is.na(magnitude)) %>%
  filter(!is.na(duration)) %>%
  filter(!is.na(Recovery))

#correlation matrix specific to species included in the magnitude dataframe
corMat.condenv <- read.csv("Almeida_et_al_Phylo_correlations_update_conddata_env.csv",sep=",",header=TRUE, row.names = 1)

#Input variance-covariance matrix describing covariance between treatments that share a common control
cov.cc.condenv = read.csv("Almeida_et_al_covarianceCC_conddata_env.csv",sep=",",header=FALSE)
diag(cov.cc.condenv) = env_data_cond$vi

#Create correlation matrix for assumption that all responses that come from the same treatment (multiple responses) are correlated at 0.5
mult.resp.condenv <- matrix(data = 0, nrow = length(env_data_cond$ES.ID), ncol = length(env_data_cond$ES.ID))

for (i in 1:length(env_data_cond$ES.ID)) {
  for (j in 1: length(env_data_cond$ES.ID)) {
    if (env_data_cond$same.trt[i] == env_data_cond$same.trt[j]) {
      mult.resp.condenv[i, j] = 0.5
    } 
  }
}

colnames(mult.resp.condenv) = env_data_cond$ES.ID
row.names(mult.resp.condenv) = env_data_cond$ES.ID

diag(mult.resp.condenv) <- c(rep(1, length(env_data_cond$yi)))




#----------------- Reproduction dataset ----------------
#dataframe with only values that have magnitude, duration, and recovery reported
env_data_repro <- Repro.data.ES %>%
  filter(!is.na(magnitude)) %>%
  filter(!is.na(duration)) %>%
  filter(!is.na(Recovery))

#correlation matrix specific to species included in the magnitude dataframe
corMat.repenv <- read.csv("Almeida_et_al_Phylo_correlations_update_reprodata_env.csv",sep=",",header=TRUE, row.names = 1)

#Input variance-covariance matrix describing covariance between treatments that share a common control
cov.cc.repenv = read.csv("Almeida_et_al_covarianceCC_reprodata_env.csv",sep=",",header=FALSE)
diag(cov.cc.repenv) = env_data_repro$vi

#Create correlation matrix for assumption that all responses that come from the same treatment (multiple responses) are correlated at 0.5
mult.resp.repenv <- matrix(data = 0, nrow = length(env_data_repro$ES.ID), ncol = length(env_data_repro$ES.ID))

for (i in 1:length(env_data_repro$ES.ID)) {
  for (j in 1: length(env_data_repro$ES.ID)) {
    if (env_data_repro$same.trt[i] == env_data_repro$same.trt[j]) {
      mult.resp.repenv[i, j] = 0.5
    } 
  }
}

colnames(mult.resp.repenv) = env_data_repro$ES.ID
row.names(mult.resp.repenv) = env_data_repro$ES.ID

diag(mult.resp.repenv) <- c(rep(1, length(env_data_repro$yi)))



#----------------- Insect dataset ----------------
#dataframe with only values that have magnitude, duration, and recovery reported
env_data_insect <- Insect.data.ES %>%
  filter(!is.na(magnitude)) %>%
  filter(!is.na(duration)) %>%
  filter(!is.na(Recovery))

#correlation matrix specific to species included in the magnitude dataframe
corMat.insectenv <- read.csv("Almeida_et_al_Phylo_correlations_update_insectdata_env.csv",sep=",",header=TRUE, row.names = 1)

#Input variance-covariance matrix describing covariance between treatments that share a common control
cov.cc.insectenv = read.csv("Almeida_et_al_covarianceCC_insectdata_env.csv",sep=",",header=FALSE)
diag(cov.cc.insectenv) = env_data_insect$vi

#Create correlation matrix for assumption that all responses that come from the same treatment (multiple responses) are correlated at 0.5
mult.resp.insectenv <- matrix(data = 0, nrow = length(env_data_insect$ES.ID), ncol = length(env_data_insect$ES.ID))

for (i in 1:length(env_data_insect$ES.ID)) {
  for (j in 1: length(env_data_insect$ES.ID)) {
    if (env_data_insect$same.trt[i] == env_data_insect$same.trt[j]) {
      mult.resp.insectenv[i, j] = 0.5
    } 
  }
}

colnames(mult.resp.insectenv) = env_data_insect$ES.ID
row.names(mult.resp.insectenv) = env_data_insect$ES.ID

diag(mult.resp.insectenv) <- c(rep(1, length(env_data_insect$yi)))




#####################################################################################################################################################





























#########################################################################################################################
############################################## Full dataset analyses ####################################################
#########################################################################################################################


############# Examining continuous moderators (characteristics of the stressful environment) on absolute value of latent effect ###############






library(car)


# --------- Model --------

# env_full <- rma.mv(abs(yi) ~ logit(magnitude) + logit(duration) + logit(Recovery), cov.cc.env,
#                    random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.env, ES.ID = mult.resp.env), Rscale = 0, data = env_data_full, method = "REML")
# summary(env_full)

env_full <- rma.mv(log(abs(yi) + 1) ~ magnitude + duration + Recovery, cov.cc.env,
                   random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.env, ES.ID = mult.resp.env), Rscale = 0, data = env_data_full, method = "REML")
summary(env_full)



#Magnitude is significant
# Mag: estimate = 0.44, p = 0.0003
# Dur: estimate = -0.01, p = 0.919
# Rec: estimate = -0.13, p = 0.267


#fitted vs. residual values
plot(fitted(env_full), rstandard(env_full)$z, pch=19, xlab = "Fitted values", ylab = "Standardized residuals")
abline(h = 0, lty = 2)
qqnorm(rstandard(env_full)$z)
qqline(rstandard(env_full)$z)
hist(rstandard(env_full)$z, xlab = "Standardized residuals", ylab = "Frequency", main = NULL)
























#***************************************************************************** Figure 2b, C, S1.6a ********************************************************************************************




x11()
figenv1a <- ggplot() +
  geom_point(aes(x = magnitude, y = log(abs(yi) + 1)), env_data_full, size = 1) +
  geom_smooth(aes(x = magnitude, y = unlist(predict(env_full)$pred)), env_data_full, method = "lm", se = FALSE, size = 0.5, color = "red") + 
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv1a




x11()
figenv1b <- ggplot() +
  geom_point(aes(x = duration, y = log(abs(yi) + 1)), env_data_full, size = 1) +
  #geom_smooth(aes(x = duration, y = unlist(predict(env_full)$pred)), env_data_full, method = "lm", se = FALSE, size = 0.5, color = "red") + 
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv1b





x11()
figenv1c <- ggplot() +
  geom_point(aes(x = Recovery, y = log(abs(yi) + 1)), env_data_full, size = 1) +
  #geom_smooth(aes(x = Recovery, y = unlist(predict(env_full)$pred)), env_data_full, method = "lm", se = FALSE, size = 0.5, color = "red") + 
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv1c




#************************************************************************************************************************************************************************************







#------------------------------------------ Plotting the 3 together---------------------------------------------

library(ggpubr)
library(cowplot)


x11()
figureenv1 <- ggarrange(figenv1a, figenv1b, figenv1c, nrow = 1, ncol = 3, labels = c("A", "B", "C"), hjust = -5) 

figureenv1

# 
# ggsave("Figure_env_full.tif", plot = last_plot(), device = "tiff",
#        width = 6.5, height = 3.5, units = "in", dpi = 300)


















#####################################################################################################################################################
############################################################ Within categories of responses #########################################################
#####################################################################################################################################################










####################################################### Condition-only ####################################################################











############# Examining continuous moderators (characteristics of the stressful environment) on absolute value of latent effect ###############






# -------- Model --------


env_cond <- rma.mv(log(abs(yi) + 1) ~ magnitude + duration + Recovery, cov.cc.condenv,
                   random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.condenv, ES.ID = mult.resp.condenv), Rscale = 0, data = env_data_cond, method = "REML")
summary(env_cond)


#Magnitude is significant
# Mag: estimate = 0.43, p = 0.031
# Dur: estimate = 0.17, p = 0.257
# Rec: estimate = -0.07, p = 0.753



#fitted vs. residual values
plot(fitted(env_cond), rstandard(env_cond)$z, pch=19, xlab = "Fitted values", ylab = "Standardized residuals")
abline(h = 0, lty = 2)
qqnorm(rstandard(env_cond)$z)
qqline(rstandard(env_cond)$z)
hist(rstandard(env_cond)$z, xlab = "Standardized residuals", ylab = "Frequency", main = NULL)





















#***************************************************************************** Figure 2b, C, S1.6a ********************************************************************************************



x11()
figenv2a <- ggplot() +
  geom_point(aes(x = magnitude, y = log(abs(yi) + 1)), env_data_cond, size = 1) +
  geom_smooth(aes(x = magnitude, y = unlist(predict(env_cond)$pred)), env_data_cond, method = "lm", se = FALSE, size = 0.5, color = "red") + 
  xlab("") +
  ylab("log + 1 transformed\n|Standardized Mean Difference|") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv2a




x11()
figenv2b <- ggplot() +
  geom_point(aes(x = duration, y = log(abs(yi) + 1)), env_data_cond, size = 1) +
  #geom_smooth(aes(x = duration, y = unlist(predict(env_cond)$pred)), env_data_cond, method = "lm", se = FALSE, size = 0.5, color = "red") + 
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv2b





x11()
figenv2c <- ggplot() +
  geom_point(aes(x = Recovery, y = log(abs(yi) + 1)), env_data_cond, size = 1) +
  #geom_smooth(aes(x = Recovery, y = unlist(predict(env_cond)$pred)), env_data_cond, method = "lm", se = FALSE, size = 0.5, color = "red") + 
  xlab("") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv2c




#************************************************************************************************************************************************************************************







#------------------------------------------ Plotting the 3 together---------------------------------------------

library(ggpubr)
library(cowplot)


x11()
figureenv2 <- ggarrange(figenv2a, figenv2b, figenv2c, nrow = 1, ncol = 3, labels = c("D", "E", "F"), hjust = -5) 

figureenv2


# ggsave("Figure_env_cond.tif", plot = last_plot(), device = "tiff",
#        width = 6.5, height = 3.5, units = "in", dpi = 300)





































####################################################### Reproduction-only #######################################################







############# Examining continuous moderators (characteristics of the stressful environment) on absolute value of latent effect ###############







# ------- Model ------

env_repro <- rma.mv(log(abs(yi) + 1) ~ magnitude + duration + Recovery, cov.cc.repenv,
                   random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.repenv, ES.ID = mult.resp.repenv), Rscale = 0, data = env_data_repro, method = "REML")
summary(env_repro)

#Magnitude is the only significant factor (estimate = 0.48, p = 0.011)




#fitted vs. residual values
plot(fitted(env_repro), rstandard(env_repro)$z, pch=19, xlab = "Fitted values", ylab = "Standardized residuals")
abline(h = 0, lty = 2)
qqnorm(rstandard(env_repro)$z)
qqline(rstandard(env_repro)$z)
hist(rstandard(env_repro)$z, xlab = "Standardized residuals", ylab = "Frequency", main = NULL)


























#***************************************************************************** Figure 2b, C, S1.6a ********************************************************************************************



x11()
figenv3a <- ggplot() +
  geom_point(aes(x = magnitude, y = log(abs(yi) + 1)), env_data_repro, size = 1) +
  geom_smooth(aes(x = magnitude, y = unlist(predict(env_repro)$pred)), env_data_repro, method = "lm", se = FALSE, size = 0.5, color = "red") + 
  xlab("Magnitude") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv3a




x11()
figenv3b <- ggplot() +
  geom_point(aes(x = duration, y = log(abs(yi) + 1)), env_data_repro, size = 1) +
  #geom_smooth(aes(x = duration, y = unlist(predict(env_repro)$pred)), env_data_repro, method = "lm", size = 0.5) + 
  xlab("Duration") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv3b





x11()
figenv3c <- ggplot() +
  geom_point(aes(x = Recovery, y = log(abs(yi) + 1)), env_data_repro, size = 1) +
  #geom_smooth(aes(x = Recovery, y = unlist(predict(env_repro)$pred)), env_data_repro, method = "lm", size = 0.5) + 
  xlab("Recovery") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv3c





#************************************************************************************************************************************************************************************







#------------------------------------------ Plotting the 3 together---------------------------------------------

library(ggpubr)
library(cowplot)


x11()
figureenv3 <- ggarrange(figenv3a, figenv3b, figenv3c, nrow = 1, ncol = 3, labels = c("G", "H", "I"), hjust = c(-4, -4, -10) )

figureenv3


# ggsave("Figure_env_repro.tif", plot = last_plot(), device = "tiff",
#        width = 6.5, height = 3.5, units = "in", dpi = 300)





















































####################################################### Insect-only ####################################################################





############# Examining continuous moderators (characteristics of the stressful environment) on absolute value of latent effect ###############





# -------- Model --------


env_insect <- rma.mv(log(abs(yi) + 1) ~ magnitude + duration + Recovery, cov.cc.insectenv,
                   random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insectenv, ES.ID = mult.resp.insectenv), Rscale = 0, data = env_data_insect, method = "REML")
summary(env_insect)

#Magnitude is the only significant factor (estimate = 1.06, p = 0.02)




#fitted vs. residual values
plot(fitted(env_insect), rstandard(env_insect)$z, pch=19, xlab = "Fitted values", ylab = "Standardized residuals")
abline(h = 0, lty = 2)
qqnorm(rstandard(env_insect)$z)
qqline(rstandard(env_insect)$z)
hist(rstandard(env_insect)$z, xlab = "Standardized residuals", ylab = "Frequency", main = NULL)
















#***************************************************************************** Figure 2b, C, S1.6a ********************************************************************************************



x11()
figenv4a <- ggplot() +
  geom_point(aes(x = magnitude, y = log(abs(yi) + 1)), env_data_insect, size = 1) +
  geom_smooth(aes(x = magnitude, y = unlist(predict(env_insect)$pred)), env_data_insect, method = "lm", size = 0.5) + 
  xlab("Magnitude") +
  ylab("log + 1 transformed\n|Standardized Mean Difference|") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv4a




x11()
figenv4b <- ggplot() +
  geom_point(aes(x = duration, y = log(abs(yi) + 1)), env_data_insect, size = 1) +
  #geom_smooth(aes(x = duration, y = unlist(predict(env_insect)$pred)), env_data_insect, method = "lm", size = 0.5) + 
  xlab("Duration") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv4b





x11()
figenv4c <- ggplot() +
  geom_point(aes(x = Recovery, y = log(abs(yi) + 1)), env_data_insect, size = 1) +
  #geom_smooth(aes(x = Recovery, y = unlist(predict(env_insect)$pred)), env_data_insect, method = "lm", size = 0.5) + 
  xlab("Recovery") +
  ylab("") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10))

figenv4c




#************************************************************************************************************************************************************************************







#------------------------------------------ Plotting the 3 together---------------------------------------------


x11()
figureenv4 <- ggarrange(figenv4a, figenv4b, figenv4c, nrow = 1, ncol = 3, labels = c("A", "B", "C"), hjust = c(-6, -5, -6) )

figureenv4


# ggsave("Figure_env_insect2.tif", plot = last_plot(), device = "tiff",
#        width = 6.5, height = 3.5, units = "in", dpi = 300)




























#################################### Plotting all together #####################################


x11()
figureenv <- ggarrange(figureenv1, figureenv2, figureenv3, nrow = 3, ncol = 1)

figureenv




# ggsave("Figure_env_main3.tif", plot = last_plot(), device = "tiff",
#        width = 165, height = 171, units = "mm", dpi = 300)


 














