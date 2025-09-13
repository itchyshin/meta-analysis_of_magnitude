#######################################################################
######## Meta-Analysis of Legacies of Nutritional Restriction ##########
#######################################################################


### Effects on direction and magnitude ###

#R version 4.0.3




#  29 March 2021
#The following is the R file used to run meta-regressions for the full dataset as well as for the Condition-only, Reproduction-only, and Insect-only analyses.
#This includes evaluating sources of heterogeneity, comparing models with likelihood ratio tests, model diagnositcs, and examining publication bias.
#This does not include evaluating the effect of environmental conditions on absolute value of effect size - See: Almeida_et_al_Meta-Analysis_environmental-conditions_code_Ecosphere_final.R


#Data inputed throughout code includes:
#  Data extracted from literature:   Almeida_et_al_Data_Ecosphere_final.csv
#  Correlations between species for full dataset:   Almeida_et_al_Phylo_correlations_update_fulldataset.csv
#  Covariance matrix for shared controls for full dataset:   Almeida_et_al_covarianceCC_fulldataset.csv

#  Correlations between species for condition-only dataset:   Almeida_et_al_Phylo_correlations_update_conddata.csv
#  Covariance matrix for shared controls for condition-only dataset:   Almeida_et_al_covarianceCC_conddata.csv

#  Correlations between species for reproduction-only dataset:   Almeida_et_al_Phylo_correlations_update_reprodata.csv
#  Covariance matrix for shared controls for reproduction-only dataset:   Almeida_et_al_covarianceCC_reprodata.csv

#  Correlations between species for insect-only dataset:   Almeida_et_al_Phylo_correlations_update_insectdata.csv
#  Covariance matrix for shared controls for insect-only dataset:   Almeida_et_al_covarianceCC_insectdata.csv







#Libraries needed
library(tidyverse) #used to arranged and rearrange data
library(rstudioapi) #used to set working directory
library(metafor) #used for meta-regressions
library(orchaRd) #used for plotting effects of moderators in orchard plots (see https://github.com/itchyshin/orchard_plot for downloading)
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

#Input correlation matrix describing relationships between species
corMat = read.csv("Almeida_et_al_Phylo_correlations_update_fulldataset.csv",sep=",",header=TRUE, row.names = 1)

#Input variance-covariance matrix describing covariance between treatments that share a common control
cov.cc = read.csv("Almeida_et_al_covarianceCC_fulldataset.csv",sep=",",header=FALSE)
diag(cov.cc) = All.data.ES$vi


#Create correlation matrix for assumption that all responses that come from the same treatment (multiple responses) are correlated at 0.5
mult.resp <- matrix(data = 0, nrow = length(All.data.ES$ES.ID), ncol = length(All.data.ES$ES.ID))

for (i in 1:length(All.data.ES$ES.ID)) {
  for (j in 1: length(All.data.ES$ES.ID)) {
    if (All.data.ES$same.trt[i] == All.data.ES$same.trt[j]) {
      mult.resp[i, j] = 0.5
    } 
  }
}

colnames(mult.resp) = All.data.ES$ES.ID
row.names(mult.resp) = All.data.ES$ES.ID

diag(mult.resp) <- c(rep(1, length(All.data.ES$yi)))



#----------------- Condition dataset ----------------
#Input correlation matrix describing relationships between species
corMat.cond = read.csv("Almeida_et_al_Phylo_correlations_update_conddata.csv",sep=",",header=TRUE, row.names = 1)

#Input variance-covariance matrix describing covariance between treatments that share a common control
cov.cc.cond = read.csv("Almeida_et_al_covarianceCC_conddata.csv",sep=",",header=FALSE)
diag(cov.cc.cond) = Cond.data.ES$vi


#Create correlation matrix for assumption that all responses that come from the same treatment (multiple responses) are correlated at 0.5
mult.resp.cond <- matrix(data = 0, nrow = length(Cond.data.ES$ES.ID), ncol = length(Cond.data.ES$ES.ID))

for (i in 1:length(Cond.data.ES$ES.ID)) {
  for (j in 1: length(Cond.data.ES$ES.ID)) {
    if (Cond.data.ES$same.trt[i] == Cond.data.ES$same.trt[j]) {
      mult.resp.cond[i, j] = 0.5
    } 
  }
}

colnames(mult.resp.cond) = Cond.data.ES$ES.ID
row.names(mult.resp.cond) = Cond.data.ES$ES.ID

diag(mult.resp.cond) <- c(rep(1, length(Cond.data.ES$yi)))



#----------------- Reproduction dataset ----------------
#Input correlation matrix describing relationships between species
corMat.rep = read.csv("Almeida_et_al_Phylo_correlations_update_reprodata.csv",sep=",",header=TRUE, row.names = 1)

#Input variance-covariance matrix describing covariance between treatments that share a common control
cov.cc.rep = read.csv("Almeida_et_al_covarianceCC_reprodata.csv",sep=",",header=FALSE)
diag(cov.cc.rep) = Repro.data.ES$vi


#Create correlation matrix for assumption that all responses that come from the same treatment (multiple responses) are correlated at 0.5
mult.resp.rep <- matrix(data = 0, nrow = length(Repro.data.ES$ES.ID), ncol = length(Repro.data.ES$ES.ID))

for (i in 1:length(Repro.data.ES$ES.ID)) {
  for (j in 1: length(Repro.data.ES$ES.ID)) {
    if (Repro.data.ES$same.trt[i] == Repro.data.ES$same.trt[j]) {
      mult.resp.rep[i, j] = 0.5
    } 
  }
}

colnames(mult.resp.rep) = Repro.data.ES$ES.ID
row.names(mult.resp.rep) = Repro.data.ES$ES.ID

diag(mult.resp.rep) <- c(rep(1, length(Repro.data.ES$yi)))



#----------------- Insect dataset ----------------
#Input correlation matrix describing relationships between species
corMat.insect = read.csv("Almeida_et_al_Phylo_correlations_update_insectdata.csv",sep=",",header=TRUE, row.names = 1)

#Input variance-covariance matrix describing covariance between treatments that share a common control
cov.cc.insect = read.csv("Almeida_et_al_covarianceCC_insectdata.csv",sep=",",header=FALSE)
diag(cov.cc.insect) = Insect.data.ES$vi


#Create correlation matrix for assumption that all responses that come from the same treatment (multiple responses) are correlated at 0.5
mult.resp.insect <- matrix(data = 0, nrow = length(Insect.data.ES$ES.ID), ncol = length(Insect.data.ES$ES.ID))

for (i in 1:length(Insect.data.ES$ES.ID)) {
  for (j in 1: length(Insect.data.ES$ES.ID)) {
    if (Insect.data.ES$same.trt[i] == Insect.data.ES$same.trt[j]) {
      mult.resp.insect[i, j] = 0.5
    } 
  }
}

colnames(mult.resp.insect) = Insect.data.ES$ES.ID
row.names(mult.resp.insect) = Insect.data.ES$ES.ID

diag(mult.resp.insect) <- c(rep(1, length(Insect.data.ES$yi)))




#####################################################################################################################################################





























#########################################################################################################################
############################################## Full dataset analyses ####################################################
#########################################################################################################################


#////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#                                           Likelihood Ratio Tests                                              ####
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////////////////////////////



#Step 1 - Univariate analyses with categorical moderator variables - comparing model with one fixed effect to a null model with no fixed effects

category_LRT_full.mr <- rma.mv(yi ~ Category, cov.cc,
                            random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "ML", data = All.data.ES)
null_LRT_full.mr <- rma.mv(yi ~ 1, cov.cc,
                        random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "ML", data = All.data.ES)
cat.v.null_full.mr <- anova(category_LRT_full.mr, null_LRT_full.mr)
cat.v.null_full.mr
#         df       AIC       BIC      AICc     logLik     LRT   pval        QE 
# Full     9 2307.1636 2347.9020 2307.4310 -1144.5818                3678.1470 
# Reduced  5 2362.9759 2385.6084 2363.0645 -1176.4880 63.8124 <.0001 3972.5880 
p.all.mr <- vector(mode = "numeric", length = 4)
p.all.mr[1] <- cat.v.null_full.mr$pval



age_LRT_full.mr <- rma.mv(yi ~ age, cov.cc,
                       random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "ML", data = All.data.ES)
age.v.null_full.mr <- anova(age_LRT_full.mr, null_LRT_full.mr)
age.v.null_full.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full     6 2364.8002 2391.9592 2364.9244 -1176.4001               3962.5932 
# Reduced  5 2362.9759 2385.6084 2363.0645 -1176.4880 0.1757 0.6751 3972.5880  
p.all.mr[2] <- age.v.null_full.mr$pval



sex_LRT_full.mr <- rma.mv(yi ~ sex, cov.cc,
                       random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "ML", data = All.data.ES)
sex.v.null_full.mr <- anova(sex_LRT_full.mr, null_LRT_full.mr)
sex.v.null_full.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full     7 2366.7032 2398.3887 2366.8692 -1176.3516               3936.7148 
# Reduced  5 2362.9759 2385.6084 2363.0645 -1176.4880 0.2727 0.8725 3972.5880 
p.all.mr[3] <- sex.v.null_full.mr$pval



taxa_LRT_full.mr <- rma.mv(yi ~ Sp.group, cov.cc,
                        random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "ML", data = All.data.ES)
taxa.v.null_full.mr <- anova(taxa_LRT_full.mr, null_LRT_full.mr)
taxa.v.null_full.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full    14 2375.8570 2439.2280 2376.4858 -1173.9285               3909.1300 
# Reduced  5 2362.9759 2385.6084 2363.0645 -1176.4880 5.1189 0.8238 3972.5880
p.all.mr[4] <- taxa.v.null_full.mr$pval


#Adjusting p-value based on Holm's correction

p.new.mr <- p.adjust(p.all.mr, method = "holm", n = length(p.all.mr))
p.new.mr
# 1.83087e-12 1.00000e+00 1.00000e+00 1.00000e+00

#Category is significant, but none of the other categorical moderators









#####  Step 2 - Examining potential for interactions among categorical responses by comparing a model with an interaction with the univariate model with only category that was significant in step 1


##---- Checking for when multiple responses are correlated at 0.5
cat_x_age_LRT_full.mr <- rma.mv(yi ~ Category + Category:age, cov.cc,
                             random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "ML", data = All.data.ES)
catage.v.cat_full.mr <- anova(cat_x_age_LRT_full.mr, category_LRT_full.mr)
catage.v.cat_full.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full    14 2310.7752 2374.1461 2311.4040 -1141.3876               3588.2662 
# Reduced  9 2307.1636 2347.9020 2307.4310 -1144.5818 6.3883 0.2702 3678.1470
p.all.mr <- vector(mode = "numeric", length = 3)
p.all.mr[1] <- catage.v.cat_full.mr$pval


cat_x_sex_LRT_full.mr <- rma.mv(yi ~ Category + Category:sex, cov.cc,
                             random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "ML", data = All.data.ES)
catsex.v.cat_full.mr <- anova(cat_x_sex_LRT_full.mr, category_LRT_full.mr)
catsex.v.cat_full.mr
#         df       AIC       BIC      AICc    logLik     LRT   pval        QE
# Full    19 2305.4330 2391.4364 2306.5793 -1133.7165                3493.4672 
# Reduced  9 2307.1636 2347.9020 2307.4310 -1144.5818 21.7306 0.0165 3678.1470 
p.all.mr[2] <- catsex.v.cat_full.mr$pval


cat_x_taxa_LRT_full.mr <- rma.mv(yi ~ Category + Category:Sp.group, cov.cc,
                              random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "ML", data = All.data.ES)
cattaxa.v.cat_full.mr <- anova(cat_x_taxa_LRT_full.mr, category_LRT_full.mr)
cattaxa.v.cat_full.mr
#         df       AIC       BIC      AICc    logLik     LRT   pval        QE
# Full    29 2306.4934 2437.7618 2309.1580 -1124.2467                3512.6280 
# Reduced  9 2307.1636 2347.9020 2307.4310 -1144.5818 40.6701 0.0041 3678.1470 
p.all.mr[3] <- cattaxa.v.cat_full.mr$pval




#Adjusting p-value based on Holm's correction

p.new.mr <- p.adjust(p.all.mr, method = "holm", n = length(p.all.mr))
p.new.mr
#  0.27024417 0.03307506 0.01231757



#Interactions with sex and taxa improve the model...






























#******************************************************************************************************************************************************
#                                       Quantifying Heterogeneity                                        ####
#******************************************************************************************************************************************************

#Null model - no moderators
hetero_randonly <- rma.mv(yi ~ 1, cov.cc, 
                      random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "REML", data = All.data.ES)
summary(hetero_randonly)


sum(hetero_randonly$sigma2)/(sum(hetero_randonly$sigma2)+(sum(1/hg$vi)*(hetero_randonly$k-1)/(sum(1/hg$vi)^2-sum((1/hg$vi)^2))))*100 # total heterogeneity 95.48%
s2t.randonly <- sum(hetero_randonly$sigma2) + sum(1/hetero_randonly$vi) * (hetero_randonly$k-1) / (sum(1/hetero_randonly$vi)^2 - sum((1/hetero_randonly$vi)^2))
hetero_randonly$sigma2[1]/s2t.randonly*100 # I^2 between study (id) < 1.0%
hetero_randonly$sigma2[2]/s2t.randonly*100 # I^2 phylogeny = 66.7%
hetero_randonly$sigma2[3]/s2t.randonly*100 # I^2 species without phylogeny < 1.0%
hetero_randonly$sigma2[4]/s2t.randonly*100 # I^2 residual = 28.7


#Comparing with model with important moderators from LRTs
hetero_full <- rma.mv(yi ~ Category + Category:sex + Category:Sp.group, cov.cc, 
                      random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "REML", data = All.data.ES)
summary(hetero_full)


sum(hetero_full$sigma2)/(sum(hetero_full$sigma2)+(sum(1/hg$vi)*(hetero_full$k-1)/(sum(1/hg$vi)^2-sum((1/hg$vi)^2))))*100 # total heterogeneity 95.91%
s2t <- sum(hetero_full$sigma2) + sum(1/hetero_full$vi) * (hetero_full$k-1) / (sum(1/hetero_full$vi)^2 - sum((1/hetero_full$vi)^2))
hetero_full$sigma2[1]/s2t*100 # I^2 between study (id) = 5.0%
hetero_full$sigma2[2]/s2t*100 # I^2 phylogeny = 70.3%
hetero_full$sigma2[3]/s2t*100 # I^2 species without phylogeny < 1.0%
hetero_full$sigma2[4]/s2t*100 # I^2 residual = 20.6%




#Differences in heterogeneity between model runs

#Total
TotalH.diff = (sum(hetero_randonly$sigma2)/(sum(hetero_randonly$sigma2)+(sum(1/hg$vi)*(hetero_randonly$k-1)/(sum(1/hg$vi)^2-sum((1/hg$vi)^2))))*100) - (sum(hetero_full$sigma2)/(sum(hetero_full$sigma2)+(sum(1/hg$vi)*(hetero_full$k-1)/(sum(1/hg$vi)^2-sum((1/hg$vi)^2))))*100) #-0.43
StudyH.diff = (hetero_randonly$sigma2[1]/s2t.randonly*100) - (hetero_full$sigma2[1]/s2t*100) # -4.8
PhyloH.diff = (hetero_randonly$sigma2[2]/s2t.randonly*100) - (hetero_full$sigma2[2]/s2t*100) # -3.7
SppH.diff = (hetero_randonly$sigma2[3]/s2t.randonly*100) - (hetero_full$sigma2[3]/s2t*100)   # < 1.0
ESH.diff = (hetero_randonly$sigma2[4]/s2t.randonly*100) - (hetero_full$sigma2[4]/s2t*100)    # 8.1






# *******************************************************************************************************************************









#------------------------------------------------------------------ Final model ------------------------------------------------------------------
final_model_full <- rma.mv(yi ~ Category + Category:sex + Category:Sp.group, cov.cc,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "REML", data = All.data.ES)
summary(final_model_full)

#-------------------------------------------------------------------------------------------------------------------------------------------------



















#Random effects
Rand.Study.int.1 <- ranef(final_model_full)$Study[1]$intrcpt
Rand.Study.se.1 <- ranef(final_model_full)$Study[2]$se
Rand.Study.min <- Rand.Study.int.1 - Rand.Study.se.1
Rand.Study.max <- Rand.Study.int.1 + Rand.Study.se.1
Rand.Study.lab <- rownames(ranef(final_model_full)$Study[1])
nn <- All.data.ES %>% group_by(Study) %>% tally
faceting.Study <- c(rep(1, 21), rep(2, 20), rep(3, 20), rep(4, 20))
Rand.Study.1 <- data.frame(Rand.Study.lab, Rand.Study.int.1, Rand.Study.min, Rand.Study.max, nn, faceting.Study)
names(Rand.Study.1)[names(Rand.Study.1) == 'n'] <- "sample.size"
As <- filter(Rand.Study.1, faceting.Study == '1')
Gs <- filter(Rand.Study.1, faceting.Study == '2')
Rs <- filter(Rand.Study.1, faceting.Study == '3')
Zs <- filter(Rand.Study.1, faceting.Study == '4')



p_random = function(dat, x, y, ymin, ymax, nsample) {
  
  x11()
  
  # Convert x-variable to a factor
  dat[,x] = as.factor(dat[, x])
  
  # Plot points
  p = ggplot(dat, aes_string(x, y)) +
    geom_point(color = 'black') + 
    geom_hline(yintercept = 0) +
    geom_pointrange(aes_string(ymax = ymax, ymin = ymin), dat) +
    coord_flip() +
    ylab("Estimated Random Intercept ± SE") + xlab("") +
    theme_classic() +
    theme(axis.line.y = element_line(color = "white"), axis.ticks = element_blank(), axis.title.x = element_text(size = 32), axis.text = element_text(size = 28), legend.title = element_text(size = 32))
  
  
  # # Summarise data to get counts by x-variable
  # nn = dat2 %>% group_by(nvar) %>% tally
  
  
  # Add counts as text labels
  p + geom_text(data = dat, aes_string(label = paste0("n = ", nsample)),
                y = min(dat[, ymin]) + 0.15 * min(dat[, ymin]), colour = "grey20", size = 7) +
    scale_y_continuous(limits = c(min(dat[, ymin]) + 0.15 * min(dat[, ymin]), max(dat[, ymax])))
  
}


p_random(As, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")
p_random(Gs, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")
p_random(Rs, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")
p_random(Zs, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")







Rand.Species.int.1 <- ranef(final_model_full)$Species[1]$intrcpt
Rand.Species.se.1 <- ranef(final_model_full)$Species[2]$se
Rand.Species.min <- Rand.Species.int.1 - Rand.Species.se.1
Rand.Species.max <- Rand.Species.int.1 + Rand.Species.se.1
Rand.Species.lab <- rownames(ranef(final_model_full)$Species[1])
nn <- All.data.ES %>% group_by(Species) %>% tally
faceting.Species <- c(rep(1, 22), rep(2, 22), rep(3, 21))
Rand.Species.1 <- data.frame(Rand.Species.lab, Rand.Species.int.1, Rand.Species.min, Rand.Species.max, nn, faceting.Species)
names(Rand.Species.1)[names(Rand.Species.1) == 'n'] <- "sample.size"
Aspp <- filter(Rand.Species.1, faceting.Species == '1')
Lspp <- filter(Rand.Species.1, faceting.Species == '2')
Zspp <- filter(Rand.Species.1, faceting.Species == '3')

p_random(Aspp, "Rand.Species.lab", "Rand.Species.int.1", "Rand.Species.min", "Rand.Species.max", "sample.size")
p_random(Lspp, "Rand.Species.lab", "Rand.Species.int.1", "Rand.Species.min", "Rand.Species.max", "sample.size")
p_random(Zspp, "Rand.Species.lab", "Rand.Species.int.1", "Rand.Species.min", "Rand.Species.max", "sample.size")




































################################################### Model diagnostics #########################################################
#fitted vs. residual values
plot(fitted(final_model_full), rstandard(final_model_full)$z, pch=19, xlab = "Fitted values", ylab = "Standardized residuals")
abline(h = 0, lty = 2)
#Some potential concerning narrowing with higher fitted values
qqnorm(rstandard(final_model_full)$z)
qqline(rstandard(final_model_full)$z)
#potentially skewed
hist(rstandard(final_model_full)$z, xlab = "Standardized residuals", ylab = "Frequency", main = NULL)


profile.rma.mv(final_model_full, sigma2 = 1)
abline(h = logLik(final_model_full) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_full, sigma2 = 2)
abline(h = logLik(final_model_full) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_full, sigma2 = 3)
abline(h = logLik(final_model_full) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_full, sigma2 = 4)
abline(h = logLik(final_model_full) - qchisq(.95, df=1)/2, lty="dashed")

################################################################################################################################


























#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Plotting figures ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^





#----------------------------- Plotting fixed effects -----------------------------


final_model_full_plot <- rma.mv(yi ~ Category + Category:sex + Category:Sp.group - 1, cov.cc,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "REML", data = All.data.ES)
summary(final_model_full_plot)

Full_table <- mod_results(final_model_full_plot, mod = "Category")
data.extract.full <- as.data.frame(Full_table$data)
mod.extract.full <- as.data.frame(Full_table$mod_table)



category_effects_full <- rma.mv(yi ~ Category - 1, cov.cc,
                                random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "REML", data = All.data.ES)
table_catsolo <- mod_results(category_effects_full, mod = "Category")
data.extract.catsolo <- as.data.frame(table_catsolo$data)
mod.extract.catsolo <- as.data.frame(table_catsolo$mod_table)




mods_of_interest.catfull <- mod.extract.full %>%
  filter(name == "Condition" | name == "Learning" | name == "Longevity" | name == "Mobility" | name == "Reproduction")

data_of_interest.catfull <- data.extract.catsolo


relist.catfull <- list(mods_of_interest.catfull, data_of_interest.catfull)
names(relist.catfull) <- c("mod_table", "data")


colors.consistent <- c("salmon", "orange", "seagreen3", "steelblue1", "plum3")


#Panel A with just category
#x11()
orchard_catfull <- orchard_plot(relist.catfull, mod = "Category", xlab = "Standardized mean difference", alpha = 0.5, 
             transfm = "none", angle = 0, cb = FALSE) + 
  theme(legend.position = c(0.45, 0.87)) +
  scale_color_manual(values = colors.consistent) +
  scale_fill_manual(values = c(rep("black", 5)))

orchard_catfull





# ggsave("orchard_catful.tif", plot = last_plot(), device = "tiff",
#       width = 6.5, height = 3, units = "in", dpi = 300)








#++++++++ Interaction between category and sex ++++++++


final_model_full.sex <- rma.mv(yi ~ Category:sex + Category:Sp.group + Category - 1, cov.cc,
                               random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "REML", data = All.data.ES)
summary(final_model_full.sex)


Catsex_table <- mod_results(final_model_full.sex, mod = "Category:sex")
data.extract.catsex <- as.data.frame(Catsex_table$data)
mod.extract.catsex <- as.data.frame(Catsex_table$mod_table)



final_model_full.sex2 <- rma.mv(yi ~ Category:sex - 1, cov.cc,
                                random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "REML", data = All.data.ES)
summary(final_model_full.sex2)


Catsex2_table <- mod_results(final_model_full.sex2, mod = "Category:sex")
data.extract.catsex2 <- as.data.frame(Catsex2_table$data)
mod.extract.catsex2 <- as.data.frame(Catsex2_table$mod_table)




mods_of_interest.sex <- mod.extract.catsex %>%
  filter(name == "CategoryCondition:sexMale" | name == "CategoryCondition:sexMixed" |
           name == "CategoryLearning:sexMale" | name == "CategoryLearning:sexMixed" |
           name == "CategoryLongevity:sexMale" | name == "CategoryLongevity:sexMixed" |
           name == "CategoryMobility:sexMale" | name == "CategoryMobility:sexMixed" |
           name == "CategoryReproduction:sexMale" | name == "CategoryReproduction:sexMixed")

data_of_interest.sex <- data.extract.catsex2 %>%
  filter(moderator == "CategoryCondition:sexMale" | moderator == "CategoryCondition:sexMixed" |
           moderator == "CategoryLearning:sexMale" | moderator == "CategoryLearning:sexMixed" |
           moderator == "CategoryLongevity:sexMale" | moderator == "CategoryLongevity:sexMixed" |
           moderator == "CategoryMobility:sexMale" | moderator == "CategoryMobility:sexMixed" |
           moderator == "CategoryReproduction:sexMale" | moderator == "CategoryReproduction:sexMixed")


relist.catsex <- list(mods_of_interest.sex, data_of_interest.sex)
names(relist.catsex) <- c("mod_table", "data")


colors.consistent.catsex <- c("salmon", "orange", "seagreen3", "steelblue1", "plum3", 
                              "salmon", "orange", "seagreen3", "steelblue1", "plum3")



#x11()
orchard_catsex <- orchard_plot(relist.catsex, mod = "Category:sex", xlab = "Difference from\nreference\n(females)", alpha = 0.5, transfm = "none", angle = 0, cb = FALSE) + 
  #xlim(-15, 15) +
  coord_cartesian(xlim = c(-15, 15)) +
  scale_y_discrete(labels = c("Condition x Male", "Learning x Male", "Longevity x Male", "Mobility x Male", "Reproduction x Male", 
                              "Condition x Mixed", "Learning x  Mixed", "Longevity x  Mixed", "Mobility x  Mixed", "Reproduction x  Mixed")) +
  theme(legend.position = "none") +
  scale_color_manual(values = colors.consistent.catsex) +
  scale_fill_manual(values = c(rep("black", 10))) +
  annotate("text", x = -5, y = 4, label = "*", size = 10)

orchard_catsex


# ggsave("orchard_catsex.tif", plot = last_plot(), device = "tiff",
#        width = 3.75, height = 6.5, units = "in", dpi = 300)










#++++++++ Interaction between category and taxa ++++++++

All.data.ES$Sp.group <- relevel(as.factor(All.data.ES$Sp.group), ref="Insect")


final_model_full.taxa <- rma.mv(yi ~ Category:Sp.group + Category + Category:sex - 1, cov.cc,
                               random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "REML", data = All.data.ES)
summary(final_model_full.taxa)


Cattaxa_table <- mod_results(final_model_full.taxa, mod = "Category:taxa")
data.extract.cattaxa <- as.data.frame(Cattaxa_table$data)
mod.extract.cattaxa <- as.data.frame(Cattaxa_table$mod_table)



final_model_full.taxa2 <- rma.mv(yi ~ Category:Sp.group - 1, cov.cc,
                                random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat, ES.ID = mult.resp), Rscale = 0, method = "REML", data = All.data.ES)
summary(final_model_full.taxa2)


Cattaxa2_table <- mod_results(final_model_full.taxa2, mod = "Category:Sp.group")
data.extract.cattaxa2 <- as.data.frame(Cattaxa2_table$data)
mod.extract.cattaxa2 <- as.data.frame(Cattaxa2_table$mod_table)




mods_of_interest.taxa <- mod.extract.cattaxa %>%
  filter(name ==      "CategoryCondition:Sp.groupAmphibian" | name ==      "CategoryCondition:Sp.groupBird" | name ==    "CategoryCondition:Sp.groupCrustacean" | name ==    "CategoryCondition:Sp.groupEchinoderm" | name ==    "CategoryCondition:Sp.groupFish" | name ==    "CategoryCondition:Sp.groupGastropod" | name ==    "CategoryCondition:Sp.groupMammal" | name ==    "CategoryCondition:Sp.groupReptile" |  
           
           name ==      "CategoryLongevity:Sp.groupAmphibian" | name ==    "CategoryLongevity:Sp.groupCrustacean" | name ==    "CategoryLongevity:Sp.groupFish" | name ==    "CategoryLongevity:Sp.groupGastropod" | name ==    "CategoryLongevity:Sp.groupReptile" | name ==    "CategoryLongevity:Sp.groupRotifer" | 
           name ==     "CategoryMobility:Sp.groupAmphibian" | name ==     "CategoryMobility:Sp.groupBird" | name ==     "CategoryMobility:Sp.groupFish" | name ==     "CategoryMobility:Sp.groupReptile" |  
           name == "CategoryReproduction:Sp.groupBird" | name == "CategoryReproduction:Sp.groupFish")

data_of_interest.taxa <- data.extract.cattaxa2 %>%
  filter(moderator ==      "CategoryCondition:Sp.groupAmphibian" | moderator ==      "CategoryCondition:Sp.groupBird" | moderator ==    "CategoryCondition:Sp.groupCrustacean" | moderator ==    "CategoryCondition:Sp.groupEchinoderm" | moderator ==    "CategoryCondition:Sp.groupFish" | moderator ==    "CategoryCondition:Sp.groupGastropod" | moderator ==    "CategoryCondition:Sp.groupMammal" | moderator ==    "CategoryCondition:Sp.groupReptile" |  
           
           moderator ==      "CategoryLongevity:Sp.groupAmphibian" |moderator ==    "CategoryLongevity:Sp.groupCrustacean" | moderator ==    "CategoryLongevity:Sp.groupFish" | moderator ==    "CategoryLongevity:Sp.groupGastropod" | moderator ==    "CategoryLongevity:Sp.groupReptile" | moderator ==    "CategoryLongevity:Sp.groupRotifer" | 
           moderator ==     "CategoryMobility:Sp.groupAmphibian" | moderator ==     "CategoryMobility:Sp.groupBird" | moderator ==     "CategoryMobility:Sp.groupFish" | moderator ==     "CategoryMobility:Sp.groupReptile" |  
           moderator == "CategoryReproduction:Sp.groupBird" | moderator == "CategoryReproduction:Sp.groupFish")



relist.cattaxa <- list(mods_of_interest.taxa, data_of_interest.taxa)
names(relist.cattaxa) <- c("mod_table", "data")


colors.consistent.cattaxa <- c("salmon", "seagreen3", "steelblue1",
                               "salmon", "steelblue1", "plum3", 
                               "salmon", "seagreen3", 
                               "salmon", 
                               "salmon", "seagreen3", "steelblue1", "plum3", 
                               "salmon", "seagreen3", 
                               "salmon",
                               "salmon", "steelblue1",
                               "seagreen3")


x11()
orchard_cattaxa <- orchard_plot(relist.cattaxa, mod = "Category:Sp.group", xlab = "Difference from\nreference\n(insects)", alpha = 0.5, transfm = "none", angle = 0, cb = FALSE) + 
  #xlim(-10, 10) +
  coord_cartesian(xlim = c(-10, 10)) +
  scale_y_discrete(labels = c("Condition x Amphibian", "Longevity x Amphibian", "Mobility x Amphibian",
                              "Condition x Bird", "Mobility x Bird", "Reproduction x Bird",
                              "Condition x Crustacean", "Longevity x Crustacean",
                              "Condition x Echinoderm",
                              "Condition x Fish", "Longevity x Fish", "Mobility x Fish", "Reproduction x Fish",
                              "Condition x Gastropod", "Longevity x Gastropod",
                              "Condition x Mammal",
                              "Condition x Reptile", "Mobility x Reptile",
                              "Longevity x Rotifer")) +
  theme(legend.position = "none") +
  scale_color_manual(values = colors.consistent.cattaxa) +
  scale_fill_manual(values = c(rep("black", 20)))


orchard_cattaxa


#ggsave("orchard_cattaxa.tif", plot = last_plot(), device = "tiff",
#       width = 3.25, height = 9, units = "in", dpi = 300)









#------------------------------------------ Plotting the 3 together---------------------------------------------

library(ggpubr)
library(cowplot)


x11()
figure1 <- ggarrange(ggarrange(orchard_catfull, orchard_catsex, nrow = 2, labels = c("A", "B"), common.legend = TRUE, legend = "bottom"), orchard_cattaxa, ncol = 2, labels = c("", "C"))

figure1


ggsave("Full_model_LRT3.tif", plot = last_plot(), device = "tiff",
       width = 165, height = 228, units = "mm", dpi = 300)

























#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Pub. bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Need to calculate residuals manually. Using 'predict.rma' to get predicted values produces residuals that are identical to the residuals I get with the function 'residuals'

predicted.vals <- predict.rma(final_model_full)$pred
original.vals <- All.data.ES$yi


x11()
Residuals_full <- original.vals - predicted.vals
Precision_full <- sqrt(1/final_model_full$vi)
plot(Residuals_full, Precision_full, xlab="Residuals", ylab="Precision [1/SE]")
abline(v=0,lty=3)

Residuals_test <- residuals(final_model_full)
plot(Residuals_test, Precision_full, xlab="Residuals", ylab="Precision [1/SE]")
abline(v=0,lty=3)


funnelmodel_full <- rma(yi = Residuals_full, sei = 1/Precision_full) 
summary(funnelmodel_full)   



#Egger's regression test
ERT <- regtest(funnelmodel_full, model="rma", predictor = "sei") #test for funnel plot asymmetry: z = -13.2464, p < 0.0001
ERT

x11()
funnela = ggplot() +
  geom_point(aes(x = Residuals_full, y = Precision_full), color = 'black', size = 0.5) +
  geom_vline(xintercept = mean(Residuals_full), linetype = 2) +
  xlim(-25, 15) + ylim(0, 8) +
  ylab(expression(paste("Precision (1/SE)"))) + xlab("") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 12), legend.title = element_text(size = 12))

funnela

ggsave("funnel_A.tif", plot = last_plot(), device = "tiff",
       width = 3, height = 3, units = "in", dpi = 300)





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



















#####################################################################################################################################################
############################################################ Within categories of responses #########################################################
#####################################################################################################################################################










####################################################### Condition-only ####################################################################





#////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#                                           Likelihood Ratio Tests                                              ####
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////////////////////////////



#Step 1 - Univariate analyses with categorical moderator variables - comparing model with one fixed effect to a null model with no fixed effects


category_LRT_cond.mr <- rma.mv(yi ~ Sub.cat, cov.cc.cond,
                               random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "ML", data = Cond.data.ES)
null_LRT_cond.mr <- rma.mv(yi ~ 1, cov.cc.cond,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "ML", data = Cond.data.ES)
cat.v.null_cond.mr <- anova(category_LRT_cond.mr, null_LRT_cond.mr)
cat.v.null_cond.mr

#         df       AIC       BIC      AICc     logLik     LRT   pval        QE 
# Full    10 1305.8740 1344.0451 1306.5509 -642.9370                2263.5422 
# Reduced  5 1322.1850 1341.2705 1322.3668 -656.0925 26.3110 <.0001 2607.2197
p.all.mr <- vector(mode = "numeric", length = 4)
p.all.mr[1] <- cat.v.null_cond.mr$pval



age_LRT_cond.mr <- rma.mv(yi ~ age, cov.cc.cond,
                          random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "ML", data = Cond.data.ES)
age.v.null_cond.mr <- anova(age_LRT_cond.mr, null_LRT_cond.mr)
age.v.null_cond.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full     6 1323.2660 1346.1686 1323.5213 -655.6330               2559.5026 
# Reduced  5 1322.1850 1341.2705 1322.3668 -656.0925 0.9190 0.3377 2607.2197 
p.all.mr[2] <- age.v.null_cond.mr$pval



sex_LRT_cond.mr <- rma.mv(yi ~ sex, cov.cc.cond,
                          random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "ML", data = Cond.data.ES)
sex.v.null_cond.mr <- anova(sex_LRT_cond.mr, null_LRT_cond.mr)
sex.v.null_cond.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full     7 1324.8228 1351.5425 1325.1642 -655.4114               2482.2313 
# Reduced  5 1322.1850 1341.2705 1322.3668 -656.0925 1.3622 0.5061 2607.2197
p.all.mr[3] <- sex.v.null_cond.mr$pval



taxa_LRT_cond.mr <- rma.mv(yi ~ Sp.group, cov.cc.cond,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "ML", data = Cond.data.ES)
taxa.v.null_cond.mr <- anova(taxa_LRT_cond.mr, null_LRT_cond.mr)
taxa.v.null_cond.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full    13 1331.0566 1380.6791 1332.1870 -652.5283               2540.1821 
# Reduced  5 1322.1850 1341.2705 1322.3668 -656.0925 7.1283 0.5229 2607.2197  
p.all.mr[4] <- taxa.v.null_cond.mr$pval


#Adjusting p-value based on Holm's correction

p.new.mr <- p.adjust(p.all.mr, method = "holm", n = length(p.all.mr))
p.new.mr
# 0.0003106349 1.0000000000 1.0000000000 1.0000000000

#Sub.cat is significant, but none of the other categorical moderators










#####  Step 2 - Examining potential for interactions among categorical responses by comparing a model with an interaction with the univariate model with only category that was significant in step 1


##---- Checking for when multiple responses are correlated at 0.5
cat_x_age_LRT_cond.mr <- rma.mv(yi ~ Sub.cat + Sub.cat:age, cov.cc.cond,
                                random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "ML", data = Cond.data.ES)
catage.v.cat_cond.mr <- anova(cat_x_age_LRT_cond.mr, category_LRT_cond.mr)
catage.v.cat_cond.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full    14 1308.6200 1362.0596 1309.9285 -640.3100               2228.8299 
# Reduced 10 1305.8740 1344.0451 1306.5509 -642.9370 5.2539 0.2622 2263.5422 
p.all.mr <- vector(mode = "numeric", length = 3)
p.all.mr[1] <- catage.v.cat_cond.mr$pval


cat_x_sex_LRT_cond.mr <- rma.mv(yi ~ Sub.cat + Sub.cat:sex, cov.cc.cond,
                                random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "ML", data = Cond.data.ES)
catsex.v.cat_cond.mr <- anova(cat_x_sex_LRT_cond.mr, category_LRT_cond.mr)
catsex.v.cat_cond.mr
#         df       AIC       BIC      AICc    logLik     LRT   pval        QE
# Full    18 1309.0930 1377.8010 1311.2507 -636.5465                2010.6891 
# Reduced 10 1305.8740 1344.0451 1306.5509 -642.9370 12.7810 0.1196 2263.5422 
p.all.mr[2] <- catsex.v.cat_cond.mr$pval


cat_x_taxa_LRT_cond.mr <- rma.mv(yi ~ Sub.cat + Sub.cat:Sp.group, cov.cc.cond,
                                 random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "ML", data = Cond.data.ES)
cattaxa.v.cat_cond.mr <- anova(cat_x_taxa_LRT_cond.mr, category_LRT_cond.mr)
cattaxa.v.cat_cond.mr
#         df       AIC       BIC      AICc    logLik     LRT   pval        QE
# Full    29 1315.8375 1426.5338 1321.5238 -628.9188                2229.0337 
# Reduced 10 1305.8740 1344.0451 1306.5509 -642.9370 28.0364 0.0827 2263.5422 
p.all.mr[3] <- cattaxa.v.cat_cond.mr$pval



p.all.mr
# 0.26222050 0.11961129 0.08273001


# No interactions improve the model.



































#******************************************************************************************************************************************************
#                                       Quantifying Heterogeneity                                        ####
#******************************************************************************************************************************************************
# Null model

hetero_randonly_cond <- rma.mv(yi ~ 1, cov.cc.cond, 
                      random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, data = Cond.data.ES, method = "REML")
summary(hetero_randonly_cond)


sum(hetero_randonly_cond$sigma2)/(sum(hetero_randonly_cond$sigma2)+(sum(1/hg.cond$vi)*(hetero_randonly_cond$k-1)/(sum(1/hg.cond$vi)^2-sum((1/hg.cond$vi)^2))))*100 # total heterogeneity 97.55%
s2t.randonly <- sum(hetero_randonly_cond$sigma2) + sum(1/hetero_randonly_cond$vi) * (hetero_randonly_cond$k-1) / (sum(1/hetero_randonly_cond$vi)^2 - sum((1/hetero_randonly_cond$vi)^2))
hetero_randonly_cond$sigma2[1]/s2t.randonly*100 # I^2 between study (id) = 2.0%
hetero_randonly_cond$sigma2[2]/s2t.randonly*100 # I^2 phylogeny = 75.1%
hetero_randonly_cond$sigma2[3]/s2t.randonly*100 # I^2 species no phylo < 1.0%
hetero_randonly_cond$sigma2[4]/s2t.randonly*100 # I^2 residual = 20.4%



# Compare to model with important moderators from LRTs
hetero_cond <- rma.mv(yi ~ Sub.cat, cov.cc.cond, 
                      random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, data = Cond.data.ES, method = "REML")
summary(hetero_cond)


sum(hetero_cond$sigma2)/(sum(hetero_cond$sigma2)+(sum(1/hg.cond$vi)*(hetero_cond$k-1)/(sum(1/hg.cond$vi)^2-sum((1/hg.cond$vi)^2))))*100 # total heterogeneity 97.24%
s2t <- sum(hetero_cond$sigma2) + sum(1/hetero_cond$vi) * (hetero_cond$k-1) / (sum(1/hetero_cond$vi)^2 - sum((1/hetero_cond$vi)^2))
hetero_cond$sigma2[1]/s2t*100 # I^2 between study (id) = 3.0%
hetero_cond$sigma2[2]/s2t*100 # I^2 phylogeny = 73.0%
hetero_cond$sigma2[3]/s2t*100 # I^2 species no phylo < 1.0%
hetero_cond$sigma2[4]/s2t*100 # I^2 residual = 21.3%




#Differences in heterogeneity between model runs

#Total
TotalH.diff.cond = (sum(hetero_randonly_cond$sigma2)/(sum(hetero_randonly_cond$sigma2)+(sum(1/hg.cond$vi)*(hetero_randonly_cond$k-1)/(sum(1/hg.cond$vi)^2-sum((1/hg.cond$vi)^2))))*100) - (sum(hetero_cond$sigma2)/(sum(hetero_cond$sigma2)+(sum(1/hg.cond$vi)*(hetero_cond$k-1)/(sum(1/hg.cond$vi)^2-sum((1/hg.cond$vi)^2))))*100) # 0.31
StudyH.diff.cond = (hetero_randonly_cond$sigma2[1]/s2t.randonly*100) - (hetero_cond$sigma2[1]/s2t*100) # -0.94
PhyloH.diff.cond = (hetero_randonly_cond$sigma2[2]/s2t.randonly*100) - (hetero_cond$sigma2[2]/s2t*100) # 2.1
SppH.diff.cond = (hetero_randonly_cond$sigma2[3]/s2t.randonly*100) - (hetero_cond$sigma2[3]/s2t*100)   # < 1.0
ESH.diff.cond = (hetero_randonly_cond$sigma2[4]/s2t.randonly*100) - (hetero_cond$sigma2[4]/s2t*100)    # -0.87







# *******************************************************************************************************************************

















#------------------------------------------------------------------------- Final model ------------------------------------------------------------------------------

final_model_cond <- rma.mv(yi ~ Sub.cat, cov.cc.cond,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "REML", data = Cond.data.ES)
summary(final_model_cond)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------







#Random effects####
Rand.Study.int.1 <- ranef(final_model_cond)$Study[1]$intrcpt
Rand.Study.se.1 <- ranef(final_model_cond)$Study[2]$se
Rand.Study.min <- Rand.Study.int.1 - Rand.Study.se.1
Rand.Study.max <- Rand.Study.int.1 + Rand.Study.se.1
Rand.Study.lab <- rownames(ranef(final_model_cond)$Study[1])
nn <- Cond.data.ES %>% group_by(Study) %>% tally
faceting.Study <- c(rep(1, 27), rep(2, 26))
Rand.Study.1 <- data.frame(Rand.Study.lab, Rand.Study.int.1, Rand.Study.min, Rand.Study.max, nn, faceting.Study)
names(Rand.Study.1)[names(Rand.Study.1) == 'n'] <- "sample.size"
As <- filter(Rand.Study.1, faceting.Study == '1')
Zs <- filter(Rand.Study.1, faceting.Study == '2')



p_random = function(dat, x, y, ymin, ymax, nsample) {
  
  x11()
  
  # Convert x-variable to a factor
  dat[,x] = as.factor(dat[, x])
  
  # Plot points
  p = ggplot(dat, aes_string(x, y)) +
    geom_point(color = 'black') + 
    geom_hline(yintercept = 0) +
    geom_pointrange(aes_string(ymax = ymax, ymin = ymin), dat) +
    coord_flip() +
    ylab("Estimated Random Intercept ± SE") + xlab("") +
    theme_classic() +
    theme(axis.line.y = element_line(color = "white"), axis.ticks = element_blank(), axis.title.x = element_text(size = 32), axis.text = element_text(size = 28), legend.title = element_text(size = 32))
  
  
  # # Summarise data to get counts by x-variable
  # nn = dat2 %>% group_by(nvar) %>% tally
  
  
  # Add counts as text labels
  p + geom_text(data = dat, aes_string(label = paste0("n = ", nsample)),
                y = min(dat[, ymin]) + 0.15 * min(dat[, ymin]), colour = "grey20", size = 7) +
    scale_y_continuous(limits = c(min(dat[, ymin]) + 0.15 * min(dat[, ymin]), max(dat[, ymax])))
  
}


p_random(As, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")
p_random(Zs, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")







Rand.Species.int.1 <- ranef(final_model_cond)$Species[1]$intrcpt
Rand.Species.se.1 <- ranef(final_model_cond)$Species[2]$se
Rand.Species.min <- Rand.Species.int.1 - Rand.Species.se.1
Rand.Species.max <- Rand.Species.int.1 + Rand.Species.se.1
Rand.Species.lab <- rownames(ranef(final_model_cond)$Species[1])
nn <- Cond.data.ES %>% group_by(Species) %>% tally
faceting.Species <- c(rep(1, 25), rep(2, 24))
Rand.Species.1 <- data.frame(Rand.Species.lab, Rand.Species.int.1, Rand.Species.min, Rand.Species.max, nn, faceting.Species)
names(Rand.Species.1)[names(Rand.Species.1) == 'n'] <- "sample.size"
Aspp <- filter(Rand.Species.1, faceting.Species == '1')
Zspp <- filter(Rand.Species.1, faceting.Species == '2')

p_random(Aspp, "Rand.Species.lab", "Rand.Species.int.1", "Rand.Species.min", "Rand.Species.max", "sample.size")
p_random(Zspp, "Rand.Species.lab", "Rand.Species.int.1", "Rand.Species.min", "Rand.Species.max", "sample.size")






















################################################### Model diagnostics #########################################################
#fitted vs. residual values
plot(fitted(final_model_cond), rstandard(final_model_cond)$z, pch=19, xlab = "Fitted values", ylab = "Standardized residuals")
abline(h = 0, lty = 2)
#Some potential concerning narrowing with higher fitted values
qqnorm(rstandard(final_model_cond)$z)
qqline(rstandard(final_model_cond)$z)
hist(rstandard(final_model_cond)$z, xlab = "Standardized residuals", ylab = "Frequency", main = NULL)


profile.rma.mv(final_model_cond, sigma2 = 1)
abline(h = logLik(final_model_cond) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_cond, sigma2 = 2)
abline(h = logLik(final_model_cond) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_cond, sigma2 = 3)
abline(h = logLik(final_model_cond) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_cond, sigma2 = 4)
abline(h = logLik(final_model_cond) - qchisq(.95, df=1)/2, lty="dashed")

################################################################################################################################











#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Plotting figures ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#----------------------------- Plotting fixed effects -----------------------------

Subcategory_effects_cond <- rma.mv(yi ~ Sub.cat - 1, cov.cc.cond, 
                                random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.cond, ES.ID = mult.resp.cond), Rscale = 0, method = "REML", data = Cond.data.ES)
summary(Subcategory_effects_cond)



table_test <- mod_results(Subcategory_effects_cond, mod = "Sub.cat")
print(table_test)


colors.consistent.cond <- c("salmon", "orange", "seagreen3", "steelblue1", "plum3", "grey")


x11()
orchard_plot(Subcategory_effects_cond, mod = "Sub-Category", xlab = "Standardized mean\ndifference", alpha = 0.5, 
             transfm = "none", angle = 0, cb = FALSE) +
  #xlim(-10, 7) +
  coord_cartesian(xlim = c(-10, 7)) +
  scale_y_discrete(labels = c("Body Composition", "Development Rate", "Feeding Rate", "Immunity", "Size & Growth", "Stress Tolerance")) +
  theme(legend.position = "top", legend.justification = c("left", "top"), legend.title = element_text(size = 8), legend.text = element_text(size = 6)) +
  scale_color_manual(values = colors.consistent.cond) +
  scale_fill_manual(values = c(rep("black", 6)))






ggsave("orchard_cond3.tif", plot = last_plot(), device = "tiff",
       width = 110, height = 110, units = "mm", dpi = 300)












































#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pub. bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test1_model_cond <- rma.mv(yi ~ Sub.cat, cov.cc.cond,
                           method = "REML", data = Cond.data.ES)
test2_model_cond <- rma.mv(yi ~ Sub.cat, cov.cc.cond,
                          random = c(list(~ 1 | Study)), method = "REML", data = Cond.data.ES)

residuals(final_model_cond)
test1.resid <- residuals(test1_model_cond)
test2.resid <- residuals(test2_model_cond)


x11()
Residuals_cond <- residuals(final_model_cond)
Precision_cond <- sqrt(1/final_model_cond$vi)
plot(Residuals_cond, Precision_cond, xlab="Residuals", ylab="Precision [1/SE]")
abline(v=0,lty=3)


funnelmodel_cond <- rma(yi = Residuals_cond, sei = 1/Precision_cond) 
summary(funnelmodel_cond)  


#Egger's regression test
regtest(funnelmodel_cond, model="rma") #test for funnel plot asymmetry: z = -13.5582, p < 0.0001



funnelb = ggplot() +
  geom_point(aes(x = Residuals_cond, y = Precision_cond), color = 'black', size = 0.5) +
  geom_vline(xintercept = mean(Residuals_cond), linetype = 2) +
  xlim(-25, 15) + ylim(0, 8) +
  ylab(expression(paste(""))) + xlab("") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 12), legend.title = element_text(size = 12))

ggsave("funnel_B_update.tif", plot = last_plot(), device = "tiff",
       width = 3, height = 3, units = "in", dpi = 300)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















####################################################### Reproduction-only #######################################################





#////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#                                           Likelihood Ratio Tests                                              ####
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////////////////////////////



#Step 1 - Univariate analyses with categorical moderator variables - comparing model with one fixed effect to a null model with no fixed effects

cat_LRT_repro.mr <- rma.mv(yi ~ Sub.cat, cov.cc.rep,
                               random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, method = "ML", data = Repro.data.ES)
null_LRT_repro.mr <- rma.mv(yi ~ 1, cov.cc.rep,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, method = "ML", data = Repro.data.ES)
cat.v.null_repro.mr <- anova(cat_LRT_repro.mr, null_LRT_repro.mr)
cat.v.null_repro.mr

#         df       AIC       BIC      AICc     logLik     LRT   pval        QE 
# Full     7 500.9562 524.7434 501.4820 -243.4781                594.4027 
# Reduced  5 508.1242 525.1151 508.4033 -249.0621 11.1680 0.0038 628.2771 
p.all.mr <- vector(mode = "numeric", length = 4)
p.all.mr[1] <- cat.v.null_repro.mr$pval



age_LRT_repro.mr <- rma.mv(yi ~ age, cov.cc.rep,
                          random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, method = "ML", data = Repro.data.ES)
age.v.null_repro.mr <- anova(age_LRT_repro.mr, null_LRT_repro.mr)
age.v.null_repro.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full     6 509.5629 529.9519 509.9554 -248.7815               627.0822 
# Reduced  5 508.1242 525.1151 508.4033 -249.0621 0.5613 0.4537 628.2771
p.all.mr[2] <- age.v.null_repro.mr$pval



sex_LRT_repro.mr <- rma.mv(yi ~ sex, cov.cc.rep,
                          random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, method = "ML", data = Repro.data.ES)
sex.v.null_repro.mr <- anova(sex_LRT_repro.mr, null_LRT_repro.mr)
sex.v.null_repro.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full     7 507.3210 531.1081 507.8468 -246.6605               617.9582 
# Reduced  5 508.1242 525.1151 508.4033 -249.0621 4.8033 0.0906 628.2771
p.all.mr[3] <- sex.v.null_repro.mr$pval



taxa_LRT_repro.mr <- rma.mv(yi ~ Sp.group, cov.cc.rep,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, method = "ML", data = Repro.data.ES)
taxa.v.null_repro.mr <- anova(taxa_LRT_repro.mr, null_LRT_repro.mr)
taxa.v.null_repro.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full     8 511.9351 539.1204 512.6143 -247.9675               624.5940 
# Reduced  5 508.1242 525.1151 508.4033 -249.0621 2.1892 0.5341 628.2771
p.all.mr[4] <- taxa.v.null_repro.mr$pval


#Adjusting p-value based on Holm's correction

p.new.mr <- p.adjust(p.all.mr, method = "holm", n = length(p.all.mr))
p.new.mr
# 0.0150298 0.9074441 0.2717118 0.9074441

# Sub-category is significant








#####  Step 2 - Examining potential for interactions among categorical responses by comparing a model with an interaction with the univariate model with only category that was significant in step 1

cat_x_age_LRT_repro.mr <- rma.mv(yi ~ Sub.cat + Sub.cat:age, cov.cc.rep,
                                random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, method = "ML", data = Repro.data.ES)
catage.v.cat_repro.mr <- anova(cat_x_age_LRT_repro.mr, cat_LRT_repro.mr)
catage.v.cat_repro.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full    10 503.9015 537.8831 504.9491 -241.9508               591.2984 
# Reduced  7 500.9562 524.7434 501.4820 -243.4781 3.0547 0.3833 594.4027  
p.all.mr <- vector(mode = "numeric", length = 3)
p.all.mr[1] <- catage.v.cat_repro.mr$pval


cat_x_sex_LRT_repro.mr <- rma.mv(yi ~ Sub.cat + Sub.cat:sex, cov.cc.rep,
                                random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, method = "ML", data = Repro.data.ES)
catsex.v.cat_repro.mr <- anova(cat_x_sex_LRT_repro.mr, cat_LRT_repro.mr)
catsex.v.cat_repro.mr
#         df       AIC       BIC      AICc    logLik     LRT   pval        QE
# Full    12 504.1005 544.8785 505.6005 -240.0503               585.9062 
# Reduced  7 500.9562 524.7434 501.4820 -243.4781 6.8557 0.2316 594.4027  
p.all.mr[2] <- catsex.v.cat_repro.mr$pval


cat_x_taxa_LRT_repro.mr <- rma.mv(yi ~ Sub.cat + Sub.cat:Sp.group, cov.cc.rep,
                                 random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, method = "ML", data = Repro.data.ES)
cattaxa.v.cat_repro.mr <- anova(cat_x_taxa_LRT_repro.mr, cat_LRT_repro.mr)
cattaxa.v.cat_repro.mr
#        df       AIC       BIC      AICc    logLik     LRT   pval        QE
# Full    14 504.0366 551.6109 506.0754 -238.0183                577.6315 
# Reduced  7 500.9562 524.7434 501.4820 -243.4781 10.9196 0.1422 594.4027   
p.all.mr[3] <- cattaxa.v.cat_repro.mr$pval



p.all.mr
# 0.3832677 0.2315945 0.1421643


# No interactions improve the model...


































#******************************************************************************************************************************************************
#                                       Quantifying Heterogeneity                                        ####
#******************************************************************************************************************************************************
# Null model

hetero_randonly_repro <- rma.mv(yi ~ 1, cov.cc.rep, 
                                random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, data = Repro.data.ES, method = "REML")
summary(hetero_randonly_repro)


sum(hetero_randonly_repro$sigma2)/(sum(hetero_randonly_repro$sigma2)+(sum(1/hg.repro$vi)*(hetero_randonly_repro$k-1)/(sum(1/hg.repro$vi)^2-sum((1/hg.repro$vi)^2))))*100 # total heterogeneity 73.18%
s2t.randonly <- sum(hetero_randonly_repro$sigma2) + sum(1/hetero_randonly_repro$vi) * (hetero_randonly_repro$k-1) / (sum(1/hetero_randonly_repro$vi)^2 - sum((1/hetero_randonly_repro$vi)^2))
hetero_randonly_repro$sigma2[1]/s2t.randonly*100 # I^2 between study (id) = 41.9%
hetero_randonly_repro$sigma2[2]/s2t.randonly*100 # I^2 phylogeny = 3.9%
hetero_randonly_repro$sigma2[3]/s2t.randonly*100 # I^2 species no phylo < 1.0%
hetero_randonly_repro$sigma2[4]/s2t.randonly*100 # I^2 residual = 27.4%



# Compare to model with important moderators from LRTs
hetero_repro <- rma.mv(yi ~ Sub.cat, cov.cc.rep, 
                       random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, data = Repro.data.ES, method = "REML")
summary(hetero_repro)


sum(hetero_repro$sigma2)/(sum(hetero_repro$sigma2)+(sum(1/hg.repro$vi)*(hetero_repro$k-1)/(sum(1/hg.repro$vi)^2-sum((1/hg.repro$vi)^2))))*100 # total heterogeneity 71.78%
s2t <- sum(hetero_repro$sigma2) + sum(1/hetero_repro$vi) * (hetero_repro$k-1) / (sum(1/hetero_repro$vi)^2 - sum((1/hetero_repro$vi)^2))
hetero_repro$sigma2[1]/s2t*100 # I^2 between study (id) = 44.1%
hetero_repro$sigma2[2]/s2t*100 # I^2 phylogeny < 1.0%
hetero_repro$sigma2[3]/s2t*100 # I^2 species no phylo < 1.0%
hetero_repro$sigma2[4]/s2t*100 # I^2 residual = 27.7%




#Differences in heterogeneity between model runs

#Total
TotalH.diff.rep = (sum(hetero_randonly_repro$sigma2)/(sum(hetero_randonly_repro$sigma2)+(sum(1/hg.repro$vi)*(hetero_randonly_repro$k-1)/(sum(1/hg.repro$vi)^2-sum((1/hg.repro$vi)^2))))*100) - (sum(hetero_repro$sigma2)/(sum(hetero_repro$sigma2)+(sum(1/hg.repro$vi)*(hetero_repro$k-1)/(sum(1/hg.repro$vi)^2-sum((1/hg.repro$vi)^2))))*100) # 1.4
StudyH.diff.rep = (hetero_randonly_repro$sigma2[1]/s2t.randonly*100) - (hetero_repro$sigma2[1]/s2t*100) # -2.2
PhyloH.diff.rep = (hetero_randonly_repro$sigma2[2]/s2t.randonly*100) - (hetero_repro$sigma2[2]/s2t*100) # 3.9
SppH.diff.rep = (hetero_randonly_repro$sigma2[3]/s2t.randonly*100) - (hetero_repro$sigma2[3]/s2t*100)   # < 1.0
ESH.diff.rep = (hetero_randonly_repro$sigma2[4]/s2t.randonly*100) - (hetero_repro$sigma2[4]/s2t*100)    # -0.37



# *******************************************************************************************************************************























#------------------------------------------------------------------------- Final model ------------------------------------------------------------------------------

final_model_repro <- rma.mv(yi ~ Sub.cat, cov.cc.rep,
                 random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, data = Repro.data.ES, method = "REML")
summary(final_model_repro)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------















#Random effects####
Rand.Study.int.1 <- ranef(final_model_repro)$Study[1]$intrcpt
Rand.Study.se.1 <- ranef(final_model_repro)$Study[2]$se
Rand.Study.min <- Rand.Study.int.1 - Rand.Study.se.1
Rand.Study.max <- Rand.Study.int.1 + Rand.Study.se.1
Rand.Study.lab <- rownames(ranef(final_model_repro)$Study[1])
nn <- Repro.data.ES %>% group_by(Study) %>% tally
faceting.Study <- c(rep(1, 18), rep(2, 18))
Rand.Study.1 <- data.frame(Rand.Study.lab, Rand.Study.int.1, Rand.Study.min, Rand.Study.max, nn, faceting.Study)
names(Rand.Study.1)[names(Rand.Study.1) == 'n'] <- "sample.size"
As <- filter(Rand.Study.1, faceting.Study == '1')
Zs <- filter(Rand.Study.1, faceting.Study == '2')



p_random = function(dat, x, y, ymin, ymax, nsample) {
  
  x11()
  
  # Convert x-variable to a factor
  dat[,x] = as.factor(dat[, x])
  
  # Plot points
  p = ggplot(dat, aes_string(x, y)) +
    geom_point(color = 'black') + 
    geom_hline(yintercept = 0) +
    geom_pointrange(aes_string(ymax = ymax, ymin = ymin), dat) +
    coord_flip() +
    ylab("Estimated Random Intercept ± SE") + xlab("") +
    theme_classic() +
    theme(axis.line.y = element_line(color = "white"), axis.ticks = element_blank(), axis.title.x = element_text(size = 32), axis.text = element_text(size = 28), legend.title = element_text(size = 32))
  
  
  # # Summarise data to get counts by x-variable
  # nn = dat2 %>% group_by(nvar) %>% tally
  
  
  # Add counts as text labels
  p + geom_text(data = dat, aes_string(label = paste0("n = ", nsample)),
                y = min(dat[, ymin]) + 0.15 * min(dat[, ymin]), colour = "grey20", size = 7) +
    scale_y_continuous(limits = c(min(dat[, ymin]) + 0.15 * min(dat[, ymin]), max(dat[, ymax])))
  
}


p_random(As, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")
p_random(Zs, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")







Rand.Species.int.1 <- ranef(final_model_repro)$Species[1]$intrcpt
Rand.Species.se.1 <- ranef(final_model_repro)$Species[2]$se
Rand.Species.min <- Rand.Species.int.1 - Rand.Species.se.1
Rand.Species.max <- Rand.Species.int.1 + Rand.Species.se.1
Rand.Species.lab <- rownames(ranef(final_model_repro)$Species[1])
nn <- Repro.data.ES %>% group_by(Species) %>% tally
faceting.Species <- c(rep(1, 14), rep(2, 14))
Rand.Species.1 <- data.frame(Rand.Species.lab, Rand.Species.int.1, Rand.Species.min, Rand.Species.max, nn, faceting.Species)
names(Rand.Species.1)[names(Rand.Species.1) == 'n'] <- "sample.size"
Aspp <- filter(Rand.Species.1, faceting.Species == '1')
Zspp <- filter(Rand.Species.1, faceting.Species == '2')


p_random(Aspp, "Rand.Species.lab", "Rand.Species.int.1", "Rand.Species.min", "Rand.Species.max", "sample.size")
p_random(Zspp, "Rand.Species.lab", "Rand.Species.int.1", "Rand.Species.min", "Rand.Species.max", "sample.size")

















































 



################################################### Model diagnostics #########################################################
#fitted vs. residual values
plot(fitted(final_model_repro), rstandard(final_model_repro)$z, pch=19, xlab = "Fitted values", ylab = "Standardized residuals")
abline(h = 0, lty = 2)
qqnorm(rstandard(final_model_repro)$z)
qqline(rstandard(final_model_repro)$z)
hist(rstandard(final_model_repro)$z, xlab = "Standardized residuals", ylab = "Frequency", main = NULL)



profile.rma.mv(final_model_repro, sigma2 = 1)
abline(h = logLik(final_model_repro) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_repro, sigma2 = 2)
abline(h = logLik(final_model_repro) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_repro, sigma2 = 3)
abline(h = logLik(final_model_repro) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_repro, sigma2 = 4)
abline(h = logLik(final_model_repro) - qchisq(.95, df=1)/2, lty="dashed")


################################################################################################################################










#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Plotting figures ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




#----------------------------- Plotting fixed effects of categorical moderator variables -----------------------------

Subcategory_effects_repro <- rma.mv(yi ~ Sub.cat - 1, cov.cc.rep, 
                                   random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.rep, ES.ID = mult.resp.rep), Rscale = 0, method = "REML", data = Repro.data.ES)
summary(Subcategory_effects_repro)



table_test <- mod_results(Subcategory_effects_repro, mod = "Sub.cat")
print(table_test)


colors.consistent.repro <- c("salmon", "orange", "steelblue1")


# x11()
orchard_plot(Subcategory_effects_repro, mod = "Sub-Category", xlab = "Standardized mean difference", alpha = 0.5, transfm = "none", angle = 0, cb = FALSE) +
  #xlim(-3, 1.5) +
  coord_cartesian(xlim = c(-3, 1.5)) +
  scale_y_discrete(labels = c("Egg & Offspring Quality", "Egg & Offspring Quantity", "Mating")) +
  theme(legend.position = "top", legend.justification = c("left", "top"), legend.title = element_text(size = 8), legend.text = element_text(size = 6)) +
  scale_color_manual(values = colors.consistent.repro) +
  scale_fill_manual(values = c(rep("black", 3))) +
  annotate("text", x = -1, y = 1:2, label = "*", size = 10)



ggsave("orchard_repro3.tif", plot = last_plot(), device = "tiff",
    width = 110, height = 110, units = "mm", dpi = 300)
























#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pub. bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



x11()
Residuals_repro <- residuals(final_model_repro)
Precision_repro <- sqrt(1/final_model_repro$vi)
plot(Residuals_repro, Precision_repro, xlab="Residuals", ylab="Precision [1/SE]")
abline(v=0,lty=3)


funnelmodel_repro <- rma(yi = Residuals_repro, sei = 1/Precision_repro) 
summary(funnelmodel_repro) 



#Egger's regression test
regtest(funnelmodel_repro, model="rma") #test for funnel plot asymmetry: z = -6.1622, p < .0001




funnelc = ggplot() +
  geom_point(aes(x = Residuals_repro, y = Precision_repro), color = 'black', size = 0.5) +
  geom_vline(xintercept = mean(Residuals_repro), linetype = 2) +
  xlim(-25, 15) + ylim(0, 8) +
  ylab(expression(paste("Precision (1/SE)"))) + xlab(expression(paste("Residuals"))) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 12), legend.title = element_text(size = 12))
funnelc

ggsave("funnel_C_update.tif", plot = last_plot(), device = "tiff",
       width = 3, height = 3, units = "in", dpi = 300)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
























####################################################### Insect-only ####################################################################





#////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#                                           Likelihood Ratio Tests                                              ####
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////////////////////////////



#Step 1 - Univariate analyses with categorical moderator variables - comparing model with one fixed effect to a null model with no fixed effects
cat_LRT_insect.mr <- rma.mv(yi ~ Category, cov.cc.insect,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, method = "ML", data = Insect.data.ES)
null_LRT_insect.mr <- rma.mv(yi ~ 1, cov.cc.insect,
                            random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, method = "ML", data = Insect.data.ES)
cat.v.null_insect.mr <- anova(cat_LRT_insect.mr, null_LRT_insect.mr)
cat.v.null_insect.mr

#         df       AIC       BIC      AICc     logLik     LRT   pval        QE 
# Full     8 1088.7050 1120.0561 1089.1016 -536.3525                2290.1720 
# Reduced  5 1157.4201 1177.0145 1157.5840 -573.7100 74.7151 <.0001 2644.4223 
p.all.mr <- vector(mode = "numeric", length = 3)
p.all.mr[1] <- cat.v.null_insect.mr$pval



age_LRT_insect.mr <- rma.mv(yi ~ age, cov.cc.insect,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, method = "ML", data = Insect.data.ES)
age.v.null_insect.mr <- anova(age_LRT_insect.mr, null_LRT_insect.mr)
age.v.null_insect.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full     6 1159.4081 1182.9215 1159.6383 -573.7041               2642.9144 
# Reduced  5 1157.4201 1177.0145 1157.5840 -573.7100 0.0119 0.9131 2644.4223
p.all.mr[2] <- age.v.null_insect.mr$pval



sex_LRT_insect.mr <- rma.mv(yi ~ sex, cov.cc.insect,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, method = "ML", data = Insect.data.ES)
sex.v.null_insect.mr <- anova(sex_LRT_insect.mr, null_LRT_insect.mr)
sex.v.null_insect.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full     7 1161.0453 1188.4775 1161.3530 -573.5226               2586.7191 
# Reduced  5 1157.4201 1177.0145 1157.5840 -573.7100 0.3748 0.8291 2644.4223
p.all.mr[3] <- sex.v.null_insect.mr$pval



#Adjusting p-value based on Holm's correction

p.new.mr <- p.adjust(p.all.mr, method = "holm", n = length(p.all.mr))
p.new.mr
# 1.25109e-15 1.00000e+00 1.00000e+00

#Category is significant










#####  Step 2 - Examining potential for interactions among categorical responses by comparing a model with an interaction with the univariate model with only category that was significant in step 1


cat_x_age_LRT_insect.mr <- rma.mv(yi ~ Category + Category:age, cov.cc.insect,
                                 random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, method = "ML", data = Insect.data.ES)
catage.v.cat_insect.mr <- anova(cat_x_age_LRT_insect.mr, cat_LRT_insect.mr)
catage.v.cat_insect.mr
#         df       AIC       BIC      AICc    logLik    LRT   pval        QE
# Full    11 1087.4644 1130.5723 1088.1978 -532.7322               2270.5003 
# Reduced  8 1088.7050 1120.0561 1089.1016 -536.3525 7.2405 0.0646 2290.1720 
p.all.mr <- vector(mode = "numeric", length = 2)
p.all.mr[1] <- catage.v.cat_insect.mr$pval


cat_x_sex_LRT_insect.mr <- rma.mv(yi ~ Category + Category:sex, cov.cc.insect,
                                 random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, method = "ML", data = Insect.data.ES)
catsex.v.cat_insect.mr <- anova(cat_x_sex_LRT_insect.mr, cat_LRT_insect.mr)
catsex.v.cat_insect.mr
#         df       AIC       BIC      AICc    logLik     LRT   pval        QE
# Full    14 1090.0601 1144.9246 1091.2365 -531.0300                2089.2803 
# Reduced  8 1088.7050 1120.0561 1089.1016 -536.3525 10.6449 0.1000 2290.1720 
p.all.mr[2] <- catsex.v.cat_insect.mr$pval




p.all.mr
# 0.06461462 0.09999146


# No interactions improve the model...



#Potential final model?
Cat_LRT_final <- rma.mv(yi ~ Category, cov.cc.insect,
                           random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, method = "REML", data = Insect.data.ES)
summary(Cat_LRT_final)





























#******************************************************************************************************************************************************
#                                       Quantifying Heterogeneity                                        ####
#******************************************************************************************************************************************************
# Null model

hetero_randonly_insect <- rma.mv(yi ~ 1, cov.cc.insect, 
                               random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, data = Insect.data.ES, method = "REML")
summary(hetero_randonly_insect)


sum(hetero_randonly_insect$sigma2)/(sum(hetero_randonly_insect$sigma2)+(sum(1/hg.insect$vi)*(hetero_randonly_insect$k-1)/(sum(1/hg.insect$vi)^2-sum((1/hg.insect$vi)^2))))*100 # total heterogeneity 90.64%
s2t.randonly <- sum(hetero_randonly_insect$sigma2) + sum(1/hetero_randonly_insect$vi) * (hetero_randonly_insect$k-1) / (sum(1/hetero_randonly_insect$vi)^2 - sum((1/hetero_randonly_insect$vi)^2))
hetero_randonly_insect$sigma2[1]/s2t.randonly*100 # I^2 between study (id) = 12.6%
hetero_randonly_insect$sigma2[2]/s2t.randonly*100 # I^2 phylogeny < 1.0%%
hetero_randonly_insect$sigma2[3]/s2t.randonly*100 # I^2 species no phylo < 1.0%
hetero_randonly_insect$sigma2[4]/s2t.randonly*100 # I^2 residual = 78.0%



# Compare to model with important moderators from LRTs
hetero_insect <- rma.mv(yi ~ Category, cov.cc.insect, 
                      random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, data = Insect.data.ES, method = "REML")
summary(hetero_insect)


sum(hetero_insect$sigma2)/(sum(hetero_insect$sigma2)+(sum(1/hg.insect$vi)*(hetero_insect$k-1)/(sum(1/hg.insect$vi)^2-sum((1/hg.insect$vi)^2))))*100 # total heterogeneity 92.5%
s2t <- sum(hetero_insect$sigma2) + sum(1/hetero_insect$vi) * (hetero_insect$k-1) / (sum(1/hetero_insect$vi)^2 - sum((1/hetero_insect$vi)^2))
hetero_insect$sigma2[1]/s2t*100 # I^2 between study (id) = 29.7%
hetero_insect$sigma2[2]/s2t*100 # I^2 phylogeny = 21.4%
hetero_insect$sigma2[3]/s2t*100 # I^2 species no phylo < 1.0%
hetero_insect$sigma2[4]/s2t*100 # I^2 residual = 41.4%




#Differences in heterogeneity between model runs

#Total
TotalH.diff.insect = (sum(hetero_randonly_insect$sigma2)/(sum(hetero_randonly_insect$sigma2)+(sum(1/hg.insect$vi)*(hetero_randonly_insect$k-1)/(sum(1/hg.insect$vi)^2-sum((1/hg.insect$vi)^2))))*100) - (sum(hetero_insect$sigma2)/(sum(hetero_insect$sigma2)+(sum(1/hg.insect$vi)*(hetero_insect$k-1)/(sum(1/hg.insect$vi)^2-sum((1/hg.insect$vi)^2))))*100) # -1.88
StudyH.diff.insect = (hetero_randonly_insect$sigma2[1]/s2t.randonly*100) - (hetero_insect$sigma2[1]/s2t*100) # -17.1
PhyloH.diff.insect = (hetero_randonly_insect$sigma2[2]/s2t.randonly*100) - (hetero_insect$sigma2[2]/s2t*100) # -21.4
SppH.diff.insect = (hetero_randonly_insect$sigma2[3]/s2t.randonly*100) - (hetero_insect$sigma2[3]/s2t*100)   # < 1.0
ESH.diff.insect = (hetero_randonly_insect$sigma2[4]/s2t.randonly*100) - (hetero_insect$sigma2[4]/s2t*100)    # 36.6







# *******************************************************************************************************************************

































#------------------------------------------------------------------------- Final model ------------------------------------------------------------------------------

final_model_insect <- rma.mv(yi ~ Category, cov.cc.insect,
                             random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, method = "REML", data = Insect.data.ES)
summary(final_model_insect)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------









Rand.Study.int.1 <- ranef(final_model_insect)$Study[1]$intrcpt
Rand.Study.se.1 <- ranef(final_model_insect)$Study[2]$se
Rand.Study.min <- Rand.Study.int.1 - Rand.Study.se.1
Rand.Study.max <- Rand.Study.int.1 + Rand.Study.se.1
Rand.Study.lab <- rownames(ranef(final_model_insect)$Study[1])
nn <- Insect.data.ES %>% group_by(Study) %>% tally
faceting.Study <- c(rep(1, 20), rep(2, 19))
Rand.Study.1 <- data.frame(Rand.Study.lab, Rand.Study.int.1, Rand.Study.min, Rand.Study.max, nn, faceting.Study)
names(Rand.Study.1)[names(Rand.Study.1) == 'n'] <- "sample.size"
As <- filter(Rand.Study.1, faceting.Study == '1')
Zs <- filter(Rand.Study.1, faceting.Study == '2')



p_random = function(dat, x, y, ymin, ymax, nsample) {
  
  x11()
  
  # Convert x-variable to a factor
  dat[,x] = as.factor(dat[, x])
  
  # Plot points
  p = ggplot(dat, aes_string(x, y)) +
    geom_point(color = 'black') + 
    geom_hline(yintercept = 0) +
    geom_pointrange(aes_string(ymax = ymax, ymin = ymin), dat) +
    coord_flip() +
    ylab("Estimated Random Intercept ± SE") + xlab("") +
    theme_classic() +
    theme(axis.line.y = element_line(color = "white"), axis.ticks = element_blank(), axis.title.x = element_text(size = 32), axis.text = element_text(size = 28), legend.title = element_text(size = 32))
  
  
  # # Summarise data to get counts by x-variable
  # nn = dat2 %>% group_by(nvar) %>% tally
  
  
  # Add counts as text labels
  p + geom_text(data = dat, aes_string(label = paste0("n = ", nsample)),
                y = min(dat[, ymin]) + 0.15 * min(dat[, ymin]), colour = "grey20", size = 7) +
    scale_y_continuous(limits = c(min(dat[, ymin]) + 0.15 * min(dat[, ymin]), max(dat[, ymax])))
  
}


p_random(As, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")
p_random(Zs, "Rand.Study.lab", "Rand.Study.int.1", "Rand.Study.min", "Rand.Study.max", "sample.size")







Rand.Species.int.1 <- ranef(final_model_insect)$Species[1]$intrcpt
Rand.Species.se.1 <- ranef(final_model_insect)$Species[2]$se
Rand.Species.min <- Rand.Species.int.1 - Rand.Species.se.1
Rand.Species.max <- Rand.Species.int.1 + Rand.Species.se.1
Rand.Species.lab <- rownames(ranef(final_model_insect)$Species[1])
nn <- Insect.data.ES %>% group_by(Species) %>% tally
faceting.Species <- c(rep(1, 16), rep(2, 16))
Rand.Species.1 <- data.frame(Rand.Species.lab, Rand.Species.int.1, Rand.Species.min, Rand.Species.max, nn, faceting.Species)
names(Rand.Species.1)[names(Rand.Species.1) == 'n'] <- "sample.size"
Aspp <- filter(Rand.Species.1, faceting.Species == '1')
Zspp <- filter(Rand.Species.1, faceting.Species == '2')


p_random(Aspp, "Rand.Species.lab", "Rand.Species.int.1", "Rand.Species.min", "Rand.Species.max", "sample.size")
p_random(Zspp, "Rand.Species.lab", "Rand.Species.int.1", "Rand.Species.min", "Rand.Species.max", "sample.size")










################################################### Model diagnostics #########################################################
#fitted vs. residual values
plot(fitted(final_model_insect), rstandard(final_model_insect)$z, pch=19, xlab = "Fitted values", ylab = "Standardized residuals")
abline(h = 0, lty = 2)
#Concerning funnel shape with increasing fitted values
qqnorm(rstandard(final_model_insect)$z)
qqline(rstandard(final_model_insect)$z)
hist(rstandard(final_model_insect)$z, xlab = "Standardized residuals", ylab = "Frequency", main = NULL)


profile.rma.mv(final_model_insect, sigma2 = 1)
abline(h = logLik(final_model_insect) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_insect, sigma2 = 2)
abline(h = logLik(final_model_insect) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_insect, sigma2 = 3)
abline(h = logLik(final_model_insect) - qchisq(.95, df=1)/2, lty="dashed")

profile.rma.mv(final_model_insect, sigma2 = 4)
abline(h = logLik(final_model_insect) - qchisq(.95, df=1)/2, lty="dashed")


################################################################################################################################













#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Plotting figures ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^





#----------------------------- Plotting fixed effects -----------------------------

Cat_effects_insect <- rma.mv(yi ~ Category - 1, cov.cc.insect, 
                                   random = c(list(~ 1 | Study, ~ 1 | Species, ~ 1 | Species.no.phylo, ~ 1 | ES.ID)), R = list(Species = corMat.insect, ES.ID = mult.resp.insect), Rscale = 0, method = "REML", data = Insect.data.ES)
summary(Cat_effects_insect)



table_test <- mod_results(Cat_effects_insect, mod = "Category")
print(table_test)


colors.consistent.insect <- c("salmon", "seagreen3", "steelblue1", "plum3")

x11()
orchard_plot(Cat_effects_insect, mod = "Category", xlab = "Standardized mean difference", alpha = 0.5, transfm = "none", angle = 0, cb = FALSE) +
  #xlim(-7.5, 5) +
  coord_cartesian(xlim = c(-7.5, 5)) +
  theme(legend.position = "top", legend.justification = c("left", "top"), legend.title = element_text(size = 8), legend.text = element_text(size = 6)) +
  scale_color_manual(values = colors.consistent.insect) +
  scale_fill_manual(values = c(rep("black", 4))) +
  annotate("text", x = -2, y = 1.1, label = "*", size = 10) +
  annotate("text", x = -1, y = 3.1, label = "*", size = 10)





ggsave("orchard_insect2.tif", plot = last_plot(), device = "tiff",
       width = 110, height = 110, units = "mm", dpi = 300)















#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pub. bias %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x11()
Residuals_insect <- residuals(final_model_insect)
Precision_insect <- sqrt(1/final_model_insect$vi)
plot(Residuals_insect, Precision_insect, xlab="Residuals", ylab="Precision [1/SE]")
abline(v=0,lty=3)


funnelmodel_insect <- rma(yi = Residuals_insect, sei = 1/Precision_insect) 
summary(funnelmodel_insect) 


#Egger's regression test
regtest(funnelmodel_insect, model="rma") #test for funnel plot asymmetry: z = -9.0774, p < .0001



funneld = ggplot() +
  geom_point(aes(x = Residuals_insect, y = Precision_insect), color = 'black', size = 0.5) +
  geom_vline(xintercept = mean(Residuals_insect), linetype = 2) +
  xlim(-25, 15) + ylim(0, 8) +
  ylab(expression(paste(""))) + xlab(expression(paste("Residuals"))) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 12), legend.title = element_text(size = 12))
funneld

ggsave("funnel_D.tif", plot = last_plot(), device = "tiff",
       width = 3, height = 3, units = "in", dpi = 300)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















