#### Library and file loading ####
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(metafor)
library(forestplot)

source("0-Function_MA_TpH.R")

Inverts <- data.table(read_csv("file_Inverts.csv"))
Fish <- data.table(read_csv("file_Fish.csv"))


#### 1.1. Testing effect of scenario for each metric ####
## "Studies with NAs omitted from model fitting" correspond to TpH treatments 
  # for which no scenario could be assigned (because pH and T scenarios differed)
## Inverts, Metabolism
InvertMetabolism <- Inverts[Category %in% c("Routine respiration","Aerobic scope respiration")]
InvertMetabolism$Category <- c("Metabolism")
InvertMet <- MA_TpH_moderator("inverts","Metabolism",InvertMetabolism)
Abs_InvertMet <- MA_TpH_abs_moderator("inverts","Metabolism",InvertMetabolism)

## Fish, Metabolism
FishMetabolism <- Fish[Category %in% c("Routine respiration","Aerobic scope respiration")]
FishMetabolism$Category <- c("Metabolism")
FishMet <- MA_TpH_moderator("fish","Metabolism",FishMetabolism)
Abs_FishMet <- MA_TpH_abs_moderator("fish","Metabolism",FishMetabolism)

## Invertebrates, Calcification
InvertCalci <- MA_TpH_moderator("inverts","Calcification", Inverts)
Abs_InvertCalci <- MA_TpH_abs_moderator("inverts","Calcification", Inverts)

## Fish, Calcification 
# not enough data to test effect of scenarios 

## Invertebrate, Biodiversity
InvertBiodiv <- MA_TpH_moderator("inverts","Biodiversity", Inverts)
Abs_InvertBiodiv <- MA_TpH_abs_moderator("inverts","Biodiversity", Inverts)

## Invertebrate, Survival
InvertSurvi <- MA_TpH_moderator("inverts","Survival", Inverts)
Abs_InvertSurvi <- MA_TpH_abs_moderator("inverts","Survival", Inverts)

## Fish, Survival
FishSurvi <- MA_TpH_moderator("fish","Survival", Fish)
Abs_FishSurvi <- MA_TpH_abs_moderator("fish","Survival", Fish)

## Invertebrate, Reproduction 
InvertRepro <- MA_TpH_moderator("inverts","Reproduction", Inverts)
Abs_InvertRepro <- MA_TpH_abs_moderator("inverts","Reproduction", Inverts)

## Fish, Reproduction 
FishRepro <- MA_TpH_moderator("fish","Reproduction", Fish)
Abs_FishRepro <- MA_TpH_abs_moderator("fish","Reproduction", Fish)

## Invertebrate, Behaviour
InvertBehav <- MA_TpH_moderator("inverts","Behavior", Inverts)
Abs_InvertBehav <- MA_TpH_abs_moderator("inverts","Behavior", Inverts)

## Fish, Behaviour 
FishBehav  <- MA_TpH_moderator("fish","Behavior", Fish)
Abs_FishBehav  <- MA_TpH_abs_moderator("fish","Behavior", Fish)

## Invertebrate, Biomechanics
InvertBio <- MA_TpH_moderator("inverts","Biomechanics", Inverts)
Abs_InvertBio <- MA_TpH_abs_moderator("inverts","Biomechanics", Inverts)

## Invertebrate, Development
InvertDev <- MA_TpH_moderator("inverts","Development", Inverts)
Abs_InvertDev <- MA_TpH_abs_moderator("inverts","Development", Inverts)

## Fish, Development
FishDev  <- MA_TpH_moderator("fish","Development", Fish)
Abs_FishDev  <- MA_TpH_abs_moderator("fish","Development", Fish)

## Invertebrate, Growth 
InvertGrowth <- MA_TpH_moderator("inverts","Growth", Inverts)
Abs_InvertGrowth <- MA_TpH_abs_moderator("inverts","Growth", Inverts)

##Fish, Growth
FishGrowth <- MA_TpH_moderator("fish","Growth", Fish)
Abs_FishGrowth <- MA_TpH_abs_moderator("fish","Growth", Fish)

## Invertebrate, Physiology
InvertPhys <- MA_TpH_moderator("inverts","Physiology", Inverts)
Abs_InvertPhys <- MA_TpH_abs_moderator("inverts","Physiology", Inverts)

## Fish, Physiology
FishPhys <-  MA_TpH_moderator("fish","Physiology", Fish)
Abs_FishPhys <-  MA_TpH_abs_moderator("fish","Physiology", Fish)

#### 1.2. Summary table with significant effects of scenarios ####
Abs_results <- rbind(Abs_FishBehav,  Abs_FishDev, Abs_FishGrowth, Abs_FishMet,
                     Abs_FishPhys, Abs_FishRepro, Abs_FishSurvi,
                     Abs_InvertBehav,Abs_InvertBio,Abs_InvertBiodiv, Abs_InvertCalci,Abs_InvertDev, Abs_InvertGrowth, 
                     Abs_InvertMet,Abs_InvertPhys, Abs_InvertRepro, Abs_InvertSurvi)
Abs_results$QM <- formatC(as.numeric(Abs_results$QM),digits=2,format="f",flag="#")
Abs_results$QMp <-formatC(as.numeric(Abs_results$QMp),digits=2,format="f",flag="#")
Abs_results$QE <- formatC(as.numeric(Abs_results$QE),digits=0,flag="#",format="f")
Abs_results$QEp <- formatC(as.numeric(Abs_results$QEp),digits=2,flag="#",format="e")
write_csv(Abs_results,"file_ScenarioEffectAbs.csv")

Rel_results<- rbind( FishBehav,   FishDev,  FishGrowth,  FishMet,
                     FishPhys,  FishRepro,  FishSurvi,
                     InvertBehav, InvertBio, InvertBiodiv,  InvertCalci, InvertDev,  InvertGrowth, 
                     InvertMet, InvertPhys,  InvertRepro,  InvertSurvi)
Rel_results$QM <- formatC(as.numeric(Rel_results$QM),digits=2,format="f",flag="#")
Rel_results$QMp <-formatC(as.numeric(Rel_results$QMp),digits=2,format="f",flag="#")
Rel_results$QE <- formatC(as.numeric(Rel_results$QE),digits=0,flag="#",format="f")
Rel_results$QEp <- formatC(as.numeric(Rel_results$QEp),digits=2,flag="#",format="e")
write_csv(Rel_results,"file_ScenarioEffectRel.csv")