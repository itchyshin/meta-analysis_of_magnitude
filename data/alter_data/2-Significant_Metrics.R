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



#### 1. MAs for each metric & scenarios ####
  ## table containing sensitivity analysis data
sensitivity <- tibble(taxa=character(), stressor=character(), metric=character(),
                      Rosenthal=numeric(),threshold=numeric(),pval=numeric(),trimfill=numeric())
## Inverts, Metabolism
InvertMetabolism <- Inverts[Category %in% c("Routine respiration","Aerobic scope respiration")]
InvertMetabolism$Category <- c("Metabolism")
InvertMet <- MA_TpH("inverts","Metabolism",InvertMetabolism,sensitivity)
Abs_InvertMet <- MA_TpH_abs("inverts","Metabolism",InvertMetabolism)

## Fish, Metabolism
FishMetabolism <- Fish[Category %in% c("Routine respiration","Aerobic scope respiration")]
FishMetabolism$Category <- c("Metabolism")
FishMet <- MA_TpH("fish","Metabolism",FishMetabolism,sensitivity)
Abs_FishMet <- MA_TpH_abs("fish","Metabolism",FishMetabolism)

## Invertebrates, Calcification
InvertCalci <- MA_TpH("inverts","Calcification", Inverts,sensitivity)
Abs_InvertCalci <- MA_TpH_abs("inverts","Calcification", Inverts)
## Fish, Calcification 
FishCalci <- MA_TpH("fish","Calcification", Fish,sensitivity)
Abs_FishCalci <- MA_TpH_abs("fish","Calcification", Fish)

## Invertebrate, Biodiversity
InvertBiodiv <- MA_TpH("inverts","Biodiversity", Inverts,sensitivity)
Abs_InvertBiodiv <- MA_TpH_abs("inverts","Biodiversity", Inverts)

## Invertebrate, Survival
InvertSurvi <- MA_TpH("inverts","Survival", Inverts,sensitivity)
Abs_InvertSurvi <- MA_TpH_abs("inverts","Survival", Inverts)

## Fish, Survival
FishSurvi <- MA_TpH("fish","Survival", Fish,sensitivity)
Abs_FishSurvi <- MA_TpH_abs("fish","Survival", Fish)

## Invertebrate, Reproduction 
InvertRepro <- MA_TpH("inverts","Reproduction", Inverts,sensitivity)
Abs_InvertRepro <- MA_TpH_abs("inverts","Reproduction", Inverts)

## Fish, Reproduction 
FishRepro <- MA_TpH("fish","Reproduction", Fish,sensitivity)
Abs_FishRepro <- MA_TpH_abs("fish","Reproduction", Fish)

## Invertebrate, Behavior
InvertBehav <- MA_TpH("inverts","Behavior", Inverts,sensitivity)
Abs_InvertBehav <- MA_TpH_abs("inverts","Behavior", Inverts)

## Fish, Behavior 
FishBehav  <- MA_TpH("fish","Behavior", Fish,sensitivity)
Abs_FishBehav  <- MA_TpH_abs("fish","Behavior", Fish)

## Invertebrate, Biomechanics
InvertBio <- MA_TpH("inverts","Biomechanics", Inverts,sensitivity)
Abs_InvertBio <- MA_TpH_abs("inverts","Biomechanics", Inverts)

## Invertebrate, Development
InvertDev <- MA_TpH("inverts","Development", Inverts,sensitivity)
Abs_InvertDev <- MA_TpH_abs("inverts","Development", Inverts)

## Fish, Development
FishDev  <- MA_TpH("fish","Development", Fish,sensitivity)
Abs_FishDev  <- MA_TpH_abs("fish","Development", Fish)

## Invertebrate, Growth 
InvertGrowth <- MA_TpH("inverts","Growth", Inverts,sensitivity)
Abs_InvertGrowth <- MA_TpH_abs("inverts","Growth", Inverts)

##Fish, Growth
FishGrowth <- MA_TpH("fish","Growth", Fish,sensitivity)
Abs_FishGrowth <- MA_TpH_abs("fish","Growth", Fish)

## Invertebrate, Physiology
InvertPhys <- MA_TpH("inverts","Physiology", Inverts,sensitivity)
Abs_InvertPhys <- MA_TpH_abs("inverts","Physiology", Inverts)

## Fish, Physiology
FishPhys <-  MA_TpH("fish","Physiology", Fish,sensitivity)
Abs_FishPhys <-  MA_TpH_abs("fish","Physiology", Fish)

#### 2.1. Summary table of absolute meta-analyses results ####
Abs_results <- rbind(Abs_FishBehav, Abs_FishCalci, Abs_FishDev, Abs_FishGrowth, Abs_FishMet,
                     Abs_FishPhys, Abs_FishRepro, Abs_FishSurvi,
                     Abs_InvertBehav,Abs_InvertBio,Abs_InvertBiodiv, Abs_InvertCalci,Abs_InvertDev, Abs_InvertGrowth, 
                     Abs_InvertMet,Abs_InvertPhys, Abs_InvertRepro, Abs_InvertSurvi)

Abs_results$significant[Abs_results$lowerCI<=0]<- c("no")
Abs_results$significant[Abs_results$lowerCI>0]<- c("yes")
Abs_summary <- Abs_results[!Abs_results$Scenario=="NA",]
Abs_summary <- Abs_summary[!Abs_summary$N==1,] # phasing out categories for which n=1
write_csv(Abs_summary,"file_forestplot_absolute.csv")

#### 2.2. Summary table of relative meta-analyses results ####
Rel_results <- rbind(FishBehav, FishCalci, FishDev,FishGrowth, FishMet, FishPhys, FishRepro, FishSurvi,
                     InvertBehav,InvertBio, InvertBiodiv, InvertCalci,InvertDev, InvertGrowth,
                     InvertMet, InvertPhys, InvertRepro, InvertSurvi)
Rel_results$lowerCI <- as.numeric(Rel_results$lowerCI)
Rel_results$upperCI <- as.numeric(Rel_results$upperCI)
Rel_results$ei <- as.numeric(Rel_results$ei)

Rel_results$significant[(Rel_results$lowerCI*Rel_results$upperCI)>0]<- c("yes")
Rel_results$significant[(Rel_results$lowerCI*Rel_results$upperCI)<=0]<- c("no")
Rel_summary <- Rel_results[!Rel_results$Scenario=="NA",]
Rel_summary <- Rel_summary[!Rel_summary$N==1,] # phasing out categories for which n=1
write_csv(Rel_summary,"file_forestplot_relative.csv")
write_csv(sensitivity,"file_sensitivity_relative.csv")

#### 3.0. Data prep for graphs ####
## Directional changes 
Rel_summary <- read_csv("file_forestplot_relative.csv")
  # attributing a direction to meta-analysis results 
Rel_summary$direction <- c("NA")
for (k in 1:nrow(Rel_summary)){
    if (Rel_summary$significant[k]=="yes" & Rel_summary$ei[k]>0){
        Rel_summary$direction[k] <- "positive"
    }
    if (Rel_summary$significant[k]=="yes" & Rel_summary$ei[k]<0){
        Rel_summary$direction[k] <- "negative"
    }
}
  # creating summary table of nb of categories positive and negative 
signif_dir <- unique(Rel_summary[,c(1,2,4)]) # table with two lines per scenario,stressor and taxa
signif_dir <- rbind(signif_dir,signif_dir) # duplicate for positive and negative
signif_dir$direction <-c(rep("positive",nrow(unique(Rel_summary[,c(1,2,4)]))),
                         rep("negative",nrow(unique(Rel_summary[,c(1,2,4)]))))

for (k in 1:nrow(signif_dir)){
    total <- Rel_summary[Rel_summary$Stressor==as.character(signif_dir[k,"Stressor"]) & 
                             Rel_summary$Scenario==as.character(signif_dir[k,"Scenario"]) &
                             Rel_summary$taxa==as.character(signif_dir[k,"taxa"]),]
    signif_dir$total[k] <- nrow(total) # nb of metrics evaluated for that scenario x stressor
    signif_dir$significant[k] <- nrow(total[total$direction==as.character(signif_dir[k,"direction"]),]) # nb of significant metrics
}
signif_dir$StressorScenario <- paste(signif_dir$Stressor,signif_dir$Scenario)
signif_dir$total[signif_dir$direction=="negative"] <- - signif_dir$total[signif_dir$direction=="negative"]
signif_dir$significant[signif_dir$direction=="negative"] <- - signif_dir$significant[signif_dir$direction=="negative"]



signif_dir$Stressor  <- factor(signif_dir$Stressor,levels=c("T","pH","TpH"))
signif_dir$Scenario  <- factor(signif_dir$Scenario,levels=c("all","RCP6","RCP8","extreme"))
signif_dir$StressorScenario <- as.factor(signif_dir$StressorScenario)
signif_dir$StressorScenario  <- factor(signif_dir$StressorScenario,levels=c("T all","T RCP6","T RCP8", "T extreme",
                                                                            "pH all","pH RCP6", "pH RCP8", "pH extreme",
                                                                            "TpH all","TpH RCP6","TpH RCP8","TpH extreme"))
signif_dir$taxa  <- factor(signif_dir$taxa,levels=c("inverts","fish"))

## Absolute deviation
Abs_summary <- read_csv("file_forestplot_absolute.csv")
signif_dev <- unique(Abs_summary[,c(1,2,4)]) # table with one line per scenario,stressor,taxa
for (k in 1:nrow(signif_dev)){
    total <- Abs_summary[Abs_summary$Stressor==as.character(signif_dev[k,"Stressor"]) & 
                             Abs_summary$Scenario==as.character(signif_dev[k,"Scenario"]) &
                             Abs_summary$taxa==as.character(signif_dev[k,"taxa"]),]
    signif_dev$total[k] <- nrow(total) # nb of metrics evaluated for that scenario x stressor
    signif_dev$significant[k] <- nrow(total[total$significant=="yes",]) # nb of significant metrics
    signif_dev$signif_direction[k] <- sum(abs(signif_dir$significant[signif_dir$Stressor==as.character(signif_dev[k,"Stressor"]) &
                                                                         signif_dir$Scenario==as.character(signif_dev[k,"Scenario"]) &
                                                                         signif_dir$taxa==as.character(signif_dev[k,"taxa"])]))
}

signif_dev$StressorScenario <- paste(signif_dev$Stressor,signif_dev$Scenario)
signif_dev$Stressor  <- factor(signif_dev$Stressor,levels=c("T","pH","TpH"))
signif_dev$Scenario  <- factor(signif_dev$Scenario,levels=c("all","RCP6","RCP8","extreme"))
signif_dev$taxa <- factor(signif_dev$taxa, levels=c("inverts","fish"))
signif_dev$StressorScenario  <- factor(signif_dev$StressorScenario,levels=c("T all","T RCP6","T RCP8", "T extreme",
                                                                            "pH all","pH RCP6", "pH RCP8", "pH extreme",
                                                                            "TpH all","TpH RCP6","TpH RCP8","TpH extreme"))
signif_dev <- signif_dev %>% 
    arrange(Stressor,Scenario,taxa) %>% 
    add_column(ID2=rep(seq(1,nrow(signif_dev),1)))
signif_dev$IDend[signif_dev$taxa=="fish"] <- signif_dev$ID2[signif_dev$taxa=="fish"]+0.9
signif_dev$IDend[signif_dev$taxa=="inverts"] <- signif_dev$ID2[signif_dev$taxa=="inverts"]+1
signif_dev$ID2[signif_dev$taxa=="fish"] <- signif_dev$ID2[signif_dev$taxa=="fish"]-0.1

#### Barplot directional change #### 
## non stacked barplot

ggplot() +
    # nb of significant directional changes 
    geom_bar(data=signif_dir,position="dodge", stat="identity",
             mapping=aes(x=StressorScenario, y=significant, fill=taxa, col=direction), size=0)+
    scale_fill_manual(values=c("#E69F00","#999999"))+
    scale_color_manual(values=c("black","black"), guide="none")+
    scale_y_continuous(limits=c(-11,11), breaks = seq(-10,10,2),
                       labels=c("10","8","6","4","2","0","2","4","6","8","10"))+
    ylab("number of metric categories affected")+
    xlab("")+
    scale_x_discrete(labels=c(rep(c("all", "RCP6","RCP8", "extreme"),2), "all","RCP6", "extreme"))+
    theme_classic()+
    # outlines of total metrics tested
    geom_bar(data=signif_dir,position="dodge", stat="identity",
             mapping=aes(x=StressorScenario, y=total, fill=taxa, col=direction), alpha=0,
             size=0.05)+
    geom_vline(xintercept = c(4.5,8.5), color = "black", size=0.4)+
    geom_vline(xintercept = c(1.5,5.5,9.5), linetype="dashed", color = "black", size=0.2)+
    theme_classic()+
    annotate(geom="text", x=3, y=11, label="Temperature",color="black")+
    annotate(geom="text", x=7, y=11, label="pH",color="black")+
    annotate(geom="text", x=11, y=11, label="T x pH",color="black")

## stacked barplot, inverts and fish together
signif_dir$total <- abs(signif_dir$total)
signif_dir$significant <- abs(signif_dir$significant)
signif_dir$xlab <- paste(signif_dir$StressorScenario, signif_dir$taxa)
signif_dir$direction <- factor(signif_dir$direction,levels=c("positive","negative"))
signif_dir$xlab <- factor(signif_dir$xlab,levels=c("T all inverts", "T all fish", "T RCP6 inverts",
                                                   "T RCP6 fish", "T RCP8 inverts", "T RCP8 fish",
                                                   "T extreme inverts", "T extreme fish","pH all inverts",
                                                   "pH all fish", "pH RCP6 inverts", "pH RCP6 fish", "pH RCP8 inverts",
                                                   "pH RCP8 fish", "pH extreme inverts", "pH extreme fish",
                                                   "TpH all inverts", "TpH all fish", "TpH RCP6 inverts", "TpH RCP6 fish",
                                                   "TpH extreme inverts", "TpH extreme fish"))
signif_dir <- signif_dir %>% 
    arrange(Stressor,Scenario,taxa) %>% 
    add_column(ID2=rep(seq(1,nrow(signif_dir),1)))

ggplot() +
    # nb of significant directional changes 
    geom_bar(data=signif_dir,position="stack", stat="identity",
             mapping=aes(x=xlab, y=significant, fill=direction), size=0)+
    scale_fill_manual(values=c("darkseagreen1","darksalmon"))+
    scale_y_continuous(limits=c(-1,11), breaks = seq(0,10,2),
                       labels=c("0","2","4","6","8","10"))+
    scale_x_discrete(labels=rep("",22))+
    ylab("number of metric categories affected")+
    xlab("")+
    theme_classic()+
    # outlines of total metrics tested
    geom_bar(data=signif_dir,position="dodge", stat="identity",
             mapping=aes(x=xlab, y=total, col=taxa), alpha=0, size=0.3)+
    scale_color_manual(values=c("chocolate","aquamarine4"), guide="none")+
    geom_vline(xintercept = c(4.25,8.25)*2, color = "black", size=0.4)+
    geom_vline(xintercept = c(1.25,5.25,9.25)*2, linetype="dashed", color = "black", size=0.2)+
    theme_classic()+
    annotate(geom="text", x=5, y=11, label="Temperature",color="black")+
    annotate(geom="text", x=12, y=11, label="pH",color="black")+
    annotate(geom="text", x=19.5, y=11, label="T x pH",color="black")+
    annotate(geom="text", x=c(seq(0.5,21.5,2)+1),y=-1, label=c(rep(c("all", "RCP6","RCP8", "extreme"),2),
                                                               c("all","RCP6", "extreme")))

#### 4. Barplot with both relative and absolute results ####
## version 1
ggplot() +
     geom_bar(data=signif_dev,position="dodge", stat="identity",
                           mapping=aes(x=StressorScenario, y=significant, fill=taxa))+
     scale_fill_manual(values=c("#E69F00","aquamarine4"))+
     scale_y_continuous(limits=c(0,11), breaks = seq(0,12,2))+
     ylab("categories with a significant deviation")+
     xlab("")+
   scale_x_discrete(labels=c(rep(c("all", "RCP6","RCP8", "extreme"),3)))+
     geom_bar(data=signif_dev,position="dodge", stat="identity",
                           mapping=aes(x=StressorScenario, y=total, col=taxa), alpha=0, size=0.1)+
     scale_color_manual(values=c("black", "black"), labels=c("fish","inverts"))+
   geom_vline(xintercept = c(4.5,8.5), color = "black", size=0.4)+
   geom_vline(xintercept = c(1.5,5.5,9.5), linetype="dashed", color = "black", size=0.2)+
    geom_segment(data=signif_dev,mapping=aes(y = signif_direction, yend=signif_direction,
                                             x=(ID2+0.2)/2,xend=(IDend)/2), colour="black",size=0.7)+
    #geom_segment(mapping=aes(y = 8, yend=8, x=12 ,xend=13), colour="deepskyblue4",size=0.7, linetype="dotted")+
     theme_classic()+
   annotate(geom="text", x=3, y=11, label="Temperature",color="black")+
   annotate(geom="text", x=7, y=11, label="pH",color="black")+
   annotate(geom="text", x=10, y=11, label="T x pH",color="black")

## version 2 (fish and invers separate, lighter colouring for absolute deviation)
signif_dev <- signif_dev %>% 
    arrange(taxa,Stressor,Scenario)%>% 
add_column(ID3=rep(seq(1,nrow(signif_dev),1)))
signif_dev$IDend[signif_dev$taxa=="fish"] <- signif_dev$ID2[signif_dev$taxa=="fish"]+0.9
signif_dev$IDend[signif_dev$taxa=="inverts"] <- signif_dev$ID2[signif_dev$taxa=="inverts"]+1
signif_dev$ID2[signif_dev$taxa=="fish"] <- signif_dev$ID2[signif_dev$taxa=="fish"]-0.1


ggplot() +
    geom_bar(data=signif_dev,position="dodge", stat="identity",
             mapping=aes(x=ID3, y=significant, fill=taxa))+
    scale_fill_manual(values=c("#E69F00","aquamarine4"))+
    scale_y_continuous(limits=c(0,11), breaks = seq(0,12,2))+
    ylab("categories with a significant deviation")+
    xlab("")+
    scale_x_continuous(breaks=seq(1,24,1),labels=c(rep(c("all", "R6","R8.5", "ex"),6)),size=10)+
    geom_bar(data=signif_dev,position="dodge", stat="identity",
             mapping=aes(x=ID3, y=total, col=taxa), alpha=0, size=0.1)+
    theme_classic()+
    geom_bar(data=signif_dev,position="dodge", stat="identity",
             mapping=aes(x=ID3, y=signif_direction, col=taxa),
             fill="white",alpha=0.6,size=0.001)+
    scale_color_manual(values=c("black", "black"), labels=c("fish","inverts"))+
    
    geom_vline(xintercept = c(12.5), color = "black", size=0.4)+
    geom_vline(xintercept = seq(4.5,24,4), linetype="dashed", color = "black", size=0.2)+
    #geom_segment(data=signif_dev,mapping=aes(y = signif_direction, yend=signif_direction,
     #                                        x=(ID2+0.1)/2,xend=(IDend)/2), colour="black",size=0.5)+
    #geom_segment(mapping=aes(y = 8, yend=8, x=12 ,xend=13), colour="deepskyblue4",size=0.7, linetype="dotted")+
    annotate(geom="text", x=c(2.5,14.5), y=11, label="T",color="black")+
    annotate(geom="text", x=c(6.5,18.5), y=11, label="pH",color="black")+
    annotate(geom="text", x=c(10.5,22.5), y=11, label="T x pH",color="black")

#### FIGURE 3 - BOTTOM PANNELS - version 3 (version 2 but in %) ####
signif_dev <- signif_dev[signif_dev$Scenario!="all",]
signif_dev <- signif_dev %>% 
  arrange(taxa,Stressor,Scenario)%>% 
  add_column(ID3=rep(seq(1,nrow(signif_dev),1)))
signif_dev$IDend[signif_dev$taxa=="fish"] <- signif_dev$ID2[signif_dev$taxa=="fish"]+0.9
signif_dev$IDend[signif_dev$taxa=="inverts"] <- signif_dev$ID2[signif_dev$taxa=="inverts"]+1
signif_dev$ID2[signif_dev$taxa=="fish"] <- signif_dev$ID2[signif_dev$taxa=="fish"]-0.1
  
  # calculating % of significant deviations and direction changes 
signif_dev$signif_direction_prop <- signif_dev$signif_direction/signif_dev$total*100
signif_dev$signif_deviation_prop <- signif_dev$significant/signif_dev$total*100

  
  ## plot
ggplot() +
  # significant deviation 
  geom_bar(data=signif_dev,position="dodge", stat="identity",
           mapping=aes(x=ID3, y=signif_deviation_prop,fill=taxa),alpha=0.5)+
  scale_fill_manual(values=c("#E69F00","aquamarine4"))+
  ylab("% metric")+
  coord_cartesian(ylim = c(0, 100),clip = 'off')+ # This focuses the x-axis on the range of interest
  xlab("")+
  scale_x_continuous(breaks=seq(1,18,1),labels=c(rep(c("R6","R8.5", "ex"),6)),position="top")+
  # outline of each barplot
  geom_bar(data=signif_dev,position="dodge", stat="identity",
           mapping=aes(x=ID3, y=100, col=taxa), alpha=0, size=0.1)+
  theme_classic()+
  # directional change only 
  geom_bar(data=signif_dev,position="dodge", stat="identity",
           mapping=aes(x=ID3, y=signif_direction_prop, col=taxa, fill=taxa),
           size=0.001)+
  scale_color_manual(values=c("black", "black"), labels=c("fish","inverts"))+
  
  #geom_vline(xintercept = c(9.5), color = "black", size=0.4)+
  geom_vline(xintercept = seq(3.5,18,3), linetype="dashed", color = "black", size=0.2)+
  annotate(geom="text", x=c(1.5,14.5), y=-13, label="T",color="black")+
  annotate(geom="text", x=c(4.5,18.5), y=-13, label="pCO2",color="black")+
  annotate(geom="text", x=c(7.5,22.5), y=-13, label="T x pCO2",color="black")+
  annotate(geom="text", x=signif_dev$ID3, y=103, label=signif_dev$total,color="black",
          size=3)+
  scale_y_reverse()


#### 5. FIGURE 4 - BOTTOM PANEL ####
 signif_dir$total <- abs(signif_dir$total)
 signif_dir$significant <- abs(signif_dir$significant)
 signif_dir$xlab <- paste(signif_dir$StressorScenario, signif_dir$taxa)
 signif_dir$direction <- factor(signif_dir$direction,levels=c("positive","negative"))
 signif_dir <- signif_dir[signif_dir$Scenario!="all",]
 signif_dir$StressorScenario  <- factor(signif_dir$StressorScenario,levels=c("T RCP6","T RCP8", "T extreme",
                                                                             "pH RCP6", "pH RCP8", "pH extreme",
                                                                             "TpH RCP6","TpH RCP8","TpH extreme"))
 # calculating % of significant changes 
 signif_dir$signif_prop <- signif_dir$significant/signif_dir$total*100
 
 ## For fish 
 fish_signif_dir <- signif_dir[signif_dir$taxa=="fish",]
 
 ggplot() +
    # nb of significant directional changes 
    geom_bar(data=fish_signif_dir,position="stack", stat="identity",
             mapping=aes(x=StressorScenario, y=signif_prop, fill=direction), size=0)+
    scale_fill_manual(values=c("#638FFF","#CC6577"))+
    scale_y_continuous(limits=c(0,8), breaks = seq(0,8,2),
                       labels=c("0","2","4","6","8"))+
    ylab("effect on biological responses")+
    xlab("")+
    theme_classic()+
    scale_x_discrete(position="top",labels=c(rep(c("RCP6","RCP8", "extreme"),3)))+
    scale_y_reverse()+
    # outlines of total metrics tested
    geom_bar(data=fish_signif_dir,position="dodge", stat="identity",
             mapping=aes(x=StressorScenario, y=100), alpha=0, size=0.3,color="grey")
    
 ## For inverts
inverts_signif_dir <- signif_dir[signif_dir$taxa=="inverts",]

ggplot() +
    # nb of significant directional changes 
    geom_bar(data=inverts_signif_dir,position="stack", stat="identity",
             mapping=aes(x=StressorScenario, y=signif_prop, fill=direction), size=0)+
  scale_fill_manual(values=c("#638FFF","#CC6577"))+
    ylab("effect on biological responses")+
    xlab("")+
    theme_classic()+
    scale_x_discrete(position="top",labels=c(rep(c("RCP6","RCP8", "extreme"),3)))+
    scale_y_reverse()+
    # outlines of total metrics tested
    geom_bar(data=inverts_signif_dir,position="dodge", stat="identity",
             mapping=aes(x=StressorScenario, y=100), alpha=0, size=0.3,color="grey")
 