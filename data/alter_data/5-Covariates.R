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
Fish <- data.table(read_csv("file_Fish.csv"))%>%
  mutate(LifeStage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(LifeStage)

stat.lifestage <- tibble(taxa=c(),stressor=c(),QM=c(),QMdf=c(),QMp=c(),QE=c(),QEdf=c(),QEp=c())
signif.diff.lifestage <- tibble(lifestage=c(),N=c(),taxa=c(),stressor=c(),
                                estimate=c(),se=c(),pval=c(),zval=c())

#### 1. Inverts ####
Inverts[Inverts$LifeStage=="JUVENILE/ADULT",]$LifeStage <- "JUVENILE"
Inverts[Inverts$LifeStage=="POSTLARVA",]$LifeStage <- "JUVENILE"
Inverts <- Inverts %>%
  mutate(LifeStage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(LifeStage)

## relative
MAlifestage <- rma(ei,vei,data=Inverts,mods=~LifeStage)
MAlifestage2 <- rma(ei,vei,data=Inverts,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2$ci.lb,upperCI=MAlifestage2$ci.ub,
                          ei=MAlifestage2$b, lifestage=rownames(MAlifestage2$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

# sample size
samplesize <- Inverts %>%
  group_by(LifeStage) %>%
  summarize(n=n())%>%
  mutate(lifestage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(lifestage)


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(-0.25,0,0.25))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  geom_text(aes(x=seq(1.2,4.2,1),y=ei), size=2.5,label=samplesize$n)+
  coord_flip(ylim=c(-0.25,0.25),clip='off')


## absolute
MAlifestage.abs <- rma(abs(ei),vei,data=Inverts,mods=~LifeStage)
MAlifestage2.abs <- rma(abs(ei),vei,data=Inverts,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2.abs$ci.lb,upperCI=MAlifestage2.abs$ci.ub,
                          ei=MAlifestage2.abs$b, lifestage=rownames(MAlifestage2.abs$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT")))%>%
  arrange(lifestage)

ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (abs(lnRR))") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,7),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                  labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(0,0.3,0.6))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  coord_flip(ylim=c(-0.1,0.6),clip='off')+
  geom_text(aes(x=rev(c(seq(1.1,4.1,1))),y=ei-0.04), size=2.5,label=samplesize$n)

#### 1.1. pH ####
Inverts[Inverts$LifeStage=="JUVENILE/ADULT",]$LifeStage <- "JUVENILE"
Inverts[Inverts$LifeStage=="POSTLARVA",]$LifeStage <- "JUVENILE"

InvertspH <- Inverts[Inverts$Stressor=="pH",]
MAlifestage <- rma(ei,vei,data=InvertspH,mods=~LifeStage)
MAlifestage2 <- rma(ei,vei,data=InvertspH,mods=~LifeStage-1)

# sample size
samplesize <- InvertspH %>%
  group_by(LifeStage) %>%
  summarize(n=n())%>%
  mutate(LifeStage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(LifeStage)

## adding to statistic summary table
summary.inv.pH <-tibble(taxa=c("invertebrate"),stressor=c("OA"),QM=c(MAlifestage[["QM"]]),QMdf=c(MAlifestage[["QMdf"]][1]),
                 QMp=c(MAlifestage[["QMp"]]),QE=c(MAlifestage[["QE"]]),QEdf=c(MAlifestage[["k"]]-3),QEp=c(MAlifestage[["QEp"]]))

lifestage.diff.inv.pH <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("invertebrates"),stressor=c("OA"),
                         estimate=c(MAlifestage[["b"]]),se=c(MAlifestage[["se"]]),pval=c(MAlifestage[["pval"]]),zval=c(MAlifestage[["zval"]]))

## effect size 
effect_size <- data.frame(lowerCI=MAlifestage2$ci.lb,upperCI=MAlifestage2$ci.ub,
                          ei=MAlifestage2$b, lifestage=rownames(MAlifestage2$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(-0.25,0,0.25))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  geom_text(aes(x=seq(1.2,4.2,1),y=ei), size=2.5,label=samplesize$n)+
  coord_flip(ylim=c(-0.25,0.25),clip='off')

## add acclimatation variable 
acclim <- lm(data=InvertspH,ei~sqrt(AcclimationDays))
summary(acclim)

acclim2 <- lm(data=InvertspH,ei~sqrt(AcclimationDays)+LifeStage)
summary(acclim2)

ggplot(data=InvertspH, aes(x=sqrt(AcclimationDays), y=ei,color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OA")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=InvertspH, aes(x=sqrt(AcclimationDays), y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 25, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))


## absolute
MAlifestage.abs <- rma(abs(ei),vei,data=InvertspH,mods=~LifeStage)
MAlifestage2.abs <- rma(abs(ei),vei,data=InvertspH,mods=~LifeStage-1)

summary.inv.pH.abs <-tibble(taxa=c("invertebrate"),stressor=c("OA"),QM=c(MAlifestage.abs[["QM"]]),QMdf=c(MAlifestage.abs[["QMdf"]][1]),
                        QMp=c(MAlifestage.abs[["QMp"]]),QE=c(MAlifestage.abs[["QE"]]),QEdf=c(MAlifestage.abs[["k"]]-3),QEp=c(MAlifestage.abs[["QEp"]]))

lifestage.diff.inv.pH.abs <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("invertebrates"),stressor=c("OA"),
                                estimate=c(MAlifestage.abs[["b"]]),se=c(MAlifestage.abs[["se"]]),pval=c(MAlifestage.abs[["pval"]]),zval=c(MAlifestage.abs[["zval"]]))


effect_size <- data.frame(lowerCI=MAlifestage2.abs$ci.lb,upperCI=MAlifestage2.abs$ci.ub,
                          ei=MAlifestage2.abs$b, lifestage=rownames(MAlifestage2.abs$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (abs(lnRR))") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(0,0.3,0.6))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  coord_flip(ylim=c(-0.1,0.6),clip='off')+
  geom_text(aes(x=c(seq(1.1,4.1,1)),y=ei-0.04), size=2.5,label=samplesize$n)

## add acclimation variable 
acclim <- lm(data=InvertspH,abs(ei)~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=InvertspH, aes(x=sqrt(AcclimationDays), y=abs(ei),color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OA")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=InvertspH, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(color="black",size=2,aes(x = 25, y = 2.5, label = paste("y = ", round(acclim$coefficients[2], 2), "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2)), hjust = 0))

#### 1.1.1. Calcification x pH ####
inv.cal.pH <- InvertspH[InvertspH$Category=="Calcification",]

## acclimation variable x relative ei
acclim <- lm(data=inv.cal.pH,ei~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=inv.cal.pH, aes(x=sqrt(AcclimationDays), y=ei,color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OA (calcification)")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=inv.cal.pH, aes(x=sqrt(AcclimationDays), y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))

## acclimation variable x absolute ei
acclim <- lm(data=inv.cal.pH,abs(ei)~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=inv.cal.pH, aes(x=sqrt(AcclimationDays), y=abs(ei),color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OA")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=inv.cal.pH, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))


#### 1.2. T ####
InvertsT <- Inverts[Inverts$Stressor=="T",]
MAlifestage <- rma(ei,vei,data=InvertsT,mods=~LifeStage)
MAlifestage2 <- rma(ei,vei,data=InvertsT,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2$ci.lb,upperCI=MAlifestage2$ci.ub,
                          ei=MAlifestage2$b, lifestage=rownames(MAlifestage2$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

# sample size
samplesize <- InvertsT %>%
  group_by(LifeStage) %>%
  summarize(n=n())%>%
  mutate(lifestage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(lifestage)

## adding to statistic summary table
summary.inv.T <-tibble(taxa=c("invertebrate"),stressor=c("OW"),QM=c(MAlifestage[["QM"]]),QMdf=c(MAlifestage[["QMdf"]][1]),
                 QMp=c(MAlifestage[["QMp"]]),QE=c(MAlifestage[["QE"]]),QEdf=c(MAlifestage[["k"]]-3),QEp=c(MAlifestage[["QEp"]]))

lifestage.diff.inv.T <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("invertebrates"),stressor=c("OW"),
                         estimate=c(MAlifestage[["b"]]),se=c(MAlifestage[["se"]]),pval=c(MAlifestage[["pval"]]),zval=c(MAlifestage[["zval"]]))

ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(-0.25,0,0.25))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  geom_text(aes(x=seq(1.2,4.2,1),y=ei), size=2.5,label=samplesize$n)+
  coord_flip(ylim=c(-0.25,0.25),clip='off')

## add acclimation variable 
acclim <- lm(data=InvertsT,ei~sqrt(AcclimationDays))
summary(acclim)

acclim2 <- lm(data=InvertsT,ei~sqrt(AcclimationDays)+LifeStage)
summary(acclim2)

ggplot(data=InvertsT, aes(x=sqrt(AcclimationDays), y=ei,color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OW")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=InvertsT, aes(x=sqrt(AcclimationDays), y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 25, y = 2.5,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))

## absolute
MAlifestage.abs <- rma(abs(ei),vei,data=InvertsT,mods=~LifeStage)
MAlifestage2.abs <- rma(abs(ei),vei,data=InvertsT,mods=~LifeStage-1)

summary.inv.T.abs <-tibble(taxa=c("invertebrate"),stressor=c("OW"),QM=c(MAlifestage.abs[["QM"]]),QMdf=c(MAlifestage.abs[["QMdf"]][1]),
                       QMp=c(MAlifestage.abs[["QMp"]]),QE=c(MAlifestage.abs[["QE"]]),QEdf=c(MAlifestage.abs[["k"]]-3),QEp=c(MAlifestage.abs[["QEp"]]))

lifestage.diff.inv.T.abs <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("invertebrates"),stressor=c("OW"),
                               estimate=c(MAlifestage.abs[["b"]]),se=c(MAlifestage.abs[["se"]]),pval=c(MAlifestage.abs[["pval"]]),zval=c(MAlifestage.abs[["zval"]]))



effect_size <- data.frame(lowerCI=MAlifestage2.abs$ci.lb,upperCI=MAlifestage2.abs$ci.ub,
                          ei=MAlifestage2.abs$b, lifestage=rownames(MAlifestage2.abs$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (abs(lnRR))") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(0,0.3,0.6))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  coord_flip(ylim=c(-0.1,0.6),clip='off')+
  geom_text(aes(x=c(seq(1.1,4.1,1)),y=ei-0.04), size=2.5,label=samplesize$n)

## add acclimation variable 
acclim <- lm(data=InvertsT,abs(ei)~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=InvertsT, aes(x=sqrt(AcclimationDays), y=abs(ei),color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OW")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=InvertsT, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 25, y = 2.5,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))

#### 1.2.3. TpH ####
InvertsTpH <- Inverts[Inverts$Stressor=="TpH",]%>%
  mutate(LifeStage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(LifeStage)

MAlifestage <- rma(ei,vei,data=InvertsTpH,mods=~LifeStage)
MAlifestage2 <- rma(ei,vei,data=InvertsTpH,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2$ci.lb,upperCI=MAlifestage2$ci.ub,
                          ei=MAlifestage2$b, lifestage=rownames(MAlifestage2$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

# sample size
samplesize <- InvertsTpH %>%
  group_by(LifeStage) %>%
  summarize(n=n())%>%
  mutate(lifestage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(lifestage)

summary.inv.TpH <-tibble(taxa=c("invertebrate"),stressor=c("TpH"),QM=c(MAlifestage[["QM"]]),QMdf=c(MAlifestage[["QMdf"]][1]),
                       QMp=c(MAlifestage[["QMp"]]),QE=c(MAlifestage[["QE"]]),QEdf=c(MAlifestage[["k"]]-3),QEp=c(MAlifestage[["QEp"]]))

lifestage.diff.inv.TpH <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("invertebrates"),stressor=c("TpH"),
                               estimate=c(MAlifestage[["b"]]),se=c(MAlifestage[["se"]]),pval=c(MAlifestage[["pval"]]),zval=c(MAlifestage[["zval"]]))


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(-0.25,0,0.25))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  geom_text(aes(x=seq(1.2,4.2,1),y=ei), size=2.5,label=samplesize$n)+
  coord_flip(ylim=c(-0.25,0.25),clip='off')

InvertsTpH.adult <- InvertsTpH[InvertsTpH$LifeStage=="ADULT",][c(0:50),]
MA.InvertsTpH.adult <-  rma(ei,vei,data=InvertsTpH.adult)
forest(MA.InvertsTpH.adult,slab = paste(InvertsTpH.adult$REF,InvertsTpH.adult$PaperYear,
                       substr(InvertsTpH.adult$PaperAuthor,1,20), InvertsTpH.adult$LifeStage),
       showweights = TRUE)

## add acclimatation variable 
acclim <- lm(data=InvertsTpH,ei~sqrt(AcclimationDays))
summary(acclim)

acclim2 <- lm(data=InvertsTpH,ei~sqrt(AcclimationDays)+LifeStage)
summary(acclim2)

ggplot(data=InvertsTpH, aes(x=sqrt(AcclimationDays), y=ei,color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OA x OW")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=InvertsTpH, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 3,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0)) 
## absolute
MAlifestage.abs <- rma(abs(ei),vei,data=InvertsTpH,mods=~LifeStage)
MAlifestage2.abs <- rma(abs(ei),vei,data=InvertsTpH,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2.abs$ci.lb,upperCI=MAlifestage2.abs$ci.ub,
                          ei=MAlifestage2.abs$b, lifestage=rownames(MAlifestage2.abs$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

summary.inv.TpH.abs <-tibble(taxa=c("invertebrate"),stressor=c("TpH"),QM=c(MAlifestage.abs[["QM"]]),QMdf=c(MAlifestage.abs[["QMdf"]][1]),
                       QMp=c(MAlifestage.abs[["QMp"]]),QE=c(MAlifestage.abs[["QE"]]),QEdf=c(MAlifestage.abs[["k"]]-3),QEp=c(MAlifestage.abs[["QEp"]]))

lifestage.diff.inv.TpH.abs <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("invertebrates"),stressor=c("TpH"),
                               estimate=c(MAlifestage.abs[["b"]]),se=c(MAlifestage.abs[["se"]]),pval=c(MAlifestage.abs[["pval"]]),zval=c(MAlifestage.abs[["zval"]]))


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (abs(lnRR))") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(0,0.3,0.6))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  coord_flip(ylim=c(-0.1,0.6),clip='off')+
  geom_text(aes(x=c(seq(1.1,4.1,1)),y=ei-0.04), size=2.5,label=samplesize$n)

## add acclimatation variable 
acclim <- lm(data=InvertsTpH,abs(ei)~sqrt(AcclimationDays))
summary(acclim)

acclim2 <- lm(data=InvertsTpH,abs(ei)~sqrt(AcclimationDays)+LifeStage)
summary(acclim2)

ggplot(data=InvertsTpH, aes(x=sqrt(AcclimationDays), y=abs(ei),color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OA x OW")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=InvertsTpH, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 3,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0)) 


#### 1.1.2. Survival x TpH ####
inv.surv.TpH <- InvertsTpH[InvertsTpH$Category=="Survival",]

## acclimation variable x relative ei
acclim <- lm(data=inv.surv.TpH,ei~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=inv.surv.TpH, aes(x=sqrt(AcclimationDays), y=ei,color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OA + OW (survival)")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=inv.surv.TpH, aes(x=sqrt(AcclimationDays), y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 10, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))

## acclimation variable x absolute ei
acclim <- lm(data=inv.surv.TpH,abs(ei)~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=inv.surv.TpH, aes(x=sqrt(AcclimationDays), y=abs(ei),color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("invertebrates' response to OA")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=inv.surv.TpH, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 10, y = 3,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))


#### 2.1.Fish ####
## relative
MAlifestage <- rma(ei,vei,data=Fish,mods=~LifeStage)
MAlifestage2 <- rma(ei,vei,data=Fish,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2$ci.lb,upperCI=MAlifestage2$ci.ub,
                          ei=MAlifestage2$b, lifestage=rownames(MAlifestage2$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

# sample size
samplesize <- Fish %>%
  group_by(LifeStage) %>%
  summarize(n=n())%>%
  mutate(lifestage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(lifestage)


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(-0.25,0,0.25))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  geom_text(aes(x=seq(1.2,4.2,1),y=ei), size=2.5,label=samplesize$n)+
  coord_flip(ylim=c(-0.25,0.25),clip='off')


## absolute
MAlifestage.abs <- rma(abs(ei),vei,data=Fish,mods=~LifeStage)
MAlifestage2.abs <- rma(abs(ei),vei,data=Fish,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2.abs$ci.lb,upperCI=MAlifestage2.abs$ci.ub,
                          ei=MAlifestage2.abs$b, lifestage=rownames(MAlifestage2.abs$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (abs(lnRR))") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(0,0.3,0.6))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  coord_flip(ylim=c(-0.1,0.6),clip='off')+
  geom_text(aes(x=c(seq(1.1,4.1,1)),y=ei-0.04), size=2.5,label=samplesize$n)

#### 2.2.1. pH ####
FishpH <- Fish[Fish$Stressor=="pH",]
MAlifestage <- rma(ei,vei,data=FishpH,mods=~LifeStage)
MAlifestage2 <- rma(ei,vei,data=FishpH,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2$ci.lb,upperCI=MAlifestage2$ci.ub,
                          ei=MAlifestage2$b, lifestage=rownames(MAlifestage2$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

# sample size
samplesize <- FishpH %>%
  group_by(LifeStage) %>%
  summarize(n=n())%>%
  mutate(lifestage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(lifestage)

summary.fish.pH <-tibble(taxa=c("fish"),stressor=c("OA"),QM=c(MAlifestage[["QM"]]),QMdf=c(MAlifestage[["QMdf"]][1]),
                       QMp=c(MAlifestage[["QMp"]]),QE=c(MAlifestage[["QE"]]),QEdf=c(MAlifestage[["k"]]-3),QEp=c(MAlifestage[["QEp"]]))

lifestage.diff.fish.pH <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("fish"),stressor=c("OA"),
                               estimate=c(MAlifestage[["b"]]),se=c(MAlifestage[["se"]]),pval=c(MAlifestage[["pval"]]),zval=c(MAlifestage[["zval"]]))


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(-0.25,0,0.25))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  geom_text(aes(x=seq(1.2,4.2,1),y=ei), size=2.5,label=samplesize$n)+
  coord_flip(ylim=c(-0.25,0.25),clip='off')

## add acclimation variable 
acclim <- lm(data=FishpH,ei~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=FishpH, aes(x=sqrt(AcclimationDays), y=ei,color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("fish response to OA")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=FishpH, aes(x=sqrt(AcclimationDays), y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))

## absolute
MAlifestage.abs <- rma(abs(ei),vei,data=FishpH,mods=~LifeStage)
MAlifestage2.abs <- rma(abs(ei),vei,data=FishpH,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2.abs$ci.lb,upperCI=MAlifestage2.abs$ci.ub,
                          ei=MAlifestage2.abs$b, lifestage=rownames(MAlifestage2.abs$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

summary.fish.pH.abs <-tibble(taxa=c("fish"),stressor=c("OA"),QM=c(MAlifestage.abs[["QM"]]),QMdf=c(MAlifestage.abs[["QMdf"]][1]),
                         QMp=c(MAlifestage.abs[["QMp"]]),QE=c(MAlifestage.abs[["QE"]]),QEdf=c(MAlifestage.abs[["k"]]-3),QEp=c(MAlifestage.abs[["QEp"]]))

lifestage.diff.fish.pH.abs <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("fish"),stressor=c("OA"),
                                 estimate=c(MAlifestage.abs[["b"]]),se=c(MAlifestage.abs[["se"]]),pval=c(MAlifestage.abs[["pval"]]),zval=c(MAlifestage.abs[["zval"]]))


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (abs(lnRR))") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(0,0.3,0.6))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  coord_flip(ylim=c(-0.1,0.6),clip='off')+
  geom_text(aes(x=c(seq(1.1,4.1,1)),y=ei-0.04), size=2.5,label=samplesize$n)

## add acclimation variable 
acclim <- lm(data=FishpH,abs(ei)~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=FishpH, aes(x=sqrt(AcclimationDays), y=abs(ei),color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("fish response to OA")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=FishpH, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))

#### 2.2.2. T ####
FishT <- Fish[Fish$Stressor=="T",]
MAlifestage <- rma(ei,vei,data=FishT,mods=~LifeStage)
MAlifestage2 <- rma(ei,vei,data=FishT,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2$ci.lb,upperCI=MAlifestage2$ci.ub,
                          ei=MAlifestage2$b, lifestage=rownames(MAlifestage2$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

# sample size
samplesize <- FishpH %>%
  group_by(LifeStage) %>%
  summarize(n=n())%>%
  mutate(lifestage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(lifestage)

summary.fish.T <-tibble(taxa=c("fish"),stressor=c("OW"),QM=c(MAlifestage[["QM"]]),QMdf=c(MAlifestage[["QMdf"]][1]),
                         QMp=c(MAlifestage[["QMp"]]),QE=c(MAlifestage[["QE"]]),QEdf=c(MAlifestage[["k"]]-3),QEp=c(MAlifestage[["QEp"]]))

lifestage.diff.fish.T <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("fish"),stressor=c("OW"),
                                 estimate=c(MAlifestage[["b"]]),se=c(MAlifestage[["se"]]),pval=c(MAlifestage[["pval"]]),zval=c(MAlifestage[["zval"]]))


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(-0.25,0,0.25))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  geom_text(aes(x=seq(1.2,4.2,1),y=ei), size=2.5,label=samplesize$n)+
  coord_flip(ylim=c(-0.25,0.25),clip='off')

## add acclimation variable 
acclim <- lm(data=FishT,ei~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=FishT, aes(x=sqrt(AcclimationDays), y=ei,color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("fish response to OW")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=FishT, aes(x=sqrt(AcclimationDays), y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))


## absolute
MAlifestage.abs <- rma(abs(ei),vei,data=FishT,mods=~LifeStage)
MAlifestage2.abs <- rma(abs(ei),vei,data=FishT,mods=~LifeStage-1)

summary.fish.T.abs <-tibble(taxa=c("fish"),stressor=c("OW"),QM=c(MAlifestage.abs[["QM"]]),QMdf=c(MAlifestage.abs[["QMdf"]][1]),
                         QMp=c(MAlifestage.abs[["QMp"]]),QE=c(MAlifestage.abs[["QE"]]),QEdf=c(MAlifestage.abs[["k"]]-3),QEp=c(MAlifestage.abs[["QEp"]]))

lifestage.diff.fish.T.abs <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("fish"),stressor=c("OW"),
                                 estimate=c(MAlifestage.abs[["b"]]),se=c(MAlifestage.abs[["se"]]),pval=c(MAlifestage.abs[["pval"]]),zval=c(MAlifestage.abs[["zval"]]))

effect_size <- data.frame(lowerCI=MAlifestage2.abs$ci.lb,upperCI=MAlifestage2.abs$ci.ub,
                          ei=MAlifestage2.abs$b, lifestage=rownames(MAlifestage2.abs$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (abs(lnRR))") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(0,0.3,0.6))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  coord_flip(ylim=c(-0.1,0.6),clip='off')+
  geom_text(aes(x=c(seq(1.1,4.1,1)),y=ei-0.04), size=2.5,label=samplesize$n)

## add acclimation variable 
acclim <- lm(data=FishT,abs(ei)~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=FishT, aes(x=sqrt(AcclimationDays), y=abs(ei),color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("fish response to OW")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=FishT, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))


#### 2.3. TpH ####
FishTpH <- Fish[Fish$Stressor=="TpH",]
MAlifestage <- rma(ei,vei,data=FishTpH,mods=~LifeStage)
MAlifestage2 <- rma(ei,vei,data=FishTpH,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2$ci.lb,upperCI=MAlifestage2$ci.ub,
                          ei=MAlifestage2$b, lifestage=rownames(MAlifestage2$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

# sample size
samplesize <- FishTpH %>%
  group_by(LifeStage) %>%
  summarize(n=n())%>%
  mutate(lifestage=factor(LifeStage, levels=c("EMBRYO","LARVA","JUVENILE","ADULT")))%>%
  arrange(lifestage)

summary.fish.TpH <-tibble(taxa=c("fish"),stressor=c("TpH"),QM=c(MAlifestage[["QM"]]),QMdf=c(MAlifestage[["QMdf"]][1]),
                         QMp=c(MAlifestage[["QMp"]]),QE=c(MAlifestage[["QE"]]),QEdf=c(MAlifestage[["k"]]-3),QEp=c(MAlifestage[["QEp"]]))

lifestage.diff.fish.TpH <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("fish"),stressor=c("TpH"),
                                 estimate=c(MAlifestage[["b"]]),se=c(MAlifestage[["se"]]),pval=c(MAlifestage[["pval"]]),zval=c(MAlifestage[["zval"]]))


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(-0.25,0,0.25))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  geom_text(aes(x=seq(1.2,4.2,1),y=ei), size=2.5,label=samplesize$n)+
  coord_flip(ylim=c(-0.25,0.25),clip='off')

## add acclimation variable 
acclim <- lm(data=FishTpH,ei~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=FishTpH, aes(x=sqrt(AcclimationDays), y=ei,color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("fish response to OA + OW")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=FishTpH, aes(x=sqrt(AcclimationDays), y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 1,
            label = paste("y = ", round(acclim$coefficients[2], 2),
            "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
            "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))


## absolute
MAlifestage.abs <- rma(abs(ei),vei,data=FishTpH,mods=~LifeStage)
MAlifestage2.abs <- rma(abs(ei),vei,data=FishTpH,mods=~LifeStage-1)

effect_size <- data.frame(lowerCI=MAlifestage2.abs$ci.lb,upperCI=MAlifestage2.abs$ci.ub,
                          ei=MAlifestage2.abs$b, lifestage=rownames(MAlifestage2.abs$beta)) %>%
  mutate(lifestage=factor(lifestage, levels=c("LifeStageEMBRYO","LifeStageLARVA","LifeStageJUVENILE","LifeStageADULT"))) %>%
  arrange(lifestage)

summary.fish.TpH.abs <-tibble(taxa=c("fish"),stressor=c("TpH"),QM=c(MAlifestage.abs[["QM"]]),QMdf=c(MAlifestage.abs[["QMdf"]][1]),
                          QMp=c(MAlifestage.abs[["QMp"]]),QE=c(MAlifestage.abs[["QE"]]),QEdf=c(MAlifestage.abs[["k"]]-3),QEp=c(MAlifestage.abs[["QEp"]]))

lifestage.diff.fish.TpH.abs <- tibble(lifestage=samplesize$LifeStage,N=samplesize$n,taxa=c("fish"),stressor=c("TpH"),
                                  estimate=c(MAlifestage.abs[["b"]]),se=c(MAlifestage.abs[["se"]]),pval=c(MAlifestage.abs[["pval"]]),zval=c(MAlifestage.abs[["zval"]]))


ggplot(data=effect_size, aes(x=lifestage, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6,stroke=0.5) + 
  xlab("") + ylab("Effect size (abs(lnRR))") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_discrete(expand = c(0.02,0.01),
                   labels=c("embryo","larvae","juvenile","adult"))+
  scale_y_continuous(position="right",breaks = c(0,0.3,0.6))+
  geom_segment(aes(y=0,yend=0,x=1,xend=4), lty=2, size=0.2)+   # adds a dotted line at x=0 after flip
  coord_flip(ylim=c(-0.1,0.6),clip='off')+
  geom_text(aes(x=c(seq(1.1,4.1,1)),y=ei-0.04), size=2.5,label=samplesize$n)

## add acclimation variable 
acclim <- lm(data=FishTpH,abs(ei)~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=FishTpH, aes(x=sqrt(AcclimationDays), y=abs(ei),color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("fish response to OA + OW")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=FishTpH, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 15, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))

#### 2.3.1. Survival x TpH ####
fish.surv.TpH <- FishTpH[FishTpH$Category=="Survival",]

## acclimation variable x relative ei
acclim <- lm(data=fish.surv.TpH,ei~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=fish.surv.TpH, aes(x=sqrt(AcclimationDays), y=ei,color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("fish response to OA + OW (survival)")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=fish.surv.TpH, aes(x=sqrt(AcclimationDays), y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 10, y = 1,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))

## acclimation variable x absolute ei
acclim <- lm(data=fish.surv.TpH,abs(ei)~sqrt(AcclimationDays))
summary(acclim)

ggplot(data=fish.surv.TpH, aes(x=sqrt(AcclimationDays), y=abs(ei),color=LifeStage)) +
  geom_point(size= 0.6)+
  #theme_classic()+
  ggtitle("fish response to OA + OW (survival)")+
  xlab("Days of acclimation (sqrt)")+
  scale_color_viridis_d(begin=0.9,end=0)+
  geom_smooth(data=fish.surv.TpH, aes(x=sqrt(AcclimationDays), y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 10, y = 2,
                                     label = paste("y = ", round(acclim$coefficients[2], 2),
                                                   "x + ", round(acclim$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(acclim)[9]),2),
                                                   "\np =", round(summary(acclim)$coefficients[2,4], 3)), 
                                     hjust = 0))
## summary statistics 
stat.summary.rel <- rbind(summary.fish.pH,summary.fish.T,summary.fish.TpH,
                        summary.inv.pH,summary.inv.T,summary.inv.TpH)
write.csv(stat.summary.rel,"lifestages.relative.stats.csv")

stat.summary.abs <- rbind(summary.fish.pH.abs,summary.fish.T.abs,summary.fish.TpH.abs,
                          summary.inv.pH.abs,summary.inv.T.abs,summary.inv.TpH.abs)
write.csv(stat.summary.abs,"lifestages.absolute.stats.csv")

diff.summary.rel <- rbind(lifestage.diff.fish.pH,lifestage.diff.fish.T,lifestage.diff.fish.TpH,
                          lifestage.diff.inv.pH,lifestage.diff.inv.T,lifestage.diff.inv.TpH)
write.csv(diff.summary.rel,"lifestages.relative.diffs.csv")

diff.summary.abs <- rbind(lifestage.diff.fish.pH.abs,lifestage.diff.fish.T.abs,lifestage.diff.fish.TpH.abs,
                          lifestage.diff.inv.pH.abs,lifestage.diff.inv.T.abs,lifestage.diff.inv.TpH.abs)
write.csv(diff.summary.abs,"lifestages.absolute.diffs.csv")


## nb of metrics across lifestage 
Fish.lifestage.metrics <- Fish %>%
  select(LifeStage,Metric,Category,MetricType)%>%
  group_by(MetricType,LifeStage) %>%
  summarize(nb.metric=n())

Invert.lifestage.metrics <- Inverts %>%
  mutate(MetricType='Metric Type') %>%
  select(LifeStage,'Metric Type',Category)%>%
  group_by('Metric Type',LifeStage) %>%
  summarize(nb.metric=n())