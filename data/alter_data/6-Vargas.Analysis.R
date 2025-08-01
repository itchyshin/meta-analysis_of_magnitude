library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(metafor)
library(forestplot)
library(ggstar)

vargas.data <- read_csv("DataMetaClimar_Nov23Vargas.csv") %>%
  filter(status=="in") %>%
  filter(Stressor=="pH") %>%
  mutate(ei=as.numeric(ei),vei=as.numeric(vei)) %>%
  mutate(vei=case_when(vei<0.0001 ~ 0.0001,
                       TRUE~vei))

### adding scenarios
# classic approach
vargas.data$Scenario[vargas.data$Stressor=="pH"] <- c("RCP8") 
vargas.data$Scenario[vargas.data$Stressor=="pH" & vargas.data$DiffpCO2<350] <- c("RCP6")
vargas.data$Scenario[vargas.data$Stressor=="pH" & vargas.data$DiffpCO2>=750] <- c("extreme")

## Vargas approach
# Adding vargas.Scenarios for pH treatment
vargas.data$vargas.Scenario[vargas.data$Stressor=="pH"] <- c("RCP8") 
vargas.data$vargas.Scenario[vargas.data$Stressor=="pH" & vargas.data$delta_pCO2_index<350] <- c("RCP6")
vargas.data$vargas.Scenario[vargas.data$Stressor=="pH" & vargas.data$delta_pCO2_index>=750] <- c("extreme")

## nb of studies kept
vargas.data %>%
  group_by(REF) %>%
  summarize (n=n())
length(unique(vargas.data$REF)) 

## distribution of data by metric category
vargas.data %>%
  group_by(Category) %>%
  summarize (datapoints=n(), studies=length(unique(REF)))

## distribution of data by metric category and scenario
vargas.data %>%
  group_by(Category,vargas.Scenario) %>%
  summarize (datapoints=n(), studies=length(unique(REF)))

#### Comparison of results ####
## distribution among climate scenarios 
vargas.data %>%
  group_by(Scenario) %>%
  summarize (n=n())

vargas.data %>%
  group_by(vargas.Scenario) %>%
  summarize (n=n())

vargas.data %>%
  mutate(unchanged.scenario=Scenario==vargas.Scenario)%>%
  group_by(unchanged.scenario) %>%
  summarize (n=n())

#### overall correlation between ei and delta pH ####
## classic method - relative
variab <- lm(data=vargas.data,ei~DiffpCO2)
summary(variab)

ggplot(data=vargas.data, aes(x=DiffpCO2, y=ei)) +
  geom_point()+
  geom_smooth(data=vargas.data, aes(x=DiffpCO2, y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 10, y = -2.5,
                                     label = paste("y = ", round(variab$coefficients[2], 4),
                                                   "x + ", round(variab$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(variab)[9]),2),
                                                   "\np =", round(summary(variab)$coefficients[2,4], 3)), 
                                     hjust = 0))+
  xlab("")

## classic method - absolute
variab <- lm(data=vargas.data,abs(ei)~DiffpCO2)
summary(variab)

ggplot(data=vargas.data, aes(x=DiffpCO2, y=abs(ei))) +
  geom_point()+
  geom_smooth(data=vargas.data, aes(x=DiffpCO2, y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 8, y = 2,
                                     label = paste("y = ", round(variab$coefficients[2], 4),
                                                   "x + ", round(variab$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(variab)[9]),2),
                                                   "\np =", round(summary(variab)$coefficients[2,4], 3)), 
                                     hjust = 0))+
  xlab("")

## vargas approach -relative
variab.vargas <- lm(data=vargas.data,ei~delta_pCO2_index)
summary(variab.vargas)

ggplot(data=vargas.data, aes(x=delta_pCO2_index, y=ei)) +
  geom_point()+
  geom_smooth(data=vargas.data, aes(x=delta_pCO2_index, y=ei),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 10, y = -2.5,
                                     label = paste("y = ", round(variab.vargas$coefficients[2], 4),
                                                   "x + ", round(variab.vargas$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(variab.vargas)[9]),2),
                                                   "\np =", round(summary(variab)$coefficients[2,4], 3)), 
                                     hjust = 0))+
  xlab("")

## vargas approach - absolute
variab.vargas <- lm(data=vargas.data,abs(ei)~delta_pCO2_index)
summary(variab.vargas)

ggplot(data=vargas.data, aes(x=delta_pCO2_index, y=abs(ei))) +
  geom_point()+
  geom_smooth(data=vargas.data, aes(x=delta_pCO2_index, y=abs(ei)),color="black",method = "lm", se = FALSE)+
  geom_text(size=2,color="black",aes(x = 10, y = 2,
                                     label = paste("y = ", round(variab.vargas$coefficients[2], 4),
                                                   "x + ", round(variab.vargas$coefficients[1], 2), "\nR^2 = ", round(as.numeric(summary(variab.vargas)[9]),2),
                                                   "\np =", round(summary(variab)$coefficients[2,4], 3)), 
                                     hjust = 0))+
  xlab("")

#### Forest plots- pCO2 index ####
MA_TpH_vargas <- function(metric,data){
  data <- data[Category==metric]
  #Sample size
  samplesize <- data %>%
    group_by(vargas.Scenario) %>%
    summarize(n=n())%>%
    mutate(vargas.Scenario=factor(vargas.Scenario, levels=c("extreme","RCP6","RCP8")))%>%
    arrange(vargas.Scenario)
  
    if (nrow(samplesize)>1){
      MAscenarioInt <- rma(ei,vei,data=data,mods=~vargas.Scenario)
      print(MAscenarioInt)
      MAscenario <- rma(ei,vei,data=data,mods=~vargas.Scenario-1)
      Stats <- data.frame(method=c("vargas"),N=sum(samplesize$n),metric=c(metric),
                          QM=c(MAscenarioInt[["QM"]]),QMdf=c(MAscenarioInt[["QMdf"]][1]),
                          QMp=c(MAscenarioInt[["QMp"]]),QE=c(MAscenarioInt[["QE"]]),QEdf=c(MAscenarioInt[["k"]]-1),
                          QEp=c(MAscenarioInt[["QEp"]]))
      Forest <- data.frame(method=c("vargas"),metric=c(metric),scenario=samplesize$vargas.Scenario,lowerCI=MAscenario$ci.lb,upperCI=MAscenario$ci.ub,
                           ei=MAscenario$b, n=samplesize$n)
      StatModerator <- bind_rows(stat.vargas, Stats)
      forestplot <- bind_rows(forestplot.vargas.rel,Forest)
      list <- list(stat.vargas=StatModerator,forestplot.vargas=forestplot)
    }
  return(list)
}
MA_TpH_classic <- function(metric,data){
  data <- data[Category==metric]
  #Sample size
  samplesize <- data %>%
    group_by(Scenario) %>%
    summarize(n=n())%>%
    mutate(Scenario=factor(Scenario, levels=c("extreme","RCP6","RCP8")))%>%
    arrange(Scenario)
  print(samplesize)
  
  if (nrow(samplesize)>1){
    MAscenarioInt <- rma(ei,vei,data=data,mods=~Scenario)
    MAscenario <- rma(ei,vei,data=data,mods=~Scenario-1)
    print(MAscenario)
    Stats <- data.frame(method=c("classic"),N=sum(samplesize$n),metric=c(metric),
                        QM=c(MAscenarioInt[["QM"]]),QMdf=c(MAscenarioInt[["QMdf"]][1]),
                        QMp=c(MAscenarioInt[["QMp"]]),QE=c(MAscenarioInt[["QE"]]),QEdf=c(MAscenarioInt[["k"]]-1),
                        QEp=c(MAscenarioInt[["QEp"]]))
    Forest <- tibble(method=c("classic"),metric=c(metric),scenario=samplesize$Scenario,lowerCI=MAscenario$ci.lb,upperCI=MAscenario$ci.ub,
                         ei=MAscenario$b, n=samplesize$n)
    StatModerator <- bind_rows(stat.classic, Stats)
    forestplot <- bind_rows(forestplot.classic.rel,Forest)
    list <- list(stat.classic=StatModerator,forestplot.classic=forestplot)
  }
  return(list)
}

MA_TpH_vargas_abs <- function(metric,data){
  data <- data[Category==metric]
  #Sample size
  samplesize <- data %>%
    group_by(vargas.Scenario) %>%
    summarize(n=n())%>%
    mutate(vargas.Scenario=factor(vargas.Scenario, levels=c("extreme","RCP6","RCP8")))%>%
    arrange(vargas.Scenario)
  
  if (nrow(samplesize)>1){
    MAscenarioInt <- rma(abs(ei),vei,data=data,mods=~vargas.Scenario)
    MAscenario <- rma(abs(ei),vei,data=data,mods=~vargas.Scenario-1)
    Stats <- data.frame(method=c("vargas"),N=sum(samplesize$n),metric=c(metric),
                        QM=c(MAscenarioInt[["QM"]]),QMdf=c(MAscenarioInt[["QMdf"]][1]),
                        QMp=c(MAscenarioInt[["QMp"]]),QE=c(MAscenarioInt[["QE"]]),QEdf=c(MAscenarioInt[["k"]]-1),
                        QEp=c(MAscenarioInt[["QEp"]]))
    Forest <- data.frame(method=c("vargas"),metric=c(metric),scenario=samplesize$vargas.Scenario,lowerCI=MAscenario$ci.lb,upperCI=MAscenario$ci.ub,
                         ei=MAscenario$b, n=samplesize$n)
    StatModerator <- bind_rows(stat.vargas.abs, Stats)
    forestplot <- bind_rows(forestplot.vargas.abs,Forest)
    list.abs <- list(stat.vargas=StatModerator,forestplot.vargas=forestplot)
  }
  return(list.abs)
}
MA_TpH_classic_abs <- function(metric,data){
  data <- data[Category==metric]
  #Sample size
  samplesize <- data %>%
    group_by(Scenario) %>%
    summarize(n=n())%>%
    mutate(Scenario=factor(Scenario, levels=c("extreme","RCP6","RCP8")))%>%
    arrange(Scenario)
  print(samplesize)
  
  if (nrow(samplesize)>1){
    MAscenarioInt <- rma(abs(ei),vei,data=data,mods=~Scenario)
    MAscenario <- rma(abs(ei),vei,data=data,mods=~Scenario-1)
    print(MAscenario)
    Stats <- data.frame(method=c("classic"),N=sum(samplesize$n),metric=c(metric),
                        QM=c(MAscenarioInt[["QM"]]),QMdf=c(MAscenarioInt[["QMdf"]][1]),
                        QMp=c(MAscenarioInt[["QMp"]]),QE=c(MAscenarioInt[["QE"]]),QEdf=c(MAscenarioInt[["k"]]-1),
                        QEp=c(MAscenarioInt[["QEp"]]))
    Forest <- tibble(method=c("classic"),metric=c(metric),scenario=samplesize$Scenario,lowerCI=MAscenario$ci.lb,upperCI=MAscenario$ci.ub,
                     ei=MAscenario$b, n=samplesize$n)
    StatModerator <- bind_rows(stat.classic.abs, Stats)
    forestplot <- bind_rows(forestplot.classic.abs,Forest)
    list.abs <- list(stat.classic=StatModerator,forestplot.classic=forestplot)
  }
  return(list.abs)
}

#### running for each metric category - vargas  ####
# initializing
stat.vargas <- data.frame(N=as.integer(),metric=c(),
                          QM=c(),QMdf=c(),QMp=c(),QE=c(),QEdf=c(),QEp=c())
forestplot.vargas.rel <- data.frame(metric=c(),lowerCI=c(),upperCI=c(),ei=c(), scenario=c())

#### running relative ####
list.vargas <- MA_TpH_vargas("Survival",data.table(vargas.data))
forestplot.vargas.rel <- list.vargas$forestplot.vargas
stat.vargas <- list.vargas$stat.vargas

list.vargas <- MA_TpH_vargas("Biomechanics",data.table(vargas.data))
forestplot.vargas.rel <- list.vargas$forestplot.vargas
stat.vargas <- list.vargas$stat

list.vargas <- MA_TpH_vargas("Behaviour",data.table(vargas.data))
forestplot.vargas.rel <- list.vargas$forestplot.vargas
stat.vargas <- list.vargas$stat

list.vargas <- MA_TpH_vargas("Calcification",data.table(vargas.data))
forestplot.vargas.rel <- list.vargas$forestplot.vargas
stat.vargas <- list.vargas$stat

list.vargas <- MA_TpH_vargas("Development",data.table(vargas.data))
forestplot.vargas.rel <- list.vargas$forestplot.vargas
stat.vargas <- list.vargas$stat

list.vargas <- MA_TpH_vargas("Growth",data.table(vargas.data))
forestplot.vargas.rel <- list.vargas$forestplot.vargas
stat.vargas <- list.vargas$stat

list.vargas <- MA_TpH_vargas("Physiology",data.table(vargas.data))
forestplot.vargas.rel <- list.vargas$forestplot.vargas
stat.vargas <- list.vargas$stat

list.vargas <- MA_TpH_vargas("Reproduction",data.table(vargas.data))
forestplot.vargas.rel <- list.vargas$forestplot.vargas
stat.vargas <- list.vargas$stat

list.vargas <- MA_TpH_vargas("Routine respiration",data.table(vargas.data))
forestplot.vargas.rel <- list.vargas$forestplot.vargas
stat.vargas <- list.vargas$stat

#### running absolute ####
stat.vargas.abs <- data.frame(N=as.integer(),metric=c(),
                              QM=c(),QMdf=c(),QMp=c(),QE=c(),QEdf=c(),QEp=c())
forestplot.vargas.abs <- data.frame(metric=c(),lowerCI=c(),upperCI=c(),ei=c(), scenario=c())

list.vargas.abs <- MA_TpH_vargas_abs("Survival",data.table(vargas.data))
forestplot.vargas.abs <- list.vargas.abs$forestplot.vargas
stat.vargas.abs <- list.vargas.abs$stat.vargas.vargas

list.vargas.abs <- MA_TpH_vargas_abs("Biomechanics",data.table(vargas.data))
forestplot.vargas.abs <- list.vargas.abs$forestplot.vargas
stat.vargas.abs <- list.vargas.abs$stat.vargas

list.vargas.abs <- MA_TpH_vargas_abs("Behaviour",data.table(vargas.data))
forestplot.vargas.abs <- list.vargas.abs$forestplot.vargas
stat.vargas.abs <- list.vargas.abs$stat.vargas

list.vargas.abs <- MA_TpH_vargas_abs("Calcification",data.table(vargas.data))
forestplot.vargas.abs <- list.vargas.abs$forestplot.vargas
stat.vargas.abs <- list.vargas.abs$stat.vargas

list.vargas.abs <- MA_TpH_vargas_abs("Development",data.table(vargas.data))
forestplot.vargas.abs <- list.vargas.abs$forestplot.vargas
stat.vargas.abs <- list.vargas.abs$stat.vargas

list.vargas.abs <- MA_TpH_vargas_abs("Growth",data.table(vargas.data))
forestplot.vargas.abs <- list.vargas.abs$forestplot.vargas
stat.vargas.abs <- list.vargas.abs$stat.vargas

list.vargas.abs <- MA_TpH_vargas_abs("Physiology",data.table(vargas.data))
forestplot.vargas.abs <- list.vargas.abs$forestplot.vargas
stat.vargas.abs <- list.vargas.abs$stat.vargas

list.vargas.abs <- MA_TpH_vargas_abs("Reproduction",data.table(vargas.data))
forestplot.vargas.abs <- list.vargas.abs$forestplot.vargas
stat.vargas.abs <- list.vargas.abs$stat.vargas

list.vargas.abs <- MA_TpH_vargas_abs("Routine respiration",data.table(vargas.data))
forestplot.vargas.abs <- list.vargas.abs$forestplot.vargas
stat.vargas.abs <- list.vargas.abs$stat.vargas

#### running for each metric category - classic ####
# initializing
stat.classic <- data.frame(method=c(),N=as.integer(),metric=c(),
                           QM=c(),QMdf=c(),QMp=c(),QE=c(),QEdf=c(),QEp=c())
forestplot.classic.rel <- data.frame(method=c(),metric=c(),scenario=c(), lowerCI=c(),upperCI=c(),ei=c(), n=c())

# running 
list.classic <- MA_TpH_classic("Survival",data.table(vargas.data))
forestplot.classic.rel <- list.classic$forestplot.classic
stat.classic <- list.classic$stat

list.classic <- MA_TpH_classic("Biomechanics",data.table(vargas.data))
forestplot.classic.rel <- list.classic$forestplot.classic
stat.classic <- list.classic$stat

list.classic <- MA_TpH_classic("Behaviour",data.table(vargas.data))
forestplot.classic.rel <- list.classic$forestplot.classic
stat.classic <- list.classic$stat

list.classic <- MA_TpH_classic("Calcification",data.table(vargas.data))
forestplot.classic.rel <- list.classic$forestplot.classic
stat.classic <- list.classic$stat

list.classic <- MA_TpH_classic("Development",data.table(vargas.data))
forestplot.classic.rel <- list.classic$forestplot.classic
stat.classic <- list.classic$stat

list.classic <- MA_TpH_classic("Growth",data.table(vargas.data))
forestplot.classic.rel <- list.classic$forestplot.classic
stat.classic <- list.classic$stat

list.classic <- MA_TpH_classic("Physiology",data.table(vargas.data))
forestplot.classic.rel <- list.classic$forestplot.classic
stat.classic <- list.classic$stat

list.classic <- MA_TpH_classic("Reproduction",data.table(vargas.data))
forestplot.classic.rel <- list.classic$forestplot.classic
stat.classic <- list.classic$stat

list.classic <- MA_TpH_classic("Routine respiration",data.table(vargas.data))
forestplot.classic.rel <- list.classic$forestplot.classic
stat.classic <- list.classic$stat

# initializing absolute
stat.classic.abs <- data.frame(method=c(),N=as.integer(),metric=c(),
                           QM=c(),QMdf=c(),QMp=c(),QE=c(),QEdf=c(),QEp=c())
forestplot.classic.abs <- data.frame(method=c(),metric=c(),scenario=c(), lowerCI=c(),upperCI=c(),ei=c(), n=c())

# running 
list.classic.abs <- MA_TpH_classic_abs("Survival",data.table(vargas.data))
forestplot.classic.abs <- list.classic.abs$forestplot
stat.classic.abs <- list.classic.abs$stat

list.classic.abs <- MA_TpH_classic_abs("Biomechanics",data.table(vargas.data))
forestplot.classic.abs <- list.classic.abs$forestplot
stat.classic.abs <- list.classic.abs$stat

list.classic.abs <- MA_TpH_classic_abs("Behaviour",data.table(vargas.data))
forestplot.classic.abs <- list.classic.abs$forestplot
stat.classic.abs <- list.classic.abs$stat

list.classic.abs <- MA_TpH_classic_abs("Calcification",data.table(vargas.data))
forestplot.classic.abs <- list.classic.abs$forestplot
stat.classic.abs <- list.classic.abs$stat

list.classic.abs <- MA_TpH_classic_abs("Development",data.table(vargas.data))
forestplot.classic.abs <- list.classic.abs$forestplot
stat.classic.abs <- list.classic.abs$stat

list.classic.abs <- MA_TpH_classic_abs("Growth",data.table(vargas.data))
forestplot.classic.abs <- list.classic.abs$forestplot
stat.classic.abs <- list.classic.abs$stat

list.classic.abs <- MA_TpH_classic_abs("Physiology",data.table(vargas.data))
forestplot.classic.abs <- list.classic.abs$forestplot
stat.classic.abs <- list.classic.abs$stat

list.classic.abs <- MA_TpH_classic_abs("Reproduction",data.table(vargas.data))
forestplot.classic.abs <- list.classic.abs$forestplot
stat.classic.abs <- list.classic.abs$stat

list.classic.abs <- MA_TpH_classic_abs("Routine respiration",data.table(vargas.data))
forestplot.classic.abs <- list.classic.abs$forestplot
stat.classic.abs <- list.classic.abs$stat

#### PLOT: Classic, relative ####
forestplot.classic.rel <- forestplot.classic.rel %>%
  mutate(scenario=factor(scenario,levels=c("extreme","RCP8","RCP6"))) %>%
  arrange(metric, scenario)%>%
  mutate(ID=seq(1:nrow(forestplot.classic.rel)),
         significant=case_when(lowerCI*upperCI>0 ~ "yes",
                                TRUE ~ "no"))

ggplot(data=forestplot.classic.rel, aes(x=ID, y=ei[,1], ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6, aes(fill=significant, color=scenario),stroke=0.5,shape=21) + 
  scale_fill_manual(values=c("white","black"))+
  xlab("") + ylab("effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(angle = 45,colour="black"),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_continuous(breaks=seq(1.5,28,3),
                  labels=c(unique(forestplot.classic.rel$metric)))+
  #coord_cartesian(xlim = c(2.5, 52),ylim=c(-0.2,0.8),clip = 'off')+
  #annotate(geom="text",x=c(c(5,11,14,17,23,29,35,41,47,54)-4),y=-0.6,label=unique(generalei$metric), size=3.5,hjust = 0)+
  geom_vline(xintercept=seq(3.5,25,3), color="grey",size=0.3)+
  geom_text(aes(label=n,x=ID,y=0.8), size=3)+
  geom_segment(aes(y=0,yend=0,x=1,xend=26), lty=2, size=0.2)

#### Classic, absolute ####
forestplot.classic.abs <- forestplot.classic.abs %>%
  mutate(scenario=factor(scenario,levels=c("extreme","RCP8","RCP6"))) %>%
  arrange(metric, scenario)%>%
  mutate(ID=seq(1:nrow(forestplot.classic.abs)),
         significant=case_when(lowerCI*upperCI>0 ~ "yes",
                               TRUE ~ "no"))

ggplot(data=forestplot.classic.abs, aes(x=ID, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6, aes(fill=significant, color=scenario),stroke=0.5,shape=21) + 
  scale_fill_manual(values=c("white","black"))+
  xlab("") + ylab("abs(lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(angle = 45,colour="black"),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_continuous(breaks=seq(1.5,28,3),
                     labels=c(unique(forestplot.classic.abs$metric)))+
  geom_vline(xintercept=seq(3.5,25,3), color="grey",size=0.3)+
  geom_text(aes(label=n,x=ID,y=0.8), size=3)+
  geom_segment(aes(y=0,yend=0,x=1,xend=26), lty=2, size=0.2)

#### Vargas, relative ####
forestplot.vargas.rel <- forestplot.vargas.rel %>%
  mutate(scenario=factor(scenario,levels=c("extreme","RCP8","RCP6"))) %>%
  arrange(metric, scenario)%>%
  mutate(ID=seq(1:nrow(forestplot.vargas.rel)),
         significant=case_when(lowerCI*upperCI>0 ~ "yes",
                               TRUE ~ "no"))

ggplot(data=forestplot.vargas.rel, aes(x=ID, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6, aes(fill=significant, color=scenario),stroke=0.5,shape=21) + 
  scale_fill_manual(values=c("white","black"))+
  xlab("") + ylab("effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(angle = 45,colour="black"),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_continuous(breaks=seq(1.5,28,3),
                     labels=c(unique(forestplot.vargas.rel$metric)))+
  #coord_cartesian(xlim = c(2.5, 52),ylim=c(-0.2,0.8),clip = 'off')+
  #annotate(geom="text",x=c(c(5,11,14,17,23,29,35,41,47,54)-4),y=-0.6,label=unique(generalei$metric), size=3.5,hjust = 0)+
  geom_vline(xintercept=seq(3.5,25,3), color="grey",size=0.3)+
  geom_text(aes(label=n,x=ID,y=0.8), size=3)+
  geom_segment(aes(y=0,yend=0,x=1,xend=28), lty=2, size=0.2)

#### Vargas, absolute ####
forestplot.vargas.abs <- forestplot.vargas.abs %>%
  mutate(scenario=factor(scenario,levels=c("extreme","RCP8","RCP6"))) %>%
  arrange(metric, scenario)%>%
  mutate(ID=seq(1:nrow(forestplot.vargas.abs)),
         significant=case_when(lowerCI*upperCI>0 ~ "yes",
                               TRUE ~ "no"))

ggplot(data=forestplot.vargas.abs, aes(x=ID, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6, aes(fill=significant, color=scenario),stroke=0.5,shape=21) + 
  scale_fill_manual(values=c("white","black"))+
  xlab("") + ylab("effect size - abs(lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(angle = 45,colour="black"),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,1),"lines"))+
  scale_x_continuous(breaks=seq(1.5,28,3),
                     labels=c(unique(forestplot.vargas.abs$metric)))+
  #coord_cartesian(xlim = c(2.5, 52),ylim=c(-0.2,0.8),clip = 'off')+
  #annotate(geom="text",x=c(c(5,11,14,17,23,29,35,41,47,54)-4),y=-0.6,label=unique(generalei$metric), size=3.5,hjust = 0)+
  geom_vline(xintercept=seq(3.5,25,3), color="grey",size=0.3)+
  geom_text(aes(label=n,x=ID,y=0.8), size=3)+
  geom_segment(aes(y=0,yend=0,x=1,xend=28), lty=2, size=0.2)
