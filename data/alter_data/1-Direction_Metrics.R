#### This script:
# 1. attribute climate scenarios to each experimental study and creates file to
    # be analyzed in following scripts
# 2. graph of nb of metrics and directionality for each category

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(metafor)
library(forestplot)

#### 0.1 Creating Invert analysis file (selecting "in" status + scenarios) ####
Inverts <- read_csv("DataInverts.csv")
Inverts <- data.table(Inverts)
Inverts <- Inverts[Inverts$status=="in",]
Inverts$ei <- as.numeric(Inverts$ei)
Inverts$vei <- as.numeric(Inverts$vei)
Inverts[Inverts$vei<0.0001,]$vei <- 0.0001
# Adding Scenarios for temperature treatment
Inverts$Scenario[Inverts$Stressor=="T"] <- c("RCP8") 
Inverts$Scenario[Inverts$Stressor=="T" & Inverts$DiffT<2] <- c("RCP6")
Inverts$Scenario[Inverts$Stressor=="T" & Inverts$DiffT>=4] <- c("extreme")
# Adding Scenarios for pH treatment
Inverts$Scenario[Inverts$Stressor=="pH"] <- c("RCP8") 
Inverts$Scenario[Inverts$Stressor=="pH" & Inverts$DiffpCO2<350] <- c("RCP6")
Inverts$Scenario[Inverts$Stressor=="pH" & Inverts$DiffpCO2>=750] <- c("extreme")
# Adding Scenarios for pH x T treatment
Inverts$Scenario[Inverts$Stressor=="TpH"] <- c(NA) 
Inverts$Scenario[Inverts$Stressor=="TpH" & Inverts$DiffpCO2<350 & Inverts$DiffT<2] <- c("RCP6")
Inverts$Scenario[Inverts$Stressor=="TpH" & Inverts$DiffpCO2>=750 & Inverts$DiffT>=4] <- c("extreme")
Inverts$Scenario[Inverts$Stressor=="TpH" & Inverts$DiffpCO2<750 & Inverts$DiffpCO2>=350 &
                   Inverts$DiffT>=2 & Inverts$DiffT<4] <- c("RCP8")
write_csv(Inverts,"file_Inverts.csv")

#### 0.2 Creating Fish analysis file (selecting "in" status + scenarios) ####
Fish <- read.csv("DataFish.csv")
Fish <- data.table(Fish)
Fish <- Fish[Fish$status=="in",]
# Adding Scenarios for temperature treatment
Fish$Scenario[Fish$Stressor=="T"] <- c("RCP8") 
Fish$Scenario[Fish$Stressor=="T" & Fish$DiffT<2] <- c("RCP6")
Fish$Scenario[Fish$Stressor=="T" & Fish$DiffT>=4] <- c("extreme")
# Adding Scenarios for pH treatment
Fish$Scenario[Fish$Stressor=="pH"] <- c("RCP8") 
Fish$Scenario[Fish$Stressor=="pH" & Fish$DiffpCO2<350] <- c("RCP6")
Fish$Scenario[Fish$Stressor=="pH" & Fish$DiffpCO2>=750] <- c("extreme")
# Adding Scenarios for pH x T treatment
Fish$Scenario[Fish$Stressor=="TpH"] <- c(NA) 
Fish$Scenario[Fish$Stressor=="TpH" & Fish$DiffpCO2<350 & Fish$DiffT<2] <- c("RCP6")
Fish$Scenario[Fish$Stressor=="TpH" & Fish$DiffpCO2>=750 & Fish$DiffT>=4] <- c("extreme")
Fish$Scenario[Fish$Stressor=="TpH" & Fish$DiffpCO2<750 & Fish$DiffpCO2>=350 &
                Fish$DiffT>=2 & Fish$DiffT<4] <- c("RCP8")
write_csv(Fish,"file_Fish.csv")

#### 1.1 Directional metrics Fish #### 
## overall
Fish <- data.table(read_csv("file_Fish.csv"))

## creating data file to review papers based on Vargas 2022
#vargas.fish <- tibble(Fish) %>% 
#  group_by(REF,Species,PaperYear,PaperAuthor) %>%
#  summarize(nb=n(),Taxa="Fish")


direction_metrics <- Fish[which(!duplicated(Fish$TransformedMetric)),"MetricType"]
nb_directions_fish <- direction_metrics[,.(.N),by=.(MetricType)]

## within each metric category
direction_category_fish <- Fish[which(!duplicated(Fish$TransformedMetric)),c("MetricType","Category")]
direction_category_fish[direction_category_fish$Category %in% c("Routine respiration","Aerobic scope respiration"),
                   "Category"] <- "Metabolism"
nb_directions_cat_fish <- direction_category_fish[,.(.N),by=.(MetricType,Category)]
colnames(nb_directions_cat_fish) <- c("direction","Category","N")
nb_directions_cat_fish$direction <- factor(nb_directions_cat_fish$direction, 
                                            levels=c("ambiguous","positive","negative"))
nb_metrics <- direction_category_fish[,.(.N),by=.(Category)]
  

ggplot(data=nb_directions_cat_fish) +
  geom_bar(position="stack", stat="identity",
           mapping=aes(x=Category, y=N, fill=direction))+
  theme_classic()+
  scale_fill_manual(values=c("seashell2","greenyellow","orangered"))+
  annotate(geom="text", x=nb_metrics$Category, y=nb_metrics$N+1,
           label=nb_metrics$N,color="black")+
  ylab("# metrics included")
  
#### 1.2 Directional metrics Invert #### 
Inverts <- data.table(read_csv("file_Inverts.csv"))

vargas.inverts <- tibble(Inverts) %>% 
  group_by(REF, Species, PaperAuthor, PaperYear) %>%
  summarize(nb=n(),Taxa="Inverts")

#vargas.joined <- rbind(vargas.fish,vargas.inverts) %>%
#  arrange(REF)
#write.csv(vargas.joined,"vargas.dataset.csv")


direction_metrics <- Inverts[which(!duplicated(Inverts$TransformedMetric)),"Metric Type"]
colnames(direction_metrics) <- "MetricType"
nb_directions_inverts <- direction_metrics[,.(.N),by=.(MetricType)]

## within each metric category
direction_category_invert <- Inverts[which(!duplicated(Inverts$TransformedMetric)),c("Metric Type","Category")]
direction_category_invert[direction_category_invert$Category %in% c("Routine respiration","Aerobic scope respiration"),
                   "Category"] <- "Metabolism"
colnames(direction_category_invert) <- c("direction","Category")
nb_directions_cat_invert <- direction_category_invert[,.(.N),by=.(direction,Category)]
nb_directions_cat_invert$direction <- factor(nb_directions_cat_invert$direction, 
                                           levels=c("ambiguous","positive","negative"))
nb_metrics <- direction_category_invert[,.(.N),by=.(Category)]


ggplot(data=nb_directions_cat_invert) +
  geom_bar(position="stack", stat="identity",
           mapping=aes(x=Category, y=N, fill=direction))+
  theme_classic()+
  scale_fill_manual(values=c("seashell2","greenyellow","orangered"))+
  annotate(geom="text", x=nb_metrics$Category, y=nb_metrics$N+1,
           label=nb_metrics$N,color="black")+
  ylab("# metrics included")

#### 1.3 - FIGURE 2 ####
overall_direction_category <- rbind(Inverts[which(!duplicated(Inverts$TransformedMetric)),c("TransformedMetric","Category","Metric Type")],
                                    Fish[which(!duplicated(Fish$TransformedMetric)),c("TransformedMetric","Category","MetricType")],
                                    use.names=F)
  # all metrics to lower case to avoid duplicates from upper case
overall_direction_category$TransformedMetric <- tolower(overall_direction_category$TransformedMetric)
  # remove duplicates and NA
overall_direction_category <- overall_direction_category[!duplicated(overall_direction_category$TransformedMetric),]
overall_direction_category[overall_direction_category$Category %in% c("Routine respiration","Aerobic scope respiration"),
                          "Category"] <- "Metabolism"
colnames(overall_direction_category) <- c("metric","Category","direction")
overall_direction_category <- overall_direction_category[!is.na(overall_direction_category$direction),]
nb_overall_direction_category <- overall_direction_category[,.(.N),by=.(direction,Category)]
nb_overall_direction_category$direction  <- factor(nb_overall_direction_category$direction, 
                                             levels=c("ambiguous","positive","negative"))
nb_metrics <- overall_direction_category %>%
  mutate(Category=case_when(Category=="Behaviour"~"Behavior",
                            TRUE~as.character(Category))) %>%
  group_by(Category) %>%
  summarise(N=n()) %>%
  mutate(Category=factor(Category)) %>%
  arrange(desc(Category))


nb_overall_direction_category <- nb_overall_direction_category %>%
  mutate(Category=case_when(Category=="Behaviour"~"Behavior",
                   TRUE~as.character(Category))) %>%
  mutate(Category=factor(Category)) %>%
  arrange(desc(Category))

ggplot(data=nb_overall_direction_category) +
  geom_bar(position="stack", stat="identity",
           mapping=aes(x=Category, y=N, fill=direction),width=0.8)+
  theme(panel.background = element_rect(fill ="white"),
    axis.text.y = element_text(size = 10,colour="black"),
    axis.text.x = element_text(size = 10,colour="black"))+
  scale_fill_manual(values=c("honeydew3","#638FFF","#CC6577"))+
  annotate(geom="text", x=nb_metrics$Category, y=nb_metrics$N+1,
           label=nb_metrics$N,color="black",size=3)+
  ylab("Number of metrics per biological response category")+
  xlab(" ")+
  coord_flip()+
  labs(fill="Direction of the metric")+
  scale_y_continuous(position="right", labels=NULL,breaks=NULL)+
  scale_x_discrete(limits = rev(levels(nb_overall_direction_category$Category)))


#### 3. Verifying effect size and variance #### 
Inverts$ei_cal <- log(Inverts$TreatmentMeanTransformed/Inverts$ControlMeanTransformed)
Inverts[`Metric Type`=="negative","ei_cal"] <- log(Inverts[`Metric Type`=="negative","ControlMeanTransformed"]/Inverts[`Metric Type`=="negative","TreatmentMeanTransformed"])
which(abs(Inverts$ei-Inverts$ei_cal)>0.1) #the only differences are lines where ratios where inverted because double negative

Inverts$vei_cal <- Inverts$ControlSE^2/Inverts$ControlMeanTransformed^2+Inverts$TreatmentSE^2/Inverts$TreatmentMeanTransformed^2
which(abs(Inverts$vei-Inverts$vei_cal)>0.1) #no differences

Fish$ei_cal <- log(Fish$TreatmentMeanTransformed/Fish$ControlMeanTransformed)
Fish[MetricType=="negative","ei_cal"] <- log(Fish[MetricType=="negative","ControlMeanTransformed"]/Fish[MetricType=="negative","TreatmentMeanTransformed"])
which(abs(Fish$ei-Fish$ei_cal)>0.1) #the only differences are lines where ratios where inverted because double negative

Fish$vei_cal <- Fish$ControlSE^2/Fish$ControlMeanTransformed^2+Fish$TreatmentSE^2/Fish$TreatmentMeanTransformed^2
which(abs(Fish$vei-Fish$vei_cal)>0.1) #no differences
