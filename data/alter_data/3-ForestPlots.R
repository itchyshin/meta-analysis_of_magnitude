library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(metafor)
library(forestplot)
library(ggstar)
forest_absolute <- read_csv("file_forestplot_absolute.csv")
forest_relative <- read_csv("file_forestplot_relative.csv")

#### 0. Inverts, relative, across scenarios (first attempts) ####
inverts_general <- forest_relative[forest_relative$taxa=="inverts" & forest_relative$Scenario=="all",]
inverts_general$Stressor <- factor(inverts_general$Stressor, levels=c("T","pH","TpH"))
inverts_general %>% 
  arrange(metric,Stressor) %>% 
  forestplot(labeltext = c(metric, Stressor, N), 
             mean=ei, lower=lowerCI,
             upper=upperCI ,
             pos=list("topright","inset" = .1,"align"="horizontal"),
             txt_gp = fpTxtGp(xlab=gpar(fontfamily="",fontsize=12),
                              label=gpar(fontfamily="",fontsize=8)))

## BETTER!! using ggplot
inverts_general$color[inverts_general$significant=="no"] <- "white"
inverts_general$color[inverts_general$significant=="yes"] <- "black"

inverts_general$Stressor <- factor(inverts_general$Stressor, levels=c("T","pH","TpH"))
inverts_general$ID <- seq(1,nrow(inverts_general),1)

ggplot(data=inverts_general, aes(x=ID, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(width=2,size= 0.5, fill=inverts_general$color, color="grey",shape=rep(c(21,22,23),10)) + 
  scale_fill_discrete(c("white","black"))+
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  #coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("lnRR (95% CI)") +
  theme_classic()+
  scale_x_continuous(breaks=seq(1,30,1), labels=inverts_general$Stressor)+
  annotate(geom="text",x=seq(2,31,3),y=-0.5,label=unique(inverts_general$metric), size=3)+
  geom_vline(xintercept=seq(3.5,32.5,3), color="grey",size=0.3)
  

#### 1.1. FIGURE 2 : forest plot of relative ei ####
generalei <- forest_relative[forest_relative$Scenario=="all",]
generalei$color[generalei$taxa=="fish"] <- "aquamarine4"
generalei$color[generalei$taxa=="inverts"] <- "chocolate"
generalei$fill[generalei$significant=="no"] <- "white"
generalei$fill[generalei$significant=="yes"] <- generalei$color[generalei$significant=="yes"]
generalei$shape[generalei$Stressor=="T"] <- 21
generalei$shape[generalei$Stressor=="pH"] <- 22
generalei$shape[generalei$Stressor=="TpH"] <- 23

generalei$Stressor <- factor(generalei$Stressor, levels=c("T","pH","TpH"))
generalei$taxa <- factor(generalei$taxa, levels=c("inverts","fish"))
generalei <- generalei %>% 
  arrange(metric,Stressor,taxa) %>% 
  add_column(ID=seq(1,nrow(generalei),1))

#pdf(file="plot_general_ei.pdf",width=10, height=6) 
ggplot(data=generalei, aes(x=ID, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6, fill=generalei$fill, color=generalei$color,shape=generalei$shape,stroke=0.4) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,7),"lines"))+
  scale_x_reverse(expand = c(0.02,0.01),breaks=generalei$ID[generalei$taxa=="inverts"],
                     labels=rep(c("OW","OA","OW + OA"),10))+
 #coord_cartesian(xlim = c(2.5, 52),ylim=c(-0.2,0.8),clip = 'off')+
  annotate(geom="text",x=c(c(5,11,14,17,23,29,35,41,47,54)-4),y=-1.35,label=unique(generalei$metric), size=3.5,hjust = 0)+
  geom_vline(xintercept=c(c(3,6,9,15,21,27,33,39,45)+3.5), color="grey",size=0.3)+
  labs(fill="taxa")+
  geom_text(aes(label=N,x=ID+0.4,y=ei-0.04), size=2.5)+
  scale_y_continuous(position="right",breaks = c(-0.4,0,0.4,0.8))+
  geom_segment(aes(y=0,yend=0,x=1,xend=54), lty=2, size=0.2)+   # adds a dotted line at x=1 after flip
  coord_flip(ylim=c(-0.7,0.8),clip='off')

#dev.off()

#### 1.2. FIGURE 4: forest plot of absolute ei ####
generalei <- forest_absolute[forest_relative$Scenario=="all",]
generalei$color[generalei$taxa=="fish"] <- "aquamarine4"
generalei$color[generalei$taxa=="inverts"] <- "chocolate"
generalei$fill[generalei$significant=="no"] <- "white"
generalei$fill[generalei$significant=="yes"] <- generalei$color[generalei$significant=="yes"]
generalei$shape[generalei$Stressor=="T"] <- 21
generalei$shape[generalei$Stressor=="pH"] <- 22
generalei$shape[generalei$Stressor=="TpH"] <- 23

generalei$Stressor <- factor(generalei$Stressor, levels=c("T","pH","TpH"))
generalei$taxa <- factor(generalei$taxa, levels=c("inverts","fish"))
generalei <- generalei %>% 
  arrange(metric,Stressor,taxa) %>% 
  add_column(ID=seq(1,nrow(generalei),1))

ggplot(data=generalei, aes(x=ID, y=ei, ymin=lowerCI, ymax=upperCI)) +
  geom_pointrange(size= 0.6, fill=generalei$fill, color=generalei$color,shape=generalei$shape,stroke=0.5) + 
  xlab("") + ylab("Effect size (lnRR)") + 
  theme(axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        panel.background = element_rect(fill ="white"),
        plot.margin = unit(c(1,1,1,7),"lines"))+
  scale_x_reverse(expand = c(0.02,0.01),breaks=generalei$ID[generalei$taxa=="inverts"],
                  labels=rep(c("T",expression(paste(italic("p"),"CO2")),
                               expression(paste("T x ",italic("p"),"CO2"))),10))+
  #coord_cartesian(xlim = c(2.5, 52),ylim=c(-0.2,0.8),clip = 'off')+
  annotate(geom="text",x=c(c(5,11,14,17,23,29,35,41,47,54)-4),y=-0.6,label=unique(generalei$metric), size=3.5,hjust = 0)+
  geom_vline(xintercept=c(c(3,6,9,15,21,27,33,39,45)+3.5), color="grey",size=0.3)+
  labs(fill="taxa")+
  geom_text(aes(label=N,x=ID+0.4,y=ei-0.04), size=2.5)+
  scale_y_continuous(position="right",breaks = c(-0.2,0,0.4,0.8))+
  geom_segment(aes(y=0,yend=0,x=1,xend=54), lty=2, size=0.2)+   # adds a dotted line at x=1 after flip
  coord_flip(ylim=c(-0.2,0.8),clip='off')


#### Inverts, relative, per scenario ####
## all effect size
relative_inverts <- forest_relative[forest_relative$taxa=="inverts",]

relative_inverts$Stressor <- factor(relative_inverts$Stressor, 
                                    levels=c("T","pH","TpH"))
relative_inverts$Scenario <- factor(relative_inverts$Scenario, 
                                    levels=c("all","RCP6","RCP8","extreme"))

pdf(file="inverts_ei.pdf",width=6, height=6)
relative_inverts %>% 
  arrange(metric,Stressor,Scenario) %>% 
  forestplot(labeltext = c(metric, Stressor, Scenario,N), 
             mean=ei, lower=lowerCI,
             upper=upperCI,
             txt_gp = fpTxtGp(xlab=gpar(fontfamily="",cex=0.5, fontsize=2),
                              label=gpar(fontfamily="serif",cex=0.5, fontsize=8)))
dev.off()

## significant effect size
relative_inverts_sig <- relative_inverts[relative_inverts$significant=="yes",]

pdf(file="inverts_ei_significant.pdf",width=6, height=6)
relative_inverts_sig %>% 
  arrange(metric,Stressor,Scenario) %>% 
  forestplot(labeltext = c(metric, Stressor, Scenario,N), 
             mean=ei, lower=lowerCI,
             upper=upperCI,
             txt_gp = fpTxtGp(xlab=gpar(fontfamily="",cex=0.5, fontsize=2),
                              label=gpar(fontfamily="serif",cex=0.5, fontsize=15)))
dev.off()

#### Fish, relative  ####
## all effect size
relative_fish <- forest_relative[forest_relative$taxa=="fish",]

relative_fish$Stressor <- factor(relative_fish$Stressor, 
                                    levels=c("T","pH","TpH"))
relative_fish$Scenario <- factor(relative_fish$Scenario, 
                                    levels=c("all","RCP6","RCP8","extreme"))

pdf(file="fish_ei.pdf",width=6, height=6)
relative_fish %>% 
  arrange(metric,Stressor,Scenario) %>% 
  forestplot(labeltext = c(metric, Stressor, Scenario,N), 
             mean=ei, lower=lowerCI,
             upper=upperCI,
             txt_gp = fpTxtGp(xlab=gpar(fontfamily="",cex=0.5, fontsize=2),
                              label=gpar(fontfamily="serif",cex=0.5, fontsize=8)))
dev.off()

## significant effect size
relative_fish_sig <- relative_fish[relative_fish$significant=="yes",]

pdf(file="fish_ei_significant.pdf",width=6, height=6)
relative_fish_sig %>% 
  arrange(metric,Stressor,Scenario) %>% 
  forestplot(labeltext = c(metric, Stressor, Scenario,N), 
             mean=ei, lower=lowerCI,
             upper=upperCI,
             txt_gp = fpTxtGp(xlab=gpar(fontfamily="",cex=0.5, fontsize=2),
                              label=gpar(fontfamily="serif",cex=0.5, fontsize=15)))
dev.off()



#### Inverts, absolute ####
## all effect size
absolute_inverts <- forest_absolute[forest_absolute$taxa=="inverts",]

absolute_inverts$Stressor <- factor(absolute_inverts$Stressor, 
                                    levels=c("T","pH","TpH"))
absolute_inverts$Scenario <- factor(absolute_inverts$Scenario, 
                                    levels=c("all","RCP6","RCP8","extreme"))

pdf(file="inverts_ei_absolute.pdf",width=6, height=6)
absolute_inverts %>% 
  arrange(metric,Stressor,Scenario) %>% 
  forestplot(labeltext = c(metric, Stressor, Scenario,N), 
             mean=ei, lower=lowerCI,
             upper=upperCI,
             txt_gp = fpTxtGp(xlab=gpar(fontfamily="",cex=0.5, fontsize=2),
                              label=gpar(fontfamily="serif",cex=0.5, fontsize=6)))
dev.off()

## significant effect size
absolute_inverts_sig <- absolute_inverts[absolute_inverts$significant=="yes",]

pdf(file="inverts_ei_absolute_significant.pdf",width=6, height=6)
absolute_inverts_sig %>% 
  arrange(metric,Stressor,Scenario) %>% 
  forestplot(labeltext = c(metric, Stressor, Scenario,N), 
             mean=ei, lower=lowerCI,
             upper=upperCI,
             txt_gp = fpTxtGp(xlab=gpar(fontfamily="",cex=0.5, fontsize=2),
                              label=gpar(fontfamily="serif",cex=0.5, fontsize=12)))
dev.off()

#### Fish, absolute  ####
## all effect size
absolute_fish <- forest_absolute[forest_absolute$taxa=="fish",]

absolute_fish$Stressor <- factor(absolute_fish$Stressor, 
                                 levels=c("T","pH","TpH"))
absolute_fish$Scenario <- factor(absolute_fish$Scenario, 
                                 levels=c("all","RCP6","RCP8","extreme"))

pdf(file="fish_ei_absolute.pdf",width=6, height=6)
absolute_fish %>% 
  arrange(metric,Stressor,Scenario) %>% 
  forestplot(labeltext = c(metric, Stressor, Scenario,N), 
             mean=ei, lower=lowerCI,
             upper=upperCI,
             txt_gp = fpTxtGp(xlab=gpar(fontfamily="",cex=0.5, fontsize=2),
                              label=gpar(fontfamily="serif",cex=0.5, fontsize=8)))
dev.off()

## significant effect size
absolute_fish_sig <- absolute_fish[absolute_fish$significant=="yes",]

pdf(file="plot_ei_absolute_significant.pdf",width=6, height=6)
absolute_fish_sig %>% 
  arrange(metric,Stressor,Scenario) %>% 
  forestplot(labeltext = c(metric, Stressor, Scenario,N), 
             mean=ei, lower=lowerCI,
             upper=upperCI,
             txt_gp = fpTxtGp(xlab=gpar(fontfamily="",cex=0.5, fontsize=2),
                              label=gpar(fontfamily="serif",cex=0.5, fontsize=12)))
dev.off()

#### 2.1.1 Heatmap - relative inverts ####
#### 2.1.1.0 - Past versions ####
relative_inverts <- forest_relative[forest_relative$taxa=="inverts",]
relative_inverts$StressorScenario <- paste(relative_inverts$Stressor,relative_inverts$Scenario)
relative_inverts$StressorScenario <- factor(relative_inverts$StressorScenario, 
                                 levels=c("T all","T RCP6","T RCP8","T extreme",
                                          "pH all","pH RCP6","pH RCP8","pH extreme",
                                          "TpH all","TpH RCP6","TpH RCP8","TpH extreme"))
relative_inverts[relative_inverts$significant=="no","significant2"] <- c(0)
relative_inverts[relative_inverts$significant=="yes","significant2"] <- c(1)
relative_inverts[relative_inverts$significant=="no","ei2"] <- c(0)
relative_inverts[relative_inverts$significant=="yes","ei2"] <- relative_inverts[relative_inverts$significant=="yes","ei"]

## v1. grey when not significant (using ei2)
colour_breaks <- c(min(relative_inverts$ei), -0.01,0,0.01, max(relative_inverts$ei))
colours <- c("#A51122","#FDF6B6", "grey","#F6F9B1", "#006228")  
#pdf(file="plot_heatmap_invert.pdf",width=10, height=6)
ggplot()+
  geom_tile(data=relative_inverts, aes(x=StressorScenario, y=metric,fill=ei2,alpha=factor(significant2)),
            color="white")+
  scale_fill_gradientn(limits=range(relative_inverts$ei),
                       colours=colours[c(1,seq_along(colours),length(colours))],
                       values=c(0, scales::rescale(colour_breaks, from = range(relative_inverts$ei)), 1))+
  scale_alpha_manual(values=c(0.3,1))+  
  ylab("")+
  xlab("")+
  theme_classic()+
  geom_vline(xintercept = c(4.5,8.5), color = "black", size=0.4)+
  scale_x_discrete(position = "top",labels=c(rep(c("all", "RCP6","RCP8", "extreme"),3)))+
  annotate(geom="text", x=2.5, y=11.1, label="T",color="black")+
  annotate(geom="text", x=6.5, y=11.1, label="pH",color="black")+
  annotate(geom="text", x=10.5, y=11.1, label="T x pH",color="black")+
  coord_cartesian(ylim = c(0, 10), clip = "off")



## v2. not transparency based on significance
colour_breaks <- c(min(relative_inverts$ei), -0.01,0.01, max(relative_inverts$ei))
colours <- c("#A51122","#FDF6B6","#F6F9B1", "#006228")  

ggplot(data=relative_inverts, aes(x=StressorScenario, y=metric,fill=ei))+
  geom_tile()+
  scale_fill_gradientn(limits=range(relative_inverts$ei),
                       colours=colours[c(1,seq_along(colours),length(colours))],
                       values=c(0, scales::rescale(colour_breaks, from = range(relative_inverts$ei)), 1))+
  ylab("metric category")+
  xlab("stressor x scenario")+
  theme_classic()+
  geom_vline(xintercept = c(4.5,8.5), color = "black", size=0.4)+
  scale_x_discrete(position = "top")

## v3. with transparent squares when not significant 
ggplot(data=relative_inverts, aes(x=StressorScenario, y=metric,fill=ei, alpha=significant2))+
  geom_tile()+
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"))+
  ylab("metric category")+
  xlab("stressor x scenario")+
  theme_classic()+
  geom_vline(xintercept = c(4.5,8.5), color = "black", size=0.4)+
  scale_x_discrete(position = "top")

#### 2.1.1.1. FIGURE 3 TOP LEFT - Final version  ####
relative_inverts <- forest_relative[forest_relative$taxa=="inverts",]
relative_inverts <- relative_inverts[relative_inverts$Scenario!="all",]
relative_inverts$StressorScenario <- paste(relative_inverts$Stressor,relative_inverts$Scenario)
relative_inverts$StressorScenario <- factor(relative_inverts$StressorScenario, 
                                            levels=c("T RCP6","T RCP8","T extreme",
                                                    "pH RCP6","pH RCP8","pH extreme",
                                                     "TpH RCP6","TpH RCP8","TpH extreme"))
relative_inverts$metric <- as.factor(relative_inverts$metric)
relative_inverts$metric <- factor(relative_inverts$metric,levels=rev(levels(relative_inverts$metric)))

colour_breaks <- c(min(forest_relative$ei), -0.02,0.02,0.15, max(forest_relative$ei))
colours <- c("#A51122","#FDF6B6","#F6F9B1", "#006228")  
colours2 <- c("#CC6577","#FDF6B6","#F6F9B1","#77C9F7", "#638FFF") 

ggplot()+
  geom_tile(data=relative_inverts, mapping=aes(x=StressorScenario, y=metric,fill=ei),colour="white",size=0.7)+
  scale_fill_gradientn(limits=range(forest_relative$ei),
                       colours=colours2[c(1,seq_along(colours2),length(colours2))],
                       values=c(0, scales::rescale(colour_breaks, from = range(forest_relative$ei)), 1))+
  ylab("")+
  xlab("")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 14,colour="black"),
        axis.text.x = element_text(size = 10,colour="black"),
        plot.margin = unit(c(2,1,1,1),"lines"))+
  geom_segment(aes(x = c(3.5), xend=3.5, y=0.5,yend=11.2), color = "black", size=0.4)+
  geom_segment(aes(x = c(6.5), xend=6.5, y=0.5,yend=11.2),color = "black", size=0.4)+
  scale_x_discrete(expand = c(0,0),position = "top",labels=c(rep(c("R6","R8", "ex"),3)))+
  scale_y_discrete(expand = c(0.055,0))+
  coord_cartesian(ylim=c(1,10),clip = "off")+
  geom_point(shape=8,mapping=aes(x=as.numeric(relative_inverts$StressorScenario[relative_inverts$significant=="yes"]),
           y=as.numeric(relative_inverts$metric[relative_inverts$significant=="yes"])), size=2)+
  annotate(geom="text", x=2, y=11.2, label="OW",color="black")+
  annotate(geom="text", x=5, y=11.2, label="OA",color="black")+
  annotate(geom="text", x=8, y=11.2, label="OW + OA",color="black")

#### 2.1.2. Heatmap - relative fish ####
#### 2.1.2.0 Past versions ####
relative_fish <- forest_relative[forest_relative$taxa=="fish",]
relative_fish$StressorScenario <- paste(relative_fish$Stressor,relative_fish$Scenario)
relative_fish$StressorScenario <- factor(relative_fish$StressorScenario, 
                                            levels=c("T all","T RCP6","T RCP8","T extreme",
                                                     "pH all","pH RCP6","pH RCP8","pH extreme",
                                                     "TpH all","TpH RCP6","TpH RCP8","TpH extreme"))
relative_fish[relative_fish$significant=="no","significant2"] <- c(0)
relative_fish[relative_fish$significant=="yes","significant2"] <- c(1)
relative_fish[relative_fish$significant=="no","ei2"] <- c(0)
relative_fish[relative_fish$significant=="yes","ei2"] <- relative_fish[relative_fish$significant=="yes","ei"]

## grey when not significant (using ei2)
#pdf(file="plot_heatmap_fish.pdf",width=10, height=6)
colour_breaks <- c(min(relative_fish$ei), -0.01,0,0.01, max(relative_fish$ei))
colours <- c("#A51122","#FDF6B6", "grey","#F6F9B1", "#006228")  
ggplot()+
    geom_tile(data=relative_fish, aes(x=StressorScenario, y=metric,fill=ei2,alpha=factor(significant2)),
              color="white")+
    scale_fill_gradientn(limits=range(relative_fish$ei),
                         colours=colours[c(1,seq_along(colours),length(colours))],
                         values=c(0, scales::rescale(colour_breaks, from = range(relative_fish$ei)), 1))+
  scale_alpha_manual(values=c(0.4,1))+  
  ylab("")+
    xlab("")+
    theme_classic()+
    geom_vline(xintercept = c(4.5,8.5), color = "black", size=0.4)+
    scale_x_discrete(position = "top",labels=c(rep(c("all", "RCP6","RCP8", "extreme"),3)))+
  annotate(geom="text", x=2.5, y=9.1, label="T",color="black")+
  annotate(geom="text", x=6.5, y=9.1, label="pH",color="black")+
  annotate(geom="text", x=10.5, y=9.1, label="T x pH",color="black")+
  coord_cartesian(ylim = c(0, 8), clip = "off")

  
  
## all ei shown regardless of significance
ggplot(data=relative_fish, aes(x=StressorScenario, y=metric,fill=ei))+
  geom_tile()+
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"))+
  ylab("metric category")+
  xlab("stressor x scenario")+
  theme_classic()+
  geom_vline(xintercept = c(4.5,8.5), color = "black", size=0.4)+
  scale_x_discrete(position = "top")

#### 2.1.2.1 FIGURE 3 TOP RIGHT - Final version ####
relative_fish <- forest_relative[forest_relative$taxa=="fish",]
relative_fish <- relative_fish[relative_fish$Scenario!="all",]
biomechanics <- c("T","RCP6",0,"fish","Biomechanics",0,0,0,0,"no")
biodiversity <- c("T","RCP6",0,"fish","Biodiversity",0,0,0,0,"no")
relative_fish <- rbind(relative_fish,biomechanics,biodiversity) %>%
  mutate(ei=as.numeric(ei))

relative_fish$StressorScenario <- paste(relative_fish$Stressor,relative_fish$Scenario)
relative_fish$StressorScenario <- factor(relative_fish$StressorScenario, 
                                         levels=c("T RCP6","T RCP8","T extreme",
                                                  "pH RCP6","pH RCP8","pH extreme",
                                              "TpH RCP6","TpH RCP8","TpH extreme"))
relative_fish$metric <- as.factor(relative_fish$metric)
relative_fish$metric <- factor(relative_fish$metric,levels=rev(levels(relative_fish$metric)))

colour_breaks <- c(min(relative_fish$ei), -0.02,0.02,0.15, max(relative_fish$ei))
colours <- c("#A51122","#FDF6B6","#F6F9B1", "#006228")  
colours2 <- c("#CC6577","#FDF6B6","#F6F9B1","#77C9F7", "#638FFF") 


ggplot()+
  geom_tile(data=relative_fish, mapping=aes(x=StressorScenario, y=metric,fill=ei),colour="white",size=0.7)+
  scale_fill_gradientn(limits=range(forest_relative$ei),
                       colours=colours2[c(1,seq_along(colours2),length(colours2))],
                       values=c(0, scales::rescale(colour_breaks, from = range(forest_relative$ei)), 1))+
  ylab("")+
  xlab("")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 14,colour="black"),
        axis.text.x = element_text(size = 10,colour="black"),
        plot.margin = unit(c(2,1,1,1),"lines"))+
  geom_segment(aes(x = c(3.5), xend=3.5, y=0.5,yend=11.2), color = "black", size=0.4)+
  geom_segment(aes(x = c(6.5), xend=6.5, y=0.5,yend=11.2),color = "black", size=0.4)+
  scale_x_discrete(expand = c(0,0),position = "top",labels=c(rep(c("R6","R8", "ex"),3)))+
  scale_y_discrete(expand = c(0.055,0))+
  coord_cartesian(ylim=c(1,10),clip = "off")+
  geom_point(shape=8,mapping=aes(x=as.numeric(relative_fish$StressorScenario[relative_fish$significant=="yes"]),
                                 y=as.numeric(relative_fish$metric[relative_fish$significant=="yes"])), size=2)+
  annotate(geom="text", x=2, y=11.2, label="OW",color="black")+
  annotate(geom="text", x=5, y=11.2, label="OA",color="black")+
  annotate(geom="text", x=8, y=11.2, label="OW + OA",color="black")




#### 2.2.1. Heatmap - absolute fish ####
absolute_fish <- forest_absolute[forest_absolute$taxa=="fish",]
absolute_fish$StressorScenario <- paste(absolute_fish$Stressor,absolute_fish$Scenario)
absolute_fish$StressorScenario <- factor(absolute_fish$StressorScenario, 
                                         levels=c("T all","T RCP6","T RCP8","T extreme",
                                                  "pH all","pH RCP6","pH RCP8","pH extreme",
                                                  "TpH all","TpH RCP6","TpH RCP8","TpH extreme"))
absolute_fish[absolute_fish$significant=="no","significant2"] <- c(0)
absolute_fish[absolute_fish$significant=="yes","significant2"] <- c(1)
absolute_fish[absolute_fish$significant=="no","ei2"] <- c(0)
absolute_fish[absolute_fish$significant=="yes","ei2"] <- absolute_fish[absolute_fish$significant=="yes","ei"]

colour_breaks <- c(0,0.01,max(absolute_fish$ei2))
colours <- c("gray88","lightskyblue", "mediumblue")  


#pdf(file="plot_heatmap_fish_absolute.pdf",width=10, height=6)
ggplot()+
  geom_tile(data=absolute_fish, aes(x=StressorScenario, y=metric,fill=ei2),
            color="white")+
  scale_fill_gradientn(limits=range(absolute_fish$ei2),
                       colours=c("gray88","lightskyblue", "mediumblue"),
                       values=c(0, scales::rescale(colour_breaks, from = range(absolute_fish$ei2)), 1))+
  ylab("")+
  xlab("")+
  theme_classic()+
  geom_vline(xintercept = c(4.5,8.5), color = "black", size=0.4)+
  scale_x_discrete(position = "top",labels=c(rep(c("all", "RCP6","RCP8", "extreme"),3)))+
  labs(fill="ei (lnRR)")+
  annotate(geom="text", x=3, y=9.1, label="T",color="black")+
  annotate(geom="text", x=7, y=9.1, label="pH",color="black")+
  annotate(geom="text", x=10.5, y=9.1, label="T x pH",color="black")+
  coord_cartesian(ylim = c(0, 8), clip = "off")
#dev.off()





#### 2.2.2. Heatmap - absolute invert ####
absolute_invert <- forest_absolute[forest_absolute$taxa=="inverts",]
absolute_invert$StressorScenario <- paste(absolute_invert$Stressor,absolute_invert$Scenario)
absolute_invert$StressorScenario <- factor(absolute_invert$StressorScenario, 
                                         levels=c("T all","T RCP6","T RCP8","T extreme",
                                                  "pH all","pH RCP6","pH RCP8","pH extreme",
                                                  "TpH all","TpH RCP6","TpH RCP8","TpH extreme"))
absolute_invert[absolute_invert$significant=="no","significant2"] <- c(0)
absolute_invert[absolute_invert$significant=="yes","significant2"] <- c(1)
absolute_invert[absolute_invert$significant=="no","ei2"] <- c(0)
absolute_invert[absolute_invert$significant=="yes","ei2"] <- absolute_invert[absolute_invert$significant=="yes","ei"]

colour_breaks <- c(0,0.01,max(absolute_invert$ei2))
colours <- c("gray98","lightskyblue", "mediumblue")  


#pdf(file="plot_heatmap_invert_absolute.pdf",width=10, height=6)
ggplot()+
  geom_tile(data=absolute_invert, aes(x=StressorScenario, y=metric,fill=ei2),
            color="white")+
  scale_fill_gradientn(limits=range(absolute_invert$ei2),
                       colours=c("gray88","lightskyblue", "mediumblue"),
                       values=c(0, scales::rescale(colour_breaks, from = range(absolute_invert$ei2)), 1))+
  ylab("")+
  xlab("")+
  theme_classic()+
  geom_vline(xintercept = c(4.5,8.5), color = "black", size=0.4)+
  scale_x_discrete(position = "top",labels=c(rep(c("all", "R6","R8", "ex"),3)))+
  labs(fill="ei (lnRR)")+
  annotate(geom="text", x=3, y=11.1, label="T",color="black")+
  annotate(geom="text", x=7, y=11.1, label="pH",color="black")+
  annotate(geom="text", x=10.5, y=11.1, label="T x pH",color="black")+
  coord_cartesian(ylim = c(0, 10), clip = "off")

dev.off()

#### 2.3. FIGURE 5 - TOP PANELS: Heatmap - difference absolute / relative ####
forest_absolute <- read_csv("file_forestplot_absolute.csv")
forest_relative <- read_csv("file_forestplot_relative.csv")

forest_absolute <- forest_absolute[forest_absolute$Scenario!="all",] 
forest_relative <- forest_relative[forest_relative$Scenario!="all",] 
forest_absolute$ID <- paste(forest_absolute$taxa,forest_absolute$metric,forest_absolute$Stressor,forest_absolute$Scenario)
forest_relative$ID <- paste(forest_relative$taxa,forest_relative$metric,forest_relative$Stressor,forest_relative$Scenario)
forest <- left_join(forest_absolute,forest_relative, by="ID",suffix=c(".abs",".rel"))
#forest <- forest[,-c(12,13,14,19,20)]

forest$diff[forest$significant.abs=="yes"&forest$significant.rel=="yes"] <- "significant change"
forest$diff[forest$significant.abs=="yes"&forest$significant.rel=="no"] <- "deviation only"
forest$diff[forest$significant.abs=="no"] <- "no effect"
forest$diff <- factor(forest$diff, levels=c("significant change", "deviation only", "no effect"))
forest$StressorScenario <- paste(forest$Stressor.abs,forest$Scenario.abs)
forest$StressorScenario <- factor(forest$StressorScenario, 
                                           levels=c("T RCP6","T RCP8","T extreme",
                                                    "pH RCP6","pH RCP8","pH extreme",
                                                    "TpH RCP6","TpH RCP8","TpH extreme"))

## fish
diff_fish <- forest[forest$taxa.abs=="fish",c("metric.abs","diff","StressorScenario")]
biomechanics <- c("Biomechanics","no effect","T RCP6")
biodiversity <- c("Biodiversity","no effect","T RCP6")
diff_fish<- rbind(diff_fish,biomechanics,biodiversity)

diff_fish$metric.abs <- as.factor(diff_fish$metric.abs)
diff_fish$metric.abs <- factor(diff_fish$metric.abs,levels=rev(levels(diff_fish$metric.abs)))

ggplot()+
  geom_tile(data=diff_fish, aes(x=StressorScenario, y=metric.abs,fill=diff),
            color="white",size=0.7)+
  scale_fill_manual(values=c("#268a73","#a3c9ba","white"))+
  ylab("")+
  xlab("")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 14,colour="black"),
        axis.text.x = element_text(size = 10,colour="black"),
        plot.margin = unit(c(2,1,1,1),"lines"))+
  geom_segment(aes(x = c(3.5), xend=3.5, y=0.5,yend=11.2), color = "black", size=0.4)+
  geom_segment(aes(x = c(6.5), xend=6.5, y=0.5,yend=11.2),color = "black", size=0.4)+
  scale_x_discrete(expand = c(0,0),position = "top",labels=c(rep(c("R6","R8", "ex"),3)))+
  scale_y_discrete(expand = c(0.055,0))+
  coord_cartesian(ylim=c(1,10),clip = "off")+
  annotate(geom="text", x=2, y=11.2, label="T",color="black",size = 4)+
  annotate(geom="text", x=5, y=11.2, label=expression(paste(italic("p"),"CO2")),color="black",size = 4)+
  annotate(geom="text", x=8, y=11.2, label=expression(paste("T x ",italic("p"),"CO2")),color="black",size = 4)+
  coord_cartesian(ylim = c(1, 10), clip = "off")

## invert
diff_invert <- forest[forest$taxa.abs=="inverts",]

diff_invert$metric.abs <- as.factor(diff_invert$metric.abs)
diff_invert$metric.abs <- factor(diff_invert$metric.abs,levels=rev(levels(diff_invert$metric.abs)))

ggplot()+
  geom_tile(data=diff_invert, aes(x=StressorScenario, y=metric.abs,fill=diff),
            color="white",size=0.7)+
  scale_fill_manual(values=c("#f39b04","#f9cb8c","white"))+
  ylab("")+
  xlab("")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size = 14,colour="black"),
        axis.text.x = element_text(size = 10,colour="black"),
        plot.margin = unit(c(2,1,1,1),"lines"))+
  geom_segment(aes(x = c(3.5), xend=3.5, y=0.5,yend=11.2), color = "black", size=0.4)+
  geom_segment(aes(x = c(6.5), xend=6.5, y=0.5,yend=11.2),color = "black", size=0.4)+
  scale_x_discrete(expand = c(0,0),position = "top",labels=c(rep(c("R6","R8", "ex"),3)))+
  scale_y_discrete(expand = c(0.055,0))+
  annotate(geom="text", x=2, y=11.2, label="T",color="black",size = 4)+
  annotate(geom="text", x=5, y=11.2, label=expression(paste(italic("p"),"CO2")),color="black",size = 4)+
  annotate(geom="text", x=8, y=11.2, label=expression(paste("T x ",italic("p"),"CO2")),color="black",size = 4)+
  coord_cartesian(ylim = c(1, 10), clip = "off")


