MA_TpH <- function(taxa,metric,data,sensitivity){
  data <- data[Category==metric]
  #Sample size
  sample1 <- data[,.(.N),by=.(Stressor)]
  StatModerator <- data.frame(Stressor=sample1$Stressor,N=sample1$N,taxa=c(taxa),metric=c(metric),
                              QM=c("NA"),QMdf=c("NA"),QMd=c("NA"),QE=c("NA"),QEdf=c("NA"),QEp=c("NA"))
  sample2 <- data[,.(.N),by=.(Stressor,Scenario)]
  sample2 <- sample2[order(Stressor,Scenario),]
  sample2 <- sample2[!which(is.na(sample2$Scenario)),]
  StatSimple <- data.frame(Stressor=sample2$Stressor,Scenario=sample2$Scenario,N=sample2$N,
                           taxa=c(taxa),metric=c(metric),pvalue=c("NA"),lowerCI=c("NA"),upperCI=c("NA"),
                           ei=c("NA"),Qt=c(""),df=c("NA"),pQt=c("NA"))
  for (k in (1:nrow(sample2))){
    MAdata <- data[Stressor==sample2[k,"Stressor"] & Scenario==sample2[k,"Scenario"]]
    MA <- rma(ei,vei,data=MAdata)
    StatSimple[k,"pvalue"] <- MA[["pval"]]
    StatSimple[k,"lowerCI"] <- MA[["ci.lb"]]
    StatSimple[k,"upperCI"] <- MA[["ci.ub"]]
    StatSimple[k,"ei"] <- MA[["b"]]
    StatSimple[k,"Qt"] <- MA[["QE"]]
    StatSimple[k,"df"] <- MA[["k"]]-1
    StatSimple[k,"pQt"] <- MA[["QEp"]]
  }
  for (k in c("T","pH","TpH")){
    MAdata <- data[Stressor==k]
    MA <- rma(ei,vei,data=MAdata)
    print(MA)
    row_values <- data.frame(Stressor=k,Scenario="all",N=MA[["k"]],taxa=taxa,
                             metric=metric,pvalue=MA[["pval"]],lowerCI=MA[["ci.lb"]],
                             upperCI=MA[["ci.ub"]],ei=MA[["b"]],Qt=MA[["QE"]],df=MA[["k"]]-1,pQt=MA[["QEp"]])
    StatSimple <- rbind(StatSimple,row_values)
    funnel(MA,xlab=paste(taxa, metric, k))
    forest(MA,slab = paste(MAdata$REF,MAdata$PaperYear,
                           substr(MAdata$PaperAuthor,1,20), MAdata$Scenario),
           showweights = TRUE, xlab=paste(taxa, metric, k))
    rosenthal <- fsn(ei,vei, data=MAdata)
    threshold <- nrow(MAdata)*5+10
    trimfill <- trimfill(rma(data=MAdata,yi=ei,vi=vei))
    sensitivity <- rbind(sensitivity,c(taxa,k,metric,rosenthal$fsnum,threshold,MA[["pval"]],trimfill$pval))
  }
  sensitivity <<-sensitivity 
  return(StatSimple)
}

MA_TpH_abs <- function(taxa,metric,data){
  data <- data[Category==metric]
  #Sample size
  sample1 <- data[,.(.N),by=.(Stressor)]
  sample2 <- data[,.(.N),by=.(Stressor,Scenario)]
  sample2 <- sample2[order(Stressor,Scenario),]
  sample2 <- sample2[!which(is.na(sample2$Scenario)),]
  StatSimple <- data.frame(Stressor=sample2$Stressor,Scenario=sample2$Scenario,N=sample2$N,
                           taxa=c(taxa),metric=c(metric),pvalue=c("NA"),lowerCI=c("NA"),upperCI=c("NA"),
                           ei=c("NA"),Qt=c(""),df=c("NA"),pQt=c("NA"))
  for (k in (1:nrow(sample2))){
    MAdata <- data[Stressor==sample2[k,"Stressor"] & Scenario==sample2[k,"Scenario"]]
    MA <- rma(abs(ei),vei,data=MAdata)
    StatSimple[k,"pvalue"] <- MA[["pval"]]
    StatSimple[k,"lowerCI"] <- MA[["ci.lb"]]
    StatSimple[k,"upperCI"] <- MA[["ci.ub"]]
    StatSimple[k,"ei"] <- MA[["b"]]
    StatSimple[k,"Qt"] <- MA[["QE"]]
    StatSimple[k,"df"] <- MA[["k"]]-1
    StatSimple[k,"pQt"] <- MA[["QEp"]]
  }
  for (k in c("T","pH","TpH")){
    MAdata <- data[Stressor==k]
    MA <- rma(abs(ei),vei,data=MAdata)
    row_values <- data.frame(Stressor=k,Scenario="all",N=MA[["k"]],taxa=taxa,
                             metric=metric,pvalue=MA[["pval"]],lowerCI=MA[["ci.lb"]],
                             upperCI=MA[["ci.ub"]],ei=MA[["b"]],Qt=MA[["QE"]],df=MA[["k"]]-1,pQt=MA[["QEp"]])
    StatSimple <- rbind(StatSimple,row_values)
  }
  return(StatSimple)
} 

MA_TpH_moderator <- function(taxa,metric,data){
  data <- data[Category==metric]
  #Sample size
  sample1 <- data[,.(.N),by=.(Stressor)]
  sample2 <- data[,.(.N),by=.(Stressor,Scenario)]
  sample2 <- sample2[!which(is.na(sample2$Scenario)),]
  StatModerator <- data.frame(Stressor=sample1$Stressor,N=sample1$N,taxa=c(taxa),metric=c(metric),
                              QM=c("NA"),QMdf=c("NA"),QMp=c("NA"),QE=c("NA"),QEdf=c("NA"),QEp=c("NA"))

  for (k in c("T","pH","TpH")){
    samplek <- sample2[sample2$N>1 & sample2$Stressor==k,]
    if (nrow(samplek)>1){
       MAdata <- data[Stressor==k]
    MAscenarioInt <- rma(ei,vei,data=MAdata,mods=~Scenario)
    MAscenario <- rma(ei,vei,data=MAdata,mods=~Scenario-1)
    StatModerator$QM[StatModerator$Stressor==k] <- MAscenarioInt[["QM"]]
    StatModerator$QMdf[StatModerator$Stressor==k] <- MAscenarioInt[["QMdf"]][1]
    StatModerator$QMp[StatModerator$Stressor==k] <- MAscenarioInt[["QMp"]]
    StatModerator$QE[StatModerator$Stressor==k] <- MAscenarioInt[["QE"]]
    StatModerator$QEdf[StatModerator$Stressor==k] <- MAscenarioInt[["k"]]-3
    StatModerator$QEp[StatModerator$Stressor==k] <- MAscenarioInt[["QEp"]]
    }
  }
  return(StatModerator)
}

MA_TpH_abs_moderator <- function(taxa,metric,data){
  data <- data[Category==metric]
  #Sample size
  sample1 <- data[,.(.N),by=.(Stressor)]
  sample2 <- data[,.(.N),by=.(Stressor,Scenario)]
  sample2 <- sample2[!which(is.na(sample2$Scenario)),]
  sample2 <- sample2[sample2$N>1,]
  StatModerator <- data.frame(Stressor=sample1$Stressor,N=sample1$N,taxa=c(taxa),metric=c(metric),
                              QM=c("NA"),QMdf=c("NA"),QMp=c("NA"),QE=c("NA"),QEdf=c("NA"),QEp=c("NA"))
  for (k in c("T","pH","TpH")){
    samplek <- sample2[sample2$Stressor==k,]
    if (nrow(samplek)>1){
      MAdata <- data[Stressor==k]
      MAscenarioInt <- rma(abs(ei),vei,data=MAdata,mods=~Scenario)
      MAscenario <- rma(abs(ei),vei,data=MAdata,mods=~Scenario-1)
      StatModerator$QM[StatModerator$Stressor==k] <- MAscenarioInt[["QM"]]
      StatModerator$QMdf[StatModerator$Stressor==k] <- MAscenarioInt[["QMdf"]][1]
      StatModerator$QMp[StatModerator$Stressor==k] <- MAscenarioInt[["QMp"]]
      StatModerator$QE[StatModerator$Stressor==k] <- MAscenarioInt[["QE"]]
      StatModerator$QEdf[StatModerator$Stressor==k] <- MAscenarioInt[["k"]]-3
      StatModerator$QEp[StatModerator$Stressor==k] <- MAscenarioInt[["QEp"]]
    }
  }
  return(StatModerator)
}