### R script for analyses on PUBLICATION BIAS in the article 
### "Species sensitivities to a global pollutant: a meta-analysis on acoustic signals in response to anthropogenic noise"
### Hansjoerg P. Kunc and Rouven Schmidt
### contact: kunc@gmx.at
### R version 3.5.2
### RStudio Version 1.1.463


### _______________________________________________________ ###
###             PUBLICATION BIAS                            ###
###                     and                                 ###
###             TIME-LAG BIAS                               ###
### _______________________________________________________ ###

# this removes everything from R and gives a clean slate
rm(list=ls()) 



# FUNCTIONS NEEDED LATER ON TO SAVE PLOTS AND DATAFRAMES ------------------------

# FUNCTION TO RETURN P-VALUES which are below 0.001 as "<0.001"  
# pvalr : returns p-values >0.001 as numeric, and those < 0.001 as "<0.001" 
pvalr <- function(pvals, sig.limit = .001, digits = 4, html = FALSE) {
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }
  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit)
      if (html)
        return(sprintf('&lt; %s', format(sig.limit))) else
          return(sprintf('<%s', format(sig.limit)))
    if (x > .1)
      return(roundr(x, digits = 4)) else
        return(roundr(x, digits = digits))
  }, sig.limit = sig.limit)
}

# function to paste only the first two words of a string into a new column;
# used to get the appropriate species names given by the tree of life to be used in forest plots later on
string_fun <- function(x) {
  ul = unlist(strsplit(x, split = "\\s+"))[1:2]
  paste(ul,collapse=" ")
}





# load packages -----------------------------------------------------------

library("xlsx")
library("metafor")
library("rotl")
library("ape")
library("corpcor") 
library("dplyr")
library("ggplot2")
library("ggtree")


# set working directory---------------------------------------
setwd("C:\\x/y")
getwd()


# read xlsx file
data<-read.xlsx("data_file.xlsx",
                1,	as.data.frame=TRUE, header=TRUE)
length(data$case.nr)
# data file consists of 122 cases
length(unique(data$study.nr))
length(unique(data$species.latin))
# 23 studies on 31 species are included



# GET all DIFFERENT SPECIES in the data file ---------------------------------
# species.number: vector of all species names
species.number<-unique(as.character(data$species.latin))
# replace "." between genus and species names by " "
# to search for species in Open Tree of Life later on
species.number<-gsub("\\.", "\\ ", species.number)
species.number
# species.number is a vector containing 31 dfferent species
# (this corresponds with the number of species in the data file)




# GET INFORMATION ON SPECIES' PHYLOGENY USING THE OPEN TREE OF LIFE (OTL)-------------------------

# SEARCH for open phylogeny SPECIES ID (=ott_id) 
# taxa.data.file:  data frame summarizing the matched taxonomic names and species ID (ott_id)
taxa.data.file<-tnrs_match_names(names= species.number, context_name = "Animals")
taxa.data.file


### WRITE SPECIES_OTT ID TO DATA FRAME -----------------------------------------

# replace dots by spaces characters in species.latin column
# important to literally match the unique names given by tree of life
data$species.latin<-gsub("\\.", "\\ ", data$species.latin)

# add species_ott column to data file for ott IDs given by OTL
data$species_ott<- NA 
# add unique_name column for species names given by OTL to data file
data$unique.name<- NA

# write species_otts to column "species_ott"
for(i in 1:length(data$species.latin)){
  data$species_ott[i]<-paste(
    taxa.data.file$unique_name[which(taxa.data.file$search_string==tolower(data$species.latin[i]))],
    taxa.data.file$ott_id[which(taxa.data.file$search_string==tolower(data$species.latin[i]))], 
    sep="_ott")
  
  # for a potential figure without "ott id"
  data$unique.name[i]<-unlist(
    lapply(taxa.data.file$unique_name[which(taxa.data.file$search_string==tolower(data$species.latin[i]))]
           , string_fun))
}

# replace space characters by underscores in species_ott column
# important to match the trees' tip labels
data$species_ott<-gsub("\\ ", "\\_", data$species_ott)



# CALCULATE EFFECT SIZES and add to data file --------------------------------------------------
data<- escalc(measure="SMDH", 
              m2i=mean.control, sd2i=sd.control, n2i=sample.size.control, 
              m1i=mean.noise, sd1i=sd.noise,  n1i=sample.size.noise.1,
              data=data)


# check for NA in effect sizes data.file
if(any(is.na(data$yi))) {
  cat("NA found in case.nr ", data[is.na(data$yi),]$case.nr) 
} else {print("No NA found in column")}
# NA found in case.nr 403, i.e. effect size cannot not be calculated.
# The reason: sd of control and treatment group was 0
# this case is to be excluded
# data.comm: data frame of all cases on communication that will be analysed
data<-data%>% filter(yi !="NA")
length(data$yi)
# consequently, the remaining 121 cases will be analysed
length(unique(data$study.nr))
length(unique(data$species.latin))
# we have 23 studies on 31 species
# the excluded case.nr 403 was on a species (Zonotrichia leucophrys)
# for which several measures in one study were reported,
# that is why we still have the same number of studies and species.



# IDENTIFY REPONSE PARAMETER GROUPS THAT CAN BE ANALYSED  -----------------

table(factor(data$parameter.rough.standardised))
# sufficient data are available to analyse
# 1. amplitude (20 ES)
# 2. complexity (7 ES)
# 3. dominant.frequency (21 ES)
# 4. duration (28 ES)
# 5. minimum.frequency (13 ES)
# 6. rate (32 ES)


# make data frames for results of analyses on each category 
# C: Egger's regression
# D: analyses of time lag bias
# these data frames will be saved to xlsx files

# component: contains the categories to be analysed (1.-6., see main R code)
component<-c("amplitude", "complexity", "dominant.frequency", "duration", "minimum.frequency", "rate")

# define the results from analyses in R that are to be saved 
values<-c("sample.size",   # number of effect sizes
          "studies",       # number of studies
          "species",       # number of species
          "mean",          # estimated coefficient of the model
          "se",            # standard error of the coefficient
          "z",             # test statistics of the coefficient
          "lower.CI",      # lower bound of the confidence interval for the coefficient
          "upper.CI",      # upper bound of the confidence interval for the coefficient
          "p")             # p-value for the test statistics
          

# generate data frames to save analyses on subgroup means
# C: Egger's regression
# D: analyses of time lag bias

# C 
# data.Egger: data frame categories X results of Egger's regression values
data.Egger<-data.frame(matrix(nrow=length(component), ncol=length(values)))
colnames(data.Egger)<- values
rownames(data.Egger)<- component
data.Egger

# D 
# data.time.lag.bias: data frame categories X results of analyses of time lag bias
data.time.lag.bias<-data.frame(matrix(nrow=length(component), ncol=length(values)))
colnames(data.time.lag.bias)<- values
rownames(data.time.lag.bias)<- component
data.time.lag.bias



for (j in component){
  
  
  # get data subset for category
  # data.pub.bias
  data.pub.bias<-data[which(data$parameter.rough.standardised==j),]
  cat("______________________________________________\n\n")
  cat(toupper(j), "\n\n")
  cat("EGGER'S REGRESSION and TIME LAG BIAS", toupper(j), "\n")
  cat("______________________________________________\n\n")
  cat("- data set consists of", length(data.pub.bias$case.nr), "cases\n")  
# SPECIES AND PHYLOGENY
# GET all SPECIES and replace "." between genus and species names by " " 
unique.species<-unique(as.character(gsub("\\.", "\\ ", data.pub.bias$unique.name)))
cat("- we have", length(unique(as.character(data.pub.bias$study))), "different studies\n")
cat("- we have", length(unique.species), "different species\n\n")
taxa.unique.species<-tnrs_match_names(names= unique.species, context_name = "Animals")
# make PHYLOGENETIC TREE
tree1.species <- suppressWarnings(tol_induced_subtree(ott_ids = ott_id(taxa.unique.species)))
# compute BRANCH LENGTH and PLOT new tree
tree2.species<-compute.brlen(tree1.species, method = "Grafen", power = 1)


# TEST whether all tip.labels of tree1 are used in species_ott
if(any(isFALSE( tree1.species$tip.label %in% data.pub.bias$species_ott))) {
  cat("CHECK! Not all tiplabel are in column species_ott!\n")
  }
# vice versa: TEST whether all species_ott are part of tip.label of tree1 
if(any(isFALSE( data.pub.bias$species_ott %in% tree1.species$tip.label))) {
  cat("CHECK! Not all species_ott are in column tiplabel!\n")
  }
# check for NA in  column species_ott
if(any(is.na(data.pub.bias$species_ott))) {
  cat("NA found in case.nr",data[is.na(data$species_ott),]$case.nr,".\n\n") 
  } 

# create correlation matrix
corr_matrix<-vcv(tree2.species, corr=TRUE)



##############################
###   Egger's regression   ###
###   Time lag bias        ###
##############################

precision<-sqrt(1/data.pub.bias$vi)
###--------------------------------------------------------------------###
# Egger's regression
data.egger<-rma.mv(yi,vi, mod = precision, 
                        random = list(~1 | case.nr,
                                      ~1 | study,
                                      ~1 | species_ott), 
                        R = list(species_ott = corr_matrix),
                        data=data.pub.bias, 
                        control=list(optimizer="optim", optmethod="Nelder-Mead"),
                        method="REML")

cat("______________________________________________\n")
cat("\nSUMMARY OF EGGER'S REGRESSION\n\n")
cat("______________________________________________\n\n")
print(summary(data.egger))

###--------------------------------------------------------------------###
# time-lag bias model
meta.data.year<-rma.mv(yi=yi, V=vi, mod=year, 
                        random = list(~1 | study,
                                      ~1 | species_ott, 
                                      ~1 | case.nr), 
                        data=data.pub.bias, 
                        control=list(optimizer="optim", optmethod="Nelder-Mead"),
                        method="REML")

cat("______________________________________________\n")
cat("\nSUMMARY OF TIME LAG BIAS\n\n")
cat("______________________________________________\n\n")
print(summary(meta.data.year))

#extract residuals
resid.1<-rstandard.rma.mv(meta.data.year, type="rstandard")

###--------------------------------------------------------------------###
#funnel plot
gg.funnel <- ggplot(
  data.pub.bias, aes(x=precision, y=resid.1$resid)) + 
  geom_point(aes(col=study), size=3.5) +
  theme_bw()+  
  theme(legend.position="None") +
  geom_hline(yintercept = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ 
  xlim(c(0, 14)) + 
  ylim(c(-15, 15)) + 
  labs( 
    y="meta-analytic residuals", 
    x="precision" 
  )

plot(gg.funnel)
ggsave.name.funnel<-paste(Sys.Date(), toupper(j),"funnel.plot.png", sep="_")
cat(paste(ggsave.name.funnel, ": ", sep=""))
ggsave(ggsave.name.funnel, plot=gg.funnel)

# copy results of Egger regression to data frame data.Egger
# the correct row [j,] and column [,x] of data frame is chosen  
# copy the number of cases in data subset
data.Egger[j,1]<-data.egger$s.nlevels[1]
# copy the number of different studies
data.Egger[j,2]<-data.egger$s.nlevels[2]
# copy the number of different species
data.Egger[j,3]<-data.egger$s.nlevels[3]
# copy the estimated coefficient (intercept)
data.Egger[j,4]<-round(data.egger$b[1], digits=2)
# copy the standard error of the coefficient (intercept)
data.Egger[j,5]<-round(data.egger$se[1], digits=2)
# copy the test statistics of the coefficient (intercept)
data.Egger[j,6]<-round(data.egger$zval[1], digits=2)
# copy the lower bound of the confidence interval for the coefficient (intercept)
data.Egger[j,7]<-round(data.egger$ci.lb[1], digits=2)
# copy the upper bound of the confidence interval for the coefficient (intercept)
data.Egger[j,8]<-round(data.egger$ci.ub[1], digits=2)
# copy the p-value for the test statistics (intercept)
data.Egger[j,9]<-pvalr(round(data.egger$pval[1], digits=4))

# copy results of time lag model to data frame data.time.lag.bias
# the correct row [j,] and column [,x] of data frame is chosen  
# copy the number of cases in data subset
data.time.lag.bias[j,1]<-meta.data.year$s.nlevels[3]
# copy the number of different studies
data.time.lag.bias[j,2]<-meta.data.year$s.nlevels[1]
# copy the number of different species
data.time.lag.bias[j,3]<-meta.data.year$s.nlevels[2]
# copy the estimated coefficient (moderator)
data.time.lag.bias[j,4]<-round(meta.data.year$b[2], digits=2)
# copy the standard error of the coefficient (moderator)
data.time.lag.bias[j,5]<-round(meta.data.year$se[2], digits=2)
# copy the test statistics of the coefficient (moderator)
data.time.lag.bias[j,6]<-round(meta.data.year$zval[2], digits=2)
# copy the lower bound of the confidence interval for the coefficient (moderator)
data.time.lag.bias[j,7]<-round(meta.data.year$ci.lb[2], digits=2)
# copy the upper bound of the confidence interval for the coefficient (moderator)
data.time.lag.bias[j,8]<-round(meta.data.year$ci.ub[2], digits=2)
# copy the p-value for the test statistics (moderator)
data.time.lag.bias[j,9]<-pvalr(round(meta.data.year$pval[2], digits=4))

###--------------------------------------------------------------------###

### plot yi over years with vi as bubble size and bubble colour as size factor
# divide the sample size variable in a categorical variable with 4 levels
data.pub.bias$category <- cut(data.pub.bias$sample.size.control, 
                             breaks=c(-Inf, 20, 30, Inf), 
                             labels=c("1","2","3"))


### plot yi over years with vi as bubble size and bubble colour as size factor

plot.yi.vi.sample <- ggplot(data=data.pub.bias, aes(x =year , y = yi)) +
  geom_point(aes(size=vi), 
             show.legend=FALSE,
             colour="black",
             stroke=1,
             fill=factor(data.pub.bias$category),
             shape=21,
             alpha=0.2)+
  scale_size_continuous(range=c(4,16))+
  theme_bw()+
  theme(legend.position="None") +
  geom_hline(yintercept = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ 
  scale_x_continuous(breaks = seq(2000, 2022, 2))+
  labs(	y="effect size", 
        x="year")#+
  #geom_smooth(method="auto",  
  #            linetype="dashed", 
  #            span = 1, 
  #            level = 0.95,
  #            color="black", 
  #            fill="grey69")

plot(plot.yi.vi.sample)

ggsave.name.bubble<-paste(Sys.Date(), toupper(j),"bubble.plot.png", sep="_")
cat(paste(ggsave.name.bubble,": ", sep=""))
ggsave(ggsave.name.bubble, plot=plot.yi.vi.sample)

#############################################
cat("\nEnd of analyses on", toupper(j), "\n")

line <- readline(prompt="Press [enter] to continue...")
}



# WRITE RESULTS OF SUBGROUP ANALYSES TO XLSX FILES  ---------------------------

# A: RAW values
data.Egger
write.xlsx(data.Egger, file=paste(Sys.Date(), "Eggers.regression.xlsx", sep="_"), 
           col.names=TRUE, row.names=TRUE, append=FALSE, password=NULL)


# B: FOLDED values
data.time.lag.bias
write.xlsx(data.time.lag.bias, file=paste(Sys.Date(), "time.lag.bias.xlsx", sep="_"), 
           col.names=TRUE, row.names=TRUE, append=FALSE, password=NULL)


# END #


### ------------------------------------------------
# sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1
#
# Matrix products: default
#
# locale:
# [1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252  LC_MONETARY=German_Germany.1252 LC_NUMERIC=C  LC_TIME=German_Germany.1252    
#
# attached base packages:
#  [1] stats    graphics    grDevices   utils   datasets    methods   base     
#
# other attached packages:
#  [1] ggtree_1.14.4    ggplot2_3.2.1   dplyr_0.8.4   corpcor_1.6.9   ape_5.2   rotl_3.0.10   metafor_2.0-0   Matrix_1.2-14   xlsx_0.6.1   
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.3          BiocManager_1.30.10   compiler_3.5.1    pillar_1.4.3      prettyunits_1.1.1   tools_3.5.1       progress_1.2.2     
# [8] digest_0.6.25       tidytree_0.3.1        lifecycle_0.1.0   gtable_0.3.0      jsonlite_1.6.1      tibble_2.1.3      nlme_3.1-137       
# [15] lattice_0.20-35    pkgconfig_2.0.3       rlang_0.4.5       rstudioapi_0.11   rvcheck_0.1.7       curl_4.3          yaml_2.2.1         
# [22] parallel_3.5.1     treeio_1.6.1          rJava_0.9-11      withr_2.1.2       httr_1.4.1          xlsxjars_0.6.1    vctrs_0.2.3        
# [29] hms_0.5.3          grid_3.5.1            tidyselect_1.0.0  glue_1.3.1        R6_2.4.1            rentrez_1.2.2     XML_3.99-0.3       
# [36] farver_2.0.3       tidyr_1.0.2           purrr_0.3.3       magrittr_1.5      scales_1.1.0        assertthat_0.2.1  colorspace_1.4-1   
# [43] labeling_0.3       lazyeval_0.2.2        munsell_0.5.0     rncl_0.8.4        crayon_1.3.4       