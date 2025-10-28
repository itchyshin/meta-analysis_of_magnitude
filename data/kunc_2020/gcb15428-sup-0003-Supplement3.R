### R script for analyses in the article 
### "Species sensitivities to a global pollutant: a meta-analysis on acoustic signals in response to anthropogenic noise"
### Hansjoerg P. Kunc and Rouven Schmidt
### contact: kunc@gmx.at



### _______________________________________________________ ###
###                   META-ANALYIS                          ###
###                                                         ###
###   - analysis of song components                         ###
###   - generate data frames for forest plots               ###
### _______________________________________________________ ###


# this removes everything from R and gives a clean slate
rm(list=ls()) 




# load packages -----------------------------------------------------------
library("xlsx")
library("metafor")
library("rotl")
library("ape")
library("corpcor") 
library("dplyr")
# for versions of RStudio and packages see sessionInfo() at the end



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


# function for the transformation,
# we want the mean of a normal distribution with this mean and SD, folded about the origin.
# Code from in Morissey 2016, J Evol Biol 29, 1922-1931  (Box 1)
mu.fnorm<-function(mu,sigma){
  sigma*sqrt(2/pi)*exp((-1*mu^2)/(2*sigma^2))+mu*(1-2*pnorm(-1*mu/sigma,0,1))
}
# different function written by Mike Morrissey and taken from the thread: 
# https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q1/021684.html
#     mu.fnorm.2 <-  function(mu, sigma){
#     dnorm(mu, 0, sigma)*2*sigma^2 + mu*(2*pnorm(mu, 0, sigma) -1)
#     }
# both functions yield the same results
# first function "mu.fnorm" will be used throughout


# function to calculate the variance 
# from Dougherty & Guillette 2018,
# initially also provided by Morrissey and taken from the same thread mentioned above
var.fnorm <- function(mu, sigma){
  mu^2 + sigma^2 - (sigma*sqrt(2/pi)*exp((-1*mu^2)/(2*sigma^2)) + mu*(1-2*pnorm(-1*mu/sigma, 0, 1)))^2
}



# set working directory---------------------------------------
setwd("C:\\x/y")
getwd()


# read xlsx file
data<-read.xlsx("data_file.xlsx", 1, as.data.frame=TRUE, header=TRUE)
str(data)

#  [1] case.nr                      : case number
#  [2] study.nr                     : study number
#  [3] study                        : authors and year of study
#  [4] year                         : year of study
#  [5] journal                      : journal
#  [6] taxon.for.plot               : taxonomic group
#  [7] species.english              : species name (English)
#  [8] species.latin                : species name (Latin)
#  [9] parameter.rough.standardised : signal component
# [10] sample.size.control          : sample size of control group
# [11] mean.control                 : mean value of control group
# [12] sd.control                   : standard deviation of mean (control group)
# [13] sample.size.noise.1          : sample size of noise exposure group
# [14] mean.noise                   : mean value of noise exposure group
# [15] sd.noise                     : standard deviation of mean (noise exposure group)


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
# species.number is a vector containing 37 dfferent species
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


# check for NA in effect sizes od data file
if(any(is.na(data$yi))) {
  cat("NA found in case.nr ", data[is.na(data$yi),]$case.nr) 
  } else {print("No NA found in column")}
# NA found in case.nr 403, i.e. effect size cannot not be calculated.
# The reason: sd of control and treatment group was 0
# this case has to be excluded
data<-data%>% filter(yi !="NA")
length(data$yi)
# consequently, the remaining 121 cases will be analysed
length(unique(data$study.nr))
length(unique(data$species.latin))
# we have 23 studies on 31 species
# the excluded case.nr 403 was on a species (Zonotrichia leucophrys)
# for which several measures in one study were reported;
# that is why we still have the same number of studies and species.



# IDENTIFY REPONSE COMPONENTS THAT CAN BE ANALYSED  -----------------

table(factor(data$parameter.rough.standardised))
# subsequently, we will analyses these components of responses
# 1. amplitude          (20 ES)
# 2. complexity         (7 ES)
# 3. dominant.frequency (21 ES)
# 4. duration           (28 ES)
# 5. minimum.frequency  (13 ES)
# 6. rate               (32 ES)


# make data frames for results of analyses on each category 
# A: analysis of RAW values
# B: analyses of FOLDED NORMAL values
# these data frames will be saved to xlsx files

# component: contains the components to be analysed (1.-6., see above)
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
          "p",             # p-value for the test statistics
          "I^2.ES",        # heterogeneity within studies 
          "I^2.study",     # heterogeneity among studies
          "I^2.species",   # heterogeneity due to phylogenetic relatedness
          "I^2.total",     # total heterogeneity
          "Q",             # test statistic for the test of heterogeneity
          "Q.df",          # degrees of freedom of test statistic for the test of heterogeneity
          "Q.p")           # p-value for the test of heterogeneity
                  
          
# generate data frames to save analyses on subgroup means, A: on raw values, B: on folded values

# A 
# data.category.RAW: data frame categories X results of analyses on RAW values
data.component.RAW<-data.frame(matrix(nrow=length(component), ncol=length(values)))
colnames(data.component.RAW)<- values
rownames(data.component.RAW)<- component
data.component.RAW

# B 
# data.category.FOLDED: data frame categories X results of analyses on FOLDED NORMAL values
data.component.FOLDED<-data.frame(matrix(nrow=length(component), ncol=length(values)))
colnames(data.component.FOLDED)<- values
rownames(data.component.FOLDED)<- component
data.component.FOLDED


# MODELS ON SINGLE PARAMETERS ------------------------------------

# CALCULATE META-ANALYSIS FOR EACH COMPONENT OF RESPONSE VARIABLES
# using a for-loop
# A: the DIRECTION of responses
#    i.e. analysis on raw values
# B: the MAGNITUDE of responses
#    i.e. analysis on folded normal values


for (j in component){

    # get data subset for each component
    # data.each.comp: data of each component
    # as subsetted by for-loop
    data.each.comp<-data[which(data$parameter.rough.standardised==j),]
    cat("_________________________________________\n\n")
    cat("ANALYSIS ON", toupper(j), "\n")
    cat("_________________________________________\n\n")
    cat("- data set consists of", length(data.each.comp$case.nr), "cases\n")
  
  
    # SPECIES AND PHYLOGENY
    # GET all SPECIES and replace "." between genus and species names by " " 
    unique.species<-unique(as.character(gsub("\\.", "\\ ", data.each.comp$unique.name)))
    cat("- we have", length(unique.species), "different species\n")
    cat("- we have", length(unique(as.character(data.each.comp$study))), "different studies\n\n")
    cat("CALCULATING TREE and saving it to", paste(Sys.Date(), toupper(j), "species_tree.tre", sep="_"), "\n")
    taxa.unique.species<-tnrs_match_names(names= unique.species, context_name = "Animals")
    # make PHYLOGENETIC TREE
    tree1.species <- suppressWarnings(tol_induced_subtree(ott_ids = ott_id(taxa.unique.species)))
    # compute BRANCH LENGTH and PLOT new tree
    tree2.species<-compute.brlen(tree1.species, method = "Grafen", power = 1)
    plot(tree2.species)
    assign(paste("tree2.", j, sep=""),tree2.species, envir = .GlobalEnv)
    write.tree(tree2.species, 
             file = paste(Sys.Date(), j, "species_tree.tre", sep="_"),
             append = FALSE, digits = 10, tree.names = FALSE)
  
  
    # TEST whether all tip.labels of tree1 are used in species_ott
    if(any(isFALSE( tree1.species$tip.label %in% data.each.comp$species_ott))) {
      cat("\nCHECK! Not all tiplabel are in column species_ott!\n")
      } else {cat("\n- tree tip.labels checked.\n")}
    # vice versa: TEST whether all species_ott are part of tip.label of tree1 
    if(any(isFALSE( data.each.comp$species_ott %in% tree1.species$tip.label))) {
      cat("CHECK! Not all species_ott are in column tiplabel!\n")
      } else {cat("- species_ott checked.\n")}
    # check for NA in  column species_ott
    if(any(is.na(data$species_ott))) {
      cat("NA found in case.nr",data.each.comp[is.na(data.each.comp$species_ott),]$case.nr,".\n\n") 
      } else {cat("- no NA in species_ott.\n\n")}
  
    # create correlation matrix
    corr_matrix<-vcv(tree2.species, corr=TRUE)
  
    # ______________________________________________________________________ 
    # META-ANALYSES 
    # ______________________________________________________________________ 
    
    # A: DIRECTION of response 
    #    i.e. on RAW VALUES
    # ______________________________________________________________________ 
    
    meta.component.RAW<-rma.mv	(yi,
                                vi,
                                random = list(~1 | case.nr, 
                                              ~1 | study, 
                                              ~1 | species_ott), 
                                R = list(species_ott = corr_matrix),
                                control=list(optimizer="optim", optmethod="Nelder-Mead"),
                                data=data.each.comp, 
                                method="REML")
    
    # copy meta-analytic results to workspace
    # make the output of the model available directly in R outside the for-loop 
    assign(paste("meta.", j,".RAW", sep=""),meta.component.RAW, envir = .GlobalEnv)
    
    # calculate heterogeneity of RAW model
    # cf.  http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
    W.RAW <- diag(1/data.each.comp$vi)
    X.RAW <- model.matrix(meta.component.RAW)
    P.RAW <- W.RAW - W.RAW %*% X.RAW %*% solve(t(X.RAW) %*% W.RAW %*% X.RAW) %*% t(X.RAW) %*% W.RAW
    # I.square.RAW.1: total heterogeneity [%]
    I.square.RAW.1<-100 * sum(meta.component.RAW$sigma2) / (sum(meta.component.RAW$sigma2) + 
                          (meta.component.RAW$k-meta.component.RAW$p)/sum(diag(P.RAW)))
    # I.square.RAW.2: partitions heterogeneity into the contribution of random factors
    # as given by the order of factors in the model above
    # factor 1: within studies (random = list(~1 | CASE.NR, ~ 1 | study, ~1 | species_ott))
    # factor 2: among studies (random = list(~1 | case.nr, ~ 1 | STUDY, ~1 | species_ott))
    # factor 3: phylogeny (random = list(~1 | case.nr, ~ 1 | study, ~1 | SPECIES_OTT))
    I.square.RAW.2<-100 * meta.component.RAW$sigma2 / (sum(meta.component.RAW$sigma2) + 
                          (meta.component.RAW$k-meta.component.RAW$p)/sum(diag(P.RAW)))
    
    
    
    # copy results of RAW model to data frame data.component.RAW
    # the correct row [j,] and column [,x] of data frame is chosen  
    # copy the number of effect sizes
    data.component.RAW[j,1]<-length(data.each.comp$species.latin)
    # copy the number of different studies
    data.component.RAW[j,2]<-length(unique(data.each.comp$study))
    # copy the number of different species
    data.component.RAW[j,3]<-length(unique(data.each.comp$species.latin))
    # copy the estimated coefficient 
    data.component.RAW[j,4]<-round(meta.component.RAW$b[1], digits=2)
    # copy the standard error of the coefficient
    data.component.RAW[j,5]<-round(meta.component.RAW$se[1], digits=2)
    # copy the test statistics of the coefficient
    data.component.RAW[j,6]<-round(meta.component.RAW$zval[1], digits=2)
    # copy the lower bound of the confidence interval for the coefficient
    data.component.RAW[j,7]<-round(meta.component.RAW$ci.lb[1], digits=2)
    # copy the upper bound of the confidence interval for the coefficient
    data.component.RAW[j,8]<-round(meta.component.RAW$ci.ub[1], digits=2)
    # copy the p-value for the test statistics
    data.component.RAW[j,9]<-pvalr(round(meta.component.RAW$pval[1], digits=4))
    # copy the heterogeneity within studies 
    data.component.RAW[j,10]<-round(I.square.RAW.2[1], digits=2)
    # copy the heterogeneity among studies
    data.component.RAW[j,11]<-round(I.square.RAW.2[2], digits=2)
    # copy heterogeneity due to phylogenetic relatedness
    data.component.RAW[j,12]<-round(I.square.RAW.2[3], digits=2)
    # copy the total heterogeneity
    data.component.RAW[j,13]<-round(I.square.RAW.1, digits=2)
    # copy the test statistic for the test of heterogeneity
    data.component.RAW[j,14]<-round(meta.component.RAW$QE, digits=1)
    # copy the degrees of freedom of test statistic for the test of heterogeneity
    data.component.RAW[j,15]<-length(data.each.comp$yi)-1
    # copy the p-value for the test of heterogeneity
    data.component.RAW[j,16]<-pvalr(round(meta.component.RAW$QEp, digits=4))
    
    
    # B: MAGNITUDE of response 
    #    i.e. on FOLDED NORMAL VALUES
    # ______________________________________________________________________ 
    
    meta.component.FOLDED<-rma.mv	(mu.fnorm(yi,sqrt(vi)),
                                    var.fnorm(yi,sqrt(vi)),
                                    random = list(~1 | case.nr, 
                                                  ~1 | study, 
                                                  ~1 | species_ott), 
                                    R = list(species_ott = corr_matrix), 
                                    control=list(optimizer="optim", optmethod="Nelder-Mead"),
                                    data=data.each.comp, 
                                    method="REML")
    
    # copy meta-analytic results to workspace
    # make the output of the model available directly in R outside the for-loop
    assign(paste("meta.", j,".FOLDED", sep=""),meta.component.FOLDED, envir = .GlobalEnv)
    
    
    # calculate heterogeneity of FOLDED model
    # cf.  http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
    W.FOLDED <- diag(1/(var.fnorm(data.each.comp$yi,sqrt(data.each.comp$vi))))
    X.FOLDED <- model.matrix(meta.component.FOLDED)
    P.FOLDED <- W.FOLDED - W.FOLDED %*% X.FOLDED %*% solve(t(X.FOLDED) %*% W.FOLDED %*% X.FOLDED) %*% t(X.FOLDED) %*% W.FOLDED
    # I.square.RAW.1: total heterogeneity [%]
    I.square.FOLDED.1<-100 * sum(meta.component.FOLDED$sigma2) / (sum(meta.component.FOLDED$sigma2) + 
                              (meta.component.FOLDED$k-meta.component.FOLDED$p)/sum(diag(P.FOLDED)))
    # I.square.RAW.2: partitions of heterogeneity into the contribution of random factors
    # as given by the order of factors in the model above
    # factor 1: within studies (1st factor in: random = list(~1 | case.nr, ~ 1 | study, ~1 | species_ott))
    # factor 2: among studies (2nd factor in:random = list(~1 | case.nr, ~ 1 | study, ~1 | species_ott))
    # factor 3: phylogeny (3rd factor in: random = list(~1 | case.nr, ~ 1 | study, ~1 | species_ott))
    I.square.FOLDED.2<-100 * meta.component.FOLDED$sigma2 / (sum(meta.component.FOLDED$sigma2) + 
                              (meta.component.FOLDED$k-meta.component.FOLDED$p)/sum(diag(P.FOLDED)))
    
    
    # copy results of FOLDED NORMAL  model to data frame data.category.FOLDED
    # the correct row [j,] and column [,x] of data frame is chosen  
    # see description for each copied value in the raw model above
    data.component.FOLDED[j,1]<-length(data.each.comp$species.latin)       
    data.component.FOLDED[j,2]<-length(unique(data.each.comp$study))
    data.component.FOLDED[j,3]<-length(unique(data.each.comp$species.latin))
    data.component.FOLDED[j,4]<-round(meta.component.FOLDED$b[1], digits=2)
    data.component.FOLDED[j,5]<-round(meta.component.FOLDED$se[1], digits=2)
    data.component.FOLDED[j,6]<-round(meta.component.FOLDED$zval[1], digits=2)
    data.component.FOLDED[j,7]<-round(meta.component.FOLDED$ci.lb[1], digits=2)
    data.component.FOLDED[j,8]<-round(meta.component.FOLDED$ci.ub[1], digits=2)
    data.component.FOLDED[j,9]<-pvalr(round(meta.component.FOLDED$pval[1], digits=4))
    data.component.FOLDED[j,10]<-round(I.square.FOLDED.2[1], digits=2)
    data.component.FOLDED[j,11]<-round(I.square.FOLDED.2[2], digits=2)
    data.component.FOLDED[j,12]<-round(I.square.FOLDED.2[3], digits=2)
    data.component.FOLDED[j,13]<-round(I.square.FOLDED.1, digits=2)
    data.component.FOLDED[j,14]<-round(meta.component.FOLDED$QE, digits=1)
    data.component.FOLDED[j,15]<-length(data.each.comp$yi)-1
    data.component.FOLDED[j,16]<-pvalr(round(meta.component.FOLDED$QEp, digits=4))
    
    
    
    
    
    # DATA SHEETS FOR FOREST PLOTS
    # generate a separate xls file for each component 
    # ______________________________________________________________________ 
  
    # get list of SPECIES NAMES
    species.data<-data.each.comp[,c("taxon.for.plot", "unique.name")]
    species.data.unique<-unique(as.character(species.data$unique.name))
    species.data.unique.2<-species.data[!duplicated(species.data), ]
    colnames(species.data.unique.2)<-c("taxon", "species.data.unique") 
  
    # combine study and journal in a new variable 						
    # to account for several different studies by the same authors in the same year
    data.each.comp$study.journal<-paste(data.each.comp$study, data.each.comp$journal, sep="_")
  
    # make data frame in which to paste estimate and CI of rma.uni or rma.mv models	#
    
    # A: data frame for analysis on RAW values
    data.raw<-data.frame(species.data.unique.2)
    data.raw$taxon<-as.character(data.raw$taxon)
  
    cat("CALCULATING DATA for FOREST PLOT and saving it to",paste(Sys.Date(), j, "RAW.xlsx", sep="_"), "\n\n")
    
    for (k in 1:length(species.data.unique)){
      data.raw$studies[k]<-length(unique(subset(data.each.comp$study.journal, data.each.comp$unique.name==species.data.unique[k])))
      data.raw$ES[k]<-length(subset(data.each.comp$yi, data.each.comp$unique.name==species.data.unique[k]))
      if(data.raw$ES[k]=="1") {
        # meta1 without any random factor
        meta1<-rma.uni(yi, vi , 
                       data=data.each.comp[ which(data.each.comp$unique.name==species.data.unique[k]),], 
                       control=list(optimizer="optim", optmethod="Nelder-Mead"), 
                       method="REML")
        data.raw$estimate[k]	<-format(round(meta1$b[1], digits=2), nsmall=2)
        data.raw$se[k]		<-format(round(meta1$se, digits=2), nsmall=2)
        data.raw$z[k]			<-format(round(meta1$z, digits=2), nsmall=2)	
        data.raw$CI.lower[k]	<-format(round(meta1$ci.lb, digits=2), nsmall=2)
        data.raw$CI.upper[k]	<-format(round(meta1$ci.ub, digits=2), nsmall=2)
        data.raw$p[k]			<-pvalr(meta1$pval, digits=4)
      }				
      else
        if(data.raw$studies[k]=="1")	{
          # meta2 with case.nr as random factor
          meta2<-suppressWarnings(rma.uni(yi, vi , 
                                          random = ~1 | case.nr, 
                                          data=data.each.comp[ which(data.each.comp$unique.name==species.data.unique[k]),], 
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"), 
                                          method="REML"))
          data.raw$estimate[k]	<-format(round(meta2$b[1], digits=2), nsmall=2)
          data.raw$se[k]		<-format(round(meta2$se, digits=2), nsmall=2)	
          data.raw$z[k]			<-format(round(meta2$z, digits=2), nsmall=2)	
          data.raw$CI.lower[k]	<-format(round(meta2$ci.lb, digits=2), nsmall=2)
          data.raw$CI.upper[k]	<-format(round(meta2$ci.ub, digits=2), nsmall=2)  
          data.raw$p[k]			<-pvalr(meta2$pval, digits=4) 
        }
      else {
        # meta3 with case.nr and study as random factors
        meta3<-rma.mv(y=yi, vi , 
                      random = list(~1 | case.nr, ~1 | study.journal), 
                      data=data.each.comp[which(data.each.comp$unique.name==species.data.unique[k]),], 
                      control=list(optimizer="optim", optmethod="Nelder-Mead"), 
                      method="REML")
        data.raw$estimate[k]	<-format(round(meta3$b[1], digits=2), nsmall=2)
        data.raw$se[k]		<-format(round(meta3$se, digits=2), nsmall=2)
        data.raw$z[k]			<-format(round(meta3$z, digits=2), nsmall=2)
        data.raw$CI.lower[k]	<-format(round(meta3$ci.lb, digits=2), nsmall=2)
        data.raw$CI.upper[k]	<-format(round(meta3$ci.ub, digits=2), nsmall=2)
        data.raw$p[k]			<-pvalr(meta3$pval, digits=4)
      }
      
    }
  
    # check for any mistakes
    if(sum(data.raw$ES)==length(data.each.comp$case.nr)) {
      cat("- sample size checked.\n\n") 
      } else {print("- check if-functions")}
    # add a line to the end of data.raw
    # and paste the relevant information on the "overall model" on direction of responses
    # as has been calculated in meta.component:RAW
    data.raw[length(data.raw$species.data.unique)+1,]<-c("all taxa",
                                                         "direction of response", 					# column1:  "overall raw model"
                                                         length(unique(data.each.comp$study.journal)), 	# column2:  number of studies
                                                         length(data.each.comp$yi), 				# column3:  number of ES
                                                         round(meta.component.RAW$b[1], digits=2),	# column4:  estimate
                                                         round(meta.component.RAW$se, digits=2), 	# column5:  se 
                                                         round(meta.component.RAW$zval, digits=2), 	# column6:  z-value 
                                                         round(meta.component.RAW$ci.lb, digits=2), 	# column7:  lower CI 
                                                         round(meta.component.RAW$ci.ub, digits=2), 	# column8:  upper CI 
                                                         pvalr(meta.component.RAW$pval))			# column9:  p-value 
   
    
    data.raw$species.data.unique<-as.factor(data.raw$species.data.unique)		# convert first variable to factor
    data.raw$taxon<-as.factor(data.raw$taxon)						# convert first variable to factor
    
    # convert estimate and CI columns to class "numeric", necessary for forestplot
    data.raw$estimate<-as.numeric(data.raw$estimate)
    data.raw$CI.lower<-as.numeric(data.raw$CI.lower)
    data.raw$CI.upper<-as.numeric(data.raw$CI.upper)
    
    # convert other columns to class "numeric", easier to write to xlsx
    data.raw$studies<-as.numeric(data.raw$studies)
    data.raw$ES<-as.numeric(data.raw$ES)
    data.raw$se<-as.numeric(data.raw$se)
    data.raw$z<-as.numeric(data.raw$z)
    
    # save data for forest plots to a xlsx file
    # file name: actual date_component_RAW.XLSX 
    write.xlsx(data.raw, file=paste(Sys.Date(), j, "RAW.xlsx", sep="_"), 
               col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE, password=NULL)
    
    # make the output of the model available directly in R outside the for-loop
    assign(paste("data",j,"raw", sep="."), data.raw, envir=.GlobalEnv)
    
    cat("- model results on direction of responses can be accessed with", paste("'meta.", j,".RAW'\n", sep=""))
    cat("- model results on magnitude of responses can be accessed with",paste("'meta.", j,".FOLDED'\n\n", sep=""))
    
    
    cat("End of analyses on", j, "\n\n")
    line <- readline(prompt="Press [enter] to continue...")
}
    


    # WRITE RESULTS OF SUBGROUP ANALYSES TO XLSX FILES  ---------------------------

    # A: RAW values
    data.component.RAW
    write.xlsx(data.component.RAW, file=paste(Sys.Date(), "component.analyses.RAW.xlsx", sep="_"), 
                col.names=TRUE, row.names=TRUE, append=FALSE)
   
  
    # B: FOLDED values
    data.component.FOLDED
    write.xlsx(data.component.FOLDED, file=paste(Sys.Date(), "component.analyses.FOLDED.xlsx", sep="_"), 
               col.names=TRUE, row.names=TRUE, append=FALSE)
    
    # list all meta-analyses variables in work space (to check)
    # variables of for-loop
    for (j in component) {
                  print(ls(pattern=j))
                
    }
   # 
   # for each component XXX we have 4 files:
   # 1.) data.XXX.raw: within component XXX, model results for each species separately
   #                   used subsequently to generate forest plot
   # 2.) meta.XXX.RAW: results of rma.mv() model on parameter XXX (direction of response)
   # 3.) meta.XXX.FOLDED: results of rma.mv() model on parameter XXX (magnitude of response)
   # 4.) tree2.XXX: species tree for component XXX
   #
   # additionally, we have 2 files summarizing the effects of noise of all components
   # - component.analyses.FOLDED: Effect of anthropogenic noise on each of 6 components.
   #                   MAGNITUDE of signal adjustment. 
   #                   - Estimates and 95% confidence intervals (CI)
   #                   - displayed in table 1a of main paper
   # - component.analyses.RAW: Effect of anthropogenic noise on each of 6 components. 
   #                   DIRECTION of signal adjustment. 
   #                   - Estimates and 95% confidence intervals (CI)
   #                   - displayed in table 1b of main paper

    
    
  ### END ###
 
 
    
### ------------------------------------------------
# sessionInfo()
#
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
#  [1] stats  graphics  grDevices  utils  datasets  methods  base     
# 
# other attached packages:
#  [1] corpcor_1.6.9  ape_5.2  metafor_2.0-0  Matrix_1.2-14  rotl_3.0.10  dplyr_0.8.4  xlsx_0.6.1   
#
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.3        pillar_1.4.3      compiler_3.5.1    prettyunits_1.1.1 tools_3.5.1       progress_1.2.2    jsonlite_1.6.1    tibble_2.1.3     
# [9] nlme_3.1-137      lattice_0.20-35   pkgconfig_2.0.3   rlang_0.4.5       rstudioapi_0.11   curl_4.3          yaml_2.2.1        parallel_3.5.1   
# [17] rJava_0.9-11     httr_1.4.1        xlsxjars_0.6.1    vctrs_0.2.3       hms_0.5.3         grid_3.5.1        tidyselect_1.0.0  glue_1.3.1       
# [25] R6_2.4.1         rentrez_1.2.2     XML_3.99-0.3      purrr_0.3.3       magrittr_1.5      assertthat_0.2.1  rncl_0.8.4        crayon_1.3.4     
