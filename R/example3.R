# Kunc??

# this removes everything from R ? gives a clean slate
rm(list=ls()) 



### load packages 

#library("xlsx")     # version 0.6.1 
library("metafor")  # version 2.0-0
library("rotl")     # version 3.0.5
library("ape")      # version 5.2
library("corpcor")  # version 1.6.9
library("ggplot2")  # version 3.0.0
library("orchaRd")
library("here")
theme_set(theme_bw())


source(here("R", "lnM_SAFE8.R"))   # defines lnM_delta1_indep()   &   safe_lnM_indep()

# lnM - SAFE
## ── 5.  SAFE ln M  (parametric bootstrap)  ─────────────────────────────────
##        B = 10 000 resamples is a good default.
get_lnM_safe <- function(m1, m2, s1, s2, n1, n2, B = 1e3) {
  out <- safe_lnM_indep(m1, m2, s1, s2, n1, n2, B = B)
  tibble(
    yi_lnM_safe = out$point,
    vi_lnM_safe = out$var,           # ← bootstrap sampling variance
    draws_kept  = out$kept,
    draws_total = out$total
  )
}

# load data

data<-read.csv(here("data", "kunc_2019.csv"), na.strings=c("","NA"))
data<-as.data.frame(data)
data.2<-data[which(data$species.latin!="NA"),]                ### exclude cases where no species were provided
str(data.2)

# GET all SPECIES 
species.overall<-unique(as.character(data.2$species.latin))
species.overall
# 109 species overall

### replace "." between genus and species names by " " ###
species.overall.2<-gsub("\\.", "\\ ", species.overall)
species.overall.2

### SEARCH for open phylogeny SPECIES ID											
taxa.overall<-tnrs_match_names(names= species.overall.2, context_name = "Animals")
taxa.overall



### For eight species, the OTT_ID (Tree of life) could not be given unambiguously.
### Therefore, these species had to be excluded to allow phylogenetic analysis.
data.3<-data.2[ which(data.2$species.latin!="Pomacentrus.amboinensis"),]		  ### exclude Pomacentrus amboinensis
data.4<-data.3[ which(data.3$species.latin!="Pseudochromis.fuscus"),]		      ### exclude Pseudochromis fuscus
data.5<-data.4[ which(data.4$species.latin!="Dascyllus.trimaculatus"),]		    ### exclude Dascyllus trimaculatus
data.6<-data.5[ which(data.5$species.latin!="Acanthochromis.polyacanthus"),]	### exclude Acanthochromis polyacanthus
data.7<-data.6[ which(data.6$species.latin!="Heterostichus.rostratus"),]	    ### exclude Heterostichus rostratus
data.8<-data.7[ which(data.7$species.latin!="Regulus.calendula"),]		        ### exclude Regulus calendula
data.9<-data.8[ which(data.8$species.latin!="Corothrella.spp."),]		          ### exclude Corothrella spp (only genus provided)
data.10<-data.9[ which(data.9$species.latin!="Chromis.chromis"),]	            ### exclude Chromis chromis


### add species ott_ids to dataframe so that models can match correlation matrix phylo IDs to species in dataframe

#Arthropoda
data.10$species_ott[data.10$species.latin=="Palaemon.serratus"]<-"Palaemon_serratus_ott727114"
data.10$species_ott[data.10$species.latin=="Pagurus.bernhardus"]<-"Pagurus_bernhardus_ott929440"
data.10$species_ott[data.10$species.latin=="Coenobita.clypeatus"]<-"Coenobita_clypeatus_ott2974729"
data.10$species_ott[data.10$species.latin=="Carcinus.maenas"]<-"Carcinus_maenas_ott1089892" 
data.10$species_ott[data.10$species.latin=="Gryllus.bimaculatus"]<-"Gryllus_bimaculatus_ott1077624"
data.10$species_ott[data.10$species.latin=="Schizocrosa.ocreata"]<-"Schizocosa_ocreata_ott914339"
data.10$species_ott[data.10$species.latin=="Jasus.edwardsii"]<-"Jasus_edwardsii_ott557395"
data.10$species_ott[data.10$species.latin=="Danaus.plexippus"]<-"Danaus_plexippus_ott190091"
data.10$species_ott[data.10$species.latin=="Teleogryllus.oceanicus"]<-"Teleogryllus_oceanicus_ott385033"
data.10$species_ott[data.10$species.latin=="Palinurus.elephas"]<-"Palinurus_elephas_ott535217"

#Mollusca
data.10$species_ott[data.10$species.latin=="Pecten.novaezelandiae"]<-"Pecten_novaezelandiae_ott1002246"
data.10$species_ott[data.10$species.latin=="Stylocheilus.striatus"]<-"Stylocheilus_striatus_ott5735083"
data.10$species_ott[data.10$species.latin=="Sepia.officinalis"]<-"Sepia_officinalis_ott528817"
data.10$species_ott[data.10$species.latin=="Sinonovacula.constricta"]<-"Sinonovacula_constricta_ott1064443"
data.10$species_ott[data.10$species.latin=="Magallana.gigas"]<-"Magallana_gigas_ott987409"

#Amphibia
data.10$species_ott[data.10$species.latin=="Dendropsophus.triangulum"]<-"Dendropsophus_triangulum_ott682688"
data.10$species_ott[data.10$species.latin=="Dendropsophus.ebraccatus"]<-"Dendropsophus_ebraccatus_ott108407"
data.10$species_ott[data.10$species.latin=="Dendropsophus.microcephalus"]<-"Dendropsophus_microcephalus_ott682674"
data.10$species_ott[data.10$species.latin=="Hyla.arborea"]<-"Hyla_arborea_ott677266"
data.10$species_ott[data.10$species.latin=="Hyla.versicolor"]<- "Dryophytes_versicolor_ott386095"
data.10$species_ott[data.10$species.latin=="Tlalocohyla.loquax"]<-"Tlalocohyla_loquax_ott632258"
data.10$species_ott[data.10$species.latin=="Tlalocohyla.picta"]<-"Tlalocohyla_picta_ott489788"
data.10$species_ott[data.10$species.latin=="Pseudacris.crucifer"]<-"Pseudacris_crucifer_ott747287" 
data.10$species_ott[data.10$species.latin=="Litoria.caerulea"]<-"Litoria_caerulea_ott386103"
data.10$species_ott[data.10$species.latin=="Agalychnis.moreletii"]<-"Agalychnis_moreletii_ott375954"
data.10$species_ott[data.10$species.latin=="Agalychnis.callidryas"]<-"Agalychnis_callidryas_ott9483"
data.10$species_ott[data.10$species.latin=="Incilius.valliceps"]<-"Incilius_valliceps_ott937819"
data.10$species_ott[data.10$species.latin=="Bufo.americanus"]<-"Anaxyrus_americanus_ott889326" 
data.10$species_ott[data.10$species.latin=="Batrachyla.leptopus"]<-"Batrachyla_leptopus_ott636141"
data.10$species_ott[data.10$species.latin=="Pelophylax.ridibundus"]<-"Pelophylax_ridibundus_ott14708"	
data.10$species_ott[data.10$species.latin=="Rana.pipiens"]<-"Rana_pipiens_ott14703"
data.10$species_ott[data.10$species.latin=="Rana.clamitans"]<-"Rana_clamitans_ott515378"
data.10$species_ott[data.10$species.latin=="Hypsiboas.bischoffi"]<-"Boana_bischoffi_ott669688"
data.10$species_ott[data.10$species.latin=="Hypsiboas.leptolineatus"]<-"Hypsiboas_leptolineatus_ott499875"
data.10$species_ott[data.10$species.latin=="Scinax.perereca"]<-"Scinax_perereca_ott3619779"
data.10$species_ott[data.10$species.latin=="Dendropsophus.minutus"]<-"Dendropsophus_minutus_ott682679"

#Fish
data.10$species_ott[data.10$species.latin=="Gadus.morhua"]<-"Gadus_morhua_ott114170"
data.10$species_ott[data.10$species.latin=="Dicentrarchus.labrax"]<-"Dicentrarchus_labrax_ott3549" 
data.10$species_ott[data.10$species.latin=="Sparus.aurata"]<-"Sparus_aurata_ott760723"
data.10$species_ott[data.10$species.latin=="Perca.fluviatilis"]<-"Perca_fluviatilis_ott186488" 
data.10$species_ott[data.10$species.latin=="Gasterosteus.aculeatus"]<-"Gasterosteus_aculeatus_ott401066"  
data.10$species_ott[data.10$species.latin=="Neolamprologus.pulcher"]<-"Neolamprologus_pulcher_ott280185" 
data.10$species_ott[data.10$species.latin=="Haplochromis.piceatus"]<-"Haplochromis_piceatus_ott135516"
data.10$species_ott[data.10$species.latin=="Amatitlania.nigrofasciata"]<-"Amatitlania_nigrofasciata_ott170115"
data.10$species_ott[data.10$species.latin=="Cyprinus.carpio"]<-"Cyprinus_carpio_ott429083" 
data.10$species_ott[data.10$species.latin=="Carassius.auratus"]<-"Carassius_auratus_ott1005907"
data.10$species_ott[data.10$species.latin=="Gobio.gobio"]<-"Gobio_gobio_ott230487"
data.10$species_ott[data.10$species.latin=="Phoxinus.phoxinus"]<-"Phoxinus_phoxinus_ott1019917"
data.10$species_ott[data.10$species.latin=="Cyprinella.venusta"]<-"Cyprinella_venusta_ott351750"
data.10$species_ott[data.10$species.latin=="Rutilus.rutilus"]<-"Rutilus_rutilus_ott421037"
data.10$species_ott[data.10$species.latin=="Danio.rerio"]<-"Danio_rerio_ott1005914"
data.10$species_ott[data.10$species.latin=="Anguilla.anguilla"]<-"Anguilla_anguilla_ott854201"
data.10$species_ott[data.10$species.latin=="Chanos.chanos"]<-"Chanos_chanos_ott575844"
data.10$species_ott[data.10$species.latin=="Halobatrachus.didactylus"]<-"Halobatrachus_didactylus_ott706472"
data.10$species_ott[data.10$species.latin=="Gobius.cruentatus"]<-"Gobius_cruentatus_ott1007315"
data.10$species_ott[data.10$species.latin=="Sciena.umbra"]<-"Sciaena_umbra_ott3634399"

#Aves
data.10$species_ott[data.10$species.latin=="Agelaius.phoeniceus"]<-"Agelaius_phoeniceus_ott226605"
data.10$species_ott[data.10$species.latin=="Baeolophus.bicolor"]<-"Baeolophus_bicolor_ott922720"
data.10$species_ott[data.10$species.latin=="Bombycilla.cedrorum"]<-"Bombycilla_cedrorum_ott606625"
data.10$species_ott[data.10$species.latin=="Carduelis.carduelis"]<-"Carduelis_carduelis_ott1083719"
data.10$species_ott[data.10$species.latin=="Carpodacus.mexicanus"]<-"Haemorhous_mexicanus_ott711865"
data.10$species_ott[data.10$species.latin=="Chloris.chloris"]<-"Chloris_chloris_ott1083715"
data.10$species_ott[data.10$species.latin=="Cyanistes.caeruleus"]<-"Cyanistes_caeruleus_ott746120"
data.10$species_ott[data.10$species.latin=="Emberiza.schoeniclus"]<-"Emberiza_schoeniclus_ott711856"
data.10$species_ott[data.10$species.latin=="Empidonax.oberholseri"]<-"Empidonax_oberholseri_ott119887"
data.10$species_ott[data.10$species.latin=="Erithacus.rubecula"]<-"Erithacus_rubecula_ott507122"
data.10$species_ott[data.10$species.latin=="Fringilla.coelebs"]<-"Fringilla_coelebs_ott828694"
data.10$species_ott[data.10$species.latin=="Parus.major"]<-"Parus_major_ott515143"
data.10$species_ott[data.10$species.latin=="Passer.domesticus"]<-"Passer_domesticus_ott745175"
data.10$species_ott[data.10$species.latin=="Pavo.cristatus"]<-"Pavo_cristatus_ott453575"
data.10$species_ott[data.10$species.latin=="Phylloscopus.collybita"]<-"Phylloscopus_collybita_ott726307"
data.10$species_ott[data.10$species.latin=="Pipilo.maculatus"]<-"Pipilo_maculatus_ott97318"
data.10$species_ott[data.10$species.latin=="Piranga.ludoviciana"]<-"Piranga_ludoviciana_ott565244"
data.10$species_ott[data.10$species.latin=="Poecile.atricapilla"]<-"Poecile_atricapillus_ott17262"
data.10$species_ott[data.10$species.latin=="Poecile.gambeli"]<-"Poecile_gambeli_ott694655"
data.10$species_ott[data.10$species.latin=="Serinus.canaria"]<-"Serinus_canaria_ott464865"
data.10$species_ott[data.10$species.latin=="Serinus.serinus"]<-"Serinus_serinus_ott1083724"
data.10$species_ott[data.10$species.latin=="Setophaga.petechia"]<-"Setophaga_petechia_ott612066"
data.10$species_ott[data.10$species.latin=="Setophaga.townsendi"]<-"Setophaga_townsendi_ott157200"
data.10$species_ott[data.10$species.latin=="Sitta.canadensis"]<-"Sitta_canadensis_ott603923"
data.10$species_ott[data.10$species.latin=="Sitta.carolinensis"]<-"Sitta_carolinensis_ott82247"
data.10$species_ott[data.10$species.latin=="Spizella.passerina"]<-"Spizella_passerina_ott989504"
data.10$species_ott[data.10$species.latin=="Streptopelia.decaocto"]<-"Streptopelia_decaocto_ott277818"
data.10$species_ott[data.10$species.latin=="Sturnus.unicolor"]<-"Sturnus_unicolor_ott366470"
data.10$species_ott[data.10$species.latin=="Tachycineta.bicolor"]<-"Tachycineta_bicolor_ott136028"
data.10$species_ott[data.10$species.latin=="Taeniopygia.guttata"]<-"Taeniopygia_guttata_ott708327"
data.10$species_ott[data.10$species.latin=="Troglodytes.troglodytes"]<-"Troglodytes_troglodytes_ott259813"
data.10$species_ott[data.10$species.latin=="Turdus.migratorius"]<-"Turdus_migratorius_ott96291"
data.10$species_ott[data.10$species.latin=="Vireo.cassinii"]<-"Vireo_cassinii_ott289074"
data.10$species_ott[data.10$species.latin=="Zosterops.lateralis"]<-"Zosterops_lateralis_ott819155"
data.10$species_ott[data.10$species.latin=="Pyrocephalus.rubinus"]<-"Pyrocephalus_rubinus_ott375103"
data.10$species_ott[data.10$species.latin=="Gallus.gallus.domesticus"]<-"Gallus_gallus_ott153563"
data.10$species_ott[data.10$species.latin=="Zonotrichia.leucophrys"]<-"Zonotrichia_leucophrys_ott265553"

#Reptilia
data.10$species_ott[data.10$species.latin=="Tiliqua.scincoides"]<-"Tiliqua_scincoides_ott1016681"

#Mammalia
data.10$species_ott[data.10$species.latin=="Antrozous.pallidus"]<-"Antrozous_pallidus_ott913941"
data.10$species_ott[data.10$species.latin=="Cynomys.ludovicianus"]<-"Cynomys_ludovicianus_ott580345"
data.10$species_ott[data.10$species.latin=="Dipodomys.stephensi"]<-"Dipodomys_stephensi_ott844928"
data.10$species_ott[data.10$species.latin=="Helogale.parvula"]<-"Helogale_parvula_ott6940"
data.10$species_ott[data.10$species.latin=="Myotis.myotis"]<-"Myotis_myotis_ott966432"
data.10$species_ott[data.10$species.latin=="Rattus.norvegicus"]<-"Rattus_norvegicus_ott271555"
data.10$species_ott[data.10$species.latin=="Megaptera.novaeangliae"]<-"Megaptera_novaeangliae_ott226198"



#############################################################################################

###########################
###						          ###
###	OVERALL ANALYSIS	  ###
###					            ###
###########################

# GET all SPECIES 
all.species<-unique(as.character(data.10$species.latin))
all.species
# 101 species in total

### replace "." between genus and species names by " " ###
all.species<-gsub("\\.", "\\ ", all.species)
all.species

### SEARCH for open phylogeny SPECIES ID											
taxa.all<-tnrs_match_names(names= all.species, context_name = "Animals")
taxa.all



# # make PHYLO TREE
# tree1.all <- tol_induced_subtree(ott_ids = ott_id(taxa.all))
# str(tree1.all)
# 
# #correct tip.labels in tree to literally match data.10$species_ott
# tree1.all$tip.label[2]		
# #is "Pecten_novaezelandiae_(species_in_Bilateria)_ott1002246"
# tree1.all$tip.label[2]<-"Pecten_novaezelandiae_ott1002246"
# #should be "Pecten_novaezelandiae_ott1002246"
# tree1.all$tip.label[2]
# #check name
# tree1.all$tip.label[100]
# #is "Gadus_morhua_(species_in_domain_Eukaryota)_ott114170"
# tree1.all$tip.label[100]<-"Gadus_morhua_ott114170" 
# #should be "Gadus_morhua_ott114170"
# tree1.all$tip.label[100]
# #check name
# 
# # compute BRANCH LENGTH and PLOT new tree
# tree2.all<-compute.brlen(tree1.all, method = "Grafen", power = 1)
# str(tree2.all)
# 
# ### create correlation matrix
# corr_matrix<-vcv(tree2.all, corr=TRUE)
# corr_matrix
# str(corr_matrix)


# ### TEST whether all tip.labels of tree1 are used in species_ott ###
# tree1.all$tip.label %in% data.10$species_ott
# ### TEST whether all species_ott are part of tip.label of tree1 ###
# data.10$species_ott %in% tree1.all$tip.label

##############################
### calculate effect sizes ###
##############################

# calculates effect sizes

#  make sd.control and sd.noise numeric

data.10$sd.control<-as.numeric(as.character(data.10$sd.control))
data.10$sd.noise<-as.numeric(as.character(data.10$sd.noise))

dat <- escalc(measure="SMDH", 
                     m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
                     m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
                     data=data.10)

# get lnM using the safe function

dat <- dat %>% mutate(
  lnM_safe = pmap_dfr(
    list(mean.control, mean.noise, sd.control, sd.noise, sample.size.control, sample.size.noise.1),
    get_lnM_safe)
) %>% unnest(lnM_safe)


#data.phylo.2<-na.omit(data.phylo)


######################################
### COVARIANCE MATRIX TO TAKE CARE ###
### 	OF NON-INDEPENDENCE		       ###
######################################
# 
# data.phylo.2$unit<-factor(1:dim(data.phylo.2)[1])
# VCV.phylo <- matrix(0,nrow = dim(data.phylo.2)[1],ncol = dim(data.phylo.2)[1])
# rownames(VCV.phylo) <- data.phylo.2$unit
# colnames(VCV.phylo) <- data.phylo.2$unit 
# 
# # find start and end coordinates for the subsets
# shared_coord <- which(data.phylo.2$independent.effect.id %in% data.phylo.2$independent.effect.id[duplicated(data.phylo.2$independent.effect.id)]==TRUE)
# # matrix of combinations of coordinates for each experiment with shared control
# combinations <- do.call("rbind", tapply(shared_coord, data.phylo.2[shared_coord,"independent.effect.id"], function(x) t(combn(x,2))))
# # calculate covariance values between  values at the positions in shared_list and place them on the matrix
# # NOTE: if control group is shared between experimental group covariance equal variance of the control group
# 
# for (i in 1:dim(combinations)[1]){
#   p1 <- combinations[i,1]
#   p2 <- combinations[i,2]
#   p1_p2_cov <-0.5*sqrt(data.phylo.2[p1,"vi"])*sqrt(data.phylo.2[p2,"vi"])
#   VCV.phylo[p1,p2] <- p1_p2_cov
#   VCV.phylo[p2,p1] <- p1_p2_cov
# }
# 
# #create variance-covariance matrix
# diag(VCV.phylo) <- data.phylo.2$vi
# is.positive.definite(VCV.phylo) #needs to be positive definite so it can be inverted in analyses (must be invertable)
# 
# filster data with no vi

dat<-dat[which(!is.na(dat$vi)),]


VCV <- vcalc(data = dat, 
                     cluster = dat$study, 
                     obs = dat$case.nr,
                     rho = 0.5, 
                     vi = dat$vi)


##################### Fit model with taxa as moderator, and study, phylogeny and case number as random variable ###########################################


###########################################

mod0 <-rma.mv	(yi=abs(yi), 
                       V=VCV, 
                       #mods= as.numeric(taxon.for.plot),
                       random = list(~1 | study, ~ 1 | species.latin, ~1 | case.nr), 
                       #R = list(species_ott = corr_matrix), 
                       #data=data.phylo.2, 
                       data = dat,
                       method="REML",
               sparse = TRUE)

summary(mod0)

orchard_plot(mod0,  xlab = "abs(SMDH)", group = "study")

# run the same model using yi_lnM_safe and vi_lnM_safe - may need to create VCV.lnM_safe first

VCV.lnM_safe <- vcalc(data = dat, 
                     cluster = dat$study, 
                     obs = dat$case.nr,
                     rho = 0.5, 
                     vi = dat$vi_lnM_safe)

mod0_lnM_safe <-rma.mv	(yi=yi_lnM_safe, 
                       V=VCV.lnM_safe, 
                       #mods= as.numeric(taxon.for.plot),
                       random = list(~1 | study, ~ 1 | species.latin, ~1 | case.nr), 
                       #R = list(species_ott = corr_matrix), 
                       #data=data.phylo.2, 
                       data = dat,
                       method="REML",
                       sparse = TRUE)

summary(mod0_lnM_safe)

# evidence aganist - they are not different or more within group variation is higher
orchard_plot(mod0_lnM_safe,  xlab = "abs(lnM SAFE)", group = "study")

# do meta-regression with taxon.for.plot

mod1 <-rma.mv	(yi=abs(yi), 
                       V=VCV, 
                       mods= ~ taxon.for.plot,
                       random = list(~1 | study, ~ 1 | species.latin, ~1 | case.nr), 
                       #R = list(species_ott = corr_matrix), 
                       #data=data.phylo.2, 
                       data = dat,
                       method="REML",
                       spares = TRUE)

summary(mod1)

orchard_plot(mod1, mod = "taxon.for.plot", xlab = "abs(SMDH)", group = "study")

# use lnM safe

mod1_lnM_safe <-rma.mv	(yi=yi_lnM_safe, 
                       V=VCV.lnM_safe, 
                       mods= ~ taxon.for.plot,
                       random = list(~1 | study, ~ 1 | species.latin, ~1 | case.nr), 
                       #R = list(species_ott = corr_matrix), 
                       #data=data.phylo.2, 
                       data = dat,
                       method="REML",
                       sparse = TRUE)

summary(mod1_lnM_safe)
i2_ml(mod1_lnM_safe)

orchard_plot(mod1_lnM_safe, mod = "taxon.for.plot", xlab = "lnM SAFE", group = "study")

#####

library(metafor)
library(multcomp)

## 1) Ensure factor and drop empty levels in the analysis data
dat$taxon.for.plot <- droplevels(factor(dat$taxon.for.plot))

## Quick sanity check
table(dat$taxon.for.plot)
stopifnot(nlevels(dat$taxon.for.plot) >= 2)

## 2) Fit with no intercept so coefs == group means
mod0 <- rma.mv(
  yi = yi_lnM_safe,
  V  = VCV.lnM_safe,
  mods   = ~ taxon.for.plot - 1,
  random = list(~ 1 | study, ~ 1 | species.latin, ~ 1 | case.nr),
  data   = dat,
  method = "REML",
  sparse = TRUE
)

summary(mod0)

## 3) Build pairwise contrast matrix from the **model's** coefficients
cn <- names(coef(mod0))          # e.g., "taxon.for.plotBird", "taxon.for.plotMammal", ...
p  <- length(cn)
if (p < 2) stop("Only one taxon level has data in the fitted model; need >= 2.")

pairs_ij <- utils::combn(seq_len(p), 2)   # indices of all pairs
K <- matrix(0, nrow = ncol(pairs_ij), ncol = p)
for (k in seq_len(ncol(pairs_ij))) {
  i <- pairs_ij[1, k]; j <- pairs_ij[2, k]
  K[k, i] <-  1
  K[k, j] <- -1
}
colnames(K) <- cn
rownames(K) <- paste0(sub("^taxon.for.plot", "", cn[pairs_ij[1, ]]),
                      " - ",
                      sub("^taxon.for.plot", "", cn[pairs_ij[2, ]]))

## 4) Simultaneous tests
cmp <- glht(mod0, linfct = K)
summary(cmp)
confint(cmp)
plot(cmp)

# build mod0 and K as you already have

cmp <- glht(mod0, linfct = K)

## Unadjusted (per-contrast) p-values:
summary(cmp, test = univariate())

## Unadjusted 95% CIs (pointwise):
confint(cmp, calpha = qnorm(0.975))

