### R script for analyses in the article "The global effects of anthropogenic noise: a meta-analysis"
### Hansjoerg P. Kunc and Rouven Schmidt
### contact: kunc@gmx.at
### R version 3.5.2
### RStudio Version 1.1.463


# this removes everything from R – gives a clean slate
	rm(list=ls()) 



### load packages 
  
library("xlsx")     # version 0.6.1 
library("metafor")  # version 2.0-0
library("rotl")     # version 3.0.5
library("ape")      # version 5.2
library("corpcor")  # version 1.6.9
library("ggplot2")  # version 3.0.0

	
# sets your working direcotry 
setwd("C:\\x/y")



	### read WHOLE sheet 1 from xlsx file 	   ###
	data<-read.xlsx("data_file.xlsx", 1,	as.data.frame=TRUE, header=TRUE, na.rm=TRUE)
	data<-as.data.frame(data)
	data.2<-data[which(data$species.latin!="NA"),]                ### exclude cases where no species were provided
	str(data.2)

	#[1]  "case.nr"			         # case number
	#[2]  "independent.effect.id # for matrix
	#[3]  "study"			           # authors and year of study
	#[4]	"year"			           # year of study
	#[5]  "journal"			         # journal
	#[6]  "taxon.for.plot"	     # taxonomic group
	#[7]  "species.english"		   # species name (English)
	#[8]  "species.latin"		     # species name (Latin)
	#[9]  "sample.size.control"	 # sample size of control group
	#[10] "mean.control"		     # mean value of control group
	#[11] "sd.control"		       # standard deviation of mean (control group)
	#[12] "sample.size.noise.1"	 # sample size of noise exposure group
	#[13] "mean.noise"		       # mean value of noise exposure group
	#[14] "sd.noise"			       # standard deviation of the mean (noise exposure group)
	
		
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



	# make PHYLO TREE
 	tree1.all <- tol_induced_subtree(ott_ids = ott_id(taxa.all))
	str(tree1.all)
	
	#correct tip.labels in tree to literally match data.10$species_ott
	tree1.all$tip.label[2]		
		#is "Pecten_novaezelandiae_(species_in_Bilateria)_ott1002246"
	tree1.all$tip.label[2]<-"Pecten_novaezelandiae_ott1002246"
		#should be "Pecten_novaezelandiae_ott1002246"
	tree1.all$tip.label[2]
		#check name
	tree1.all$tip.label[100]
		#is "Gadus_morhua_(species_in_domain_Eukaryota)_ott114170"
	tree1.all$tip.label[100]<-"Gadus_morhua_ott114170" 
		#should be "Gadus_morhua_ott114170"
	tree1.all$tip.label[100]
		#check name

	# compute BRANCH LENGTH and PLOT new tree
	tree2.all<-compute.brlen(tree1.all, method = "Grafen", power = 1)
	str(tree2.all)
	
	### create correlation matrix
	corr_matrix<-vcv(tree2.all, corr=TRUE)
	corr_matrix
	str(corr_matrix)

	
	### TEST whether all tip.labels of tree1 are used in species_ott ###
	tree1.all$tip.label %in% data.10$species_ott
	### TEST whether all species_ott are part of tip.label of tree1 ###
	data.10$species_ott %in% tree1.all$tip.label

	##############################
	### calculate effect sizes ###
	##############################

	# calculates effect sizes
	data.phylo <- escalc(measure="SMDH", 
				m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
				m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
				data=data.10)

	data.phylo.2<-na.omit(data.phylo)


	######################################
	### COVARIANCE MATRIX TO TAKE CARE ###
	### 	OF NON-INDEPENDENCE		       ###
	######################################

	data.phylo.2$unit<-factor(1:dim(data.phylo.2)[1])
	VCV.phylo <- matrix(0,nrow = dim(data.phylo.2)[1],ncol = dim(data.phylo.2)[1])
	rownames(VCV.phylo) <- data.phylo.2$unit
	colnames(VCV.phylo) <- data.phylo.2$unit 

	# find start and end coordinates for the subsets
	shared_coord <- which(data.phylo.2$independent.effect.id %in% data.phylo.2$independent.effect.id[duplicated(data.phylo.2$independent.effect.id)]==TRUE)
	# matrix of combinations of coordinates for each experiment with shared control
	combinations <- do.call("rbind", tapply(shared_coord, data.phylo.2[shared_coord,"independent.effect.id"], function(x) t(combn(x,2))))
	# calculate covariance values between  values at the positions in shared_list and place them on the matrix
	# NOTE: if control group is shared between experimental group covariance equal variance of the control group

	for (i in 1:dim(combinations)[1]){
  		p1 <- combinations[i,1]
  		p2 <- combinations[i,2]
  		p1_p2_cov <-0.5*sqrt(data.phylo.2[p1,"vi"])*sqrt(data.phylo.2[p2,"vi"])
  		VCV.phylo[p1,p2] <- p1_p2_cov
  		VCV.phylo[p2,p1] <- p1_p2_cov
	}

	#create variance-covariance matrix
	diag(VCV.phylo) <- data.phylo.2$vi
	is.positive.definite(VCV.phylo) #needs to be positive definite so it can be inverted in analyses (must be invertable)



	
	
 
	##################### Fit model with taxa as moderator, and study, phylogeny and case number as random variable ###########################################
	
	
	###########################################
	
	meta.phylo.mm<-rma.mv	(yi=abs(yi), 
					V=VCV.phylo, 
					mods= as.numeric(taxon.for.plot),
					random = list(~1 | study, ~ 1 | species_ott, ~1 | case.nr), 
					R = list(species_ott = corr_matrix), 
					data=data.phylo.2, 
					method="REML")

	meta.phylo.mm

	
	

		### HETEROGENEITY cf.  http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate ###
		
		W <- diag(1/data.phylo.2$vi)
		X <- model.matrix(meta.phylo.mm)
		P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
		I.square.all.1<-100 * sum(meta.phylo.mm$sigma2) / (sum(meta.phylo.mm$sigma2) + (meta.phylo.mm$k-meta.phylo.mm$p)/sum(diag(P)))
		round(I.square.all.1, digits=2)
		I.square.all.2<-100 * meta.phylo.mm$sigma2 / (sum(meta.phylo.mm$sigma2) + (meta.phylo.mm$k-meta.phylo.mm$p)/sum(diag(P)))
		round(I.square.all.2, digits=2)


############################################################################################################
############################################################################################################

#########################################
###							                      ###
### MODELS ON SINGLE TAXONOMIC GROUPS ###
###							                      ###
#########################################

	data.11<-data[ which(data$species.latin!="NA"),]
		data.responses <- escalc(measure="SMDH", 
		                          m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
		                          m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
		                          data=data.11)
		
		data.responses.2<-na.omit(data.responses)


###################
###	ARTHROPODA	###
###################

	test.arthropoda<-data.11[ which(data.11$taxon.for.plot=='arthropoda'),]
	head(test.arthropoda)

	##############################
	### calculate effect sizes ###
	##############################

	data.arthropoda <- escalc(measure="SMDH", 
				m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
				m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
				data=test.arthropoda)

	data.arthropoda.2<-na.omit(data.arthropoda)

	#######################################
	### COVARIANCE MATRIX TO TAKE CARE	###
	### 	OF NON-INDEPENDENCE		        ###
	#######################################

	data.arthropoda.2$unit<-factor(1:dim(data.arthropoda.2)[1])
	VCV.arthropoda <- matrix(0,nrow = dim(data.arthropoda.2)[1],ncol = dim(data.arthropoda.2)[1])
	rownames(VCV.arthropoda) <- data.arthropoda.2$species.latin
	colnames(VCV.arthropoda) <- data.arthropoda.2$unit 

	# find start and end coordinates for the subsets
	shared_coord <- which(data.arthropoda.2$independent.effect.id %in% data.arthropoda.2$independent.effect.id[duplicated(data.arthropoda.2$independent.effect.id)]==TRUE)
	# matrix of combinations of coordinates for each experiment with shared control
	combinations <- do.call("rbind", tapply(shared_coord, data.arthropoda.2[shared_coord,"independent.effect.id"], function(x) t(combn(x,2))))
	# calculate covariance values between  values at the positions in shared_list and place them on the matrix
	# NOTE: if control group is shared between experimental group covariance equal variance of the control group

	for (i in 1:dim(combinations)[1]){
  		p1 <- combinations[i,1]
  		p2 <- combinations[i,2]
  		p1_p2_cov <-0.5*sqrt(data.arthropoda.2[p1,"vi"])*sqrt(data.arthropoda.2[p2,"vi"])
  		VCV.arthropoda[p1,p2] <- p1_p2_cov
  		VCV.arthropoda[p2,p1] <- p1_p2_cov
	}

	#create variance-covariance matrix
	diag(VCV.arthropoda) <- data.arthropoda.2$vi
	is.positive.definite(VCV.arthropoda) #needs to be positive definite so it can be inverted in analyses (must be invertable)

	
	
	###################### FIT MODEL with study, species and case number as random factors ######################

	meta.arthropoda.mm<-rma.mv(yi = abs(yi), V = VCV.arthropoda, random = list(~1 | study, ~1 | case.nr, ~1 | species.latin), 
	     method = "REML", data = data.arthropoda.2)
	meta.arthropoda.mm

	### HETEROGENEITY cf.  http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate ###
		
		W <- diag(1/data.arthropoda.2$vi)
		X <- model.matrix(meta.arthropoda.mm)
		P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
		I.square.arthropoda.1<-100 * sum(meta.arthropoda.mm$sigma2) / (sum(meta.arthropoda.mm$sigma2) + (meta.arthropoda.mm$k-meta.arthropoda.mm$p)/sum(diag(P)))
		round(I.square.arthropoda.1, digits=2)
		I.square.arthropoda.2<-100 * meta.arthropoda.mm$sigma2 / (sum(meta.arthropoda.mm$sigma2) + (meta.arthropoda.mm$k-meta.arthropoda.mm$p)/sum(diag(P)))
		round(I.square.arthropoda.2, digits=2)


###############################################################################################################
###############################################################################################################
	
	

#################
###	MOLLUSCA  ###
#################

	test.mollusca<-data.11[ which(data.11$taxon.for.plot=='mollusca'),]
	head(test.mollusca)

	##############################
	### calculate effect sizes ###
	##############################

	data.mollusca <- escalc(measure="SMDH", 
				m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
				m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
				data=test.mollusca)
	str(data.mollusca)
	data.mollusca.2<-na.omit(data.mollusca)
	str(data.mollusca.2)

	
	#######################################
	### COVARIANCE MATRIX TO TAKE CARE	###
	### 	OF NON-INDEPENDENCE		    ###
	#######################################

	data.mollusca.2$unit<-factor(1:dim(data.mollusca.2)[1])
	VCV.mollusca<- matrix(0,nrow = dim(data.mollusca.2)[1],ncol = dim(data.mollusca.2)[1])
	rownames(VCV.mollusca) <- data.mollusca.2$unit
	colnames(VCV.mollusca) <- data.mollusca.2$unit 

	# find start and end coordinates for the subsets
	shared_coord <- which(data.mollusca.2$independent.effect.id %in% data.mollusca.2$independent.effect.id[duplicated(data.mollusca.2$independent.effect.id)]==TRUE)
	# matrix of combinations of coordinates for each experiment with shared control
	combinations <- do.call("rbind", tapply(shared_coord, data.mollusca.2[shared_coord,"independent.effect.id"], function(x) t(combn(x,2))))
	# calculate covariance values between  values at the positions in shared_list and place them on the matrix
	# NOTE: if control group is shared between experimental group covariance equal variance of the control group

	for (i in 1:dim(combinations)[1]){
  		p1 <- combinations[i,1]
  		p2 <- combinations[i,2]
  		p1_p2_cov <-0.5*sqrt(data.mollusca.2[p1,"vi"])*sqrt(data.mollusca.2[p2,"vi"])
  		VCV.mollusca[p1,p2] <- p1_p2_cov
  		VCV.mollusca[p2,p1] <- p1_p2_cov
	}

	#create variance-covariance matrix
	diag(VCV.mollusca) <- data.mollusca.2$vi
	is.positive.definite(VCV.mollusca) #needs to be positive definite so it can be inverted in analyses (must be invertable)

	
	###################### FIT MODEL with study, species and case number as random factors ######################

	meta.mollusca.mm<-rma.mv(yi = abs(yi), V = VCV.mollusca, random = list(~1 | study, ~1 | case.nr, ~1 | species.latin), 
	     method = "REML", data = data.mollusca)
	meta.mollusca.mm

	### HETEROGENEITY cf.  http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate ###
		
		W <- diag(1/data.mollusca.2$vi)
		X <- model.matrix(meta.mollusca.mm)
		P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
		I.square.mollusca.1<-100 * sum(meta.mollusca.mm$sigma2) / (sum(meta.mollusca.mm$sigma2) + (meta.mollusca.mm$k-meta.mollusca.mm$p)/sum(diag(P)))
		round(I.square.mollusca.1, digits=2)
		I.square.mollusca.2<-100 * meta.mollusca.mm$sigma2 / (sum(meta.mollusca.mm$sigma2) + (meta.mollusca.mm$k-meta.mollusca.mm$p)/sum(diag(P)))
		round(I.square.mollusca.2, digits=2)




###############################################################################################################
###############################################################################################################
		
		


#################
###	AMPHIBIA  ###
#################

	test.amphibia<-data.11[ which(data.11$taxon.for.plot=='amphibia'),]
	head(test.amphibia)

	##############################
	### calculate effect sizes ###
	##############################

	data.amphibia <- escalc(measure="SMDH", 
				m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
				m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
				data=test.amphibia)
	str(data.amphibia)
	data.amphibia.2<-na.omit(data.amphibia)
	str(data.amphibia.2)

	#######################################
	### COVARIANCE MATRIX TO TAKE CARE	###
	### 	OF NON-INDEPENDENCE		        ###
	#######################################

	data.amphibia.2$unit<-factor(1:dim(data.amphibia.2)[1])
	VCV.amphibia<- matrix(0,nrow = dim(data.amphibia.2)[1],ncol = dim(data.amphibia.2)[1])
	rownames(VCV.amphibia) <- data.amphibia.2$unit
	colnames(VCV.amphibia) <- data.amphibia.2$unit 

	# find start and end coordinates for the subsets
	shared_coord <- which(data.amphibia.2$independent.effect.id %in% data.amphibia.2$independent.effect.id[duplicated(data.amphibia.2$independent.effect.id)]==TRUE)
	# matrix of combinations of coordinates for each experiment with shared control
	combinations <- do.call("rbind", tapply(shared_coord, data.amphibia.2[shared_coord,"independent.effect.id"], function(x) t(combn(x,2))))
	# calculate covariance values between  values at the positions in shared_list and place them on the matrix
	# NOTE: if control group is shared between experimental group covariance equal variance of the control group

	for (i in 1:dim(combinations)[1]){
  		p1 <- combinations[i,1]
  		p2 <- combinations[i,2]
  		p1_p2_cov <-0.5*sqrt(data.amphibia.2[p1,"vi"])*sqrt(data.amphibia.2[p2,"vi"])
  		VCV.amphibia[p1,p2] <- p1_p2_cov
  		VCV.amphibia[p2,p1] <- p1_p2_cov
	}

	#create variance-covariance matrix
	diag(VCV.amphibia) <- data.amphibia.2$vi
	is.positive.definite(VCV.amphibia) #needs to be positive definite so it can be inverted in analyses (must be invertable)


	###################### FIT MODEL with study, species and case number as random factors ######################

	meta.amphibia.mm<-rma.mv(yi = abs(yi), V = VCV.amphibia, random = list(~1 | study, ~1 | case.nr, ~1 | species.latin), 
	     method = "REML", data = data.amphibia.2)
	meta.amphibia.mm

	### HETEROGENEITY cf.  http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate ###
		
		W <- diag(1/data.amphibia.2$vi)
		X <- model.matrix(meta.amphibia.mm)
		P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
		I.square.amphibia.1<-100 * sum(meta.amphibia.mm$sigma2) / (sum(meta.amphibia.mm$sigma2) + (meta.amphibia.mm$k-meta.amphibia.mm$p)/sum(diag(P)))
		round(I.square.amphibia.1, digits=2)
		I.square.amphibia.2<-100 * meta.amphibia.mm$sigma2 / (sum(meta.amphibia.mm$sigma2) + (meta.amphibia.mm$k-meta.amphibia.mm$p)/sum(diag(P)))
		round(I.square.amphibia.2, digits=2)
		
		
		
		
###############################################################################################################
###############################################################################################################
		


###############
###	FISH	  ###
###############

	test.fish<-data.11[ which(data.11$taxon.for.plot=='fish'),]
	head(test.fish)

	##############################
	### calculate effect sizes ###
	##############################

	data.fish <- escalc(measure="SMDH", 
				m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
				m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
				data=test.fish)
	str(data.fish)
	data.fish.2<-na.omit(data.fish)
	str(data.fish.2)

	#######################################
	### COVARIANCE MATRIX TO TAKE CARE	###
	### 	OF NON-INDEPENDENCE		        ###
	#######################################

	data.fish.2$unit<-factor(1:dim(data.fish.2)[1])
	VCV.fish<- matrix(0,nrow = dim(data.fish.2)[1],ncol = dim(data.fish.2)[1])
	rownames(VCV.fish) <- data.fish.2$unit
	colnames(VCV.fish) <- data.fish.2$unit 

	# find start and end coordinates for the subsets
	shared_coord <- which(data.fish.2$independent.effect.id %in% data.fish.2$independent.effect.id[duplicated(data.fish.2$independent.effect.id)]==TRUE)
	# matrix of combinations of coordinates for each experiment with shared control
	combinations <- do.call("rbind", tapply(shared_coord, data.fish.2[shared_coord,"independent.effect.id"], function(x) t(combn(x,2))))
	# calculate covariance values between  values at the positions in shared_list and place them on the matrix
	# NOTE: if control group is shared between experimental group covariance equal variance of the control group

	for (i in 1:dim(combinations)[1]){
  		p1 <- combinations[i,1]
  		p2 <- combinations[i,2]
  		p1_p2_cov <-0.5*sqrt(data.fish.2[p1,"vi"])*sqrt(data.fish.2[p2,"vi"])
  		VCV.fish[p1,p2] <- p1_p2_cov
  		VCV.fish[p2,p1] <- p1_p2_cov
	}

	#create variance-covariance matrix
	diag(VCV.fish) <- data.fish.2$vi
	is.positive.definite(VCV.fish) #needs to be positive definite so it can be inverted in analyses (must be invertable)

	###################### FIT MODEL with study, species and case number as random factors ######################

	meta.fish.mm<-rma.mv(yi = abs(yi), V = VCV.fish, random = list(~1 | study, ~1 | case.nr, ~1 | species.latin), 
	     method = "REML", data = data.fish.2)
	meta.fish.mm

	### HETEROGENEITY cf.  http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate ###
		
		W <- diag(1/data.fish.2$vi)
		X <- model.matrix(meta.fish.mm)
		P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
		I.square.fish.1<-100 * sum(meta.fish.mm$sigma2) / (sum(meta.fish.mm$sigma2) + (meta.fish.mm$k-meta.fish.mm$p)/sum(diag(P)))
		round(I.square.fish.1, digits=2)
		I.square.fish.2<-100 * meta.fish.mm$sigma2 / (sum(meta.fish.mm$sigma2) + (meta.fish.mm$k-meta.fish.mm$p)/sum(diag(P)))
		round(I.square.fish.2, digits=2)


###############################################################################################################
###############################################################################################################
		
		
	
################
###	AVES	   ###
################

	test.aves<-data.11[ which(data.11$taxon.for.plot=='aves'),]
	head(test.aves)

	##############################
	### calculate effect sizes ###
	##############################

	data.aves <- escalc(measure="SMDH", 
				m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
				m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
				data=test.aves)
	
	data.aves.2 <- na.omit(data.aves)

	#######################################
	### COVARIANCE MATRIX TO TAKE CARE	###
	### 	OF NON-INDEPENDENCE		        ###
	#######################################

	data.aves.2$unit<-factor(1:dim(data.aves.2)[1])
	VCV.aves<- matrix(0,nrow = dim(data.aves.2)[1],ncol = dim(data.aves.2)[1])
	rownames(VCV.aves) <- data.aves.2$unit
	colnames(VCV.aves) <- data.aves.2$unit 

	# find start and end coordinates for the subsets
	shared_coord <- which(data.aves.2$independent.effect.id %in% data.aves.2$independent.effect.id[duplicated(data.aves.2$independent.effect.id)]==TRUE)
	# matrix of combinations of coordinates for each experiment with shared control
	combinations <- do.call("rbind", tapply(shared_coord, data.aves.2[shared_coord,"independent.effect.id"], function(x) t(combn(x,2))))
	# calculate covariance values between  values at the positions in shared_list and place them on the matrix
	# NOTE: if control group is shared between experimental group covariance equal variance of the control group

	for (i in 1:dim(combinations)[1]){
  		p1 <- combinations[i,1]
  		p2 <- combinations[i,2]
  		p1_p2_cov <-0.5*sqrt(data.aves.2[p1,"vi"])*sqrt(data.aves.2[p2,"vi"])
  		VCV.aves[p1,p2] <- p1_p2_cov
  		VCV.aves[p2,p1] <- p1_p2_cov
	}

	#create variance-covariance matrix
	diag(VCV.aves) <- data.aves.2$vi
	is.positive.definite(VCV.aves) #needs to be positive definite so it can be inverted in analyses (must be invertable)



	###################### FIT MODEL with study, species and case number as random factors ######################

	meta.aves.mm<-rma.mv(yi = abs(yi), V = VCV.aves, random = list(~1 | study, ~1 | case.nr, ~1 | species.latin), 
	     method = "REML", data = data.aves.2)
	meta.aves.mm

	### HETEROGENEITY cf.  http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate ###
		
		W <- diag(1/data.aves.2$vi)
		X <- model.matrix(meta.aves.mm)
		P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
		I.square.aves.1<-100 * sum(meta.aves.mm$sigma2) / (sum(meta.aves.mm$sigma2) + (meta.aves.mm$k-meta.aves.mm$p)/sum(diag(P)))
		round(I.square.aves.1, digits=2)
		I.square.aves.2<-100 * meta.aves.mm$sigma2 / (sum(meta.aves.mm$sigma2) + (meta.aves.mm$k-meta.aves.mm$p)/sum(diag(P)))
		round(I.square.aves.2, digits=2)


###############################################################################################################
###############################################################################################################
		
	

#################
###	MAMMALIA  ###
#################
	
	test.mammalia<-data.11[ which(data.11$taxon.for.plot=='mammalia'),]

	head(test.mammalia)

	##############################
	### calculate effect sizes ###
	##############################

	data.mammalia <- escalc(measure="SMDH", 
				m1i=mean.control, sd1i=sd.control, n1i=sample.size.control, 
				m2i=mean.noise, sd2i=sd.noise,  n2i=sample.size.noise.1,
				data=test.mammalia)

	data.mammalia.2 <- na.omit(data.mammalia)


	#######################################
	### COVARIANCE MATRIX TO TAKE CARE	###
	### 	OF NON-INDEPENDENCE		        ###
	#######################################

	data.mammalia.2$unit<-factor(1:dim(data.mammalia.2)[1])
	VCV.mammalia<- matrix(0,nrow = dim(data.mammalia.2)[1],ncol = dim(data.mammalia.2)[1])
	rownames(VCV.mammalia) <- data.mammalia.2$unit
	colnames(VCV.mammalia) <- data.mammalia.2$unit 

	# find start and end coordinates for the subsets
	shared_coord <- which(data.mammalia.2$independent.effect.id %in% data.mammalia.2$independent.effect.id[duplicated(data.mammalia.2$independent.effect.id)]==TRUE)
	# matrix of combinations of coordinates for each experiment with shared control
	combinations <- do.call("rbind", tapply(shared_coord, data.mammalia.2[shared_coord,"independent.effect.id"], function(x) t(combn(x,2))))
	# calculate covariance values between  values at the positions in shared_list and place them on the matrix
	# NOTE: if control group is shared between experimental group covariance equal variance of the control group

	for (i in 1:dim(combinations)[1]){
  		p1 <- combinations[i,1]
  		p2 <- combinations[i,2]
  		p1_p2_cov <-0.5*sqrt(data.mammalia.2[p1,"vi"])*sqrt(data.mammalia.2[p2,"vi"])
  		VCV.mammalia[p1,p2] <- p1_p2_cov
  		VCV.mammalia[p2,p1] <- p1_p2_cov
	}

	#create variance-covariance matrix
	diag(VCV.mammalia) <- data.mammalia.2$vi
	is.positive.definite(VCV.mammalia) #needs to be positive definite so it can be inverted in analyses (must be invertable)



	###################### FIT MODEL with study, species and case number as random factors ######################

	meta.mammalia.mm<-rma.mv(yi = abs(yi), V = VCV.mammalia, random = list(~1 | study, ~1 | case.nr, ~1 | species.latin), 
	     method = "REML", data = data.mammalia.2)
	meta.mammalia.mm

	### HETEROGENEITY cf.  http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate ###
		
		W <- diag(1/data.mammalia.2$vi)
		X <- model.matrix(meta.mammalia.mm)
		P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
		I.square.mammalia.1<-100 * sum(meta.mammalia.mm$sigma2) / (sum(meta.mammalia.mm$sigma2) + (meta.mammalia.mm$k-meta.mammalia.mm$p)/sum(diag(P)))
		round(I.square.mammalia.1, digits=2)
		I.square.mammalia.2<-100 * meta.mammalia.mm$sigma2 / (sum(meta.mammalia.mm$sigma2) + (meta.mammalia.mm$k-meta.mammalia.mm$p)/sum(diag(P)))
		round(I.square.mammalia.2, digits=2)



		
####################################################################################################################
####      PUBLICATION BIAS                                                                                    ######
####################################################################################################################		
		

		precision<-sqrt(1/data.phylo.2$vi)
		str(data.phylo.2)
		
		# Egger regression
		data.phylo.2.egger<-rma.mv(yi,vi, mod = precision, 
		                           random = list(~1 | study,
		                                         ~1 | species.latin,
		                                         ~1 | case.nr), 
		                           data=data.phylo.2, 
		                           method="REML")
		
		summary(data.phylo.2.egger)
		
		
		# time-lag bias
		meta.data.phylo.2.raw<-rma.mv (yi=yi, mod=year, 
		                               V=vi, 
		                               random = list(~1 | study, ~ 1 | species.latin, ~1 | case.nr), 
		                               data=data.phylo.2, 
		                               method="REML")
		meta.data.phylo.2.raw
		
		#extract residuals
		resid.1<-rstandard.rma.mv(meta.data.phylo.2.raw, type="rstandard")
		
		#funnel plot
		gg.funnel <- ggplot(
		  data.phylo.2, aes(x=precision, y=resid.1$resid)) + 
		  geom_point(aes(col=study), size=2) +
		  theme_bw()+  
		  theme(legend.position="None") +
		  geom_hline(yintercept = 0) +
		  theme(axis.line = element_line(colour = "black"),
		        panel.grid.major = element_blank(),
		        panel.grid.minor = element_blank(),
		        panel.border = element_blank(),
		        panel.background = element_blank())+ 
		  xlim(c(0, 20)) + 
		  ylim(c(-20, 20)) + 
		  labs( 
		    y="meta-analytic residuals", 
		    x="precision" 
		  )
		
		plot(gg.funnel)
		ggsave("funnel.plot.SMDH.png", plot=gg.funnel)
		
		
		
		
		# bubble plot
		# divide the sample size variable in a categorical variable with 4 levels
		data.phylo.2$category <- cut(data.phylo.2$sample.size.control, 
		                             breaks=c(-Inf, 20, 30, Inf), 
		                             labels=c("1","2","3"))
		
		
		### plot yi over years with vi as bubble size and bubble colour as size factor
		plot.yi.vi.sample <- ggplot(data=data.phylo.2, aes(x =year , y = yi)) +
		  geom_point(aes(size=vi), 
		             show.legend=FALSE,
		             colour="black",
		             stroke=1,
		             fill=factor(data.phylo.2$category),
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
		  scale_x_continuous(breaks = seq(2000, 2017, 2))+
		  labs(	y="effect size", 
		        x="year")+
		  geom_smooth(method="auto",  
		              linetype="dashed", 
		              span = 1, 
		              level = 0.95,
		              color="black", 
		              fill="grey69")
		
		
		plot.yi.vi.sample
		
		
		
		# time-lag bias model
		
		meta.data.phylo.2.year<-rma.mv (yi=yi, mod=year, 
		                        V=vi, 
		                        random = list(~1 | study, ~ 1 | species.latin, ~1 | case.nr), 
		                        data=data.phylo.2, 
		                        method="REML")
		meta.data.phylo.2.year	
		
		