#############################################################
###      Analysis of NACP_TERRA trait and stand data
###         to assess within-species trait variation
##############################################################
# this code performs the analyses and creates the figures 
# for Anderegg et al. (2018) Within-species patterns challenge our understanding of the Leaf Economics Spectrum. Ecol Let. 21:5 https://doi.org/10.1111/ele.12945
# and writes results and figures into a new directory named in L43
# It requires:
# - setting the working directory to a folder with all the necessary data files

##### NOTE 1: All data are available online, either as previously published datasets:
# https://doi.org/10.1038/nature02403 for GLOPNET,
# dx.doi.org/10.3334/ORNLDAAC/1292 for NACP TERR-PNW,
# https:/doi.org/10.1111/1365-2435.12790 for Coffea Arabica).
# or in the Dryad repository https://doi.org/10.5061/dryad.c1dn34b for previously unpublished trait data.

##### NOTE 2: the Source data for this code are somewhat modified from these sources (e.g. family names were added to GLOPNET). 
# Thus, this code is not fully repeatable. If you would like the Source Data, please email Leander Anderegg (leanderegg@gmail.com)
# the derived dataset SlopesAndCorrelations_alltaxa_clean.cs is provided in the Github repository,
# which can be used to replicate the statistical analyses and figure plotting for Figs 2-4

##### NOTE 3: if you want to take Source Data and create the derived datasets:
# - all.data (all trait data from GLOPNET, PNW and additional intraspecific studies)
# - biomass (plot data with trait Community Weighted Means)
# - traits (trait data for the PNW plot network)
# run Lines 68-539 (BEGIN to END of "Dataset creation and cleaning")

# ** if you want to import derived datasets and just run analyses, 
# uncomment code at Line 553 "BEGIN: Variance Decompostion and Figure 1"
# to import all.data, biomass, and traits

# **** if you do not want to run the analysis to calculate Trait-Trait relationships and Null models at different taxonomic scales
# (this can take multiple hours)
# skip the section:
# ". RUN ANALYSIS: Calculate Trait-Trait Relationship & Null Models"
#instead, uncomment Lines after 1768 under "BEGIN: Taxonomic Scale Statistical Tests"




# load required packages #
require(lme4)
require(lmerTest)
require(MuMIn)
require(lmodel2)
require(dplyr)
require(reshape)
require(RColorBrewer)
require(lattice)

require(mgcv)
require(car)
require(ggplot2)
require(stringr)
require(stringi)


# standard error function
sterr <- function(x,...){
  sd(x, na.rm=T)/ length(na.omit(x))
}


#set color palette
mypal <- brewer.pal(n=9, "Set1")
palette(mypal)
# set the colors that wil denote within-species, within-genus, within family and across CWMs
colchoices <- c(1,2,4,3,6)

# set working directory
setwd("/Users/leeanderegg/Dropbox/NACP_Traits/NACP_Traits_Rcode/FinalExample")

# Create a directory to store results in
results_dirname <- "20180219_results"
dir.create(results_dirname)



#__________________________________________________________________________________
######## ** BEGIN: Dataset creation and cleaning: #####################################
#__________________________________________________________________________________

############ . Load GLOPNET data ########
# these data were downloaded from https://www.nature.com/articles/nature02403 on 1/15/17
# species names were then cleaned and family data added using the {taxize} R package, and finalized by hand
# for script detailing this cleaning, please email Leander Anderegg at leanderegg@gmail.com
# 9 measurements were unable to be identified to family and were removed from the analysis
# all measurements ID'd to genus (e.g. "Pelea sp") were treated as a unique species
#note: this results in the averaging of 3 "Rubus sp" and 2 "Nephrolepsis sp", but all other unidentified species are sp1,sp2 etc.

# import site climate data, unmodified from Wright et al. 2004 Nature
LESclim <- read.csv("SourceData/Wright2004_sitedata.csv")
# import cleaned GLOPNET trait database. For cleaning code to recreate this
# dataset, please email Leander Anderegg at leanderegg@gmail.com
LES <- read.csv("SourceData/GLOPNETtraitdata_cleaned_121917_v5.csv", row.names=1)
# still have 9 NAs for family that couldn't be taxonomically id'ed

# create unlogged trait values for averaging to higher taxonomic levels
LES$LMA <- 10^LES$log.LMA 
LES$LL <- 10^LES$log.LL
LES$Nmass <- 10^LES$log.Nmass
LES$Narea <- 10^LES$log.Narea
LES$Aarea <- 10^LES$log.Aarea
LES$Amass <- 10^LES$log.Amass
LES$SLA <- 1/LES$LMA # in GLOPNET, SLA and LMA are in different units (g/m2 vs cm2/g)

## add in climate
LES$MAT <- LESclim$MAT[match(LES$Dataset, LESclim$Dataset)]
LES$MAP <- LESclim$Rain[match(LES$Dataset, LESclim$Dataset)]
LES$VPD <- LESclim$VPD[match(LES$Dataset, LESclim$Dataset)]
LES$RAD <- LESclim$RAD[match(LES$Dataset, LESclim$Dataset)]
LES$PET <- LESclim$PET[match(LES$Dataset, LESclim$Dataset)]








######## . Load PNW data #########################
#these data were downloaded from  http://dx.doi.org/10.3334/ORNLDAAC/1292 on 1/23/16 by LDLA
# see Berner & Law 2016 for data description
# additional climate data was added to the trait dataset using PRISM and VIC
# (see Methods in main text)


## _____________Species ID lookup table to relate traits and biomass ______________________
speciesID <- read.csv("SourceData/PNW_SpeciesNames_lookup_031316.csv", header=T)
# this table relates the species abbreviations used in the NACP_TERRA_PNW forest_biomass and trait measurement datasets.


## ____________________________Soil Data______________________________
soil <- read.csv("SourceData/NACP_TERRA_PNW_soil_cleaned.csv", header=T, na.strings = "-9999")

# Cleaning some data inconsistencies
#soil[which(soil$Layer=="top" & soil$UpperDepth>5),]
# Plot 252 has two top layers. need to switch the second to 'bottom'
# Plot 86 has top and bottom layers flipped
## and something on the order of 9 plots don't have a 'top layer
soil$Layer[which(soil$PLOT_ID==252)] <- c("top", "middle","bottom")
soil$Layer[which(soil$PLOT_ID==86)] <- c("top", "bottom")
probs2 <- names(xtabs(~PLOT_ID, soil)[which(xtabs(~PLOT_ID, soil)==2)])
probs2names <- xtabs(~PLOT_ID, soil[which(soil$PLOT_ID %in% probs2 & soil$Layer != "top"),])
## 8 plots have 'middle' and 'bottom', but middle starts at 0. I'm going to just make all 'middles' into 'tops
soil$Layer[which(soil$PLOT_ID %in% names(probs2names[which(probs2names==2)]) & soil$Layer=="middle")] <- "top" # 8 plots have only 2 layers and one of them is not 'top'
# OK, that seemed to solve it. Everything has a 'top' layer, but below that is unknown. could be 'middle' could be 'bottom'...
# thus, we'll work with only the top layer for now.
soil.top <- soil[which(soil$Layer=="top"), ]



###_____________ Stand characteristics dataset______________________
# import stand characteristics (biomass, species basal area fractions, productivity)
biomass <- read.csv("SourceData/NACP_TERRA_PNW_forest_biomass_productivity_v2.csv", header= T,na.strings = "-9999" )

# merge soil characteristics and biomass data
biomass$soil_N <- soil.top$soil_N[match(biomass$PLOT_ID, soil.top$PLOT_ID)]
biomass$soil_pH <- soil.top$soil_pH[match(biomass$PLOT_ID, soil.top$PLOT_ID)]



### _______________________Traits dataset ____________________________________________

# this is the traits data as downloaded from the ORNL repo, plus climate variables pulled by L Berner (see main text)
traits <- read.csv("SourceData/NACP_TERRA_PNW_leaf_trait_v1_plusClim.csv", header=T, na.strings="-9999")
# make sure the correct columns are factors
fac.colstraits <- c(1,5:7,12:16)
for(j in fac.colstraits){ traits[,j] <- factor(traits[,j])}
# make sure the correct columns are numeric
num.colstraits <- c(8:11,17:27)
for (i in num.colstraits){traits[,i] <- as.numeric(as.character(traits[,i]))}

# add in Genus.species column, some clim columns and a column for species IDs used in biomass dataset
traits$GE.SP <- paste(traits$GENUS, traits$SPECIES, sep=".")
traits$SP.ID <- speciesID$bio.sp[match(traits$GE.SP, speciesID$traits.sp)]
#sum soil moisture for total moisture content
traits$soilmoist.all.mm <- apply(traits[, c("soilmoist.lvl1.mm", "soilmoist.lvl2.mm", "soilmoist.lvl3.mm")],MARGIN = 1, FUN=sum)
# Id the measurements that come from a dominant species,,,
traits$FOREST_TYPE <- biomass$SPP_O1_ABBREV[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$dominance <- biomass$SPP_O1_BASAL_AREA_FRACTION[match(traits$PLOT_ID, biomass$PLOT_ID)]
# note: this only gives the BA fraction of the dominant spp, not of the spp for which the trait is measured. if FOREST_TYPE!=GE.SP, this is meaningless
# Approximate Stand Age
traits$ASA <- biomass$ASA[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$ELEVATION <- biomass$ELEVATION[match(traits$PLOT_ID, biomass$PLOT_ID)]
traits$AG_TGROWTH <- biomass$AG_PROD_TREE_TOTAL_AS_CARBON[match(traits$PLOT_ID, biomass$PLOT_ID)]

## add in soil characteristics
traits$soil_N <- soil.top$soil_N[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$soil_C <- soil.top$soil_C[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$soil_pH <- soil.top$soil_pH[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$TotalSoilDepth <- soil.top$TotalDepth[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$TopSoilDepth <- soil.top$LowerDepth[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$Bulk_Density <- soil.top$Bulk_Density[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$pSAND <- soil.top$pSAND[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$pSILT <- soil.top$pSILT[match(traits$PLOT_ID, soil.top$PLOT_ID)]
traits$pCLAY <- soil.top$pCLAY[match(traits$PLOT_ID, soil.top$PLOT_ID)]


## calculate traits used in LES
traits$LMA <- with(traits, LEAF_DRY_WT/(LEAF_HSA*1/100^2)) # LES works in LMA rather than SLA
# also SLA is in gC rather than LEAF_DRY_WT, and in cm2 rather than m2 as in Wright 2004. so I need to recalculate based on leaf dry weight and hemispheric leaf area
# bunch of EPA plots are missing LEAF_DRY_WT and LEAF_HSA, so have to back calculate from LEAF_CARBON and SLA_HSA
traits$LMA[which(is.na(traits$LMA))] <- with(traits[which(is.na(traits$LMA)),], 1/(SLA_HSA * LEAF_CARBON/100 * 1/100^2))
# Note: recalculating LMA for all plots using this method yeilds LMA values w/ mean difference of 0.00044, so essentially rounding error
traits$SLA_drymass <- 1/traits$LMA # make a non-carbon SLA
traits$log.LMA <- with(traits, log(LMA,base = 10))
traits$LLmonths <- traits$LEAF_LIFE * 12 # LES works in months rather than years
traits$log.LL <- log(traits$LLmonths,base = 10)
traits$log.Nmass <- log(traits$LEAF_NITROGEN,base = 10)
traits$Narea <- traits$LEAF_NITROGEN/100 / traits$SLA_drymass #Narea = Nmass * Ml/Al (LMA)
traits$log.Narea <- log(traits$Narea,base = 10)


traits$FullSpecies <- paste(traits$GENUS, traits$SPECIES, sep=" ")

# make a PCA of climate variables and add them to the df
climpca.traits <- prcomp(traits[,c(grep("gy", colnames(traits)), which(colnames(traits) %in% c("soilmoist.lvl1.mm","soilmoist.all.mm")))],scale=T)
traits$climPC1 <- climpca.traits$x[,1]
traits$climPC2 <- climpca.traits$x[,2]
traits$climPC3 <- climpca.traits$x[,3]


## ______________ adding in family names to Traits ___________________
PNWfams <- read.csv("SourceData/PNW_Families_031717.csv")
traits$Family <- PNWfams$family[match(traits$FullSpecies, PNWfams$query)]

write.csv(traits, paste0("./",results_dirname,"/PNW_Trait_data_for_Anderegg_etal_2018_EcolLet.csv"))


############### . calculate CWMs for PNW dataset ############

# note on nomenclature:
# TRAIT1 = species mean trait value for species 1
# wTRAIT1 = BA weighted species mean trait value
# TRAITp1 = plot/site mean value for the species species 1 - originally used s, but that was horrific. changed 06.14.17
# wTRAITp1 = BA weighted plot/site mean value for species 1

## Species Means (across all sites)
# Note: mLMA = LMA_HSA, and mLMA_PSA = LMA_PSA, I've used LMA_HSA in full dataset creation
# also, mlog.Trait = mean of logged traits
# log.Trait = log of mean traits

spp.traits <- traits %>% group_by(SP.ID) %>% summarise(nsample = n(), SLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), CN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                       , LIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), nplots =length(unique(PLOT_ID))
                                                       , mLMA = mean(LMA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNmass = mean(LEAF_NITROGEN, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                       , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                       , climPC1 = mean(climPC1, na.rm=T), climPC2 = mean(climPC2, na.rm=T), climPC3 = mean(climPC3, na.rm=T)
                                                       , soil_N= mean(soil_N, na.rm=T), soil_pH=mean(soil_pH, na.rm=T), ASA = mean(ASA, na.rm=T), LAI_O=mean(LAI_O, na.rm=T), AG_TGROWTH = mean(AG_TGROWTH, na.rm=T))
# log species mean traits for between species analysis
spp.traits$log.LMA <- log(spp.traits$mLMA, base=10)
spp.traits$log.LL <- log(spp.traits$mLLmonths, base=10)
spp.traits$log.Nmass <- log(spp.traits$mNmass, base=10)
spp.traits$log.Narea <- log(spp.traits$mNarea, base=10)


# calculate Species means for each plot in which they occur
spp.plot.traits <- traits %>% group_by(SP.ID, PLOT_ID) %>% summarise(nsample = n(), mSLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), mCN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
                                                                     , mLIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), mNmass = mean(LEAF_NITROGEN, na.rm=T), nplots =length(unique(PLOT_ID))
                                                                     , mLMA = mean(LMA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNarea=mean(Narea, na.rm=T)
                                                                     , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
                                                                     , climPC1 = unique(climPC1), climPC2 = unique(climPC2), climPC3 = unique(climPC3)
)
# create unique species-plot tag
spp.plot.traits$SP.PLOT <- paste(spp.plot.traits$SP.ID, spp.plot.traits$PLOT_ID, sep="-")


### make unique identifiers for matching biomass and spp.plot.traits rows
biomass$SP1.PLOT <- paste(biomass$SPP_O1_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP2.PLOT <- paste(biomass$SPP_O2_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP3.PLOT <- paste(biomass$SPP_O3_ABBREV,biomass$PLOT_ID, sep="-")
biomass$SP4.PLOT <- paste(biomass$SPP_O4_ABBREV,biomass$PLOT_ID, sep="-")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############ ++ CWMs based on species-plot mean trait values #######################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# naming convention:
#   "wTraitp1" = the plot mean trait value for species 1, weighted by its basal area fraction at that plot
#   "wTrait1" = the species mean trait value for species 1, weighted by its basal area fraction (for infilling species with no measurements at a plot)
#   "cw_Traitp" = the sum of basal area fraction-weighted traits, i.e. the community weighted mean (CWM) for a plot based on plot level trait measurements
#    the p indicates trait values were measured for that species, at that plot. no p in variable name indicates that species mean trait value was used.
#   "if" or "_if" indicates values with missing plot-level data that have been infilled with species mean data

#____________________________________ LMA  ______________________________________
biomass$wLMAp1 <- spp.plot.traits$mLMA[match(as.character(biomass$SP1.PLOT), as.character(spp.plot.traits$SP.PLOT))] * biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLMAp2 <- spp.plot.traits$mLMA[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLMAp3 <- spp.plot.traits$mLMA[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLMAp4 <- spp.plot.traits$mLMA[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)] * biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_LMAp <- apply(biomass[,c("wLMAp1", "wLMAp2","wLMAp3","wLMAp4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_LMAp[which(is.na(biomass$wLMAp1))] <- NA # 86 plots lacking LMA data of SPP01
biomass$cw_LMAp[which(is.na(biomass$wLMAp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
#18 sites lacking LMA data for SPP02 where SPP02 makes up >30 % of BA, but only 5 have SPP01
biomass$cw_LMAp[which(biomass$cw_LMAp==0)] <- NA
# 91 plots missing substantial LMA data
# length(which(is.na(biomass$cw_LMAp)))




#____________________________________ Leaf Life  ______________________________________
biomass$wLLp1 <- spp.plot.traits$mLLmonths[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLLp2 <- spp.plot.traits$mLLmonths[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLLp3 <- spp.plot.traits$mLLmonths[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLLp4 <- spp.plot.traits$mLLmonths[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_LLp <- apply(biomass[,c("wLLp1", "wLLp2","wLLp3","wLLp4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_LLp[which(is.na(biomass$wLLp1))] <- NA # 94 sites lacking LeafLL of SPP01
biomass$cw_LLp[which(is.na(biomass$wLLp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
# 22 sites lack SPP02 leaf life, 9 of them unique
length(which(is.na(biomass$cw_LLp))) # 103 plots sans significant LL

#____________________________________ Nmass  ______________________________________
biomass$wNmassp1 <- spp.plot.traits$mNmass[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNmassp2 <- spp.plot.traits$mNmass[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNmassp3 <- spp.plot.traits$mNmass[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNmassp4 <- spp.plot.traits$mNmass[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_Nmassp <- apply(biomass[,c("wNmassp1", "wNmassp2","wNmassp3","wNmassp4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_Nmassp[which(is.na(biomass$wNmassp1))] <- NA # 98 sites lacking LeafNITROGEN of SPP01
biomass$cw_Nmassp[which(is.na(biomass$wNmassp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
biomass$cw_Nmassp[which(biomass$cw_Nmassp==0)] <- NA
# 26 sites lack SPP02 leaf NITROGEN, 3 of them unique
# length(which(is.na(biomass$cw_Nmassp))) # 102 sites sans Nmass




#____________________________________ Narea  ______________________________________
biomass$wNareap1 <- spp.plot.traits$mNarea[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNareap2 <- spp.plot.traits$mNarea[match(biomass$SP2.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNareap3 <- spp.plot.traits$mNarea[match(biomass$SP3.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNareap4 <- spp.plot.traits$mNarea[match(biomass$SP4.PLOT, spp.plot.traits$SP.PLOT)]* biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_Nareap <- apply(biomass[,c("wNareap1", "wNareap2","wNareap3","wNareap4")],MARGIN=1,FUN=sum, na.rm=T)
# now need to remove values that got smoothed over due to na.rm=T in the apply function
biomass$cw_Nareap[which(is.na(biomass$wNareap1))] <- NA # 98 sites lacking LeafNarea of SPP01
biomass$cw_Nareap[which(is.na(biomass$wNareap2) & biomass$SPP_O2_BASAL_AREA_FRACTION>30)] <- NA
biomass$cw_Nareap[which(biomass$cw_Nareap==0)] <- NA
# 26 sites lack SPP02 leaf Narea, 3 of them unique



# logging community weighted trait means based on plot level trait measurements
biomass$log.cw_LMAp <- log(biomass$cw_LMAp, base=10)
biomass$log.cw_LLp <- log(biomass$cw_LLp, base=10)
biomass$log.cw_Nmassp <- log(biomass$cw_Nmassp, base=10)
biomass$log.cw_Nareap <- log(biomass$cw_Nareap, base=10)

###### add the climPCs to biomass
biomass$climPC1 <-spp.plot.traits$climPC1[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]
biomass$climPC2 <-spp.plot.traits$climPC2[match(biomass$SP1.PLOT, spp.plot.traits$SP.PLOT)]



#_____________________________________________________________________________
###### ++ Infilling spp trait values with spp means when they're missing ###########
#_____________________________________________________________________________
### Infilling LMA of spp with <25% of BA with species mean values


###________________________________ Creating BA-weighted LMA values using spp means ________________________________
biomass$wLMA1 <- spp.traits$mLMA[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLMA2 <- spp.traits$mLMA[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLMA3 <- spp.traits$mLMA[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLMA4 <- spp.traits$mLMA[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100

biomass$cw_LMA <- apply(biomass[,c("wLMA1", "wLMA2","wLMA3","wLMA4")],MARGIN=1,FUN=sum, na.rm=T)
biomass$cw_LMA[which(biomass$cw_LMA==0)] <- NA

wLMAif1 <- biomass$wLMAp1 # pull out the plot-level basal-area weighted values so I can infill missing values with the species mean
wLMAif1[which(is.na(wLMAif1))] <- biomass$wLMA1[which(is.na(wLMAif1))] # infill species means for missing SP1 traits
wLMAif2 <- biomass$wLMAp2
wLMAif2[which(is.na(wLMAif2))] <- biomass$wLMA2[which(is.na(wLMAif2))] # infill species means for missing SP2 traits
wLMAif3 <- biomass$wLMAp3
wLMAif3[which(is.na(wLMAif3))] <- biomass$wLMA3[which(is.na(wLMAif3))] # infill species means for missing SP3 traits
wLMAif4 <- biomass$wLMAp4
wLMAif4[which(is.na(wLMAif4))] <- biomass$wLMA4[which(is.na(wLMAif4))] # infill species means for missing SP4 traits
# calculate CWM, but with all missing trait values infilled with species means. Note: this infills some missing dominant species with species mean traits
cw_LMAp_if  <- apply(data.frame(wLMAif1, wLMAif2, wLMAif3, wLMAif4),MARGIN=1, FUN=sum, na.rm=T)
# get rid of plots that still have no data (na.rm=T means empty plots return 0)
cw_LMAp_if[which(cw_LMAp_if==0)] <- NA
# now remove plots where missing species represent >25% of basal area
cw_LMAp_if[which(is.na(biomass$wLMAp1) & biomass$SPP_O1_BASAL_AREA_FRACTION>25)] <- NA # 80 plots lacking LMA data of SPP01 & SPO1 has >25% of the basal area
cw_LMAp_if[which(is.na(biomass$wLMAp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>25)] <- NA # 26 plots lacking LMA of SPP2 w>25% of BA, 6 plots that weren't killed above
cw_LMAp_if[which(is.na(biomass$wLMAp3) & biomass$SPP_O3_BASAL_AREA_FRACTION>25)] <- NA # 2 plots lacking a prevalent SPP3, only one novel
#93 plots are removed because they are missing trait data from a dominant species (species with >25% of stand basal area)

biomass$cw_LMAp_if <- cw_LMAp_if # add infilled CWM LMA values to biomass




###________________________________ Leaf Lifespan infilling with spp mean values ________________________________
biomass$wLL1 <- spp.traits$mLLmonths[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wLL2 <- spp.traits$mLLmonths[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wLL3 <- spp.traits$mLLmonths[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wLL4 <- spp.traits$mLLmonths[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_LL <- apply(biomass[,c("wLL1", "wLL2","wLL3","wLL4")],MARGIN=1,FUN=sum, na.rm=T)
biomass$cw_LL[which(biomass$cw_LL==0)] <- NA

wLLif1 <- biomass$wLLp1
wLLif1[which(is.na(wLLif1))] <- biomass$wLL1[which(is.na(wLLif1))]
wLLif2 <- biomass$wLLp2
wLLif2[which(is.na(wLLif2))] <- biomass$wLL2[which(is.na(wLLif2))]
wLLif3 <- biomass$wLLp3
wLLif3[which(is.na(wLLif3))] <- biomass$wLL3[which(is.na(wLLif3))]
wLLif4 <- biomass$wLLp4
wLLif4[which(is.na(wLLif4))] <- biomass$wLL4[which(is.na(wLLif4))]
cw_LLp_if  <- apply(data.frame(wLLif1, wLLif2, wLLif3, wLLif4),MARGIN=1, FUN=sum, na.rm=T)
# get rid of plots that still have no data
cw_LLp_if[which(cw_LLp_if==0)] <- NA
cw_LLp_if[which(is.na(biomass$wLLp1) & biomass$SPP_O1_BASAL_AREA_FRACTION>25)] <- NA # 88 plots lacking LL data of SPP01 & SPO1 has >25% of the basal area
cw_LLp_if[which(is.na(biomass$wLLp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>25)] <- NA # 31 plots lacking LL of SPP2 w>25% of BA, 11 plots that weren't killed above
cw_LLp_if[which(is.na(biomass$wLLp3) & biomass$SPP_O3_BASAL_AREA_FRACTION>25)] <- NA # 2 plots lacking a prevalent SPP3, only one novel

biomass$cw_LLp_if <- cw_LLp_if


###________________________________ Nmass infilling with spp mean values ________________________________
biomass$wNmass1 <- spp.traits$mNmass[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNmass2 <- spp.traits$mNmass[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNmass3 <- spp.traits$mNmass[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNmass4 <- spp.traits$mNmass[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_Nmass <- apply(biomass[,c("wNmass1", "wNmass2","wNmass3","wNmass4")],MARGIN=1,FUN=sum, na.rm=T)
biomass$cw_Nmass[which(biomass$cw_Nmass==0)] <- NA

wNmassif1 <- biomass$wNmassp1
wNmassif1[which(is.na(wNmassif1))] <- biomass$wNmass1[which(is.na(wNmassif1))]
wNmassif2 <- biomass$wNmassp2
wNmassif2[which(is.na(wNmassif2))] <- biomass$wNmass2[which(is.na(wNmassif2))]
wNmassif3 <- biomass$wNmassp3
wNmassif3[which(is.na(wNmassif3))] <- biomass$wNmass3[which(is.na(wNmassif3))]
wNmassif4 <- biomass$wNmassp4
wNmassif4[which(is.na(wNmassif4))] <- biomass$wNmass4[which(is.na(wNmassif4))]
cw_Nmassp_if  <- apply(data.frame(wNmassif1, wNmassif2, wNmassif3, wNmassif4),MARGIN=1, FUN=sum, na.rm=T)
# get rid of plots that still have no data
cw_Nmassp_if[which(cw_Nmassp_if==0)] <- NA
cw_Nmassp_if[which(is.na(biomass$wNmassp1) & biomass$SPP_O1_BASAL_AREA_FRACTION>25)] <- NA # 93 plots lacking Nmass data of SPP01 & SPO1 has >25% of the basal area
cw_Nmassp_if[which(is.na(biomass$wNmassp2) & biomass$SPP_O2_BASAL_AREA_FRACTION>25)] <- NA # 36 plots lacking Nmass of SPP2 w>25% of BA, 4 plots that weren't killed above
cw_Nmassp_if[which(is.na(biomass$wNmassp3) & biomass$SPP_O3_BASAL_AREA_FRACTION>25)] <- NA # 2 plots lacking a prevalent SPP3, only one novel

biomass$cw_Nmassp_if <- cw_Nmassp_if



##________________________________ Narea infilling with spp mean values ________________________________
biomass$wNarea1 <- spp.traits$mNarea[match(biomass$SPP_O1_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O1_BASAL_AREA_FRACTION/100
biomass$wNarea2 <- spp.traits$mNarea[match(biomass$SPP_O2_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O2_BASAL_AREA_FRACTION/100
biomass$wNarea3 <- spp.traits$mNarea[match(biomass$SPP_O3_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O3_BASAL_AREA_FRACTION/100
biomass$wNarea4 <- spp.traits$mNarea[match(biomass$SPP_O4_ABBREV, spp.traits$SP.ID)]*biomass$SPP_O4_BASAL_AREA_FRACTION/100
biomass$cw_Narea <- apply(biomass[,c("wNarea1", "wNarea2","wNarea3","wNarea4")],MARGIN=1,FUN=sum, na.rm=T)
biomass$cw_Narea[which(biomass$cw_Narea==0)] <- NA

wNareaif1 <- biomass$wNareap1
wNareaif1[which(is.na(wNareaif1))] <- biomass$wNarea1[which(is.na(wNareaif1))]
wNareaif2 <- biomass$wNareap2
wNareaif2[which(is.na(wNareaif2))] <- biomass$wNarea2[which(is.na(wNareaif2))]
wNareaif3 <- biomass$wNareap3
wNareaif3[which(is.na(wNareaif3))] <- biomass$wNarea3[which(is.na(wNareaif3))]
wNareaif4 <- biomass$wNareap4
wNareaif4[which(is.na(wNareaif4))] <- biomass$wNarea4[which(is.na(wNareaif4))]
cw_Nareap_if  <- apply(data.frame(wNareaif1, wNareaif2, wNareaif3, wNareaif4),MARGIN=1, FUN=sum, na.rm=T)
# get rid of plots that still have no data
cw_Nareap_if[which(cw_Nareap_if==0)] <- NA
cw_Nareap_if[which(is.na(biomass$wNareap1) & biomass$SPP_O1_BASAL_AREA_FRACTION>25)] <- NA # 93 plots lacking Narea data of SPP01 & SPO1 has >25% of the basal area
cw_Nareap_if[which(is.na(biomass$wNareap2) & biomass$SPP_O2_BASAL_AREA_FRACTION>25)] <- NA # 36 plots lacking Narea of SPP2 w>25% of BA, 4 plots that weren't killed above
cw_Nareap_if[which(is.na(biomass$wNareap3) & biomass$SPP_O3_BASAL_AREA_FRACTION>25)] <- NA # 2 plots lacking a prevalent SPP3, only one novel

biomass$cw_Nareap_if <- cw_Nareap_if


biomass$log.cw_LMAp_if <- log(biomass$cw_LMAp_if, base=10)
biomass$log.cw_LLp_if <- log(biomass$cw_LLp_if, base=10)
biomass$log.cw_Nmassp_if <- log(biomass$cw_Nmassp_if, base=10)
biomass$log.cw_Nareap_if <- log(biomass$cw_Nareap_if, base=10)

write.csv(biomass, paste0("./",results_dirname,"/PNW_Biomass_data_for_Anderegg_etal_2018_EcolLet.csv"))


#________________ End CWM trait calculation for PNW dataset _________________________







############### . Load supplemental w/in spp data ############
# This comes from Martin et al. 2015 on Coffea arabica
# plus a number of unpublished datasets collected primarily by LDL Anderegg
# these were combined into one dataframe. please email LDL Anderegg at leanderegg@gmail.com for additional metadata.


data.supp <- read.csv("SourceData/WithinSpecies_Trait_Data_040117.csv", header=T, row.names=1)
# make climate columns that align with PNW and LES
data.supp$MAT <- NA
data.supp$MAP <- NA
data.supp$VPD <- NA




############### . Load additional w/in spp data from Mt Rainier############
# this includes data from 5 perennial wildflowers, collected across an elevation gradient by by ML Sethi in 2017
Rainier <- read.csv("SourceData/WithinSpecies_MtRainier_alltraits_v1_171211.csv", row.names=1)
Rainier$MAT <- NA
Rainier$MAP <- NA
Rainier$VPD <- NA
Rainier$Project <- "Rainier"






############### . Making Master Dataset ###################


###### combine PNW and glopnet dataset. ###
data1 <- traits %>% select(FullSpecies,log.LMA, log.LL, log.Nmass, log.Narea, GENUS, Family, tmean.gy.c,ppt.gy.mm,vpd.gy.max)
colnames(data1)[c(1,2,6,8:10)]<- c("Species", "log.LMA","Genus","MAT","MAP","VPD")
data1$Project <- rep("PACNW", times=nrow(data1))
data1$log.LL[which(data1$log.LL<1.2)] <- NA # all the deciduous species have '1yr' lifespan, but really that's wrong
data2 <- LES %>% filter(!is.na(Family)) %>% select(Species, log.LMA, log.LL, log.Nmass, log.Narea,Genus, Family, MAT, MAP, VPD)
# currently removes 9 records as of 12.19.17
data2$Project <- rep("GLOPNET", times=nrow(data2))

# select relevant columns from supplemental data
data3 <- data.supp %>% select(Species,log.LMA,log.LL, log.Nmass,log.Narea,Genus,Family,MAT, MAP, VPD, Project)

# select relevant columns from Mt. Rainier data
data4 <- Rainier %>% select(Species=Species_full, log.LMA, log.LL, log.Nmass, log.Narea, Genus, Family, MAT, MAP, VPD, Project)

## make combined dataset with all the taxonomically resolved species
data.all <- rbind(data1, data2, data3, data4) # new w/ Rainier: 4267 measurements total
data.all$Species <- factor(data.all$Species)
data.all$Genus <- factor(data.all$Genus)
data.all$Family <- factor(data.all$Family)
data.all$Project <- factor(data.all$Project)

# write the full trait data to the Results Directory
write.csv(data.all, paste0("./",results_dirname,"/AllTraitData_for_Anderegg_etal_2018_EcolLet.csv"))

#__________________________________________________________________________________
################ ** END: Dataset creation and cleaning #####################
#__________________________________________________________________________________








#__________________________________________________________________________________
################ BEGIN: Variance Decomposition and Figure 1 #####################
#__________________________________________________________________________________

# _____________________________________________________________________
##### to skip creation of **full datasets**, uncomment this code: 
# biomass <- read.csv("DerivedData/PNW_Biomass_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names = 1)
# traits <- read.csv("DerivedData/PNW_Trait_data_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names=1)
# data.all <- read.csv("DerivedData/AllTraitData_for_Anderegg_etal_2018_EcolLet.csv", header=T, row.names = 1)
# # create a dataset of species mean traits just for PNW dataset:
# spp.traits <- traits %>% group_by(SP.ID) %>% summarise(nsample = n(), SLA = mean(SLA_HSA, na.rm=T), nSLA = n()- length(which(is.na(SLA_HSA))), CN = mean(LEAF_CN, na.rm=T), nCN = n()- length(which(is.na(LEAF_CN)))
#                                                        , LIFE = mean(LEAF_LIFE, na.rm=T), nLIFE=n()- length(which(is.na(LEAF_LIFE))), nplots =length(unique(PLOT_ID))
#                                                        , mLMA = mean(LMA, na.rm=T), mLLmonths= mean(LLmonths, na.rm=T), mNmass = mean(LEAF_NITROGEN, na.rm=T), mNarea=mean(Narea, na.rm=T)
#                                                        , mlog.LMA = mean(log.LMA, na.rm=T), mlog.LL = mean(log.LL, na.rm=T), mlog.Narea=mean(log.Narea, na.rm=T), mlog.Nmass=mean(log.Nmass, na.rm=T)
#                                                        , climPC1 = mean(climPC1, na.rm=T), climPC2 = mean(climPC2, na.rm=T), climPC3 = mean(climPC3, na.rm=T)
#                                                        , soil_N= mean(soil_N, na.rm=T), soil_pH=mean(soil_pH, na.rm=T), ASA = mean(ASA, na.rm=T), LAI_O=mean(LAI_O, na.rm=T), AG_TGROWTH = mean(AG_TGROWTH, na.rm=T))
# # log species mean traits for between species analysis
# spp.traits$log.LMA <- log(spp.traits$mLMA, base=10)
# spp.traits$log.LL <- log(spp.traits$mLLmonths, base=10)
# spp.traits$log.Nmass <- log(spp.traits$mNmass, base=10)
# spp.traits$log.Narea <- log(spp.traits$mNarea, base=10)

# _____________________________________________________________________





###### . Full dataset, global variation #######################

## with all spp used for the hierarchical analysis
logLMAvar <- lmer(log.LMA~ 1 + (1|Project) + (1|Family) + (1|Genus) + (1|Species), data.all)
logLLvar <- lmer(log.LL~ 1 + (1|Project)  + (1|Family) + (1|Genus) + (1|Species), data.all[-which(data.all$Project=="CO"),])
logNmassvar <- lmer(log.Nmass~ 1 + (1|Project)  + (1|Family) + (1|Genus) + (1|Species), data.all)
logNareavar <- lmer(log.Narea~ 1 + (1|Project)  + (1|Family) + (1|Genus) + (1|Species), data.all)
  # technically, Genus and Species random effects should be nested within Family, but the nested model never converged.
  # tests on subsets yeild extremely similar variance estimates for nested and non-nested effects.

# create dataframes with variance parameter estimates
LMAvariance <- data.frame(VarCorr(logLMAvar))
LLvariance <- data.frame(VarCorr(logLLvar))
Nmassvariance <- data.frame(VarCorr(logNmassvar))
Nareavariance <- data.frame(VarCorr(logNareavar))

# combine all variance estimates, leaving out "Project"
traitvars <- data.frame(LMAvariance[which(LMAvariance$grp !="Project"),4], LLvariance[which(LLvariance$grp !="Project"),4], Nmassvariance[which(Nmassvariance$grp !="Project"),4], Nareavariance[which(Nareavariance$grp !="Project"),4])
colnames(traitvars) <- c("logLMA", "logLL", "logNmass", "logNarea")
rownames(traitvars) <- c("BtwSpecies", "BtwGenera", "BtwFamilies", "WtinSpecies")

traitvars_scaled <- traitvars  
for(i in 1:ncol(traitvars)){
  traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
}
# re-order rows
traitvars_scaled2 <- traitvars_scaled[c(3,2,1,4),]

##### raw traits
# LMAvar <- lmer(10^log.LMA~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
# LLvar <- lmer(10^log.LL~ Project + (1|Family) + (1|Genus) + (1|Species), data.all[-which(data.all$Project=="CO"),])
# Nmassvar <- lmer(10^log.Nmass~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
# Nareavar <- lmer(10^log.Narea~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
# #logNmassvar <- lmer(log.Nmass~ Project + (1|Family) + (1|Genus) + (1|Species), data.all)
# # in absolute scale (months), the variation appears to be HUGE w/ind spp, then w/in genera, then between Families
# # all levels of this analysis are hugely non-normal, except Families aren't too bad

# raw traits
# rLMAvariance <- data.frame(VarCorr(LMAvar))
# rLLvariance <- data.frame(VarCorr(LLvar))
# rNmassvariance <- data.frame(VarCorr(Nmassvar))
# rNareavariance <- data.frame(VarCorr(Nareavar))

# rtraitvars <- data.frame(rLMAvariance[,4], rLLvariance[,4], rNmassvariance[,4], rNareavariance[,4])
# colnames(rtraitvars) <- c("LMA", "LL", "Nmass", "Narea")
# rownames(rtraitvars) <- c("BtwSpecies", "BtwGenera", "BtwFamilies", "WtinSpecies")

# rtraitvars_scaled <- rtraitvars  
# for(i in 1:ncol(rtraitvars)){
#   rtraitvars_scaled[,i] <- rtraitvars[,i]/sum(rtraitvars[,i])
# }
# rtraitvars_scaled2 <- rtraitvars_scaled[c(3,2,1,4),]




################ . full PNW and evergreen needle leaf PFT ####################


# make a function to perform variance decomps of log and raw traits
variance.decomp <- function(dataz, fitfam=F){
  
  #### logged trait values
  if (fitfam ==F){
    SLAvar2 <- lmer(SLA_HSA~ (1|GENUS/SPECIES), dataz)
    LMAvar2 <- lmer(log.LMA~ (1|GENUS/SPECIES), dataz)
    LeafLifevar2 <- lmer(log.LL~ (1|GENUS/SPECIES), dataz)
    LeafCNvar2 <- lmer(LEAF_CN~ (1|GENUS/SPECIES), dataz)
    LeafCarbvar2 <- lmer(LEAF_CARBON~ (1|GENUS/SPECIES), dataz)
    LeafNitvar2 <- lmer(log.Nmass~ (1|GENUS/SPECIES), dataz)
    Nareavar2 <- lmer(log.Narea~ (1|GENUS/SPECIES), dataz)
    
    SLAvariance <- data.frame(VarCorr(SLAvar2))
    LMAvariance <- data.frame(VarCorr(LMAvar2))
    LeafLifevariance <- data.frame(VarCorr(LeafLifevar2))
    LeafCNvariance <- data.frame(VarCorr(LeafCNvar2))
    LeafCarvariance <- data.frame(VarCorr(LeafCarbvar2))
    LeafNitvariance <- data.frame(VarCorr(LeafNitvar2))
    Nareavariance <- data.frame(VarCorr(Nareavar2))
    
    traitvars <- data.frame(LMAvariance[,4], LeafLifevariance[,4], LeafNitvariance[,4], Nareavariance[,4], LeafCarvariance[,4], LeafCNvariance[,4])
    colnames(traitvars) <- c("log.LMA", "log.LL","log.Nmass","log.Narea","Cmass", "LEAF_CN")
    rownames(traitvars) <- c("BtwSpecies", "BtwGenus", "WtinSpecies")
    traitvars_scaled <- traitvars  
    for(i in 1:ncol(traitvars)){
      traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
    }
    # reorder appropriately so that it plots with w.inPlots on top and btw Genera on bottom.
    traitvars_scaled1 <- traitvars_scaled[c(2,1,3),]
    return(traitvars_scaled1)
  }
  
  if(fitfam==T){
    SLAvar2 <- lmer(SLA_HSA~ (1|Family/GENUS/SPECIES), dataz)
    LMAvar2 <- lmer(log.LMA~ (1|Family/GENUS/SPECIES), dataz)
    LeafLifevar2 <- lmer(log.LL~ (1|Family/GENUS/SPECIES), dataz)
    LeafCNvar2 <- lmer(LEAF_CN~ (1|Family/GENUS/SPECIES), dataz)
    LeafCarbvar2 <- lmer(LEAF_CARBON~ (1|Family/GENUS/SPECIES), dataz)
    LeafNitvar2 <- lmer(log.Nmass~ (1|Family/GENUS/SPECIES), dataz)
    Nareavar2 <- lmer(log.Narea~ (1|Family/GENUS/SPECIES), dataz)
    
    SLAvariance <- data.frame(VarCorr(SLAvar2))
    LMAvariance <- data.frame(VarCorr(LMAvar2))
    LeafLifevariance <- data.frame(VarCorr(LeafLifevar2))
    LeafCNvariance <- data.frame(VarCorr(LeafCNvar2))
    LeafCarvariance <- data.frame(VarCorr(LeafCarbvar2))
    LeafNitvariance <- data.frame(VarCorr(LeafNitvar2))
    Nareavariance <- data.frame(VarCorr(Nareavar2))
    
    traitvars <- data.frame(LMAvariance[,4], LeafLifevariance[,4], LeafNitvariance[,4], Nareavariance[,4], LeafCarvariance[,4], LeafCNvariance[,4])
    colnames(traitvars) <- c("log.LMA", "log.LL","log.Nmass","log.Narea","Cmass", "LEAF_CN")
    rownames(traitvars) <- c("BtwSpecies", "BtwGenus","BtwFam", "WtinSpecies")
    traitvars_scaled <- traitvars  
    for(i in 1:ncol(traitvars)){
      traitvars_scaled[,i] <- traitvars[,i]/sum(traitvars[,i])
    }
    # reorder appropriately so that it plots with w.inPlots on top and btw Genera on bottom.
    traitvars_scaled1 <- traitvars_scaled[c(3,2,1,4),]
    return(traitvars_scaled1)
  }
  
  
}


####### full PNW dataset

alltraitsvars <- variance.decomp(dataz=traits, fitfam = T) 
allPNWvars <- data.frame(alltraitsvars[,c(1:4)],"Type"=c(rep("all_PNW", times=nrow(alltraitsvars))))


####### Evergreen Needleleaf PFT
# all the spp codes for evergreen conifers
ENF <- c("ABIAMA", "ABICON", "ABIGRA", "ABILAS", "ABIMAG", "ABIPRO", "CALDEC", "JUNOCC", "PICENG",
         "PICSIT", "PINCON", "PINFLE", "PINJEF", "PINLAM", "PINMON", "PINPON", "PSEMEN", "THUPLI", "TSUHET", "TSUMER")
traits.conifers <- subset(traits, subset=SP.ID %in% ENF)
traits.conifers$SP.ID <- factor(traits.conifers$SP.ID)

# do variance decomposition
domconvars <- variance.decomp(dataz=traits.conifers)
PNWconvars <- data.frame(domconvars[,c(1:4)], "Type"=rep("ENF", times=nrow(domconvars)))
PNWconvars[4,] <- rep(NA, times=5)
PNWconvars[4,5] <- "ENF"
PNWconvars <- PNWconvars[c(4,1:3),]
row.names(PNWconvars)[1] <- "BtwFam" 












######## Table S2: Variance Decompositions ##############
globalvars <- traitvars_scaled2
row.names(globalvars) <- row.names(allPNWvars)
colnames(globalvars) <- colnames(allPNWvars)[-5]
globalvars$Type <- rep("global", times=nrow(globalvars))
vardecomps <- rbind(globalvars, allPNWvars, PNWconvars)



write.csv(vardecomps, file = paste0("./",results_dirname,"/VarianceDecompositionResults_GlobalPNWcons.csv")) # for Eco Let revision w/ Rainier











#________________________________________________________________
###### .FIG 1: Variance Decomp ######
#________________________________________________________________
# panel a) log.Narea vs log.LMA for diff subsets
# panel b) global var decomp (traitvars_scaled2)
# panel c) all PNW data (alltraitsvars.comb)
# panel d) evergreen needle-leaf PFT from PNW (domconvars1.comb)
# Ecol Let width options = 3.23, 4.33 or 6.81

# first need to make dataset with only species level averages, so well sampled species don't throw of visualization of global trait space
allspp <- data.all %>% group_by(Species, Genus, Family) %>% summarise( lslog.LL = mean(log.LL, na.rm=T), lslog.LMA = mean(log.LMA, na.rm=T), lslog.Nmass = mean(log.Nmass, na.rm=T), lslog.Narea = mean(log.Narea, na.rm=T)
                                                                       ,slog.LL = log(mean(10^log.LL, na.rm=T),base=10), slog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), slog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), slog.Narea = log(mean(10^log.Narea, na.rm=T),base=10)
                                                                       ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
colnames(allspp) <- gsub("slog", "log", colnames(allspp))


#jpeg(width=4.33, height=4.73, units = "in", res=600,filename = paste0("./",results_dirname,"/Fig1_VarDecomps.jpeg"))
pdf(width=4.33, height=4.73, file = paste0("./",results_dirname,"/Fig1_VarDecomps.pdf"))
#quartz(width=4.33, height=4.73)
par(mfrow=c(2,2), mgp=c(2,.7,0), cex.lab=1.1, cex.axis=1.1, mar=c(4.5,2,1.5,2),oma=c(0,2,3.8,0))


plot(log.Narea~log.LMA, allspp, col="grey", pch=16, xlab=expression(paste(log[10](LMA))))
mtext(expression(paste(log[10](N[area]))), side=2, line=2)
points(log.Narea~log.LMA, traits)
points(log.Narea~log.LMA, traits.conifers, pch=16, cex=.8, col=mypal[3])
points(biomass$log.cw_Nareap_if~biomass$log.cw_LMAp_if, bg=mypal[5], pch=24)
mtext("a)", side=3, adj=-.1, line=.3)
legend(x=1.1, y=2.5, xpd=NA,legend = c("Global","PNW woody plants", "Evgrn Needle PFT", "PNW CWMs"), pch=c(16,1,16,24), col=c("grey","black",mypal[3],"black"), pt.bg = mypal[5], bty="n")


#par(mgp=c(3,.7,0), cex.lab=1.1, cex.axis=1.1, mfrow=c(1,3), mar=c(5,2,5,2), oma=c(0,3,0,0))
cols <- rev(c(mypal[colchoices[c(1,2,3)]],"#CCCCCC"))#brewer.pal(11, "RdBu")[c(11,9,1, 6)]
# With old colors, alpha=CC, with new colors, alpha=99
bp <-barplot(as.matrix(traitvars_scaled2),beside=F,legend.text = F,xpd = T, names.arg = c("","","",""),las=2,args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"99")
             , ylab= "Proportion of total Variance" #ylab="Proportion of total Variance\n(log traits)"
             , xlab="", mgp=c(2.5,.8,0))
#legend(xpd=T, x = 0, y=1.2, legend=rev(c("btw Fams","w/in Fam","w/in Gen","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=2, bty="n",  cex=1.2)
text(x = bp, y= par("usr")[3]-.05,labels =  c(expression(paste(log[10](LMA))), expression(paste(log[10](LL))),expression(paste(log[10](N[mass]))),expression(paste(log[10](N[area])))), srt=40, adj=1,xpd=NA, cex=1.1, font=1)
mtext("b)", side=3, adj=-.1, line=0.3)
mtext("Prop of Total Variance", side=2, line=2.3)
mtext("Global", side=3, line=0.2)

legend(xpd=NA, x = 0, y=1.95, legend=rev(c("btw Fams","w/in Fam","w/in Gen","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=1, bty="n",  cex=1)



bp <-barplot(as.matrix(allPNWvars[,1:4]),beside=F,legend.text = F,xpd = T, names.arg = c("","","",""),las=2,args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols,"99")
             , ylab= "Proportion of total Variance" #ylab="Proportion of total Variance\n(log traits)"
             , xlab="", mgp=c(2.5,.8,0), las=2)
#legend(xpd=T, x = 0, y=1.2, legend=rev(c("btw Fams","w/in Fam","w/in Gen","w/in Spp")), fill=rev(paste0(cols,"99")), ncol=2, bty="n",  cex=1.2)
text(x = bp, y= par("usr")[3]-.05,labels =  c(expression(paste(log[10](LMA))), expression(paste(log[10](LL))),expression(paste(log[10](N[mass]))),expression(paste(log[10](N[area])))), srt=40, adj=1,xpd=NA, cex=1.1, font=1)
mtext("c)", side=3, adj=-.2, line=0.3)
mtext("Prop of Total Variance", side=2, line=2.3)
mtext("PNW woody plants", side=3, line=0.3)



barplot(as.matrix(PNWconvars[2:4,1:4]),beside=F,legend.text = F,xpd = T, names.arg = c("", "","",""),args.legend = list(x=4, y=1.3, ncol=2), col = paste0(cols[-1],"99"), ylab="", mgp=c(2,1,0))
text(x = bp, y= par("usr")[3]-.05,labels =  c(expression(paste(log[10](LMA))), expression(paste(log[10](LL))),expression(paste(log[10](N[mass]))),expression(paste(log[10](N[area])))), srt=40, adj=1,xpd=NA, cex=1.1, font=1)
mtext("d)", side=3, adj=-.2, line=.3)
mtext("Evgrn Needle PFT", side=3, line=0.3)

dev.off()





########## END: Variance Decomposition and Figure 1 ################







#__________________________________________________________________________________
################ BEGIN: Trait-trait scale dependence Analysis #####################
#__________________________________________________________________________________


########### . Create w/in spp, species means, genus means, and family means datasets ###############
  # then creating a dataset of genera w/ >5 species, and families w/ >5 genera
# select species w/ > 5 measurements
commonspp <-  names(which(xtabs(~Species, data.all)>=5))
spp.data<- data.all %>% filter(Species %in% commonspp)
spp.data$Species <- factor(spp.data$Species)
spp.data$Genus <- factor(spp.data$Genus)
spp.data$Family <- factor(spp.data$Family)


########## Species level trait averages ______________________________________________

allspp <- data.all %>% group_by(Species, Genus, Family) %>% summarise( lslog.LL = mean(log.LL, na.rm=T), lslog.LMA = mean(log.LMA, na.rm=T), lslog.Nmass = mean(log.Nmass, na.rm=T), lslog.Narea = mean(log.Narea, na.rm=T)
                                                                       ,slog.LL = log(mean(10^log.LL, na.rm=T),base=10), slog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), slog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), slog.Narea = log(mean(10^log.Narea, na.rm=T),base=10)
                                                                       ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
colnames(allspp) <- gsub("slog", "log", colnames(allspp))




#### Genus means ___________________________________________________________________
allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], lglog.LL = mean(llog.LL, na.rm=T),lglog.LMA = mean(llog.LMA, na.rm=T),lglog.Nmass = mean(llog.Nmass, na.rm=T), lglog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,glog.LL = log(mean(10^log.LL, na.rm=T),base=10), glog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), glog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), glog.Narea = log(mean(10^log.Narea, na.rm=T),base=10), nspp = n() 
                                                    ,MAT = mean(MAT, na.rm=T), MAP=mean(MAP, na.rm=T), VPD=mean(VPD, na.rm=T) )
# taking mean of raw values or logged values doesn't matter all that much yet
colnames(allgen) <- gsub("glog", "log", colnames(allgen))
colnames(allgen)[2] <- "Family"
# 939 genera from 211 families
# allgen <- allspp %>% group_by(Genus)  %>% summarise(Fam = sort(unique(Family))[1], Needle.Broad = unique(na.omit(NeedleBroad))[1], glog.LL = mean(log.LL, na.rm=T), glog.LMA = mean(log.LMA, na.rm=T), glog.Nmass = mean(log.Nmass, na.rm=T),glog.Narea = mean(log.Narea, na.rm=T)
#                                                     ,rglog.LL = log(mean(10^rlog.LL, na.rm=T),base=10), rglog.LMA = log(mean(10^rlog.LMA, na.rm=T),base=10), rglog.Nmass = log(mean(10^rlog.Nmass, na.rm=T),base=10), nspp = n() )



#### Family Means _______________________________________________________________________
allfam <- allgen %>% group_by(Family) %>% summarise(lflog.LL = mean(llog.LL, na.rm=T), lflog.LMA = mean(llog.LMA, na.rm=T), lflog.Nmass = mean(llog.Nmass, na.rm=T),lflog.Narea = mean(llog.Narea, na.rm=T)
                                                    ,flog.LL = log(mean(10^log.LL, na.rm=T),base=10), flog.LMA = log(mean(10^log.LMA, na.rm=T),base=10), flog.Nmass = log(mean(10^log.Nmass, na.rm=T),base=10), flog.Narea = log(mean(10^log.Narea, na.rm=T), base=10)
                                                    , tnspp = sum(nspp), ngen=n() )
colnames(allfam) <- gsub("flog", "log", colnames(allfam))

#### Family Means for fams w/ > 3 species
fam.data <- allfam # new: 211 families (up from 189)
fam.dataclean <- allfam[which(allfam$tnspp>2),] # 101 families, up from 97 families




##### dataset of only genera w/ >5 species _____________________________________________
# for w/in genus analysis
gen.data <- allspp[which(allspp$Genus %in% names(which(xtabs(~Genus, allspp)>=5))),]
# 750 measurements of 73 genera.
gen.data$Genus <- factor(gen.data$Genus) # get rid of unused levels


##### Spp. in Family level data, for families w/ >5 species in them ________________________
sppinfam.data <- allspp[which(allspp$Family %in% names(which(xtabs(~Family,allspp)>=5))),]
sppinfam.data$Family <- factor(sppinfam.data$Family) # get rid of unused levels



#### gen in Family level data, for families w/ >5 genera in them _________________________
# 556 measurements of 40 Families , added 2 Families from last batch, and probably quite a few obs (old n_LMA.N = 459)
# 04.01.17 - 684 obs of 50 Families. again, big win from taxonomic cleaning!
geninfam.data <- allgen[which(allgen$Family %in% names(which(xtabs(~Family,allgen)>=5))),]
geninfam.data$Family <- factor(geninfam.data$Family) # get rid of unused levels
geninfam.data$Genus <- factor(geninfam.data$Genus)










#_______________________________________________________________________

############ . **** RUN ANALYSIS: Calculate Trait-Trait Relationship & Null Models ################
#_______________________________________________________________________



#++++++++++++++++++++++ Do not run to bypass fitting SMA and null models ++++++++++++++++++++++++++++++
#  This can take a number of hours on a laptop. Takes less if you set n-iters to 10 or 100 but null model results can be sensitive
# Intstead jump to Line 1722: "BEGIN: Taxonomic Scale Statistical Tests"




#### Function for fitting MARs

fit.MAR <- function(xvar, yvar, data, method="SMA") {
  if(method=="SMA") meth = 3
  if(method=="MA") meth = 2
  if(method=="OLS") meth =1
  if(length(t(as.vector(data[!(is.na(data[,xvar]) | is.na(data[,yvar])), xvar])))<3){
    return(rep(NA, times=9))
  }
  else{
    if(var(data[,yvar], na.rm=T)==0 | var(data[,xvar],na.rm=T)==0){
      return(rep(NA, times=9))
    }
    else{
      tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), data))
      intercept <- tmp.mod$regression.results$Intercept[meth]
      slope <- tmp.mod$regression.results$Slope[meth]
      slope.lci <- tmp.mod$confidence.intervals[meth,4]
      slope.uci <- tmp.mod$confidence.intervals[meth,5]
      rho <- tmp.mod$r
      r.sq <- tmp.mod$rsquare
      n <- tmp.mod$n
      varX <- var(data[,xvar], na.rm=T)
      varY <- var(data[,yvar], na.rm=T)
      results <- c(intercept, slope,slope.lci, slope.uci, rho, r.sq, n, varX, varY)
      return(results)
    }
  }
}


# from Cd-NullModel.R
test.sig <- function(x, test){
  if(x<test[1] | x>test[6]) return(0.025)
  else
    if(x<test[2] | x>test[5]) return(0.05)
  else
    if(x<test[3] | x>test[4]) return(0.1)
  else return(1)
}
# from Cd-NullModel.R - fits a null model based on the mean and range of the data
fit.null <- function(xvar, yvar, observed, nulldata, nits){
  # find the trait ranges to sample
  rangeX <- range(observed[,xvar], na.rm=T)
  difX <- rangeX[2]-rangeX[1]
  #meanX <- mean(observed[,xvar], na.rm=T)
  rangeY <- range(observed[,yvar], na.rm=T)
  difY <- rangeY[2]-rangeY[1]
  #meanY <- mean(observed[,yvar], na.rm=T)
  # cut down the null data to just have non NAs
  nulldata <- nulldata[which(!is.na(nulldata[,xvar]) & !is.na(nulldata[,yvar])),]
  #nulldata2 <- nulldata %>% filter(!is.na(xvar) & !is.na(yvar))
  restrictednull <- nulldata[which(nulldata[,xvar]> min(nulldata[,xvar], na.rm=T)+difX/2 & nulldata[,xvar]< max(nulldata[,xvar], na.rm=T)-difX/2 & nulldata[,yvar]> min(nulldata[,yvar], na.rm=T)+difY/2 & nulldata[,yvar]< max(nulldata[,yvar], na.rm=T)-difY/2),] 
  if(nrow(restrictednull)<3){restrictednull <- nulldata}
  nullcor <- c(rep(NA, times=nits))
  for(i in 1:nits){
    center <- restrictednull[sample(nrow(restrictednull),size = 1),]
    #null <- nulldata %>% filter(xvar > center[,xvar]-rangeX/2) #, xvar < center[,xvar]+rangeX/2, yvar > center[,yvar]-rangeY/2, yvar < center[,yvar]-rangeY/2)
    null <- nulldata[which(nulldata[,xvar] > as.numeric(center[,xvar])-difX/2 & nulldata[,xvar] < as.numeric(center[,xvar]+difX/2) & 
                             nulldata[,yvar] > as.numeric(center[,yvar])-difY/2 & nulldata[,yvar] < as.numeric(center[,yvar]+difY/2)),]
    if(nrow(null)<5){# var(null[,xvar])==0 | var(null[,yvar])==0){ # this is crappy as hell, but to get rid of bad centers I'm just going to try repicking them and hope the probability of getting two bad points in a row is low...
      center <- restrictednull[sample(nrow(restrictednull),size = 1),]
      #null <- nulldata %>% filter(xvar > center[,xvar]-rangeX/2) #, xvar < center[,xvar]+rangeX/2, yvar > center[,yvar]-rangeY/2, yvar < center[,yvar]-rangeY/2)
      null <- nulldata[which(nulldata[,xvar] > as.numeric(center[,xvar])-difX/2 & nulldata[,xvar] < as.numeric(center[,xvar]+difX/2) & 
                               nulldata[,yvar] > as.numeric(center[,yvar])-difY/2 & nulldata[,yvar] < as.numeric(center[,yvar]+difY/2)),]
    }
    ndist <- null[sample(nrow(null), size = nrow(observed), replace = TRUE),]
    nullcor[i] <- cor(ndist[,c(xvar,yvar)])[2,1]
  }
  nullquantiles <- quantile (nullcor, probs=c(0.025,0.05,0.1,.9,.95,.975), na.rm=T)
  names(nullquantiles) <- c("lci_2.5", "lci_5","lci_10","uci_10","uci_5","uci_2.5")
  return(nullquantiles)
}

# calculate significance of a pearson correlation coefficient
rho.sig <- function(rho, n){
  t <- rho/sqrt((1-rho^2)/(n-2))
  pt(-abs(t), n-2)*2
}




############ . Fitting SMA slopes and null models ##############################
niters <- 10000 # number of null model iterations to use
niters.spinfam <- 10 # number of null model iterations to use for species-in-family analysis (not presented in ms) 


############# + LMA vs Nmass ###########

###_______________ Species level analysis _______________________________________________


ptm <- proc.time()
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}
proc.time()-ptm





# Seems to be working!!!

###_______________ Species w/in Genus level analysis _______________________________________________
t0 <- proc.time()
gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
  
}
proc.time()-t0


###_______________ Species w./in Fam level analysis  (not in ms)_______________________________________________
t0<- proc.time()
sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = niters.spinfam)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
  
}
proc.time()-t0

#_______________ Genus w/in Family level analysis _______________________________________________
t0<-proc.time()
geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}
proc.time()-t0

###_______________ Family level analysis _______________________________________________

# # currently just working with LES until I combine the PACNW dataset into this.
# famcorr.all <- cor(x=fam.data[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")
# famcorr.clean <- cor(x=fam.dataclean[,c("log.LL","log.LMA","log.Nmass")], use = "pairwise.complete.obs")

fam.res_LMA.N <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.data),rep(NA, times=7),"Fam") 
names(fam.res_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig" ,"Type")
fam.resclean_LMA.N <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.resclean_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.N <- c("global", fit.MAR(xvar='log.LMA',yvar="log.Nmass",data=allspp),rep(NA, times=7), "global")
names(global_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")


###_______________ Combining into 1 df _______________________________________________
# first add a "Type" column
spp.results$Type <- rep("w/inSpp", times=nrow(spp.results))
gen.results$Type <- rep("w/inGen", times=nrow(gen.results))
sppinfam.results$Type <- rep("Sppw/inFam", times=nrow(sppinfam.results))
geninfam.results$Type <- rep("Genw/inFam", times=nrow(geninfam.results))
# now make the column names all match
colnames(spp.results)[1] <- "Taxo.Unit"
colnames(gen.results)[1] <- "Taxo.Unit"
colnames(sppinfam.results)[1] <- "Taxo.Unit"
colnames(geninfam.results)[1] <- "Taxo.Unit"
all.results.LMAN <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LMA.N, fam.resclean_LMA.N, global_LMA.N)
all.results.LMAN$Type <- factor(all.results.LMAN$Type)
levels(all.results.LMAN$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")

write.csv(all.results.LMAN, file = paste0("./",dirname,"/SMA_Results_LMAvNmass.csv"))





####### + LMA and LL #####################

###_______________ Species level analysis _______________________________________________

set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


###_______________ Species w/in Genus level analysis _______________________________________________


gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



###_______________ Spp w/in Family level analysis (not in ms) _______________________________________________


sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allfam, nits = niters)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}


###_______________ Genera w/in Family level analysis _______________________________________________


geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.LL",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.LL", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}


###_______________ Family level analysis _______________________________________________

fam.res_LMA.LL <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.LL",data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LMA.LL) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LMA.LL <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.LL",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.resclean_LMA.LL) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.LL <- c("global", fit.MAR(xvar='log.LMA',yvar="log.LL",data=allspp),rep(NA, times=7), "global")
names(global_LMA.LL) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")



###_______________ Combining all into one df _______________________________________________

# first add a "Type" column
spp.results$Type <- rep("w/inSpp", times=nrow(spp.results))
gen.results$Type <- rep("w/inGen", times=nrow(gen.results))
sppinfam.results$Type <- rep("Sppw/inFam", times=nrow(sppinfam.results))
geninfam.results$Type <- rep("Genw/inFam", times=nrow(geninfam.results))
# now make the column names all match
colnames(spp.results)[1] <- "Taxo.Unit"
colnames(gen.results)[1] <- "Taxo.Unit"
colnames(sppinfam.results)[1] <- "Taxo.Unit"
colnames(geninfam.results)[1] <- "Taxo.Unit"
all.results.LMALL <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LMA.LL, fam.resclean_LMA.LL, global_LMA.LL)
all.results.LMALL$Type <- factor(all.results.LMALL$Type)
levels(all.results.LMALL$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")


write.csv(all.results.LMALL, file = paste0("./",dirname,"/SMA_Results_LMAvLL.csv"))





####### + LL and Nmass #####################

###_______________ Species level analysis _______________________________________________
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


###_______________ Species w/in Genus level analysis _______________________________________________

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



###_______________ Species w/in Family level analysis (not in ms) _______________________________________________

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = niters)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}





###_______________ Genera w.in Family level analysis _______________________________________________
geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Nmass",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Nmass", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}





###_______________ Family level analysis _______________________________________________


fam.res_LL.N <- c("fam.all", fit.MAR(xvar='log.LL',yvar="log.Nmass",data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LL.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LL.N <- c("fam.clean", fit.MAR(xvar='log.LL',yvar="log.Nmass",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.res_LL.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LL.N <- c("global", fit.MAR(xvar='log.LL',yvar="log.Nmass",data=allspp),rep(NA, times=7), "global")
names(global_LL.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")





###_______________ combine into one df  _______________________________________________

# first add a "Type" column
spp.results$Type <- rep("w/inSpp", times=nrow(spp.results))
gen.results$Type <- rep("w/inGen", times=nrow(gen.results))
sppinfam.results$Type <- rep("Sppw/inFam", times=nrow(sppinfam.results))
geninfam.results$Type <- rep("Genw/inFam", times=nrow(geninfam.results))
# now make the column names all match
colnames(spp.results)[1] <- "Taxo.Unit"
colnames(gen.results)[1] <- "Taxo.Unit"
colnames(sppinfam.results)[1] <- "Taxo.Unit"
colnames(geninfam.results)[1] <- "Taxo.Unit"
all.results.LLNmass <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LL.N, fam.resclean_LL.N, global_LL.N)
all.results.LLNmass$Type <- factor(all.results.LLNmass$Type)
levels(all.results.LLNmass$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")

write.csv(all.results.LLNmass, file = paste0("./",dirname,"/SMA_Results_LLvNmass.csv"))







####### + LMA and Narea  #####################

###_______________ Species level analysis _______________________________________________
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


###_______________ Species w/in Genus level analysis _______________________________________________

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



###_______________ Species w/in Family level analysis (not in ms)_______________________________________________

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = niters)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}





###_______________ Genera w/in Family level analysis _______________________________________________

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LMA',yvar="log.Narea",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LMA', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}





###_______________ Family level analysis _______________________________________________


fam.res_LMA.Narea <- c("fam.all", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LMA.Narea <- c("fam.clean", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.res_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LMA.Narea <- c("global", fit.MAR(xvar='log.LMA',yvar="log.Narea",data=allspp),rep(NA, times=7), "global")
names(global_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")




###_______________ Combining into one df _______________________________________________

# first add a "Type" column
spp.results$Type <- rep("w/inSpp", times=nrow(spp.results))
gen.results$Type <- rep("w/inGen", times=nrow(gen.results))
sppinfam.results$Type <- rep("Sppw/inFam", times=nrow(sppinfam.results))
geninfam.results$Type <- rep("Genw/inFam", times=nrow(geninfam.results))
# now make the column names all match
colnames(spp.results)[1] <- "Taxo.Unit"
colnames(gen.results)[1] <- "Taxo.Unit"
colnames(sppinfam.results)[1] <- "Taxo.Unit"
colnames(geninfam.results)[1] <- "Taxo.Unit"
all.results.LMANarea <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LMA.Narea, fam.resclean_LMA.Narea, global_LMA.Narea)
all.results.LMANarea$Type <- factor(all.results.LMANarea$Type)
levels(all.results.LMANarea$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")


write.csv(all.results.LMANarea, file = paste0("../",dirname,"/SMA_Results_LMAvNarea.csv"))






####### + LL and Narea #####################

###_______________ Species level analysis _______________________________________________
set.seed(42)
spp.results <- data.frame(matrix(NA, nrow=length(unique(spp.data$Species)), ncol=17))
colnames(spp.results) <- c("Species", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(spp.data$Species))){
  species <- levels(spp.data$Species)[i]
  print(species)
  dataz <- spp.data[which(spp.data$Species==species),]
  res <- fit.MAR(xvar="log.LL", yvar='log.Narea',data=dataz)
  spp.results[i,1] <- species
  spp.results[i,2:10] <- res
  if (!is.na(res[1]) &res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allspp, nits = niters)  
    spp.results[i, 11:16] <- nullbounds
    spp.results[i,17] <- test.sig(x=spp.results$Rho[i], test=nullbounds)
  }
}


###_______________ Species w/in Genus level analysis _______________________________________________

gen.results <- data.frame(matrix(NA, nrow=length(unique(gen.data$Genus)), ncol=17))
colnames(gen.results) <- c("Genus", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(gen.data$Genus))){
  genus <- levels(gen.data$Genus)[i]
  print(genus)
  dataz <- gen.data[which(gen.data$Genus==genus),]
  gen.results[i,1] <- genus
  res <- fit.MAR(xvar="log.LL",yvar='log.Narea',data=dataz)
  gen.results[i,2:10] <- res 
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points, and the fit.MAR actually ran
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allgen, nits = niters)  
    gen.results[i, 11:16] <- nullbounds
    gen.results[i,17] <- test.sig(x=gen.results$Rho[i], test=nullbounds)
  }
}



###_______________ Species w/in Family level analysis (not in ms) _______________________________________________

sppinfam.results <- data.frame(matrix(NA, nrow=length(unique(sppinfam.data$Family)), ncol=17))
colnames(sppinfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(sppinfam.data$Family))){
  family <- levels(sppinfam.data$Family)[i]
  print(family)
  dataz <- sppinfam.data[which(sppinfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Narea",data=dataz)
  sppinfam.results[i,1] <- family
  sppinfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = niters)  
    sppinfam.results[i, 11:16] <- nullbounds
    sppinfam.results[i,17] <- test.sig(x=sppinfam.results$Rho[i], test=nullbounds)
  }
}





###_______________ Genera w/in Family level analysis _______________________________________________

geninfam.results <- data.frame(matrix(NA, nrow=length(unique(geninfam.data$Family)), ncol=17))
colnames(geninfam.results) <- c("Family", "Int","Slope","Slope.lci","Slope.uci","Rho","r.sq","n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig")
for(i in 1:length(unique(geninfam.data$Family))){
  family <- levels(geninfam.data$Family)[i]
  print(family)
  dataz <- geninfam.data[which(geninfam.data$Family==family),]
  res <- fit.MAR(xvar='log.LL',yvar="log.Narea",data=dataz)
  geninfam.results[i,1] <- family
  geninfam.results[i,2:10] <- res
  if (!is.na(res[1]) & res[7]>5){ # only fit null model if there are >5 data points
    nullbounds <- fit.null(xvar='log.LL', yvar="log.Narea", observed = dataz, nulldata = allfam, nits = niters)  
    geninfam.results[i, 11:16] <- nullbounds
    geninfam.results[i,17] <- test.sig(x=geninfam.results$Rho[i], test=nullbounds)
  }
}





###_______________ Family level analysis _______________________________________________


fam.res_LL.Narea <- c("fam.all", fit.MAR(xvar="log.LL",yvar='log.Narea',data=fam.data),rep(NA, times=7),"Fam")
names(fam.res_LL.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
fam.resclean_LL.Narea <- c("fam.clean", fit.MAR(xvar="log.LL",yvar='log.Narea',data=fam.dataclean),rep(NA, times=7), "Fam.clean")
names(fam.res_LL.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")
global_LL.Narea <- c("global", fit.MAR(xvar="log.LL",yvar='log.Narea',data=allspp),rep(NA, times=7), "global")
names(global_LL.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")





###_______________ Combine into one df _______________________________________________

# first add a "Type" column
spp.results$Type <- rep("w/inSpp", times=nrow(spp.results))
gen.results$Type <- rep("w/inGen", times=nrow(gen.results))
sppinfam.results$Type <- rep("Sppw/inFam", times=nrow(sppinfam.results))
geninfam.results$Type <- rep("Genw/inFam", times=nrow(geninfam.results))
# now make the column names all match
colnames(spp.results)[1] <- "Taxo.Unit"
colnames(gen.results)[1] <- "Taxo.Unit"
colnames(sppinfam.results)[1] <- "Taxo.Unit"
colnames(geninfam.results)[1] <- "Taxo.Unit"
all.results.LLNarea <-rbind(spp.results,gen.results, sppinfam.results, geninfam.results, fam.res_LL.Narea, fam.resclean_LL.Narea, global_LL.Narea)
all.results.LLNarea$Type <- factor(all.results.LLNarea$Type)
levels(all.results.LLNarea$Type) <- list(w.inSpp = "w/inSpp",   w.inGen= "w/inGen",    Sppw.inFam="Sppw/inFam", Genw.inFam="Genw/inFam", Fam = "Fam", Famclean = "Fam.clean", global="global")


write.csv(all.results.LLNarea, file = paste0("./",dirname,"/SMA_Results_LLvNarea.csv"))





############ + CMW analysis #####################


CWM_LMA.LL <- c("CWM", fit.MAR(xvar='log.cw_LMAp_if',yvar="log.cw_LLp_if",data=biomass),rep(NA, times=7), "CWM")
names(CWM_LMA.LL) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varLL","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")


CWM_LMA.N <- c("CWM", fit.MAR(xvar='log.cw_LMAp_if',yvar="log.cw_Nmassp_if",data=biomass),rep(NA, times=7), "CWM")
names(CWM_LMA.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNmass", "lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")

CWM_LL.N <- c("CWM", fit.MAR(xvar='log.cw_LLp_if',yvar="log.cw_Nmassp_if",data=biomass),rep(NA, times=7), "CWM")
names(CWM_LL.N) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNmass","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")

CWM_LMA.Narea <- c("CWM", fit.MAR(xvar='log.cw_LMAp_if',yvar="log.cw_Nareap_if",data=biomass),rep(NA, times=7), "CWM")
names(CWM_LMA.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLMA","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")

CWM_LL.Narea <- c("CWM", fit.MAR(xvar="log.cw_LLp_if",yvar='log.cw_Nareap_if',data=biomass),rep(NA, times=7), "CWM")
names(CWM_LL.Narea) <- c("Taxo.Unit","Int","Slope","Slope.lci","Slope.uci","Rho","r.sq", "n","varLL","varNarea","lci_2.5","lci_5","lci_10","uci_10","uci_5","uci_2.5", "sig","Type")












########### + Merging All Results ###########################

LMALL <- all.results.LMALL
colnames(LMALL)[c(2:8, 11:17)] <- paste(colnames(all.results.LMALL)[c(2:8, 11:17)],"LMA.LL", sep="_")
LMAN <- all.results.LMAN
colnames(LMAN)[c(2:8, 11:17)] <- paste(colnames(all.results.LMAN)[c(2:8, 11:17)],"LMA.N", sep="_")
LLNmass <- all.results.LLNmass
colnames(LLNmass)[c(2:8, 11:17)] <- paste(colnames(all.results.LLNmass)[c(2:8, 11:17)],"LL.N", sep="_")
LMANarea <- all.results.LMANarea
colnames(LMANarea)[c(2:8, 11:17)] <- paste(colnames(all.results.LMANarea)[c(2:8, 11:17)],"LMA.Narea", sep="_")
LLNarea <- all.results.LLNarea
colnames(LLNarea)[c(2:8, 11:17)] <- paste(colnames(all.results.LLNarea)[c(2:8, 11:17)],"LL.Narea", sep="_")

all.results <- cbind(LMALL, LMAN[,-c(1,9,18)], LLNmass[,-c(1,9,10,18)], LMANarea[,-c(1,9,18)], LLNarea[,-c(1,9,10,18)]) # drop the duplicate 'var' columns


### add in CWMs to dataframe post hoc
allcwms <- c(CWM_LMA.LL, CWM_LMA.N[-c(1,9,18)], CWM_LL.N[-c(1,9,10,18)], CWM_LMA.Narea[-c(1,9,18)], CWM_LL.Narea[-c(1,9,10,18)])
all.results.cwm <- all.results
all.results.cwm$Taxo.Unit <- as.character(all.results.cwm$Taxo.Unit)
all.results$Type <- as.character(all.results.cwm$Type)
  # add in a final row with CWM values 
all.results.cwm[nrow(all.results.cwm)+1,] <- allcwms



write.csv(all.results.cwm, paste0("./",results_dirname,"/SMA_Results_All_CWM.csv"))


# test <- read.csv("./20171220_results/SMA_Results_All_CWM.csv", header=T)
# all.resultsold <- read.csv("Results_SimpleMAreg_v11rawavgs_20171211_wCWM.csv")



#### a number of columns save as character but need to be numeric. easiest way is just to reload the created .csv
all.results.cwm <- read.csv(paste0(results_dirname,"/SMA_Results_All_CWM.csv"), row.names = 1)
levels(all.results.cwm$Type) <- list(w.inSpp = "w.inSpp", w.inGen = "w.inGen", Sppw.inFam= "Sppw.inFam",Genw.inFam="Genw.inFam", Fam="Fam",Famclean="Famclean", global="global", CWM="CWM")

# filter out the levels I don't care about
all.results.cwm.cl <- all.results.cwm %>% filter(Type %in% c("w.inSpp","w.inGen","Genw.inFam","Famclean","global", "CWM"))
all.results.cwm.cl$Type <- factor(all.results.cwm.cl$Type)
all.results.cwm.cl$Taxo.Unit <- factor(all.results.cwm.cl$Taxo.Unit)
# add columns of whether the correlation is significant (useful for plotting)
all.results.cwm.cl$rho.sig_LMA.N <- rho.sig(all.results.cwm.cl$Rho_LMA.N, all.results.cwm.cl$n_LMA.N)
all.results.cwm.cl$rho.sig_LMA.Narea <- rho.sig(all.results.cwm.cl$Rho_LMA.Narea, all.results.cwm.cl$n_LMA.Narea)
all.results.cwm.cl$rho.sig_LMA.LL <- rho.sig(all.results.cwm.cl$Rho_LMA.LL, all.results.cwm.cl$n_LMA.LL)
all.results.cwm.cl$rho.sig_LL.N <- rho.sig(all.results.cwm.cl$Rho_LL.N, all.results.cwm.cl$n_LL.N)
all.results.cwm.cl$rho.sig_LL.Narea <- rho.sig(all.results.cwm.cl$Rho_LL.Narea, all.results.cwm.cl$n_LL.Narea)

# note: intermediate steps are saved as SMA_Results_type...csv in case of computer crashes. Final version is "SlopesAndCorrelations_Alltaxa_clean.csv"
write.csv(all.results.cwm.cl, paste0("./",results_dirname,"/SlopesAndCorrelations_Alltaxa_clean.csv"))


##### . **** END ANALYSIS: Calculate Trait-Trait Relationship & Null Models #############################
# (skip to here if you wish to avoid the time-consuming calculations above)






#____________________________________________________________________
###### BEGIN: Taxonomic Scale Statistical Tests #####################
#____________________________________________________________________
# tests to see whether correlations or SMA slopes differ significantly across taxonomic scales
# using multiple weighting methods to make sure inferences are robust
# results are in Table S3

#_________________________________________________________________________
# **** if by-passing the fitting of trait-trait relationships, run this:

 all.results.cwm.cl <- read.csv("DerivedData/SlopesAndCorrelations_Alltaxa_clean.csv", header=T, row.name=1)
 # levels(all.results.cwm.cl$Type) <- list(w.inSpp = "w.inSpp", w.inGen = "w.inGen",Genw.inFam="Genw.inFam",Famclean="Famclean", global="global", CWM="CWM")

#_________________________________________________________________________



####### Table S3: Summary of significances with different weighting methods

# create an empty df to hold results
tableS3 <- data.frame(Comparison = rep(c("LMA v Nmass","LL v Nmass", "LMA v LL", "LMA v Narea","LL v Narea"), each=2), Type = rep(c("SMA slope","Correlation"), times=5), unweighted = NA, sample.size=NA, V1.variance=NA, V2.variance=NA)

# create function to fit linear model with various weighting schemes and report ANOVA results
run.stats <- function ( trait1 = "LMA",trait2="LL", variable="Slope"){
  output <- c(rep(NA, times=4))
  mod1 <- lm(get(paste(variable,"_",trait1,".",trait2, sep=""))~Type, all.results.cwm.cl[which(all.results.cwm.cl[,paste("n","_",trait1,".",trait2, sep="")]>5& !all.results.cwm.cl$Type %in% c("CWM","global")),])
  output[1] <- round(anova(mod1)[[5]][1], 3)
  
  # sample size weighting
  mod2 <- lm(get(paste(variable,"_",trait1,".",trait2, sep=""))~Type, all.results.cwm.cl[which(all.results.cwm.cl[,paste("n","_",trait1,".",trait2, sep="")]>5& !all.results.cwm.cl$Type %in% c("CWM","global")),], weights = get(paste("n","_",trait1,".",trait2, sep="")))
  output[2] <- round(anova(mod2)[[5]][1], 3)
  
  # Var trait1 weighting
  mod3 <- lm(get(paste(variable,"_",trait1,".",trait2, sep=""))~Type, all.results.cwm.cl[which(all.results.cwm.cl[,paste("n","_",trait1,".",trait2, sep="")]>5& !all.results.cwm.cl$Type %in% c("CWM","global")),], weights = get(paste("var",trait1, sep="")))
  output[3] <- round(anova(mod3)[[5]][1], 3)
  
  # Var trait2 weighting
    # note: because I named Nmass just N for most column names except varNmass, I have to add this caveat 
  if(nchar(trait2)<2){
    mod4 <- lm(get(paste(variable,"_",trait1,".",trait2, sep=""))~Type, all.results.cwm.cl[which(all.results.cwm.cl[,paste("n","_",trait1,".",trait2, sep="")]>5& !all.results.cwm.cl$Type %in% c("CWM","global")),], weights = get(paste("var","Nmass", sep="")))
  }
  else{
  mod4 <- lm(get(paste(variable,"_",trait1,".",trait2, sep=""))~Type, all.results.cwm.cl[which(all.results.cwm.cl[,paste("n","_",trait1,".",trait2, sep="")]>5& !all.results.cwm.cl$Type %in% c("CWM","global")),], weights = get(paste("var",trait2, sep="")))
  }
  output[4] <- round(anova(mod4)[[5]][1], 3)
  return(output)
  
}

tableS3[1,3:6] <- run.stats(trait1 = "LMA",trait2="N",variable = "Slope")
tableS3[2,3:6] <- run.stats(trait1 = "LMA",trait2="N",variable = "Rho")

tableS3[3,3:6] <- run.stats(trait1 = "LL",trait2="N",variable = "Slope")
tableS3[4,3:6] <- run.stats(trait1 = "LL",trait2="N",variable = "Rho")

tableS3[5,3:6] <- run.stats(trait1 = "LMA",trait2="LL",variable = "Slope")
tableS3[6,3:6] <- run.stats(trait1 = "LMA",trait2="LL",variable = "Rho")

tableS3[7,3:6] <- run.stats(trait1 = "LMA",trait2="Narea",variable = "Slope")
tableS3[8,3:6] <- run.stats(trait1 = "LMA",trait2="Narea",variable = "Rho")

tableS3[9,3:6] <- run.stats(trait1 = "LL",trait2="Narea",variable = "Slope")
tableS3[10,3:6] <- run.stats(trait1 = "LL",trait2="Narea",variable = "Rho")

write.csv(tableS3, file = paste0("./",results_dirname,"/TableS3_StatsTaxonomicDifferences.csv"))





#_____________________________________________________________________
######## BEGIN: Figure 2-4 Plotting ##########
#_____________________________________________________________________

# function that plots the major axis regression line for a dataset
plot.MAR <- function(xvar, yvar, data, method="SMA", linecol, lwd=1, lty=1) {
  if(method=="SMA") meth = 3
  if(method=="MA") meth = 2
  if(method=="OLS") meth =1
  if(length(t(as.vector(data[!(is.na(data[,xvar]) | is.na(data[,yvar])), xvar])))<3){
    #return(rep(NA, times=7))
    break()
  }
  else{
    if(var(data[,yvar], na.rm=T)==0){
      break()
    }
    else{
      tmp.mod <- suppressMessages(lmodel2(get(yvar)~get(xvar), data))
      intercept <- tmp.mod$regression.results$Intercept[meth]
      slope <- tmp.mod$regression.results$Slope[meth]
      xvals <- sort(tmp.mod$x)
      yhat <- xvals * slope + intercept
      lines(yhat~xvals, col=linecol, lwd=lwd, lty=lty)
      
    }
  }
}







####.FIG 2ab: LMA v Nmass #########
colchoices <- c(1,2,4,3,6)
cex.sig <- 1.1 # size of points for significant slopes and correlations
cex.ns <- .9 # size of non-sig points
palette(mypal[colchoices])
crit <- 0.05 # alpha for null model significance (in all.results.cwm.cl, $sig_XX.XX = crit means that significance is < than crit, i.e. test is significant at the alpha level reported in the $sig_XX.XX column)


mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)


#jpeg(width=7.008, height=4, units = "in", res=600,filename = paste0("./",results_dirname,"/Fig2ab_LMA_Nmass2.jpeg"))
pdf(width=7.008, height=4, file = paste0("./",results_dirname,"/Fig2ab_LMA_Nmass2.pdf"))
#quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1))

# Slope boxplots
# Slope boxplots
p <- boxplot(Slope_LMA.N~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.N>5& !all.results.cwm.cl$Type %in% c("global","Famclean","CWM")),]
             , ylim=c(-2.2,1.5),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n", plot=F)
boxplot(Slope_LMA.N~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.N>5& !all.results.cwm.cl$Type %in% c("global","Famclean","CWM")&is.na(all.results.cwm.cl$Taxo.Unit) ),]
        , ylim=c(-2.7,1.7),las=3, ylab="SMA Slope" #, main="log(LMA)~log(Nmass) strict"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n")
# add global and bewteen family points and error bars
m <- lmodel2(log.Nmass~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.Nmass~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
m <- lmodel2(log.cw_Nmassp_if~log.cw_LMAp_if, biomass)
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
points(y=m$regression.results$Slope[3],x=6, pch=17, cex=1.3, col=mypal[5])

abline(h=0, lty=2)
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= expression(paste(N[mass]," vs LMA", sep=" ")), side=3, line=.2)

### add in individual points with error bars for significant correlations
tmp <- all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.N>5 & !is.na(all.results.cwm.cl$sig_LMA.N)),]
tmp$rho.sig <- rho.sig(rho=tmp$Rho_LMA.N, n= tmp$n_LMA.N)
tmp.ns <- tmp[which(tmp$rho.sig>0.05),]
tmp.sig <- tmp[which(tmp$rho.sig<=0.05),]
#palette(paste0(mypal[colchoices], "55"))
set.seed(42)
points(Slope_LMA.N~jitter(as.numeric(Type)), tmp.ns, pch=1, col=Type, cex=cex.ns)
#palette(mypal[colchoices])
set.seed(52)
points(Slope_LMA.N~jitter(as.numeric(Type)), tmp.sig, pch=16, col=Type, cex=cex.sig)
legend('topright',legend = c("sig from 0","ns from 0"), pch=c(16,1),cex=.9, pt.cex  =c(cex.sig, cex.ns), col= mypal[colchoices[1]], bty="n")





# Rho boxplots
par(mar=c(6,4,0,1))
boxplot(Rho_LMA.N~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.N>5 & !all.results.cwm.cl$Type %in% c("global","Famclean", "CWM") & is.na(all.results.cwm.cl$Taxo.Unit) ),]
        , ylim=c(-1,1.1),las=3, ylab="Correlation"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n")
# add global and between family points
m <- cor.test(x=allspp$log.LMA, y=allspp$log.Nmass)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.Nmass)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
m <- cor.test(x=biomass$log.cw_LMAp_if, y=biomass$log.cw_Nmassp_if)
arrows(x0=6,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col=mypal[5])
points(y=m$estimate,x=6, pch=17, cex=1.3, col=mypal[5])


# layering points on top of null model
set.seed(42)
points(Rho_LMA.N~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.N>5 & all.results.cwm.cl$sig_LMA.N<=crit),], pch=16, col=Type, cex=cex.sig)
#palette(paste0(mypal[colchoices], "55"))
set.seed(42)
points(Rho_LMA.N~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.N>5 & all.results.cwm.cl$sig_LMA.N>crit),], pch=1, col=Type, cex=cex.ns)
#axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","Global"), at=c(1,2,3,4,5), las=3)
#palette(mypal[colchoies])

# add axes and sample sizes
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","Global","PNW CWM"), at=c(1,2,3,4,5,6), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2,x=c(4,5,6), labels=paste0("(",all.results.cwm.cl$n_LMA.N[c((nrow(all.results.cwm.cl)-2),nrow(all.results.cwm.cl)-1,nrow(all.results.cwm.cl))],")"), cex=0.9)
legend('topright',title = "",legend = c("sig from null","ns from null"), pch=c(16,1),cex=.9, pt.cex  =c(cex.sig, cex.ns), col= mypal[colchoices[1]], bty="n")



# scatterplot w/ lines
par(mar=c(4,4,1.5,1.5), cex.lab=1.2, mgp=c(2.5,.7,0))
ltysig <- 1
ltyns <-  2
lwdsig <-  1.8
lwdns <- 1

plot(log.Nmass~log.LMA, allspp, col="grey", pch=16, xlab=expression(paste(log[10],"(LMA)")), ylab=expression(paste(log[10],(N[mass]), sep=" ")))
tax <- "w.inGen"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.N>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.N[which(all.results.cwm.cl$Taxo.Unit==i)]<=0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lty=linety, lwd=linewd)
  
}
tax <- "Genw.inFam"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.N>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.N[which(all.results.cwm.cl$Taxo.Unit==i)]<=0.05) {linety <-ltysig ; linewd=lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lty=linety, lwd=linewd)
  
}
# across family relationship (fams w/ >3 spp)
abline(a=all.results.cwm.cl$Int_LMA.N[which(all.results.cwm.cl$Type=="Famclean")], b=all.results.cwm.cl$Slope_LMA.N[which(all.results.cwm.cl$Type=="Famclean")], lwd=3, col="black")
# global relationship across all species means
abline(a=all.results.cwm.cl$Int_LMA.N[which(all.results.cwm.cl$Type=="global")], b=all.results.cwm.cl$Slope_LMA.N[which(all.results.cwm.cl$Type=="global")], lwd=3, col="black", lty=2)
tax <- "w.inSpp"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.N>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.N[which(all.results.cwm.cl$Taxo.Unit==i)]<=0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.Nmass",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lty =linety, lwd=linewd)
}

mtext(text = "b)", side = 3, adj=0, line=.2)

#points(biomass$log.cw_Nmassp_if~biomass$log.cw_LMAp_if, pch=24, bg=mypal[5])
plot.MAR(xvar = "log.cw_LMAp_if", yvar="log.cw_Nmassp_if", data=biomass, linecol = mypal[5], lwd = 3)
legend('topright',legend = c("sig from 0","ns from 0"), lty=c(ltysig, ltyns),cex=.9, lwd  =c(lwdsig, lwdns), col= mypal[colchoices[1]], bty="n")
dev.off()








####### .FIG 2bc: LL v Nmass #########################
mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
cex.sig <- 1.1
cex.ns <- .9
palette(mypal[colchoices])

#jpeg(width=7.008, height=4, units = "in", res=600,filename = paste0("./",results_dirname,"/Fig2cd_LL_Nmass.jpeg"))
pdf(width=7.008, height=4, file = paste0("./",results_dirname,"/Fig2cd_LL_Nmass2.pdf"))

#quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1))
p <- boxplot(Slope_LL.N~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LL.N>5& !all.results.cwm.cl$Type %in% c("global","Famclean") ),]
             , ylim=c(-1.2,1.2),las=3, ylab="MA Slope" #, main="log(LMA)~log(Nmass) strict"
             , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
             ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
             , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
             , boxwex=.7, xaxt="n"
             , xlim=c(0.5,6.5), plot=F)
boxplot(Slope_LL.N~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LL.N>5& !all.results.cwm.cl$Type %in% c("global","Famclean") & is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-1.2,1.2),las=3, ylab="SMA Slope" #, main="log(LMA)~log(Nmass) strict"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(0.5,6.5))
m <- lmodel2(log.Nmass~log.LL, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.Nmass~log.LL, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
m <- lmodel2(log.cw_Nmassp_if~log.cw_LLp_if, biomass)
points(y=m$regression.results$Slope[3],x=6, pch=17, cex=1.3, col=mypal[5])
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])

abline(h=0, lty=2)
mtext(text = "c)", side = 3, adj=0, line=.2)
mtext(text= expression(paste(N[mass], " vs LL")), side=3, line=.2)



### add in individual points with error bars for significant correlations
tmp <- all.results.cwm.cl[which(all.results.cwm.cl$n_LL.N>5 & !is.na(all.results.cwm.cl$sig_LL.N)),]
tmp$rho.sig <- rho.sig(rho=tmp$Rho_LL.N, n= tmp$n_LL.N)
tmp.ns <- tmp[which(tmp$rho.sig>0.05),]
tmp.sig <- tmp[which(tmp$rho.sig<=0.05),]
#palette(paste0(mypal[colchoices], "55"))
set.seed(50)
points(Slope_LL.N~jitter(as.numeric(Type)), tmp.ns, pch=1, col=Type, cex=cex.ns)
#palette(mypal[colchoices])
set.seed(42)
points(Slope_LL.N~jitter(as.numeric(Type)), tmp.sig, pch=16, col=Type, cex=cex.sig)
set.seed(42)
#arrows(x0 = jitter(as.numeric(tmp.sig$Type)),y0=tmp.sig$Slope.lci_LL.N, y1=tmp.sig$Slope.uci_LL.N, length = 0, col=tmp.sig$Type)



# Rho boxplot
par(mar=c(6,4,0,1))
boxplot(Rho_LL.N~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LL.N>5& !all.results.cwm.cl$Type %in% c("global","Famclean") & is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-1,1.1),las=3, ylab="Correlation"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=16, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(0.5,6.5))
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","Global","PNW CWM"), at=c(1,2,3,4,5,6), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5,6), labels=paste0("(",all.results.cwm.cl$n_LL.N[c((nrow(all.results.cwm.cl)-2),(nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=.9)
m <- cor.test(x=allspp$log.LL, y=allspp$log.Nmass)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.Nmass)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
m <- cor.test(x=biomass$log.cw_LLp_if, y=biomass$log.cw_Nmassp_if)
points(y=m$estimate,x=6, pch=17, cex=1.3, col=mypal[5])
arrows(x0=6,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col=mypal[5])



# layering points on top of null model
points(Rho_LL.N~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LL.N>5 & all.results.cwm.cl$sig_LL.N<crit),], pch=16, col=Type,cex=cex.sig)
#palette(paste0(mypal[colchoices], "55"))
points(Rho_LL.N~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LL.N>5 & all.results.cwm.cl$sig_LL.N>=crit),], pch=1, col=Type, cex=cex.ns)
#palette(mypal[colchoies])




# scatterplot w/ lines
par(mar=c(4,4,1.5,1.5), cex.lab=1.2, mgp=c(2.5,.7,0))
ltysig <- 1
ltyns <-  2
lwdsig <-  1.8
lwdns <- 1


plot(log.Nmass~log.LL, allspp, col="grey", pch=16, xlab=expression(paste(log[10](LL))), ylab=expression(paste(log[10],(N[mass]))))
tax <- "w.inGen"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LL.N>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LL.N[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LL", yvar = "log.Nmass",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lty=linety, lwd=linewd)
  
}
tax <- "Genw.inFam"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LL.N>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LL.N[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <-ltysig ; linewd=lwdsig}
  plot.MAR(xvar = "log.LL", yvar = "log.Nmass",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lty=linety, lwd=linewd)
  
}
abline(a=all.results.cwm.cl$Int_LL.N[which(all.results.cwm.cl$Type=="Famclean")], b=all.results.cwm.cl$Slope_LL.N[which(all.results.cwm.cl$Type=="Famclean")], lwd=3, col="black")
abline(a=all.results.cwm.cl$Int_LL.N[which(all.results.cwm.cl$Type=="global")], b=all.results.cwm.cl$Slope_LL.N[which(all.results.cwm.cl$Type=="global")], lwd=3, col="black", lty=2)
tax <- "w.inSpp"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LL.N>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LL.N[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LL", yvar = "log.Nmass",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lty =linety, lwd=linewd)
  
}
mtext(text = "d)", side = 3, adj=0, line=.2)


#points(biomass$log.cw_Nmassp_if~biomass$log.cw_LLp_if, pch=24, bg=mypal[5])
plot.MAR(xvar = "log.cw_LLp_if", yvar="log.cw_Nmassp_if", data=biomass, linecol = mypal[5], lwd = 3)


dev.off()








###### .FIG 3:  LMA vs LL ###################
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])
cex.sig <- 1.1
cex.ns <- .9

#jpeg(width=7.008, height=4, units = "in", res=600,filename = paste0("./",results_dirname,"/Fig3_LMA_LL.jpeg"))
pdf(width=7.008, height=4, file = paste0("./",results_dirname,"/Fig3_LMA_LL3.pdf"))
#quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1))
p <- boxplot(Slope_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5& !all.results.cwm.cl$Type %in% c("global","Famclean")),], plot=F
             , xlim=c(.5,6.5), ylim=c(-2.5,4))
boxplot(Slope_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5& is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-2.5,4),las=3, ylab="SMA Slope" #, main="log(LMA)~log(Nmass) strict"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(.5,6.5))
abline(h=0, lty=2)
m <- lmodel2(log.LL~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.LL~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
m <- lmodel2(log.cw_LLp_if~log.cw_LMAp_if, biomass)
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
#points(y=m$regression.results$Slope[3],x=6, pch=24, cex=1.3, col="black", bg=mypal[5])
points(y=m$regression.results$Slope[3],x=6, pch=17, cex=1.3, col=mypal[5])
# # only evergreens
# m <- lmodel2(log.LL~log.LMA, LES[which(LES$Decid.E.green=="E"),])
# points(y=m$regression.results$Slope[3],x=5, pch=1, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
# arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")

tmp <- all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & !is.na(all.results.cwm.cl$sig_LMA.LL)),]
tmp$rho.sig <- rho.sig(rho=tmp$Rho_LMA.LL, n= tmp$n_LMA.LL)
tmp.ns <- tmp[which(tmp$rho.sig>0.1),]
tmp.sig <- tmp[which(tmp$rho.sig<=0.1),]
#palette(paste0(mypal[colchoices], "55"))
set.seed(42)
points(Slope_LMA.LL~jitter(as.numeric(Type)), tmp.ns, pch=1, col=Type, cex=cex.ns)
#palette(mypal[colchoices])
set.seed(42)
points(Slope_LMA.LL~jitter(as.numeric(Type)), tmp.sig, pch=16, col=Type,cex=cex.sig)
#set.seed(42)
#arrows(x0 = jitter(as.numeric(tmp.sig$Type)),y0=tmp.sig$Slope.lci_LMA.LL, y1=tmp.sig$Slope.uci_LMA.LL, length = 0, col=tmp.sig$Type)

#points(Slope~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig<=0.1),])),amount=.3), vib.results[which(vib.results$rho.sig<=0.1),], pch=16, col=mypal[colchoices[1]])
#points(Slope~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig>0.1),])), amount=.3), vib.results[which(vib.results$rho.sig>0.1),], pch=1, col=mypal[colchoices[1]])

#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= "LL vs LMA", side=3, line=.2)
par(mar=c(6,4,0,1))
#p <- boxplot(Rho_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5& !all.results.cwm.cl$Type %in% c("global","Famclean")),]
boxplot(Rho_LMA.LL~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 &  is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-1,1.4),las=3, ylab="Correlation"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(.5,6.5))
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","Global","PNW CWM"), at=c(1,2,3,4,5,6), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5,6), labels=paste0("(",all.results.cwm.cl$n_LMA.LL[c((nrow(all.results.cwm.cl)-2),(nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=.9)
#polygon(x = c(1,2,3,3,2,1), y=na.omit(c(meancis$m10_LMA.LL, rev(meancis$m90_LMA.LL))), border="#EEEEEE",col="#EEEEEE")#lightgrey",col = "lightgrey")
m <- cor.test(x=allspp$log.LMA, y=allspp$log.LL)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.LL)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
m <- cor.test(x=biomass$log.cw_LMAp_if, y=biomass$log.cw_LLp_if)
arrows(x0=6,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col=mypal[5])
#points(y=m$estimate,x=6, pch=24, cex=1.3, col="black", bg=mypal[5])
points(y=m$estimate,x=6, pch=17, cex=1.3, col=mypal[5])
# # only evergreens
# m <- cor.test(x=LES$log.LMA[which(LES$Decid.E.green=="E")], y=LES$log.LL[which(LES$Decid.E.green=="E")])
# points(y=m$estimate,x=5, pch=1, cex=1.3, col="darkgrey")
# arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")

# layering points on top of null model
set.seed(62)
points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$sig_LMA.LL<=crit),], pch=16, col=Type, cex=cex.sig)
#palette(paste0(mypal[colchoices], "55"))
set.seed(42)
points(Rho_LMA.LL~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.LL>5 & all.results.cwm.cl$sig_LMA.LL>crit),], pch=1, col=Type ,cex=cex.ns)
#palette(mypal[colchoices])
#axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)

#points(Rho~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig<=0.1),])),amount = .5), vib.results[which(vib.results$rho.sig<=0.1),], pch=16, col=mypal[colchoices[1]])
#points(Rho~jitter(rep(1, times=nrow(vib.results[which(vib.results$rho.sig>0.1),]))), vib.results[which(vib.results$rho.sig>0.1),], pch=1, col=paste0(mypal[colchoices[1]], "55"))


## scatterplot LMA LL _________________________
par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
ltysig <- 1
ltyns <-  2
lwdsig <-  1.8
lwdns <- 1

plot(log.LL~log.LMA, allspp, col="grey", pch=16, ylab=expression(paste(log[10](LL))), xlab=expression(paste(log[10](LMA))))

tax <- "w.inGen"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.LL[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lty=linety, lwd=linewd)
}

tax <- "Genw.inFam"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.LL[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <-ltysig ; linewd=lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lty=linety, lwd=linewd)
}


abline(a=all.results.cwm.cl$Int_LMA.LL[which(all.results.cwm.cl$Type=="Famclean")], b=all.results.cwm.cl$Slope_LMA.LL[which(all.results.cwm.cl$Type=="Famclean")], lwd=3, col="black")
abline(a=all.results.cwm.cl$Int_LMA.LL[which(all.results.cwm.cl$Type=="global")], b=all.results.cwm.cl$Slope_LMA.LL[which(all.results.cwm.cl$Type=="global")], lwd=3, col="black", lty=2)
tax <- "w.inSpp"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.LL>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.LL[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lty =linety, lwd=linewd)
}
# # add in Mt. Rainier
# for (i in c("Aster alpigenus","Castilleja parviflora", "Erythronium montanum", "Lupinus arcticus","Valeriana sitchensis","Veratrum viride")){
#   linety <- ltyns
#   linewd <- lwdns
#   if(all.results.cwm.cl$rho.sig_LMA.LL[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
#   plot.MAR(xvar = "log.LMA", yvar = "log.LL",data= Rainier[which(Rainier$Species_full==i),], linecol = mypal[colchoices[1]], lty =linety, lwd=linewd)
# }

#points(biomass$log.cw_LLp_if~biomass$log.cw_LMAp_if, pch=24, bg=mypal[5])
#all.results.cwm.cl$rho.sig_LMA.LL[170] # p of rho=0.08
plot.MAR(xvar = "log.cw_LMAp_if", yvar="log.cw_LLp_if", data=biomass, linecol = mypal[5], lwd = 3, lty=ltyns)

#plot.MAR(xvar="log.LMA", yvar="log.LL", data= fam.data, linecol = mypal[colchoices[4]], lwd=2)

#plot.MAR(xvar="log.LMA", yvar="log.LL", data= LES, linecol = "black", lwd=2)
mtext(text = "b)", side = 3, adj=0, line=.2)

dev.off()










#________________________________________________________________
###### .FIG 4ab: LMA v Narea ######
#________________________________________________________________
mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
cex.sig <- 1.1
cex.ns <- .9
palette(mypal[colchoices])
crit <- 0.05


#jpeg(width=7.008, height=4, units = "in", res=600,filename = paste0("./",results_dirname,"/Fig4ab_LMA_Narea.jpeg"))
pdf(width=7.008, height=4, file = paste0("./",results_dirname,"/Fig4ab_LMA_Narea.pdf"))

#quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1))
p <- boxplot(Slope_LMA.Narea~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.Narea>5& !all.results.cwm.cl$Type %in% c("global","Famclean")),], plot=F
             , xlim=c(.5,6.5), ylim=c(-2.5,4))
boxplot(Slope_LMA.Narea~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.Narea>5&  is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-0.5,2.5),las=3, ylab="SMA Slope" #, main="log(LMA)~log(Nmass) strict"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        ,xlim=c(.5,6.5))
# had to remove the two crazy outliers (Protea repens and genus Abies)
abline(h=0, lty=2)
abline(h=1, col="grey")
m <- lmodel2(log.Narea~log.LMA, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.Marea~log.LMA, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
m <- lmodel2(log.cw_Nareap_if~log.cw_LMAp_if, biomass)
points(y=m$regression.results$Slope[3],x=6, pch=17, cex=1.3, col=mypal[5])
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])



tmp <- all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.Narea>5 & !is.na(all.results.cwm.cl$sig_LMA.Narea)),]
tmp$rho.sig <- rho.sig(rho=tmp$Rho_LMA.Narea, n= tmp$n_LMA.Narea)
tmp.ns <- tmp[which(tmp$rho.sig>0.1),]
tmp.sig <- tmp[which(tmp$rho.sig<=0.1),]
#palette(paste0(mypal[colchoices], "55"))
set.seed(42)
points(Slope_LMA.Narea~jitter(as.numeric(Type)), tmp.ns, pch=1, col=Type,cex=cex.ns)
palette(mypal[colchoices])
set.seed(42)
points(Slope_LMA.Narea~jitter(as.numeric(Type)), tmp.sig, pch=16, col=Type,cex=cex.sig)
#set.seed(42)
#arrows(x0 = jitter(as.numeric(tmp.sig$Type)),y0=tmp.sig$Slope.lci_LMA.Narea, y1=tmp.sig$Slope.uci_LMA.Narea, length = 0, col=tmp.sig$Type)



#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "a)", side = 3, adj=0, line=.2)
mtext(text= expression(paste(N[area]," vs LMA")), side=3, line=.2)
par(mar=c(6,4,0,1))
boxplot(Rho_LMA.Narea~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.Narea>5& !all.results.cwm.cl$Type %in% c("global","Famclean")& all.results.cwm.cl$Slope_LMA.Narea>-1 & is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-1,1.4),las=3, ylab="Correlation"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.7, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim = c(0.5,6.5))
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","Global", "PNW CWM"), at=c(1,2,3,4,5,6), las=3)
abline(h=0, lty=2)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5,6), labels=paste0("(",all.results.cwm.cl$n_LMA.Narea[c((nrow(all.results.cwm.cl)-2),(nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=1)
m <- cor.test(x=allspp$log.LMA, y=allspp$log.Narea)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LMA, y=fam.dataclean$log.Narea)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
m <- cor.test(x=biomass$log.cw_LMAp_if, y=biomass$log.cw_Nareap_if)
points(y=m$estimate,x=6, pch=17, cex=1.3, col=mypal[5])
arrows(x0=6,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col=mypal[5])

# layering points on top of null model
points(Rho_LMA.Narea~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.Narea>5 & all.results.cwm.cl$sig_LMA.Narea<crit),], pch=16, col=Type,cex=cex.sig)
#palette(paste0(mypal[colchoices], "55"))
points(Rho_LMA.Narea~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LMA.Narea>5 & all.results.cwm.cl$sig_LMA.Narea>=crit),], pch=1, col=Type,cex=cex.ns)
palette(mypal[colchoices])
#axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","Global"), at=c(1,2,3,4,5), las=3)



## scatterplot
par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
ltysig <- 1
ltyns <-  2
lwdsig <-  1.8
lwdns <- 1


plot(log.Narea~log.LMA, allspp, col="grey", pch=16, ylab=expression(paste(log[10](N[area]))), xlab=expression(paste(log[10](LMA))))

tax <- "w.inGen"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.Narea>5 & all.results.cwm.cl$Slope_LMA.Narea > -1)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.Narea[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lty=linety, lwd=linewd)
  
}

tax <- "Genw.inFam"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.Narea>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.Narea[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <-ltysig ; linewd=lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lty=linety, lwd=linewd)
  
}

abline(a=all.results.cwm.cl$Int_LMA.Narea[which(all.results.cwm.cl$Type=="Famclean")], b=all.results.cwm.cl$Slope_LMA.Narea[which(all.results.cwm.cl$Type=="Famclean")], lwd=3, col="black")
abline(a=all.results.cwm.cl$Int_LMA.Narea[which(all.results.cwm.cl$Type=="global")], b=all.results.cwm.cl$Slope_LMA.Narea[which(all.results.cwm.cl$Type=="global")], lwd=3, col="black", lty=2)

# plot.MAR(xvar="log.LMA", yvar="log.Narea", data= fam.data, linecol = mypal[5], lwd=2)
# 
# plot.MAR(xvar="log.LMA", yvar="log.Narea", data= LES, linecol = "black", lwd=2)
tax <- "w.inSpp"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LMA.Narea>5& all.results.cwm.cl$Slope_LMA.Narea>-1)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LMA.Narea[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LMA", yvar = "log.Narea",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lty =linety, lwd=linewd)
  
}
mtext(text = "b)", side = 3, adj=0, line=.2)

#points(biomass$log.cw_Nareap_if~biomass$log.cw_LMAp_if, pch=24, bg=mypal[5])
plot.MAR(xvar = "log.cw_LMAp_if", yvar="log.cw_Nareap_if", data=biomass, linecol = mypal[5], lwd = 3)

dev.off()



#________________________________________________________________
###### .FIG 4cd: LL v Narea ######
#________________________________________________________________
# 2 panel figure with boxplots and scatterplot
# width = 1.5 columns -> 11.4 cm,4.89 in
#       = 2 columns -> 17.8 cm,7.008 in
# I think the easiest thing will be to make two quartez: top with boxplots and bottom with funnel
mat <- matrix(c(1,3,
                2,3), nrow=2, byrow = T)
colchoices <- c(1,2,4,3,6)
palette(mypal[colchoices])

#jpeg(width=7.008, height=4, units = "in", res=600,filename = paste0("./",results_dirname,"/Fig4cd_LL_Narea.jpeg"))
pdf(width=7.008, height=4, file = paste0("./",results_dirname,"/Fig4cd_LL_Narea2.pdf"))

#quartz(width=7.008, height=4)
layout(mat, heights=c(1,1.5))
par(mar=c(0,4,1.5,1))
p <- boxplot(Slope_LL.Narea~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LL.Narea>5 & !all.results.cwm.cl$Type %in% c("global","Famclean")),]
             , ylim=c(-2,2),las=3, ylab="SMA Slope" #, main="log(LMA)~log(Nmass) strict"
             , xlim=c(0.5,6.5), plot=F)
boxplot(Slope_LL.Narea~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LL.Narea>5 & !all.results.cwm.cl$Type %in% c("global","Famclean") & is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-2,2),las=3, ylab="SMA Slope" #, main="log(LMA)~log(Nmass) strict"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.5, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(0.5,6.5))
abline(h=0, lty=2)
#text(y=par()$usr[3]+.2, x=c(1,2,3,4,5,6), labels = p$n)
#text(y=par()$usr[4]-.2,x=.5, labels = "a)")
mtext(text = "c)", side = 3, adj=0, line=.2)
mtext(text= expression(paste(N[area]," vs LL")), side=3, line=.2)
m <- lmodel2(log.Narea~log.LL, allspp)
points(y=m$regression.results$Slope[3],x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
arrows(x0=5,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col="darkgrey")
m <- lmodel2(log.Narea~log.LL, fam.dataclean)
points(y=m$regression.results$Slope[3],x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3)
arrows(x0=4,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3)
m <- lmodel2(log.cw_Nareap_if~log.cw_LLp_if, biomass)
points(y=m$regression.results$Slope[3],x=6, pch=17, cex=1.3, col=mypal[5])
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`2.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])
arrows(x0=6,y0=m$regression.results$Slope[3], y1= m$confidence.intervals$`97.5%-Slope`[3], length = 0,lwd=3, col=mypal[5])


tmp <- all.results.cwm.cl[which(all.results.cwm.cl$n_LL.Narea>5 & !is.na(all.results.cwm.cl$sig_LL.Narea)),]
tmp$rho.sig <- rho.sig(rho=tmp$Rho_LL.Narea, n= tmp$n_LL.Narea)
tmp.ns <- tmp[which(tmp$rho.sig>0.1),]
tmp.sig <- tmp[which(tmp$rho.sig<=0.1),]
#palette(paste0(mypal[colchoices], "55"))
set.seed(42)
points(Slope_LL.Narea~jitter(as.numeric(Type)), tmp.ns, pch=1, col=Type,cex=cex.ns)
palette(mypal[colchoices])
set.seed(42)
points(Slope_LL.Narea~jitter(as.numeric(Type)), tmp.sig, pch=16, col=Type,cex=cex.sig)
#set.seed(42)
#arrows(x0 = jitter(as.numeric(tmp.sig$Type)),y0=tmp.sig$Slope.lci_LL.Narea, y1=tmp.sig$Slope.uci_LL.Narea, length = 0, col=tmp.sig$Type)


par(mar=c(6,4,0,1))
boxplot(Rho_LL.Narea~Type, all.results.cwm.cl[which(all.results.cwm.cl$n_LL.Narea>5& !all.results.cwm.cl$Type %in% c("global","Famclean") & is.na(all.results.cwm.cl$Taxo.Unit)),]
        , ylim=c(-1,1.4),las=3, ylab="Correlation"
        , col=paste0(mypal[colchoices],"66"), boxcol=paste0(mypal[colchoices],"66")
        ,whisklty=1, whisklwd=3, whiskcol=paste0(mypal[colchoices],"AA")
        , staplelwd=0, outpch=NA, outcex=.5, outcol=mypal[colchoices]
        , boxwex=.7, xaxt="n"
        , xlim=c(0.5,6.5))
abline(h=0, lty=2)
axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","Global","PNW CWM"), at=c(1,2,3,4,5,6), las=3)
text(y=par()$usr[4]-.2, x=c(1,2,3), labels = p$n[1:3])
text(y=par()$usr[4]-.2, x=c(4,5,6), labels=paste0("(",all.results.cwm.cl$n_LL.Narea[c((nrow(all.results.cwm.cl)-2),(nrow(all.results.cwm.cl)-1),nrow(all.results.cwm.cl))],")"),cex=1)
m <- cor.test(x=allspp$log.LL, y=allspp$log.Narea)
points(y=m$estimate,x=5, pch=16, cex=1.3, col="darkgrey")
arrows(x0=5,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col="darkgrey")
m <- cor.test(x=fam.dataclean$log.LL, y=fam.dataclean$log.Narea)
points(y=m$estimate,x=4, pch=16, cex=1.3, col="black")
arrows(x0=4,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3)
m <- cor.test(x=biomass$log.cw_LLp_if, y=biomass$log.cw_Nareap_if)
points(y=m$estimate,x=6, pch=17, cex=1.3, col=mypal[5])
arrows(x0=6,y0=m$conf.int[2], y1= m$conf.int[1], length = 0,lwd=3, col=mypal[5])

# layering points on top of null model
points(Rho_LL.Narea~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LL.Narea>5 & all.results.cwm.cl$sig_LL.Narea<crit),], pch=16, col=Type,cex=cex.sig)
#palette(paste0(mypal[colchoices], "55"))
points(Rho_LL.Narea~jitter(as.numeric(Type)), all.results.cwm.cl[which(all.results.cwm.cl$n_LL.Narea>5 & all.results.cwm.cl$sig_LL.Narea>=crit),], pch=1, col=Type,cex=cex.ns)
palette(mypal[colchoices])
#axis(1,labels = c("w/in Spp","w/in Gen", "w/in Fam","btw Fam","global"), at=c(1,2,3,4,5), las=3)


### scatter plot with SMA regs 
par(mar=c(4,4,1.5,1.5), mgp=c(2.5,1,0))
## scatterplot
ltysig <- 1 # line type for slopes w/ significant correlations
ltyns <-  2 # line type for slopes w/ ns corrs
lwdsig <-  1.8 # line width for slopes w/ sig corrs
lwdns <- 1 # line width for slopes w/ ns corrs

plot(log.Narea~log.LL, allspp, col="grey", pch=16,  ylab=expression(paste(log[10](N[area]))), xlab=expression(paste(log[10](LL))))

tax <- "w.inSpp"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LL.Narea>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LL.Narea[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LL", yvar = "log.Narea",data= spp.data[which(spp.data$Species==i),], linecol = mypal[colchoices[1]], lty =linety, lwd=linewd)
  
}

tax <- "w.inGen"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LL.Narea>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LL.Narea[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <- ltysig; linewd <- lwdsig}
  plot.MAR(xvar = "log.LL", yvar = "log.Narea",data= gen.data[which(gen.data$Genus==i),], linecol = mypal[colchoices[2]], lty=linety, lwd=linewd)
  
}

tax <- "Genw.inFam"
for (i in as.character(all.results.cwm.cl$Taxo.Unit[which(all.results.cwm.cl$Type==tax & all.results.cwm.cl$n_LL.Narea>5)])){
  linety <- ltyns
  linewd <- lwdns
  if(all.results.cwm.cl$rho.sig_LL.Narea[which(all.results.cwm.cl$Taxo.Unit==i)]<0.05) {linety <-ltysig ; linewd=lwdsig}
  plot.MAR(xvar = "log.LL", yvar = "log.Narea",data= geninfam.data[which(geninfam.data$Family==i),], linecol = mypal[colchoices[3]], lty=linety, lwd=linewd)
  
}

abline(a=all.results.cwm.cl$Int_LL.Narea[which(all.results.cwm.cl$Type=="Famclean")], b=all.results.cwm.cl$Slope_LL.Narea[which(all.results.cwm.cl$Type=="Famclean")], lwd=3, col="black")
abline(a=all.results.cwm.cl$Int_LL.Narea[which(all.results.cwm.cl$Type=="global")], b=all.results.cwm.cl$Slope_LL.Narea[which(all.results.cwm.cl$Type=="global")], lwd=3, col="black", lty=2)

# plot.MAR(xvar="log.Narea", yvar="log.LL", data= fam.dataclean, linecol = mypal[5], lwd=2)
# 
# plot.MAR(xvar="log.Narea", yvar="log.LL", data= LES, linecol = "black", lwd=2)
mtext(text = "d)", side = 3, adj=0, line=.2)

#points(biomass$log.cw_Nareap_if~biomass$log.cw_LLp_if, pch=24, bg=mypal[5])
plot.MAR(xvar = "log.cw_LLp_if", yvar="log.cw_Nareap_if", data=biomass, linecol = mypal[5], lwd = 3)

dev.off()








####### Table S4: Correlation and Slopes #############
# summarize the mean correlations and mean SMA slope (of things with significant correlations)
lma.n <- all.results.cwm.cl %>% group_by(Type) %>% filter(n_LMA.N>5) %>% summarise(cor.m = mean(Rho_LMA.N, na.rm=T), cor.se = sterr(Rho_LMA.N), slope.m = mean(Slope_LMA.N[which(rho.sig_LMA.N<0.1)], na.rm=T), slope.se = sterr(Slope_LMA.N[which(rho.sig_LMA.N<0.1)]))
lma.ll <- all.results.cwm.cl %>% group_by(Type) %>% filter(n_LMA.LL>5) %>% summarise(cor.m = mean(Rho_LMA.LL, na.rm=T), cor.se = sterr(Rho_LMA.LL), slope.m = mean(Slope_LMA.LL[which(rho.sig_LMA.LL<0.1)], na.rm=T), slope.se = sterr(Slope_LMA.LL[which(rho.sig_LMA.LL<0.1)]))
lma.narea <- all.results.cwm.cl %>% group_by(Type) %>% filter(n_LMA.Narea>5) %>% summarise(cor.m = mean(Rho_LMA.Narea, na.rm=T), cor.se = sterr(Rho_LMA.Narea), slope.m = mean(Slope_LMA.Narea[which(rho.sig_LMA.Narea<0.1)], na.rm=T), slope.se = sterr(Slope_LMA.Narea[which(rho.sig_LMA.Narea<0.1)]))
ll.narea <- all.results.cwm.cl %>% group_by(Type) %>% filter(n_LL.Narea>5) %>% summarise(cor.m = mean(Rho_LL.Narea, na.rm=T), cor.se = sterr(Rho_LL.Narea), slope.m = mean(Slope_LL.Narea[which(rho.sig_LL.Narea<0.1)], na.rm=T), slope.se = sterr(Slope_LL.Narea[which(rho.sig_LL.Narea<0.1)]))
n.ll <- all.results.cwm.cl %>% group_by(Type) %>% filter(n_LL.N>5) %>% summarise(cor.m = mean(Rho_LL.N, na.rm=T), cor.se = sterr(Rho_LL.N), slope.m = mean(Slope_LL.N[which(rho.sig_LL.N<0.1)], na.rm=T), slope.se = sterr(Slope_LL.N[which(rho.sig_LL.N<0.1)]))
# create function to extract mean and SE values for Table S4
extract.values <- function (dataz, taxo, variable){
  if(is.na(dataz[taxo,paste(variable,"se", sep=".")])){
    value <- as.character(round(dataz[taxo, paste(variable,"m", sep=".")],digits = 2))
  }
  else{
    value <- paste(round(dataz[taxo,paste(variable, "m", sep=".")], digits=2), round(dataz[taxo, paste(variable, "se", sep=".")], digits=2), sep="±")
  }
  return(value)
}

lmas <- list(LMA.LL=lma.ll, LMA.Nmass=lma.n, LMA.Narea=lma.narea)
lls <- list(n.ll, ll.narea)


###### Correlations ########
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=1, variable="cor")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=1, variable="cor")}))
table.cor.spp <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=4, variable="cor")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=4, variable="cor")}))
table.cor.btwfam <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=5, variable="cor")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=5, variable="cor")}))
table.cor.global <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
# new CWMs
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=6, variable="cor")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=6, variable="cor")}))
table.cor.cwm <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))


##### SMA Slopes #############
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=1, variable="slope")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=1, variable="slope")}))
table.slope.spp <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=4, variable="slope")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=4, variable="slope")}))
table.slope.btwfam <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=5, variable="slope")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=5, variable="slope")}))
table.slope.global <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))
## new for CWM
c1 <- unlist(lapply(X=lmas, FUN = function(X){extract.values(dataz=X, taxo=6, variable="slope")}))
c2 <- unlist(lapply(X=lls, FUN = function(X){extract.values(dataz=X, taxo=6, variable="slope")}))
table.slope.cwm <- data.frame(LMA = c1, LL = c("-", c2), row.names=c("LL","Nmass","Narea"))


tableS4cor <- cbind(table.cor.spp, table.cor.global, table.cor.cwm)
write.csv(tableS4cor, file = paste0("./",results_dirname,"/TableS4a_Correlations.csv"))
tableS4slope <- cbind(table.slope.spp, table.slope.global, table.slope.cwm)
# fill in the CWM slope for LL-Narea, which is removed because it has a ns correlation
tableS4slope[,6] <- as.character(tableS4slope[,6])
tableS4slope[3,6] <- round(all.results.cwm.cl$Slope_LL.Narea[which(all.results.cwm.cl$Type=="CWM")],2)
write.csv(tableS4slope, file = paste0("./",results_dirname,"/TableS4b_Slopes.csv"))



########## END: Taxonomic Scale Statistical Tests ###################







#__________________________________________________________________________________
################ BEGIN: Trait-Environment Analysis & Fig 5 #####################
#__________________________________________________________________________________





####################### .Model Selection/Ensambling ####################
# set na.action to "na.fail" so MuMIn functions work
options(na.action = "na.fail")


####### Function for fitting models without ECOREGION and with full range for Z-score ######
trait.mods.ne.z <- function(traitdata =traits.common, species, trait="log.LL", modcrit=F  ){
    # subset data for a species (only complete cases)
  dataz <- data.frame(traitdata) %>% filter(SP.ID==species) %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
  print(nrow(dataz))
  dataz$log.ASA <- log(dataz$ASA)
  datazsc <- scale(dataz[,-c(1:2)]) # scale traits and predictors by mean and sd of species subset
  colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
  datazall <- data.frame(dataz, datazsc)
  # rescale with the full mean and sd from biomass (i.e. across all plots in the network)
  datazall$climPC1sc <- (datazall$climPC1- mean(biomass$climPC1, na.rm=T))/sd(biomass$climPC1, na.rm=T)
  datazall$climPC2sc <- (datazall$climPC2- mean(biomass$climPC2, na.rm=T))/sd(biomass$climPC2, na.rm=T)
  datazall$soil_Nsc <- (datazall$soil_N- mean(biomass$soil_N, na.rm=T))/sd(biomass$soil_N, na.rm=T)
  datazall$log.ASAsc <- (datazall$log.ASA- mean(log(biomass$ASA), na.rm=T))/sd(log(biomass$ASA), na.rm=T)
  datazall$LAI_Osc <- (datazall$LAI_O- mean(biomass$LAI_O, na.rm=T))/sd(biomass$LAI_O, na.rm=T)
  datazall$AG_TGROWTHsc <- (datazall$AG_TGROWTH- mean(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T))/sd(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T)
  traitsc <- paste0(trait, "sc")
  traitmod <- lmer(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc + (1|PLOT_ID), datazall, REML=F)
  traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM))
  traitmodcalls <- dredge(traitmod, evaluate = F)
  traitmodavg <- model.avg(traitdredge,subset=delta<=4)
  
  bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
  if(modcrit==T){
    scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
    plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
    qqp(resid(bestmodobject), main="residuals")
    qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')
  }
  results <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","r.squaredGLMM.R2c",'AICc')]
                  , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")]
                  , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                  , fulldredge = traitdredge
                  , fullmodavg = traitmodavg
                  , modnum = rownames(traitdredge)[1]
                  , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                  , bestmodobject = bestmodobject
                  , n=nrow(datazall)
                  , nplots = length(ranef(bestmodobject)[[1]][,1])
                  , necos = length(unique(datazall$ECOREGION))
                  , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])
  
  return(results)
}





#___________________________________________________________________________________
############ . Within-species models for common conifers #############
#___________________________________________________________________________________


PSEMENll.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PSEMEN",trait = "log.LL", modcrit=T)
PSEMENlma.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PSEMEN",trait = "log.LMA", modcrit=T)
PSEMENnarea.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PSEMEN",trait = "log.Narea", modcrit=T)
PSEMENnmass.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PSEMEN",trait = "log.Nmass", modcrit=T)

PINPONll.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINPON",trait = "log.LL")
PINPONlma.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINPON",trait = "log.LMA")
PINPONnarea.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINPON",trait = "log.Narea")
PINPONnmass.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINPON",trait = "log.Nmass")

ABICONll.ne.z <- trait.mods.ne.z(traitdata = traits, species = "ABICON",trait = "log.LL")
ABICONlma.ne.z <- trait.mods.ne.z(traitdata = traits, species = "ABICON",trait = "log.LMA")
ABICONnarea.ne.z <- trait.mods.ne.z(traitdata = traits, species = "ABICON",trait = "log.Narea")
ABICONnmass.ne.z <- trait.mods.ne.z(traitdata = traits, species = "ABICON",trait = "log.Nmass")

TSUHETll.ne.z <- trait.mods.ne.z(traitdata = traits, species = "TSUHET",trait = "log.LL")
TSUHETlma.ne.z <- trait.mods.ne.z(traitdata = traits, species = "TSUHET",trait = "log.LMA")
TSUHETnarea.ne.z <- trait.mods.ne.z(traitdata = traits, species = "TSUHET",trait = "log.Narea") # might be one outlier that's leveraging things?
TSUHETnmass.ne.z <- trait.mods.ne.z(traitdata = traits, species = "TSUHET",trait = "log.Nmass") # might be one outlier that's leveraging things?

PINCONll.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINCON",trait = "log.LL")
PINCONlma.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINCON",trait = "log.LMA")
PINCONnarea.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINCON",trait = "log.Narea")
PINCONnmass.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINCON",trait = "log.Nmass")

PINJEFll.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINJEF",trait = "log.LL")
PINJEFlma.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINJEF",trait = "log.LMA") # two outliers
PINJEFnarea.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINJEF",trait = "log.Narea")
PINJEFnmass.ne.z <- trait.mods.ne.z(traitdata = traits, species = "PINJEF",trait = "log.Nmass")




####### . Community Weighted Traits ##########

####### + Leaf Lifespan######
trait <- "log.cw_LLp_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
colnames(dataz)[grep("AG_", colnames(dataz))] <- "AG_TGROWTH" # change this column name for compatibility
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
datazall <- data.frame(dataz, datazsc)
# z-score standardize from entire plot dataset, not just the plots with complete info (this shouldn't change much)
datazall$climPC1sc <- (datazall$climPC1- mean(biomass$climPC1, na.rm=T))/sd(biomass$climPC1, na.rm=T)
datazall$climPC2sc <- (datazall$climPC2- mean(biomass$climPC2, na.rm=T))/sd(biomass$climPC2, na.rm=T)
datazall$soil_Nsc <- (datazall$soil_N- mean(biomass$soil_N, na.rm=T))/sd(biomass$soil_N, na.rm=T)
datazall$log.ASAsc <- (datazall$log.ASA- mean(log(biomass$ASA), na.rm=T))/sd(log(biomass$ASA), na.rm=T)
datazall$LAI_Osc <- (datazall$LAI_O- mean(biomass$LAI_O, na.rm=T))/sd(biomass$LAI_O, na.rm=T)
datazall$AG_TGROWTHsc <- (datazall$AG_TGROWTH- mean(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T))/sd(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T)

tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

CWll.ne.z <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                  , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                  , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                  , fulldredge = traitdredge
                  , fullmodavg = traitmodavg
                  , modnum = rownames(traitdredge)[1]
                  , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                  , bestmodobject = bestmodobject
                  , n=nrow(datazall)
                  , nplots = nrow(datazall)
                  , necos = length(unique(datazall$ECOREGION))
                  , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


####### + CW LMA ####
trait="log.cw_LMAp_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
colnames(dataz)[grep("AG_", colnames(dataz))] <- "AG_TGROWTH" # change this column name for compatibility
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
datazall <- data.frame(dataz, datazsc)
# z-score standardize from entire plot dataset, not just the plots with complete info (this shouldn't change much)
datazall$climPC1sc <- (datazall$climPC1- mean(biomass$climPC1, na.rm=T))/sd(biomass$climPC1, na.rm=T)
datazall$climPC2sc <- (datazall$climPC2- mean(biomass$climPC2, na.rm=T))/sd(biomass$climPC2, na.rm=T)
datazall$soil_Nsc <- (datazall$soil_N- mean(biomass$soil_N, na.rm=T))/sd(biomass$soil_N, na.rm=T)
datazall$log.ASAsc <- (datazall$log.ASA- mean(log(biomass$ASA), na.rm=T))/sd(log(biomass$ASA), na.rm=T)
datazall$LAI_Osc <- (datazall$LAI_O- mean(biomass$LAI_O, na.rm=T))/sd(biomass$LAI_O, na.rm=T)
datazall$AG_TGROWTHsc <- (datazall$AG_TGROWTH- mean(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T))/sd(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T)

tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)

bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])

scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')


CWlma.ne.z <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                   , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                   , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                   , fulldredge = traitdredge
                   , fullmodavg = traitmodavg
                   , modnum = rownames(traitdredge)[1]
                   , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                   , bestmodobject = bestmodobject
                   , n=nrow(datazall)
                   , nplots = nrow(datazall)
                   , necos = length(unique(datazall$ECOREGION))
                   , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


#### + CW Nmass ####
trait="log.cw_Nmassp_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
colnames(dataz)[grep("AG_", colnames(dataz))] <- "AG_TGROWTH" # change this column name for compatibility
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
datazall <- data.frame(dataz, datazsc)
# z-score standardize from entire plot dataset, not just the plots with complete info (this shouldn't change much)
datazall$climPC1sc <- (datazall$climPC1- mean(biomass$climPC1, na.rm=T))/sd(biomass$climPC1, na.rm=T)
datazall$climPC2sc <- (datazall$climPC2- mean(biomass$climPC2, na.rm=T))/sd(biomass$climPC2, na.rm=T)
datazall$soil_Nsc <- (datazall$soil_N- mean(biomass$soil_N, na.rm=T))/sd(biomass$soil_N, na.rm=T)
datazall$log.ASAsc <- (datazall$log.ASA- mean(log(biomass$ASA), na.rm=T))/sd(log(biomass$ASA), na.rm=T)
datazall$LAI_Osc <- (datazall$LAI_O- mean(biomass$LAI_O, na.rm=T))/sd(biomass$LAI_O, na.rm=T)
datazall$AG_TGROWTHsc <- (datazall$AG_TGROWTH- mean(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T))/sd(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T)

tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)

bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])

scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')


CWnmass.ne.z <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                     , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                     , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                     , fulldredge = traitdredge
                     , fullmodavg = traitmodavg
                     , modnum = rownames(traitdredge)[1]
                     , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                     , bestmodobject = bestmodobject
                     , n=nrow(datazall)
                     , nplots = nrow(datazall)
                     , necos = length(unique(datazall$ECOREGION))
                     , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


###### + CW Narea ###########
trait <- "log.cw_Nareap_if"
dataz <- data.frame(biomass)  %>% select(PLOT_ID, ECOREGION, get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_PROD_TREE_TOTAL_AS_CARBON) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
colnames(dataz)[grep("AG_", colnames(dataz))] <- "AG_TGROWTH" # change this column name for compatibility
datazsc <- scale(dataz[,-c(1:2)]) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz[,-c(1:2)]),'sc')
datazall <- data.frame(dataz, datazsc)
# z-score standardize from entire plot dataset, not just the plots with complete info (this shouldn't change much)
datazall$climPC1sc <- (datazall$climPC1- mean(biomass$climPC1, na.rm=T))/sd(biomass$climPC1, na.rm=T)
datazall$climPC2sc <- (datazall$climPC2- mean(biomass$climPC2, na.rm=T))/sd(biomass$climPC2, na.rm=T)
datazall$soil_Nsc <- (datazall$soil_N- mean(biomass$soil_N, na.rm=T))/sd(biomass$soil_N, na.rm=T)
datazall$log.ASAsc <- (datazall$log.ASA- mean(log(biomass$ASA), na.rm=T))/sd(log(biomass$ASA), na.rm=T)
datazall$LAI_Osc <- (datazall$LAI_O- mean(biomass$LAI_O, na.rm=T))/sd(biomass$LAI_O, na.rm=T)
datazall$AG_TGROWTHsc <- (datazall$AG_TGROWTH- mean(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T))/sd(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T)

tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
#  traitmodavg <- model.avg(traitdredge,subset=cumsum(weight)<=.90)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)

bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])

scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')


CWnarea.ne.z <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                     , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                     , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                     , fulldredge = traitdredge
                     , fullmodavg = traitmodavg
                     , modnum = rownames(traitdredge)[1]
                     , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                     , bestmodobject = bestmodobject
                     , n=nrow(datazall)
                     , nplots = nrow(datazall)
                     , necos = length(unique(datazall$ECOREGION))
                     , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])




### rename R2 column calculated by dredge so it aligns with w/in species conditional R2. Will use this R2 rather than the marginal R2 from r.squaredGLMM
colnames(CWll.ne.z$best)[9] <- "r.squaredGLMM.R2c"
colnames(CWlma.ne.z$best)[9] <- "r.squaredGLMM.R2c"
colnames(CWnmass.ne.z$best)[9] <- "r.squaredGLMM.R2c"
colnames(CWnarea.ne.z$best)[9] <- "r.squaredGLMM.R2c"






###### . Species mean traits  ###########
# uses object spp.traits calculated for producing plot CWMs 

####### + SPP Leaf Lifespan  ######
trait <- "mlog.LL"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
# z-score standardize from entire plot dataset, not species means
datazall$climPC1sc <- (datazall$climPC1- mean(biomass$climPC1, na.rm=T))/sd(biomass$climPC1, na.rm=T)
datazall$climPC2sc <- (datazall$climPC2- mean(biomass$climPC2, na.rm=T))/sd(biomass$climPC2, na.rm=T)
datazall$soil_Nsc <- (datazall$soil_N- mean(biomass$soil_N, na.rm=T))/sd(biomass$soil_N, na.rm=T)
datazall$log.ASAsc <- (datazall$log.ASA- mean(log(biomass$ASA), na.rm=T))/sd(log(biomass$ASA), na.rm=T)
datazall$LAI_Osc <- (datazall$LAI_O- mean(biomass$LAI_O, na.rm=T))/sd(biomass$LAI_O, na.rm=T)
datazall$AG_TGROWTHsc <- (datazall$AG_TGROWTH- mean(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T))/sd(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T)

tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPll.ne.z <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                   , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                   , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                   , fulldredge = traitdredge
                   , fullmodavg = traitmodavg
                   , modnum = rownames(traitdredge)[1]
                   , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                   , bestmodobject = bestmodobject
                   , n=nrow(datazall)
                   #             , nplots = nrow(datazall)
                   #             , necos = length(unique(datazall$ECOREGION))
                   , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])


####### + SPP LMA  ####
trait="mlog.LMA"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
# z-score standardize from entire plot dataset, not species means
datazall$climPC1sc <- (datazall$climPC1- mean(biomass$climPC1, na.rm=T))/sd(biomass$climPC1, na.rm=T)
datazall$climPC2sc <- (datazall$climPC2- mean(biomass$climPC2, na.rm=T))/sd(biomass$climPC2, na.rm=T)
datazall$soil_Nsc <- (datazall$soil_N- mean(biomass$soil_N, na.rm=T))/sd(biomass$soil_N, na.rm=T)
datazall$log.ASAsc <- (datazall$log.ASA- mean(log(biomass$ASA), na.rm=T))/sd(log(biomass$ASA), na.rm=T)
datazall$LAI_Osc <- (datazall$LAI_O- mean(biomass$LAI_O, na.rm=T))/sd(biomass$LAI_O, na.rm=T)
datazall$AG_TGROWTHsc <- (datazall$AG_TGROWTH- mean(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T))/sd(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T)

tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPlma.ne.z <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                    , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                    , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                    , fulldredge = traitdredge
                    , fullmodavg = traitmodavg
                    , modnum = rownames(traitdredge)[1]
                    , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                    , bestmodobject = bestmodobject
                    , n=nrow(datazall)
                    #             , nplots = nrow(datazall)
                    #             , necos = length(unique(datazall$ECOREGION))
                    , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])

####  + SPP Nmass ####
trait="mlog.Nmass"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
# z-score standardize from entire plot dataset, not species means
datazall$climPC1sc <- (datazall$climPC1- mean(biomass$climPC1, na.rm=T))/sd(biomass$climPC1, na.rm=T)
datazall$climPC2sc <- (datazall$climPC2- mean(biomass$climPC2, na.rm=T))/sd(biomass$climPC2, na.rm=T)
datazall$soil_Nsc <- (datazall$soil_N- mean(biomass$soil_N, na.rm=T))/sd(biomass$soil_N, na.rm=T)
datazall$log.ASAsc <- (datazall$log.ASA- mean(log(biomass$ASA), na.rm=T))/sd(log(biomass$ASA), na.rm=T)
datazall$LAI_Osc <- (datazall$LAI_O- mean(biomass$LAI_O, na.rm=T))/sd(biomass$LAI_O, na.rm=T)
datazall$AG_TGROWTHsc <- (datazall$AG_TGROWTH- mean(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T))/sd(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T)

tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPnmass.ne.z <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                      , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                      , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                      , fulldredge = traitdredge
                      , fullmodavg = traitmodavg
                      , modnum = rownames(traitdredge)[1]
                      , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                      , bestmodobject = bestmodobject
                      , n=nrow(datazall)
                      #             , nplots = nrow(datazall)
                      #             , necos = length(unique(datazall$ECOREGION))
                      , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])

######  + SPP Narea ###########
trait <- "mlog.Narea"
dataz <- data.frame(spp.traits)  %>% select( get(trait), climPC1, climPC2, soil_N, ASA, LAI_O, AG_TGROWTH) %>% filter(complete.cases(.))
print(nrow(dataz))
dataz$log.ASA <- log(dataz$ASA)
datazsc <- scale(dataz) # originally, this just scaled the predictors, but I think I also want to scale the traits.
colnames(datazsc) <- paste0(colnames(dataz),'sc')
#colnames(datazsc)[grep("AG_", colnames(datazsc))] <- "AG_TGROWTHsc" # change this column name for compatibility
datazall <- data.frame(dataz, datazsc)
# z-score standardize from entire plot dataset, not species means
datazall$climPC1sc <- (datazall$climPC1- mean(biomass$climPC1, na.rm=T))/sd(biomass$climPC1, na.rm=T)
datazall$climPC2sc <- (datazall$climPC2- mean(biomass$climPC2, na.rm=T))/sd(biomass$climPC2, na.rm=T)
datazall$soil_Nsc <- (datazall$soil_N- mean(biomass$soil_N, na.rm=T))/sd(biomass$soil_N, na.rm=T)
datazall$log.ASAsc <- (datazall$log.ASA- mean(log(biomass$ASA), na.rm=T))/sd(log(biomass$ASA), na.rm=T)
datazall$LAI_Osc <- (datazall$LAI_O- mean(biomass$LAI_O, na.rm=T))/sd(biomass$LAI_O, na.rm=T)
datazall$AG_TGROWTHsc <- (datazall$AG_TGROWTH- mean(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T))/sd(biomass$AG_PROD_TREE_TOTAL_AS_CARBON, na.rm=T)

tn <- trait
traitsc <- paste(tn, "sc", sep="")
traitmod <- lm(get(traitsc)~climPC1sc + climPC2sc + soil_Nsc + log.ASAsc + LAI_Osc + AG_TGROWTHsc , datazall)
traitdredge <- dredge(traitmod, extra=list(r.squaredGLMM, "R^2"))
traitmodcalls <- dredge(traitmod, evaluate = F)
traitmodavg <- model.avg(traitdredge,subset=delta<=4)
bestmodobject <- eval(traitmodcalls[[rownames(traitdredge)[1]]])
scatter.smooth(resid(bestmodobject)~fitted(bestmodobject)); abline(h=0)
plot(datazall[,traitsc]~predict(bestmodobject, re.form=NA));abline(a=0,b=1)
qqp(resid(bestmodobject), main="residuals")
# qqp(ranef(bestmodobject)[[1]][,1], main='Random Effects')

SPPnarea.ne.z <- list(best = traitdredge[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc","r.squaredGLMM.R2m","R^2",'AICc')]
                      , avg = traitmodavg$coefficients[1,c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc")] # note: If I try to get ecoregion effects, this can blow up. I think I need to leave ECOREGION out of this
                      , imp = as.vector(traitmodavg$importance)[match(c("(Intercept)","climPC1sc","climPC2sc","soil_Nsc","log.ASAsc","LAI_Osc","AG_TGROWTHsc"),names(traitmodavg$importance))]
                      , fulldredge = traitdredge
                      , fullmodavg = traitmodavg
                      , modnum = rownames(traitdredge)[1]
                      , bestmodcall = traitmodcalls[[rownames(traitdredge)[1]]]
                      , bestmodobject = bestmodobject
                      , n=nrow(datazall)
                      #             , nplots = nrow(datazall)
                      #             , necos = length(unique(datazall$ECOREGION))
                      , deltaNULL = traitdredge[which(rownames(traitdredge)==1),"delta"])



###### rename R2 column for later concatenation with w/in spp models
colnames(SPPll.ne.z$best)[9] <- "r.squaredGLMM.R2c"
colnames(SPPlma.ne.z$best)[9] <- "r.squaredGLMM.R2c"
colnames(SPPnmass.ne.z$best)[9] <- "r.squaredGLMM.R2c"
colnames(SPPnarea.ne.z$best)[9] <- "r.squaredGLMM.R2c"












#________________________________________________________________
####### . Combining all environmental trait models (Table S5) ######
#________________________________________________________________


llbestmods.ne.z <- rbind(PSEMENll.ne.z$best, PINPONll.ne.z$best,PINCONll.ne.z$best,PINJEFll.ne.z$best, ABICONll.ne.z$best, TSUHETll.ne.z$best,SPPll.ne.z$best, CWll.ne.z$best)
llbestmods.ne.z$n <- rbind(PSEMENll.ne.z$n, PINPONll.ne.z$n,PINCONll.ne.z$n,PINJEFll.ne.z$n, ABICONll.ne.z$n, TSUHETll.ne.z$n,SPPll.ne.z$n,CWll.ne.z$n)
llbestmods.ne.z$nplots <- rbind(PSEMENll.ne.z$nplots, PINPONll.ne.z$nplots,PINCONll.ne.z$nplots,PINJEFll.ne.z$nplots, ABICONll.ne.z$nplots, TSUHETll.ne.z$nplots, NA, CWll.ne.z$nplots)
llbestmods.ne.z$necos <- rbind(PSEMENll.ne.z$necos, PINPONll.ne.z$necos,PINCONll.ne.z$necos,PINJEFll.ne.z$necos, ABICONll.ne.z$necos, TSUHETll.ne.z$necos,NA, CWll.ne.z$necos)
llbestmods.ne.z$deltas <-rbind(PSEMENll.ne.z$deltaNULL, PINPONll.ne.z$deltaNULL,PINCONll.ne.z$deltaNULL,PINJEFll.ne.z$deltaNULL, ABICONll.ne.z$deltaNULL, TSUHETll.ne.z$deltaNULL,SPPll.ne.z$deltaNULL,CWll.ne.z$deltaNULL)
llbestmods.ne.z$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "SPPmean","CWmean")

lmabestmods.ne.z <- rbind(PSEMENlma.ne.z$best, PINPONlma.ne.z$best,PINCONlma.ne.z$best,PINJEFlma.ne.z$best, ABICONlma.ne.z$best, TSUHETlma.ne.z$best, SPPlma.ne.z$best, CWlma.ne.z$best)
lmabestmods.ne.z$n <- rbind(PSEMENlma.ne.z$n, PINPONlma.ne.z$n,PINCONlma.ne.z$n,PINJEFlma.ne.z$n, ABICONlma.ne.z$n, TSUHETlma.ne.z$n, SPPlma.ne.z$n, CWlma.ne.z$n)
lmabestmods.ne.z$nplots <- rbind(PSEMENlma.ne.z$nplots, PINPONlma.ne.z$nplots,PINCONlma.ne.z$nplots,PINJEFlma.ne.z$nplots, ABICONlma.ne.z$nplots, TSUHETlma.ne.z$nplots, NA, CWlma.ne.z$nplots)
lmabestmods.ne.z$necos <- rbind(PSEMENlma.ne.z$necos, PINPONlma.ne.z$necos,PINCONlma.ne.z$necos,PINJEFlma.ne.z$necos, ABICONlma.ne.z$necos, TSUHETlma.ne.z$necos, NA, CWlma.ne.z$necos)
lmabestmods.ne.z$deltas <-rbind(PSEMENlma.ne.z$deltaNULL, PINPONlma.ne.z$deltaNULL,PINCONlma.ne.z$deltaNULL,PINJEFlma.ne.z$deltaNULL, ABICONlma.ne.z$deltaNULL, TSUHETlma.ne.z$deltaNULL, SPPlma.ne.z$deltaNULL, CWlma.ne.z$deltaNULL)
lmabestmods.ne.z$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean", "CWmean")


nareabestmods.ne.z <- rbind(PSEMENnarea.ne.z$best, PINPONnarea.ne.z$best,PINCONnarea.ne.z$best,PINJEFnarea.ne.z$best, ABICONnarea.ne.z$best, TSUHETnarea.ne.z$best, SPPnarea.ne.z$best, CWnarea.ne.z$best)
nareabestmods.ne.z$n <- rbind(PSEMENnarea.ne.z$n, PINPONnarea.ne.z$n,PINCONnarea.ne.z$n,PINJEFnarea.ne.z$n, ABICONnarea.ne.z$n, TSUHETnarea.ne.z$n, SPPnarea.ne.z$n, CWnarea.ne.z$n)
nareabestmods.ne.z$nplots <- rbind(PSEMENnarea.ne.z$nplots, PINPONnarea.ne.z$nplots,PINCONnarea.ne.z$nplots,PINJEFnarea.ne.z$nplots, ABICONnarea.ne.z$nplots, TSUHETnarea.ne.z$nplots, NA, CWnarea.ne.z$nplots)
nareabestmods.ne.z$necos <- rbind(PSEMENnarea.ne.z$necos, PINPONnarea.ne.z$necos,PINCONnarea.ne.z$necos,PINJEFnarea.ne.z$necos, ABICONnarea.ne.z$necos, TSUHETnarea.ne.z$necos, NA, CWnarea.ne.z$necos)
nareabestmods.ne.z$deltas <-rbind(PSEMENnarea.ne.z$deltaNULL, PINPONnarea.ne.z$deltaNULL,PINCONnarea.ne.z$deltaNULL,PINJEFnarea.ne.z$deltaNULL, ABICONnarea.ne.z$deltaNULL, TSUHETnarea.ne.z$deltaNULL,  SPPnarea.ne.z$deltaNULL,CWnarea.ne.z$deltaNULL)
nareabestmods.ne.z$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "SPPmean", "CWmean")
#nareabestmods.ne.z$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla")


nmassbestmods.ne.z <- rbind(PSEMENnmass.ne.z$best, PINPONnmass.ne.z$best,PINCONnmass.ne.z$best,PINJEFnmass.ne.z$best, ABICONnmass.ne.z$best, TSUHETnmass.ne.z$best,  SPPnmass.ne.z$best,CWnmass.ne.z$best)
nmassbestmods.ne.z$n <- rbind(PSEMENnmass.ne.z$n, PINPONnmass.ne.z$n,PINCONnmass.ne.z$n,PINJEFnmass.ne.z$n, ABICONnmass.ne.z$n, TSUHETnmass.ne.z$n,  SPPnmass.ne.z$n,CWnmass.ne.z$n)
nmassbestmods.ne.z$nplots <- rbind(PSEMENnmass.ne.z$nplots, PINPONnmass.ne.z$nplots,PINCONnmass.ne.z$nplots,PINJEFnmass.ne.z$nplots, ABICONnmass.ne.z$nplots, TSUHETnmass.ne.z$nplots, NA,CWnmass.ne.z$nplots)
nmassbestmods.ne.z$necos <- rbind(PSEMENnmass.ne.z$necos, PINPONnmass.ne.z$necos,PINCONnmass.ne.z$necos,PINJEFnmass.ne.z$necos, ABICONnmass.ne.z$necos, TSUHETnmass.ne.z$necos, NA,CWnmass.ne.z$necos)
nmassbestmods.ne.z$deltas <-rbind(PSEMENnmass.ne.z$deltaNULL, PINPONnmass.ne.z$deltaNULL,PINCONnmass.ne.z$deltaNULL,PINJEFnmass.ne.z$deltaNULL, ABICONnmass.ne.z$deltaNULL, TSUHETnmass.ne.z$deltaNULL,  SPPnmass.ne.z$deltaNULL,CWnmass.ne.z$deltaNULL)
nmassbestmods.ne.z$SP.ID <-  c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET", "SPPmean","CWmean")
#nmassbestmods.ne.z$Species <-  c("Pseudotsuga menziesii","Pinus ponderosa","Pinus contorta","Pinus jeffreyii","Abies concolor","Tsuga heterophylla")





# Write results to csv 

write.csv(llbestmods.ne.z, paste0("./",results_dirname,"/LeafLife_trait-environment_bestmodels.csv"))
write.csv(lmabestmods.ne.z, paste0("./",results_dirname,"/LMA_trait-environment_bestmodels.csv"))
write.csv(nareabestmods.ne.z, paste0("./",results_dirname,"/Nmass_trait-environment_bestmodels.csv"))
write.csv(nmassbestmods.ne.z, paste0("./",results_dirname,"/Narea_trait-environment_bestmodels.csv"))
# 
# llbestmods.ne.z <- read.csv("./Trait_models/LeafLife_bestmodels_noECOREG_Z_20171217.csv", row.names=1, header=T)
# lmabestmods.ne.z <- read.csv("./Trait_models/LMA_bestmodels_noECOREG_Z_20171217.csv", row.names=1, header=T)
# nmassbestmods.ne.z <- read.csv("./Trait_models/Nmass_bestmodels_noECOREG_Z_20171217.csv", row.names=1, header=T)
# nareabestmods.ne.z <- read.csv("./Trait_models/Narea_bestmodels_noECOREG_Z_20171217.csv", row.names=1, header=T)

## make Table S5:
tableS5 <- rbind(llbestmods.ne.z,lmabestmods.ne.z,nmassbestmods.ne.z,nareabestmods.ne.z)

tableS5 <- tableS5 %>% select (SP.ID,climPC1sc,climPC2sc,soil_Nsc,log.ASAsc,LAI_Osc,AG_TGROWTHsc,r.squaredGLMM.R2m,r.squaredGLMM.R2c,n)
tableS5[,-1] <- apply(tableS5[,-1],MARGIN = 2,FUN = round, digits=2)
tableS5$Trait <- rep(c('LL',"LMA","Nmass","Narea"), each=8)

write.csv(tableS5, paste0("./", results_dirname,"/TableS5_AllBestModels.csv"))



#________________________________________________________________________________________________________
##### . Model Averaged Effects Sizes (Figure 5) ######
#________________________________________________________________________________________________________


llavgmods.ne.z <- data.frame(rbind(PSEMENll.ne.z$avg, PINPONll.ne.z$avg,PINCONll.ne.z$avg,PINJEFll.ne.z$avg, ABICONll.ne.z$avg, TSUHETll.ne.z$avg,SPPll.ne.z$avg, CWll.ne.z$avg))
llavgmods.ne.z$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean", "CWmean")
llavgmods.ne.z$trait <- rep("LeafLife", times=nrow(llavgmods.ne.z))

lmaavgmods.ne.z <- data.frame(rbind(PSEMENlma.ne.z$avg, PINPONlma.ne.z$avg,PINCONlma.ne.z$avg,PINJEFlma.ne.z$avg, ABICONlma.ne.z$avg, TSUHETlma.ne.z$avg,SPPlma.ne.z$avg, CWlma.ne.z$avg))
lmaavgmods.ne.z$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean","CWmean")
lmaavgmods.ne.z$trait <- rep("LMA", times=nrow(lmaavgmods.ne.z))

nareaavgmods.ne.z <- data.frame(rbind(PSEMENnarea.ne.z$avg, PINPONnarea.ne.z$avg,PINCONnarea.ne.z$avg,PINJEFnarea.ne.z$avg, ABICONnarea.ne.z$avg, TSUHETnarea.ne.z$avg, SPPnarea.ne.z$avg, CWnarea.ne.z$avg))
nareaavgmods.ne.z$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean","CWmean")
nareaavgmods.ne.z$trait <- rep("Narea", times=nrow(nareaavgmods.ne.z))

nmassavgmods.ne.z <- data.frame(rbind(PSEMENnmass.ne.z$avg, PINPONnmass.ne.z$avg,PINCONnmass.ne.z$avg,PINJEFnmass.ne.z$avg, ABICONnmass.ne.z$avg, TSUHETnmass.ne.z$avg, SPPnmass.ne.z$avg, CWnmass.ne.z$avg))
nmassavgmods.ne.z$SP.ID <- c("PSEMEN","PINPON","PINCON","PINJEF","ABICON","TSUHET","SPPmean","CWmean")
nmassavgmods.ne.z$trait <- rep("Nmass", times=nrow(nmassavgmods.ne.z))

avgmods.ne.z <- rbind(llavgmods.ne.z, lmaavgmods.ne.z, nareaavgmods.ne.z, nmassavgmods.ne.z)[,-1]
colnames(avgmods.ne.z) <- c("climPC1", "climPC2","soil_N","Stand_Age","LAI","Growth", "SP.ID", "trait")


write.csv(avgmods.ne.z, paste0("./",results_dirname,"/EffectSizes_trait-environment_ensemblemeans.csv"))



# melt to long rather than wide df
avglong.ne.z <- melt(avgmods.ne.z, id.vars = c("SP.ID","trait"))


# add x value column for plotting
avglong.ne.z$xvals <- as.numeric(avglong.ne.z$variable)
# change level order to get Nmass before Narea
avglong.ne.z$trait <- factor(avglong.ne.z$trait)
levels(avglong.ne.z$trait) <- list(LeafLife= "LeafLife", LMA="LMA",Nmass="Nmass",Narea="Narea")

########### + Figure 5b Trait-Environment Effect Sizes ############## 

#set colors and sizes
fig5pal <- brewer.pal(n=5, "Set3")[c(1,3,4,5)]#paste0(brewer.pal(n=5, "Dark2"),"99")
panlab.cex <- .9
panlab.ln <- -1.2
CWbg <- "white"#mypal[6]#"white"
SPPbg <- "black"#mypal[colchoices[2]]#"black"#mypal[2]
boxcol <- fig5pal#"grey"#paste0(mypal[colchoices[1]],"66")#"grey"
ptcex=1.1
boxwex = 0.4


#________________________________________________________
#### Panel a): R2s for the different traits & scales


pdf(width=4.33, height=2, file = paste0("./",results_dirname,"/Fig5a_Trait-Env_effectsizes_color4_v3.pdf"))

#quartz(width=4.33, height=2)

par(mar=c(.2,5,3,3), mgp=c(2,.75,0))
boxplot(llbestmods.ne.z$r.squaredGLMM.R2m[-which(llbestmods.ne.z$SP.ID %in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=4, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol
        , range=0)
boxplot(lmabestmods.ne.z$r.squaredGLMM.R2m[-which(lmabestmods.ne.z$SP.ID%in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=3, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol[2], col=boxcol[2], outcol=boxcol[2]
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol[2], add=T , range=0)
boxplot(nmassbestmods.ne.z$r.squaredGLMM.R2m[-which(nmassbestmods.ne.z$SP.ID%in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=2, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol[3], col=boxcol[3], outcol=boxcol[3]
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol[3], add=T , range=0)
boxplot(nareabestmods.ne.z$r.squaredGLMM.R2m[-which(nareabestmods.ne.z$SP.ID%in% c("SPPmean","CWmean"))]~rep(1, times=6), horizontal=T, at=1, xlim=c(0,5), ylim=c(0,1)
        , yaxt="n", xaxt="n"
        ,outpch=1, border=boxcol[4], col=boxcol[4], outcol=boxcol[4]
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=.7, whiskcol=boxcol[4], add=T , range=0)


#plot(llbestmods.ne.z$r.squaredGLMM.R2m[-which(llbestmods.ne.z$SP.ID=="CWmean")]~rep(1, times=6), ylim=c(0,1), xlim=c(0,5), xaxt="n", ylab=expression(paste(R^2," or marginal ",R^2)), pch=16)
#points(lmabestmods.ne.z$r.squaredGLMM.R2m[-which(lmabestmods.ne.z$SP.ID=="CWmean")]~rep(2, times=6), pch=16)
#points(nmassbestmods.ne.z$r.squaredGLMM.R2m[-which(nmassbestmods.ne.z$SP.ID=="CWmean")]~rep(3, times=6), pch=16)
#points(nareabestmods.ne.z$r.squaredGLMM.R2m[-which(nareabestmods.ne.z$SP.ID=="CWmean")]~rep(4, times=6), pch=16)
points(c(4)~llbestmods.ne.z$r.squaredGLMM.R2c[which(llbestmods.ne.z$SP.ID=="CWmean")], pch=24,bg=CWbg) # the simple R^2 is stored in R2c for CWmean
points(c(3)~lmabestmods.ne.z$r.squaredGLMM.R2c[which(lmabestmods.ne.z$SP.ID=="CWmean")], pch=24,bg=CWbg)
points(c(2)~nmassbestmods.ne.z$r.squaredGLMM.R2c[which(nmassbestmods.ne.z$SP.ID=="CWmean")], pch=24,bg=CWbg)
points(c(1)~nareabestmods.ne.z$r.squaredGLMM.R2c[which(nareabestmods.ne.z$SP.ID=="CWmean")], pch=24,bg=CWbg)
#SPP means
points(c(4)~llbestmods.ne.z$r.squaredGLMM.R2c[which(llbestmods.ne.z$SP.ID=="SPPmean")], pch=25, bg=SPPbg) # the simple R^2 is stored in R2c for CWmean
points(c(3)~lmabestmods.ne.z$r.squaredGLMM.R2c[which(lmabestmods.ne.z$SP.ID=="SPPmean")], pch=25, bg=SPPbg)
points(c(2)~nmassbestmods.ne.z$r.squaredGLMM.R2c[which(nmassbestmods.ne.z$SP.ID=="SPPmean")], pch=25, bg=SPPbg)
points(c(1)~nareabestmods.ne.z$r.squaredGLMM.R2c[which(nareabestmods.ne.z$SP.ID=="SPPmean")], pch=25, bg=SPPbg)

# add in ecoregion only rsquareds
# points(c(4,3,2,1)~colMeans(ecoregmods[,-1]), pch=4,bg=CWbg) # the simple R^2 is stored in R2c for CWmean
# points(c(4,3,2,1)~CWecoregmods, pch=23,bg=CWbg) # the simple R^2 is stored in R2c for CWmean



axis(3)
axis(2, labels = c(expression(paste(italic("N")[area])),expression(paste(italic("N")[mass])),"LMA","Leaf Life"), at=c(1,2,3,4), las=2)
mtext(side=3,text = expression(paste(italic("R")^2," or marginal ", italic("R")^2)), line=1.7)
legend(x=.675, y=1, legend = "Ind. Spp", fill=boxcol, bty="n",col = boxcol,border = boxcol, cex=.8)
legend(x=.7, y=2.3, legend = c("CW mean","Spp mean"), pt.bg=c(CWbg, SPPbg), bty="n",pch=c(24,25), cex=.8)
# legend(x=.7-.1, y=3+2.5, legend = c("CWM Ecoreg", "Spp Ecoreg"), bty="n", pch=c(23,4), cex=.8)

# legend(x=.75, y=2.5, legend = c("CWM-eco","CWM-full","Spp-eco"), bg=CWbg, bty="n",pch=c(23,24,4), cex=.7)
# legend(x=.73, y=.9, legend = "Spp-env", fill="gray", bty="n",col = boxcol,border = boxcol, cex=.7)
mtext("a)", cex=panlab.cex, side=3, adj=0.02, line=panlab.ln)

dev.off()


#________________________________________________________
#### Panel b): model averaged std effect sizes

#jpeg(width=7.008, height=4, units = "in", res=600,filename = paste0("./",results_dirname,"/Fig4cd_LL_Narea.jpeg"))
pdf(width=3.33, height=6.5, file = paste0("./",results_dirname,"/Fig5b_Trait-Env_effectsizes_color4_v2.pdf"))

par(mfcol=c(6,1), mar=c(0,3.5,0,0), oma=c(5,0,2,1), mgp=c(2.5,.75,0), cex=.9)

### Effect Sizes 
boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="climPC1" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
abline(h=0, col="grey")
# replot boxes over grey '0' line
boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="climPC1" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0, yaxt="n", add=T)

points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="climPC1"  & avglong.ne.z$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="climPC1"  & avglong.ne.z$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "b)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "wetness",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="climPC2" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
#plot(value~xvals, avglong.ne.z[which(avglong.ne.z$trait=="LeafLife" & avglong.ne.z$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0, col="grey")
# replot over grey line
boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="climPC2" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0, yaxt="n", add=T)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="climPC2"  & avglong.ne.z$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="climPC2"  & avglong.ne.z$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "c)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "warmth",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="soil_N" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
#plot(value~xvals, avglong.ne.z[which(avglong.ne.z$trait=="LeafLife" & avglong.ne.z$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0, col="grey")
boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="soil_N" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0, yaxt="n", add=T)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="soil_N"  & avglong.ne.z$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="soil_N"  & avglong.ne.z$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "d)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "soil N",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)

boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="Stand_Age" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
#plot(value~xvals, avglong.ne.z[which(avglong.ne.z$trait=="LeafLife" & avglong.ne.z$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0, col="grey")
boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="Stand_Age" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0, yaxt="n", add=T)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="Stand_Age"  & avglong.ne.z$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="Stand_Age"  & avglong.ne.z$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "e)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )
mtext(text = "log(Age)",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
#text(x=7, y=0, pos=3, labels = "NA",cex=.7)
mtext(text="Standardized Effect Sizes",side = 2, adj=-.1, line=2.1)


boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="LAI" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
#plot(value~xvals, avglong.ne.z[which(avglong.ne.z$trait=="LeafLife" & avglong.ne.z$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0, col="grey")
boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="LAI" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0, yaxt="n", add=T)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="LAI"  & avglong.ne.z$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="LAI"  & avglong.ne.z$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "LAI",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
mtext(text = "f)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )


boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="Growth" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0)
#plot(value~xvals, avglong.ne.z[which(avglong.ne.z$trait=="LeafLife" & avglong.ne.z$SP.ID!="CWmean"),], pch=16, xaxt="n", ylab="", ylim=c(-0.7,1), xlab="", xlim=c(0.6,6.4))
abline(h=0, col="grey")
boxplot(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="Growth" & avglong.ne.z$SP.ID!="CWmean" & avglong.ne.z$SP.ID!="SPPmean" ),], at=c(1,2,3,4)
        , ylim=c(-.7,1), xaxt="n", xlim=c(0.5,4.5)
        ,outpch=1, border=boxcol, col=boxcol, outcol=boxcol
        , staplewex=0, notch=F, medcol="black" #, medlwd=0
        , whisklty=1, whisklwd=2, boxwex=boxwex, whiskcol=boxcol #, boxwex=1, medlwd=3
        ,las=2, range=0, yaxt="n",add=T)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="Growth"  & avglong.ne.z$SP.ID=="SPPmean"),], pch=25,cex=ptcex, bg=SPPbg)
points(value~trait, avglong.ne.z[which(avglong.ne.z$variable=="Growth"  & avglong.ne.z$SP.ID=="CWmean"),], pch=24, bg=CWbg, cex=ptcex)
mtext(text = "NPP",side = 3,line = panlab.ln,adj=0.9,cex=panlab.cex)
mtext(text = "g)",side = 3,line = panlab.ln,adj=0.02,cex=panlab.cex )

axis(side = 1,at = c(1,2,3,4), labels = c("Leaf Life", "LMA",expression(paste(italic("N")[mass])), expression(paste(italic("N")[area]))), xpd=NA, las=3, srt=50, xlim=c(0.6,6.4))

dev.off()






