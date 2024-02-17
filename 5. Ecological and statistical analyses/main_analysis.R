######################################################################################################################################################################################################################################
#
# AUTHORS: Elaine Shen
# DATE CREATED: 02/09/2023
# DATE MODIFIED: 02/17/24 (modified for brevity)
# DESCRIPTION: Main statistical/ecological analyses for "Preliminary characterization of coral reef diversity using environmental DNA in a hyper-diverse context"
#
# INPUTS:
#   CO1v4_rdp_raw.rds : Taxonomy-dependent eDNA data (ASVs), assigned taxonomy using 1) the CO1v4 database and RDP algorithm and 2) an iterative blast search for remaining unidentified sequences
#   LULU_curation_otu_97.rds: Taxonomy-independent eDNA data (OTUs)
#   CO1v4_RDP_otus.csv: Taxonomy table for OTUs to remove land contaminant OTUs only
#   MASTER_fish_data_061621.csv: Raw visual survey of fishes data (contact me for access)
#   MASTER_site_data_t_FINAL_SESYNC_20210720.csv: Site metadata (contact me for access)
#   uvc_unique_ab_functions_workingdoc_061621.csv: A "taxonomy" table for visual surveys of fishes that contains functional group assignments (contact me for access).
#   benthic_phylo.RDS: Visual surveys of corals (formatted for phyloseq) (contact me for access)
#
########################################################################################################################################################################################################################################
rm(list=ls())
####### 0. LOAD LIBRARIES #######
library(reshape)
library(dplyr)
library(tidyr)
library(stringr)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(rich)
library(microbiome) # accessing more of phyloseq
library(iNEXT) # rarefaction curve
library(metagMisc) # format for iNEXT, correct for singletons
library(colorspace) # for plotting colors
library(FSA) # post hoc Dunn test (nonparametric multiple comparisons)
library(multcomp) # parametric multiple comparisons
library(vegan) # permanova
library(geosphere) # geographic distances
library(ggordiplots)

####### 1. IMPORT AND PREPROCESS DATA #######
#### 1.0 Site list #### 
# List of sites by region that has all four data types (eDNA water, eDNA sediment, UVC, benthic)
lombok_sites = c("LB_GSU1", "LB_GSU2", "LB_GSU3", "LB_GLAW", "LB_GPET","LB_GKON", "LB_GNAN", "LB_BUNU", "LB_KABU", "LB_GREN", "LB_GIGE", "LB_TRA1", "LB_TRA2")
misool_sites = c("MI_BARR", "MI_WABN", "MI_WRSM", "MI_WKES", "MI_BOOE", "MI_YILL", "MI_FEME","MI_NECO", "MI_WAFS", "MI_YEKW", "MI_WAGW", "MI_DARS", "MI_CLRF") 
waigeo_sites = c("WG_YENB", "WG_CKRI", "WG_ARBO","WG_SAWA", "WG_YEB1", "WG_MISK", "WG_FRIW","WG_GGOF","WG_YANS", "WG_AREF", "WG_TIMI", "WG_RUVA","WG_PIAY", "WG_PANI","WG_YENW")

site_list = c(lombok_sites, misool_sites, waigeo_sites)

#### 1.1 Environmental DNA data ####
## 1.1.0 Import ##
merged_phyloseq = readRDS("/Users/elaineshen/Desktop/Indo_eDNA/genetic_data/CO1v4_rdp_raw.rds")

# Remove resequenced samples and Wakatobi samples
merged_phyloseq = subset_samples(merged_phyloseq, region_code%in%site_list)
merged_phyloseq = subset_samples(merged_phyloseq, sample_names(merged_phyloseq) != "EB826")
merged_phyloseq = subset_samples(merged_phyloseq, sample_names(merged_phyloseq) != "EB136")

# Only look at 0.4 um filter data 
merged_phyloseq = subset_samples(merged_phyloseq, filter == "0.4um")
merged_phyloseq = prune_taxa(taxa_sums(merged_phyloseq) > 0, merged_phyloseq)

# Take out land-based families (flies, cockroaches, humans, cows, chickens, etc) or misidentifications (sticklebacks, river stingrays, freshwater eels) from analysis
land_to_remove = c("Sciaridae", "Accipitridae", "Blattidae", "Bovidae", "Culicidae", "Formicidae","Fundulidae","Gasterosteidae","Hominidae","Liposcelididae", "Phasianidae","Phoridae","Potamotrygonidae")
genera_to_remove = c("Anguilla") # this is a bony fish, but a freshwater eel 

merged_phyloseq = subset_taxa(merged_phyloseq, !family_name %in% land_to_remove)
merged_phyloseq = subset_taxa(merged_phyloseq, !genus_name %in% genera_to_remove)
merged_phyloseq = prune_taxa(taxa_sums(merged_phyloseq) > 0, merged_phyloseq) # see how many ASVs are in samples that are not used in analysis

## 1.1.1 Subset to 97% taxonomic IDs, separate into sample type (water, sediment), only use 0.4 um ##
merged_phyloseq_97 = subset_taxa(merged_phyloseq, confidence>=0.97)

## Correct taxonomic table - iterative BLAST results returned more incomplete taxonomic assignments at various taxonomic ranks
# First modify taxonomy table so all undefined phyla are NA
merged_phyloseq_97_tax = as.data.frame(tax_table(merged_phyloseq_97))
merged_phyloseq_97_tax$phylum_name = replace_na(merged_phyloseq_97_tax$phylum_name, "undef_Eukaryota")
merged_phyloseq_97_tax$phylum_name[merged_phyloseq_97_tax$phylum_name == "" |merged_phyloseq_97_tax$phylum_name == "unknown" | merged_phyloseq_97_tax$phylum_name == "undef_undef_Eukaryota"| merged_phyloseq_97_tax$phylum_name == "undef_Eukaryota"] <- "undef_Eukaryota"

# Make sure that all unidentified Bacteria assignments have a phylum specification of "undef_Bacteria"
merged_phyloseq_97_tax$phylum_name[merged_phyloseq_97_tax$id == "undef_Bacteria"] = "undef_Bacteria"
# Phyla Apicomplexa, Bacillariophyta, Cryptista, should be in kingdom Chromista
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$phylum_name == "Apicomplexa"] = "Chromista"
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$phylum_name == "Bacillariophyta"] = "Chromista"
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$phylum_name == "Cryptista"] = "Chromista"
# Phyla Arthropoda, Chordata, Cnidaria, Mollusca should be in kingdom Metazoa
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$phylum_name == "Arthropoda"] = "Metazoa"
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$phylum_name == "Chordata"] = "Metazoa"
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$phylum_name == "Cnidaria"] = "Metazoa"
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$phylum_name == "Mollusca"] = "Metazoa"
# Phyla Chlorophyta should be in kingdom Plantae
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$phylum_name == "Chlorophyta"] = "Plantae"
# Phyla Rhodophyta should be in kingdom Protista
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$phylum_name == "Rhodophyta"] = "Protista"
# Class Dinophyceae should be in kingdom Chromista and phylum Dinoflagellata
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$class_name == "Dinophyceae"] = "Chromista"
merged_phyloseq_97_tax$phylum_name[merged_phyloseq_97_tax$class_name == "Dinophyceae"] = "Dinoflagellata"
# Class Phaeophyceae, Raphidophyceae, Oomycetes should be in Kingdom Chromista and phylum Gyrista
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$class_name == "Phaeophyceae"] = "Chromista"
merged_phyloseq_97_tax$phylum_name[merged_phyloseq_97_tax$class_name == "Phaeophyceae"] = "Gyrista"
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$class_name == "Raphidophyceae"] = "Chromista"
merged_phyloseq_97_tax$phylum_name[merged_phyloseq_97_tax$class_name == "Raphidophyceae"] = "Gyrista"
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$class_name == "Oomycetes"] = "Chromista"
merged_phyloseq_97_tax$phylum_name[merged_phyloseq_97_tax$class_name == "Oomycetes"] = "Gyrista"
# family Paramoebidae should be in and kingdom Protista and phylum Amoebozoa
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$family_name == "Paramoebidae"] = "Protista"
merged_phyloseq_97_tax$phylum_name[merged_phyloseq_97_tax$family_name == "Paramoebidae"] = "Amoebozoa"
# id Cunea_thuwala should be in kingdom Protista and phylum Amoebozoa
merged_phyloseq_97_tax$kingdom_name[merged_phyloseq_97_tax$id == "Cunea_thuwala"] = "Protista"
merged_phyloseq_97_tax$phylum_name[merged_phyloseq_97_tax$id == "Cunea_thuwala"] = "Amoebozoa"

# Input this corrected tax table as the new tax table
tax_table(merged_phyloseq_97) = tax_table(as.matrix(merged_phyloseq_97_tax))

# Animals and non-animals, these should be inverses of each other
merged_phyloseq_97_animal = subset_taxa(merged_phyloseq_97, kingdom_name == "Metazoa")
merged_phyloseq_97_nonanimal = subset_taxa(merged_phyloseq_97, kingdom_name != "Metazoa" | is.na(kingdom_name))

# Separate by sample type
merged_phyloseq_sed = subset_samples(merged_phyloseq_97, sampleType =="sediment")
merged_phyloseq_water = subset_samples(merged_phyloseq_97, sampleType =="water")

merged_phyloseq_97_animal_w = subset_samples(merged_phyloseq_97_animal, sampleType == "water")
merged_phyloseq_97_animal_s = subset_samples(merged_phyloseq_97_animal, sampleType == "sediment")

merged_phyloseq_97_nonanimal_w = subset_samples(merged_phyloseq_97_nonanimal, sampleType == "water")
merged_phyloseq_97_nonanimal_s = subset_samples(merged_phyloseq_97_nonanimal, sampleType == "sediment")

## 1.1.3 Import OTU data (for taxonomy-independent analyses) ## 
otu_tab = readRDS("/Users/elaineshen/Desktop/cox1_raw/combined/lulu_vsearch/LULU_curation_otu_97.rds")
otu_tab = otu_tab$curated_table

# Taxonomy table (to remove contaminants, land-families, etc)
otu_tax = as.matrix(read.csv("/Users/elaineshen/Desktop/cox1_raw/combined/lulu_vsearch/CO1v4_RDP_otus.csv", row.names = 1))
otu_tax = tax_table(otu_tax)

# Create phyloseq file 
otu_full = phyloseq(otu_table(otu_tab,taxa_are_rows = T), sample_data(merged_phyloseq), otu_tax)

# Remove contamination
land_to_remove = c("Sciaridae", "Accipitridae", "Blattidae", "Bovidae", "Culicidae", "Formicidae","Fundulidae","Gasterosteidae","Hominidae","Liposcelididae", "Phasianidae","Phoridae","Potamotrygonidae")
genera_to_remove = c("Anguilla") # freshwater eel
otu_full = subset_taxa(otu_full, !family_name %in% land_to_remove)
otu_full = subset_taxa(otu_full, !genus_name %in% genera_to_remove)

# Subset to filter size and sites used
otu_full = subset_samples(otu_full, region_code%in%site_list)
otu_full = subset_samples(otu_full, filter == "0.4um")

# Get rid of OTUs that have 0 reads assigned to an OTU 
otu_full = prune_taxa(taxa_sums(otu_full) > 0, otu_full)

otu_water = subset_samples(otu_full, sampleType =="water")
otu_water = subset_samples(otu_water, filter=="0.4um")

otu_sediment = subset_samples(otu_full, sampleType == "sediment")
otu_sediment = subset_samples(otu_sediment, sample_names(otu_sediment) != "EB136") # remove sample that was resequenced later

#### 1.2 Underwater Visual Census data ####

## 1.2.1 Import UVC data ## 
uvc_raw = read.csv("/Users/elaineshen/Desktop/Indo_eDNA/UVC_Data/MASTER_fish_data_061621.csv")
site_sesync = read.csv("/Users/elaineshen/Desktop/Indo_eDNA/site_data/MASTER_site_data_t_FINAL_SESYNC_20210720.csv")
tax_table = read.csv("/Users/elaineshen/Desktop/Indo_eDNA/UVC_Data/archive_uvc/uvc_unique_ab_functions_workingdoc_061621.csv")

## 1.2.2. Make into OTU table ## 
# Total (sum of all transects) - would need to make a new metadata file with site_code_t to use
uvc_site.species <- as.data.frame(reshape::cast(uvc_raw, site_code_t ~ species, value='abundance',sum))

# set site_code_1 as rownames to format into otu table
uvc_site.species <- data.frame(uvc_site.species[,-1], row.names = uvc_site.species[,1])
colnames_uvc_site.species = colnames(uvc_site.species)
colnames_uvc_site.species = gsub("."," ",colnames_uvc_site.species, fixed=TRUE)
colnames(uvc_site.species) = colnames_uvc_site.species

## 1.2.4 Create taxonomy table ## 
rownames(tax_table) = tax_table$species
tax_table = as.matrix(tax_table)
tax_table = tax_table[,-1]

## 1.2.5 Modify site metadata table ## 

rownames(site_sesync) = site_sesync$site_code_t
site_sesync$region = str_to_title(site_sesync$region)

## 1.2.6 Import into phyloseq ##
otu_table = otu_table(uvc_site.species, taxa_are_rows = FALSE)
tax_table = tax_table(tax_table)
sample_data = sample_data(site_sesync)
uvc_phylo = phyloseq(otu_table, tax_table, sample_data)

## 1.2.7 Subset to global site_list 
uvc_phylo = subset_samples(uvc_phylo, region_code%in%site_list)

## 1.2.8 Make one that's merged by region_code
uvc_phylo1 = merge_samples(uvc_phylo, "region_code")

# Put back region value
region_site_code = as.data.frame(sample_data(uvc_phylo))[,c("region", "region_code", "site_name","date")]
rownames(region_site_code) = NULL
region_site_code= unique(region_site_code, by="region_code")

rownames(region_site_code) = region_site_code$region_code
sample_data(uvc_phylo1)$region = sort(region_site_code$region)
sample_data(uvc_phylo1)$region_code = sort(region_site_code$region_code)
sample_data(uvc_phylo1)$site_name = region_site_code$site_name
sample_data(uvc_phylo1)$date = region_site_code$date

uvc_phylo = subset_samples(uvc_phylo, region != "Wakatobi")
uvc_phylo = prune_taxa(taxa_sums(uvc_phylo) > 0, uvc_phylo)
uvc_phylo1 = subset_samples(uvc_phylo1, region != "Wakatobi")
uvc_phylo1 = prune_taxa(taxa_sums(uvc_phylo1) > 0, uvc_phylo1)

#### 1.3 Benthic Survey data ####
benthic_phylo = readRDS("/Users/elaineshen/Desktop/Indo_eDNA/benthic_data/benthic_phylo.RDS")
benthic_phylo = subset_samples(benthic_phylo,site_code_1%in%site_list )
benthic_phylo = subset_samples(benthic_phylo, region != "wakatobi")

## 1.3.1 Make version that's merged by site_code_1, add region back into metadata
benthic_phylo1 = merge_samples(benthic_phylo, "site_code_1")

benthic_region_site = as.data.frame(sample_data(benthic_phylo))[,c("region","site_code_1", "site_code","site_name")]
rownames(benthic_region_site) = NULL
benthic_region_site = unique(benthic_region_site, by="site_code_1")
sample_data(benthic_phylo1)$region = sort(benthic_region_site$region)
sample_data(benthic_phylo1)$region_code = sort(benthic_region_site$site_code_1)
sample_data(benthic_phylo1)$site_code = sort(benthic_region_site$site_code)
sample_data(benthic_phylo1)$site_name = benthic_region_site$site_name

benthic_phylo1_sam = as.data.frame(as.matrix(sample_data(benthic_phylo1)))[,-c(4:7)]
benthic_phylo1_sam$region = str_to_title(benthic_phylo1_sam$region)
uvc_phylo1_sam = as.data.frame(as.matrix(sample_data(uvc_phylo1)))[,-c(4:6)]

benthic_phylo1_sam_new = left_join(benthic_phylo1_sam, uvc_phylo1_sam, by="region_code")
rownames(benthic_phylo1_sam_new) = benthic_phylo1_sam_new$region_code

sample_data(benthic_phylo1) = benthic_phylo1_sam_new


####### 2. SUMMARY TABLE  #######
#### 2.1 Environmental DNA ####

## 2.1.0 Look at raw ASV data (no need to have 97% taxon ID)
merged_phyloseq_water1 = subset_samples(merged_phyloseq, sampleType == "water")
merged_phyloseq_water1 = subset_samples(merged_phyloseq_water1, filter == "0.4um")

merged_phyloseq_sed1 = subset_samples(merged_phyloseq, sampleType == "sediment")

## 2.1.1 Water samples

# Lombok
merged_ps_water_lbk = subset_samples(merged_phyloseq_water1, region == "Lombok")
sum(sample_sums(merged_ps_water_lbk)) 
nsamples(merged_ps_water_lbk) 
unique(sample_data(merged_ps_water_lbk)$region_code)

merged_ps_water_lbk_merge <- merge_samples(merged_ps_water_lbk , "region")
richness(merged_ps_water_lbk_merge, index="observed") 

otu_water_lbk = subset_samples(otu_water, region == "Lombok")
sum(sample_sums(otu_water_lbk)) 
nsamples(otu_water_lbk)

otu_water_lbk_merge = merge_samples(otu_water_lbk, "region")
richness(otu_water_lbk_merge, index = "observed")

# Misool
merged_ps_water_mi = subset_samples(merged_phyloseq_water1, region == "Misool")
sum(sample_sums(merged_ps_water_mi)) 
nsamples(merged_ps_water_mi)
unique(sample_data(merged_ps_water_mi)$region_code)

merged_ps_water_mi_merge <- merge_samples(merged_ps_water_mi , "region")
richness(merged_ps_water_mi_merge, index="observed")

otu_water_mi = subset_samples(otu_water, region == "Misool")
sum(sample_sums(otu_water_mi)) 
nsamples(otu_water_mi)

otu_water_mi_merge = merge_samples(otu_water_mi, "region")
richness(otu_water_mi_merge, index = "observed")

# Waigeo
merged_ps_water_wg = subset_samples(merged_phyloseq_water1, region == "Waigeo")
sum(sample_sums(merged_ps_water_wg)) 
nsamples(merged_ps_water_wg)
unique(sample_data(merged_ps_water_wg)$region_code)

merged_ps_water_wg_merge <- merge_samples(merged_ps_water_wg , "region")
richness(merged_ps_water_wg_merge, index="observed")

otu_water_wg = subset_samples(otu_water, region == "Waigeo")
sum(sample_sums(otu_water_wg))
nsamples(otu_water_wg)

otu_water_wg_merge = merge_samples(otu_water_wg, "region")
richness(otu_water_wg_merge, index = "observed")

## 2.1.2 Sediment samples

# Lombok
merged_ps_sed_lbk = subset_samples(merged_phyloseq_sed1, region == "Lombok")
sum(sample_sums(merged_ps_sed_lbk)) 
nsamples(merged_ps_sed_lbk)
unique(sample_data(merged_ps_sed_lbk)$region_code)

merged_ps_sed_lbk_merge <- merge_samples(merged_ps_sed_lbk , "region")
richness(merged_ps_sed_lbk_merge, index="observed") 

otu_sed_lbk = subset_samples(otu_sediment, region == "Lombok")
sum(sample_sums(otu_sed_lbk)) 
nsamples(otu_sed_lbk)

otu_sed_lbk_merge = merge_samples(otu_sed_lbk, "region")
richness(otu_sed_lbk_merge, index = "observed") 

# Misool
merged_ps_sed_mi = subset_samples(merged_phyloseq_sed1, region == "Misool")
sum(sample_sums(merged_ps_sed_mi)) 
nsamples(merged_ps_sed_mi)
unique(sample_data(merged_ps_sed_mi)$region_code)

merged_ps_sed_mi_merge <- merge_samples(merged_ps_sed_mi , "region")
richness(merged_ps_sed_mi_merge, index="observed") 

otu_sed_mi = subset_samples(otu_sediment, region == "Misool")
sum(sample_sums(otu_sed_mi)) 
nsamples(otu_sed_mi)

otu_sed_mi_merge = merge_samples(otu_sed_mi, "region")
richness(otu_sed_mi_merge, index = "observed") 

# Waigeo
merged_ps_sed_wg = subset_samples(merged_phyloseq_sed1, region == "Waigeo")
sum(sample_sums(merged_ps_sed_wg)) 
nsamples(merged_ps_sed_wg)
unique(sample_data(merged_ps_sed_wg)$region_code)

merged_ps_sed_wg_merge <- merge_samples(merged_ps_sed_wg , "region")
richness(merged_ps_sed_wg_merge, index="observed") 

otu_sed_wg = subset_samples(otu_sediment, region == "Waigeo")
otu_sed_wg = prune_samples(sample_sums(otu_sed_wg)>0, otu_sed_wg)
sum(sample_sums(otu_sed_wg)) 
nsamples(otu_sed_wg)

otu_sed_wg_merge = merge_samples(otu_sed_wg, "region")
richness(otu_sed_wg_merge, index = "observed") 

## 2.1.3 Total OTUs ## 
otu_full_merge = merge_samples(otu_full, "primer")
richness(otu_full_merge, "observed") 

## 2.1.4 Water samples - metazoans only

# Lombok
merged_ps_water_animal_lbk = subset_samples(merged_phyloseq_97_animal_w, region == "Lombok")
sum(sample_sums(merged_ps_water_animal_lbk)) 
nsamples(merged_ps_water_animal_lbk)
unique(sample_data(merged_ps_water_animal_lbk)$region_code)

merged_ps_water_animal_lbk_merge <- merge_samples(merged_ps_water_animal_lbk , "region")
richness(merged_ps_water_animal_lbk_merge, index="observed") 

# Misool
merged_ps_water_animal_mi = subset_samples(merged_phyloseq_97_animal_w, region == "Misool")
sum(sample_sums(merged_ps_water_animal_mi)) 
nsamples(merged_ps_water_animal_mi)
unique(sample_data(merged_ps_water_animal_mi)$region_code)

merged_ps_water_animal_mi_merge <- merge_samples(merged_ps_water_animal_mi , "region")
richness(merged_ps_water_animal_mi_merge, index="observed")

# Waigeo
merged_ps_water_animal_wg = subset_samples(merged_phyloseq_97_animal_w, region == "Waigeo")
sum(sample_sums(merged_ps_water_animal_wg))
nsamples(merged_ps_water_animal_wg)
unique(sample_data(merged_ps_water_animal_wg)$region_code)

merged_ps_water_animal_wg_merge <- merge_samples(merged_ps_water_animal_wg , "region")
richness(merged_ps_water_animal_wg_merge, index="observed")

## 2.1.7 Sediment samples - metazoans only

# Lombok
merged_ps_sed_animal_lbk = subset_samples(merged_phyloseq_97_animal_s, region == "Lombok")
sum(sample_sums(merged_ps_sed_animal_lbk))
nsamples(merged_ps_sed_animal_lbk)
unique(sample_data(merged_ps_sed_animal_lbk)$region_code)

merged_ps_sed_animal_lbk_merge <- merge_samples(merged_ps_sed_animal_lbk , "region")
richness(merged_ps_sed_animal_lbk_merge, index="observed")

# Misool
merged_ps_sed_animal_mi = subset_samples(merged_phyloseq_97_animal_s, region == "Misool")
sum(sample_sums(merged_ps_sed_animal_mi))
nsamples(merged_ps_sed_animal_mi)
unique(sample_data(merged_ps_sed_animal_mi)$region_code)

merged_ps_sed_animal_mi_merge <- merge_samples(merged_ps_sed_animal_mi , "region")
richness(merged_ps_sed_animal_mi_merge, index="observed")

# Waigeo
merged_ps_sed_animal_wg = subset_samples(merged_phyloseq_97_animal_s, region == "Waigeo")
sum(sample_sums(merged_ps_sed_animal_wg)) 
nsamples(merged_ps_sed_animal_wg)
unique(sample_data(merged_ps_sed_animal_wg)$region_code)

merged_ps_sed_animal_wg_merge <- merge_samples(merged_ps_sed_animal_wg , "region")
richness(merged_ps_sed_animal_wg_merge, index="observed") 

## 2.1.8 Water samples - non-metazoans only

# Lombok
merged_ps_water_nonanimal_lbk = subset_samples(merged_phyloseq_97_nonanimal_w, region == "Lombok")
sum(sample_sums(merged_ps_water_nonanimal_lbk))
nsamples(merged_ps_water_nonanimal_lbk)
unique(sample_data(merged_ps_water_nonanimal_lbk)$region_code)

merged_ps_water_nonanimal_lbk_merge <- merge_samples(merged_ps_water_nonanimal_lbk , "region")
richness(merged_ps_water_nonanimal_lbk_merge, index="observed")

# Misool
merged_ps_water_nonanimal_mi = subset_samples(merged_phyloseq_97_nonanimal_w, region == "Misool")
sum(sample_sums(merged_ps_water_nonanimal_mi))
nsamples(merged_ps_water_nonanimal_mi)
unique(sample_data(merged_ps_water_nonanimal_mi)$region_code)

merged_ps_water_nonanimal_mi_merge <- merge_samples(merged_ps_water_nonanimal_mi , "region")
richness(merged_ps_water_nonanimal_mi_merge, index="observed")

# Waigeo
merged_ps_water_nonanimal_wg = subset_samples(merged_phyloseq_97_nonanimal_w, region == "Waigeo")
sum(sample_sums(merged_ps_water_nonanimal_wg)) 
nsamples(merged_ps_water_nonanimal_wg)
unique(sample_data(merged_ps_water_nonanimal_wg)$region_code)

merged_ps_water_nonanimal_wg_merge <- merge_samples(merged_ps_water_nonanimal_wg , "region")
richness(merged_ps_water_nonanimal_wg_merge, index="observed")

## 2.1.9 Sediment samples - nonmetazoans only

# Lombok
merged_ps_sed_nonanimal_lbk = subset_samples(merged_phyloseq_97_nonanimal_s, region == "Lombok")
sum(sample_sums(merged_ps_sed_nonanimal_lbk)) 
nsamples(merged_ps_sed_nonanimal_lbk)
unique(sample_data(merged_ps_sed_nonanimal_lbk)$region_code)

merged_ps_sed_nonanimal_lbk_merge <- merge_samples(merged_ps_sed_nonanimal_lbk , "region")
richness(merged_ps_sed_nonanimal_lbk_merge, index="observed") 

# Misool
merged_ps_sed_nonanimal_mi = subset_samples(merged_phyloseq_97_nonanimal_s, region == "Misool")
sum(sample_sums(merged_ps_sed_nonanimal_mi))
nsamples(merged_ps_sed_nonanimal_mi)
unique(sample_data(merged_ps_sed_nonanimal_mi)$region_code)

merged_ps_sed_nonanimal_mi_merge <- merge_samples(merged_ps_sed_nonanimal_mi , "region")
richness(merged_ps_sed_nonanimal_mi_merge, index="observed") 

# Waigeo
merged_ps_sed_nonanimal_wg = subset_samples(merged_phyloseq_97_nonanimal_s, region == "Waigeo")
sum(sample_sums(merged_ps_sed_nonanimal_wg)) 
nsamples(merged_ps_sed_nonanimal_wg)
unique(sample_data(merged_ps_sed_nonanimal_wg)$region_code)

merged_ps_sed_nonanimal_wg_merge <- merge_samples(merged_ps_sed_nonanimal_wg , "region")
richness(merged_ps_sed_nonanimal_wg_merge, index="observed")

#### 2.2 Underwater Visual Census ####

## 2.2.1 Fill in summary table ## 
# Lombok
uvc_phylo_lbk = subset_samples(uvc_phylo, region == "Lombok")
nsamples(uvc_phylo_lbk)
unique(sample_data(uvc_phylo_lbk)$region_code)
sum(sample_sums(uvc_phylo_lbk))

uvc_phylo_lbk_merge = merge_samples(uvc_phylo_lbk, "region")
richness(uvc_phylo_lbk_merge, index = "observed")

# Misool
uvc_phylo_mis = subset_samples(uvc_phylo, region == "Misool")
uvc_phylo_mis_otu = as.data.frame(otu_table(uvc_phylo_mis))

nsamples(uvc_phylo_mis)
unique(sample_data(uvc_phylo_mis)$region_code)
uvc_phylo_mis = prune_samples(sample_sums(uvc_phylo_mis)>0, uvc_phylo_mis)
sum(sample_sums(uvc_phylo_mis))

uvc_phylo_mis_merge = merge_samples(uvc_phylo_mis, "region")
richness(uvc_phylo_mis_merge, index = "observed")

# Waigeo
uvc_phylo_wg = subset_samples(uvc_phylo, region == "Waigeo")
nsamples(uvc_phylo_wg)
unique(sample_data(uvc_phylo_wg)$region_code)
sum(sample_sums(uvc_phylo_wg))

uvc_phylo_wg_merge = merge_samples(uvc_phylo_wg, "region")
richness(uvc_phylo_wg_merge, index = "observed")

## 2.2.2 Get total fish richness ##

sample_data(uvc_phylo1)$merge = "same"
uvc_phylo2 = prune_samples(sample_sums(uvc_phylo1)>1, uvc_phylo1)
uvc_phylo_merge = merge_samples(uvc_phylo2, "merge") # all same value
richness(uvc_phylo_merge, "observed")

#### 2.3 Benthic Data ####

## 2.3.1 Fill in summary table ## 
# Lombok
benthic_phylo_lbk = subset_samples(benthic_phylo1, region == "Lombok")
nsamples(benthic_phylo_lbk)
unique(sample_data(benthic_phylo_lbk)$region_code)
sum(sample_sums(benthic_phylo_lbk))

benthic_phylo_lbk_merge = merge_samples(benthic_phylo_lbk, "region")
richness(benthic_phylo_lbk_merge, index = "observed")

# Misool
benthic_phylo_mis = subset_samples(benthic_phylo1, region == "Misool")
benthic_phylo_mis_otu = as.data.frame(otu_table(benthic_phylo_mis))
nsamples(benthic_phylo_mis)
unique(sample_data(benthic_phylo_mis)$region_code)
sum(sample_sums(benthic_phylo_mis))

benthic_phylo_mis_merge = merge_samples(benthic_phylo_mis, "region")
richness(benthic_phylo_mis_merge, index = "observed")

# Waigeo
benthic_phylo_wg = subset_samples(benthic_phylo1, region == "Waigeo")
nsamples(benthic_phylo_wg)
unique(sample_data(benthic_phylo_wg)$region_code)
sum(sample_sums(benthic_phylo_wg))

benthic_phylo_wg_merge = merge_samples(benthic_phylo_wg, "region")
richness(benthic_phylo_wg_merge, index = "observed")

## 2.3.2 Get total genus-level richness ## 
sample_data(benthic_phylo)$merge = "same"
benthic_phylo_merge = merge_samples(benthic_phylo1, "merge") # all same value
richness(benthic_phylo_merge, "observed")

####### 3. RAREFACTION CURVES #######
#### 3.0 Set region color value #### 
color_table_region = c("#CC3311", "#009988", "#0077BB")
names(color_table_region) = c("Lombok", "Misool", "Waigeo")                       
color_scale_region = scale_color_manual(name = "Region",values = color_table_region)

#### 3.1 Environmental DNA data ####

## 3.1.1 Water samples ## 
otu_misool_water = subset_samples(otu_water, region=="Misool") %>% otu_table()
otu_misool_water = otu_misool_water@.Data
otu_misool_water = as.matrix((otu_misool_water>0)+0)

otu_waigeo_water = subset_samples(otu_water, region=="Waigeo") %>% otu_table()
otu_waigeo_water = otu_waigeo_water@.Data # get only the species x site matrix
otu_waigeo_water = as.matrix((otu_waigeo_water>0)+0) # Make into incidence matrix

otu_lombok_water = subset_samples(otu_water, region=="Lombok") %>% otu_table()
otu_lombok_water = otu_lombok_water@.Data # get only the species x site matrix
otu_lombok_water = as.matrix((otu_lombok_water>0)+0) # Make into incidence matrix

# Make a list of regional matrices to preserve sample-level data in analysis
otu_region_mat_water = list(Misool = otu_misool_water,
                            Waigeo = otu_waigeo_water,
                            Lombok = otu_lombok_water)

otu_region_mat_inext1 = iNEXT(otu_region_mat_water, datatype="incidence_raw")
otu_rarecurve = ggiNEXT(otu_region_mat_inext1, type=1, se=T, facet.var="none", color.var="site", grey=FALSE) +scale_shape_manual(values = c(5,6,7))+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB"))+
  color_scale_region  + scale_fill_manual(values=c("Lombok" = "#CC3311", "Misool" = "#009988","Waigeo"="#0077BB"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"), 
        legend.background = element_rect(fill="transparent", color= NA),
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+
  ggtitle("Water eDNA sampling effort")  

otu_rarecurve_1 = ggiNEXT(otu_region_mat_inext1, type=2, se=T, facet.var="none", color.var="site", grey=FALSE) +scale_shape_manual(values = c(5,6,7))+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB"))+
  color_scale_region  + scale_fill_manual(values=c("Lombok" = "#CC3311", "Misool" = "#009988","Waigeo"="#0077BB"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"), 
        legend.background = element_rect(fill="transparent", color= NA),
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+
  ggtitle("Water eDNA sampling effort")  

#otu_region_mat_inext1_sc = iNEXT(otu_region_mat_water, datatype="incidence_raw", knots = 200, endpoint = 200) # See when sampling coverage reaches around 100%

## 3.1.2 Sediment samples ## 
otu_misool_sed = subset_samples(otu_sediment, region=="Misool") %>% otu_table()
otu_misool_sed = otu_misool_sed@.Data
otu_misool_sed = as.matrix((otu_misool_sed>0)+0)

otu_waigeo_sed = subset_samples(otu_sediment, region=="Waigeo") %>% otu_table()
otu_waigeo_sed = otu_waigeo_sed@.Data # get only the species x site matrix
otu_waigeo_sed = as.matrix((otu_waigeo_sed>0)+0) # Make into incidence matrix

otu_lombok_sed = subset_samples(otu_sediment, region=="Lombok") %>% otu_table()
otu_lombok_sed = otu_lombok_sed@.Data # get only the species x site matrix
otu_lombok_sed = as.matrix((otu_lombok_sed>0)+0) # Make into incidence matrix

# Make a list of regional matrices to preserve sample-level data in analysis
otu_region_mat_sediment = list(Misool = otu_misool_sed,
                               Waigeo = otu_waigeo_sed,
                               Lombok = otu_lombok_sed)

otu_region_sed_inext1 = iNEXT(otu_region_mat_sediment, datatype="incidence_raw") 
edna_sed_rarecurve = ggiNEXT(otu_region_sed_inext1, type=1, se=T, facet.var="none", color.var="site", grey=FALSE)+scale_shape_manual(values = c(5,6,7))+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" ))+
  scale_fill_manual(values =c("#CC3311", "#009988", "#0077BB" ))+
  color_scale_region +  scale_fill_manual(values=c("Lombok" = "#CC3311", "Misool" = "#009988","Waigeo"="#0077BB"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA),axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+
  ggtitle("Sediment eDNA sampling effort")  

edna_sed_rarecurve_1 =  ggiNEXT(otu_region_sed_inext1, type=2, se=T, facet.var="none", color.var="site", grey=FALSE)+scale_shape_manual(values = c(5,6,7))+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" ))+
  scale_fill_manual(values =c("#CC3311", "#009988", "#0077BB" ))+
  color_scale_region +  scale_fill_manual(values=c("Lombok" = "#CC3311", "Misool" = "#009988","Waigeo"="#0077BB"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA),axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+
  ggtitle("Sediment eDNA sampling effort")

#otu_region_sed_inext1_sc = iNEXT(otu_region_mat_sediment, datatype="incidence_raw", knots = 500, endpoint = 500) # See when sampling coverage reaches around 100%

#### 3.2 Underwater Visual Census data ####
uvc_phylo = prune_samples(sample_sums(uvc_phylo)>0, uvc_phylo)

uvc_misool = subset_samples(uvc_phylo, region =="Misool") %>% otu_table()
uvc_misool = uvc_misool@.Data
uvc_misool = as.matrix((uvc_misool>0)+0) %>% t()

uvc_waigeo = subset_samples(uvc_phylo, region=="Waigeo") %>% otu_table()
uvc_waigeo = uvc_waigeo@.Data # get only the species x site matrix
uvc_waigeo = as.matrix((uvc_waigeo>0)+0) %>% t() # Make into incidence matrix

uvc_lombok = subset_samples(uvc_phylo, region=="Lombok") %>% otu_table()
uvc_lombok = uvc_lombok@.Data # get only the species x site matrix
uvc_lombok = as.matrix((uvc_lombok>0)+0) %>% t() # Make into incidence matrix

# Make a list of regional matrices to preserve sample-level data in analysis
uvc_region_mat = list(Misool = uvc_misool,
                      Waigeo = uvc_waigeo,
                      Lombok = uvc_lombok)

uvc_region_mat_inext1 = iNEXT(uvc_region_mat, datatype="incidence_raw")
uvc_rarecurve = ggiNEXT(uvc_region_mat_inext1, type=1, se=T, facet.var="none", color.var="site", grey=FALSE) +scale_shape_manual(values = c(5,6,7))+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" ))+
  color_scale_region +  scale_fill_manual(values=c("Lombok" = "#CC3311", "Misool" = "#009988","Waigeo"="#0077BB"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"),
                                                                                                                          legend.background = element_rect(fill = "transparent", color = NA),
                                                                                                                          legend.box.background = element_rect(fill = "transparent", color=NA),
                                                                                                                          legend.key = element_rect(fill = "transparent", color = NA))+
  ggtitle("UVC sampling effort")  

uvc_region_inext_sc = iNEXT(uvc_region_mat, datatype="incidence_raw", knots = 200, endpoint = 200) 

#### 3.3 Benthic Survey data #### 
benthic_misool = subset_samples(benthic_phylo, region =="misool") %>% otu_table()
benthic_misool = benthic_misool@.Data
benthic_misool = as.matrix((benthic_misool>0)+0) %>% t()

benthic_waigeo = subset_samples(benthic_phylo, region=="waigeo") %>% otu_table()
benthic_waigeo = benthic_waigeo@.Data # get only the species x site matrix
benthic_waigeo = as.matrix((benthic_waigeo>0)+0) %>% t() # Make into incidence matrix

benthic_lombok = subset_samples(benthic_phylo, region=="lombok") %>% otu_table()
benthic_lombok = benthic_lombok@.Data # get only the species x site matrix
benthic_lombok = as.matrix((benthic_lombok>0)+0) %>% t() # Make into incidence matrix

# Make a list of regional matrices to preserve sample-level data in analysis
benthic_region_mat = list(Misool = benthic_misool,
                          Waigeo = benthic_waigeo,
                          Lombok = benthic_lombok)

benthic_region_mat_inext1 = iNEXT(benthic_region_mat, datatype="incidence_raw")
benthic_rarecurve = ggiNEXT(benthic_region_mat_inext1, type=1, se=T, facet.var="none", color.var="site", grey=FALSE) +scale_shape_manual(values = c(5,6,7))+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" ))+
  color_scale_region +  scale_fill_manual(values=c("Lombok" = "#CC3311", "Misool" = "#009988","Waigeo"="#0077BB"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"), 
                                                                                                                         legend.background = element_rect(fill="transparent", color= NA),
                                                                                                                         legend.box.background = element_rect(fill = "transparent", color=NA),
                                                                                                                         legend.key = element_rect(fill = "transparent",color=NA))+
  ggtitle("Benthic sampling effort, genus-level ID's")  + scale_x_continuous("Number of transect points") +scale_y_continuous("Genus diversity")


#### 3.4 Arrange into one plot #### 
rarecurve_eDNA_combo = ggarrange(otu_rarecurve+ggtitle("Environmental DNA - water"), edna_sed_rarecurve+ggtitle("Environmental DNA - sediment"), ncol=2, labels = c("(a)","(b)"), common.legend = TRUE, legend="right")
ggsave("rarecurve_eDNA_combo_v3.png",rarecurve_eDNA_combo,bg="transparent", width=15, height = 9)

rarecurve_UVC_combo = ggarrange(uvc_rarecurve+ggtitle("Visual surveys of fishes"), benthic_rarecurve+ggtitle("Visual surveys of benthos"), ncol=2, labels = c("(a)","(b)"), common.legend=TRUE, legend="right")
ggsave("rarecurve_UVC_combo_v3.png", bg="transparent",rarecurve_UVC_combo, width=15, height = 9) # this may be useful for SI, but maybe not

####### 4. SETTING COLORS FOR TAXONOMY BARPLOTS #######
#### 4.0 Set consistent fish color palette ####
## 4.0.1 Set consistent colors for fish families in eDNA and UVC by establishing global family color key ## 
ps_97_fish = subset_taxa(merged_phyloseq_97, class_name == "Actinopteri")
merged_phyloseq_fish_tax = as.data.frame(tax_table(ps_97_fish))
taxnames1 = sort(unique(merged_phyloseq_fish_tax$family_name))
uvc_fish = as.data.frame(tax_table(uvc_phylo))
taxnames2 = sort(unique(uvc_fish$family))
taxnames_fish_full = sort(unique(c(taxnames1, taxnames2)))

color_table_fish = qualitative_hcl(56, palette = "Dark 2")
names(color_table_fish)= sort(taxnames_fish_full)

color_scale_fish = scale_fill_manual(name = "Family",values = color_table_fish)

#### 4.1 Environmental DNA data ####

## 4.1.0 Set consistent colors for phyla based on sediment data ##
merged_phyloseq_sed_tax = as.data.frame(tax_table(merged_phyloseq_sed))
color_table = c(qualitative_hcl(20, palette = "Dark 2"),"grey", "darkgrey")
names(color_table)= sort(unique(merged_phyloseq_sed_tax$phylum_name))

color_scale = scale_fill_manual(name = "Phylum",values = color_table)

#### 4.3 Benthic Survey data #### 

## 3.3.0 Set consistent colors for family ##
benthic_phylo_tax = as.data.frame(tax_table(benthic_phylo))
color_table_benthic = qualitative_hcl(21, palette = "Dark 2")
names(color_table_benthic)= sort(unique(benthic_phylo_tax$family))

color_scale_benthic = scale_fill_manual(name = "Family",values = color_table_benthic)


####### 5. TAXONOMY BARPLOTS - RELATIVE ABUNDANCE #######
#### 5.1 Environmental DNA data (full) #### 

## 5.1.1 Sediment samples ## 
sed_97_phyla <- tax_glom(merged_phyloseq_sed, "phylum_name")
sed_97_relabund_0 <- transform_sample_counts(sed_97_phyla, function(x) x / sum(x))
sed_97_relabund_0 = prune_samples(sample_sums(sed_97_relabund_0)==1, sed_97_relabund_0)
sed_97_relabund_1 <- merge_samples(sed_97_relabund_0 , "region")
sed_97_relabund_2 <- transform_sample_counts(sed_97_relabund_1, function(x) x / sum(x))

sed_97_plot1 = plot_bar(sed_97_relabund_2, fill="phylum_name") +
  geom_bar(stat="identity", position="stack")+
  color_scale + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
                      axis.text.x = element_text(angle=0, hjust=0.5, size =11),
                      legend.background = element_rect(fill="transparent", color= NA),
                      legend.box.background = element_rect(fill = "transparent", color=NA),
                      legend.key = element_rect(fill = "transparent",color=NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill=guide_legend(ncol=1))+ ylab("Relative abundance (%)")+
  ggtitle("Environmental DNA - sediment (97% identity)")

# Plot unknowns at top rather than bottom
sed_97_relabund_3 = psmelt(sed_97_relabund_2)
sed_97_relabund_3$phylum_name_1 = factor(sed_97_relabund_3$phylum_name, levels = c("undef_Bacteria", "undef_Eukaryota", "Amoebozoa",  "Annelida",  "Apicomplexa", "Arthropoda", "Ascomycota",     
                                                                                   "Bacillariophyta", "Basidiomycota",   "Chaetognatha",    "Chlorophyta",     "Chordata",       
                                                                                   "Cnidaria","Dinoflagellata", "Echinodermata",   "Gyrista", "Haptista", "Mollusca", "Mucoromycota","Porifera", "Proteobacteria",  "Rhodophyta"))

sed_97_plot = ggplot(sed_97_relabund_3[order(sed_97_relabund_3$phylum_name_1),], aes(x=Sample, y=Abundance, fill = phylum_name_1,))+
  geom_bar(stat="identity")+color_scale + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
                                                axis.text.x = element_text(angle=0, hjust=0.5, size =11),
                                                legend.background = element_rect(fill="transparent", color= NA),
                                                legend.box.background = element_rect(fill = "transparent", color=NA),
                                                legend.key = element_rect(fill = "transparent",color=NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill=guide_legend(ncol=1))+ ylab("Relative abundance (%)")+
  ggtitle("Environmental DNA - sediment (97% identity)")

## 5.1.2 Water samples ## 
ps_97_phyla <- tax_glom(merged_phyloseq_water, "phylum_name")
ps_97_relabund_0 <- transform_sample_counts(ps_97_phyla, function(x) x / sum(x))
ps_97_relabund_0 = prune_samples(sample_sums(ps_97_relabund_0)==1, ps_97_relabund_0)
ps_97_relabund_1 <- merge_samples(ps_97_relabund_0 , "region")
ps_97_relabund_2 <- transform_sample_counts(ps_97_relabund_1, function(x) x / sum(x))

water_97_plot1 = plot_bar(ps_97_relabund_2, fill="phylum_name") +
  geom_bar(stat="identity", position="stack")+
  color_scale + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
                      axis.text.x = element_text(angle=0, hjust=0.5, size =11),
                      legend.background = element_rect(fill="transparent", color= NA),
                      legend.box.background = element_rect(fill = "transparent", color=NA),
                      legend.key = element_rect(fill = "transparent",color=NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill=guide_legend(ncol=1))+ ylab("Relative abundance (%)")+
  ggtitle("Environmental DNA - water (97% identity)")

# Put unknown taxa at the top rather than the bottom by creating a factor for phyla
ps_97_relabund_3 = psmelt(ps_97_relabund_2)
ps_97_relabund_3$phylum_name_1 = factor(ps_97_relabund_3$phylum_name,levels = c("undef_Bacteria", "undef_Eukaryota", "Amoebozoa",  "Annelida",  "Apicomplexa", "Arthropoda", "Ascomycota",     
                                                                                "Bacillariophyta", "Basidiomycota",   "Chaetognatha",    "Chlorophyta",     "Chordata",       
                                                                                "Cnidaria","Dinoflagellata", "Echinodermata",   "Gyrista", "Haptista", "Mollusca", "Mucoromycota","Porifera", "Proteobacteria",  "Rhodophyta"))

water_97_plot = ggplot(ps_97_relabund_3[order(ps_97_relabund_3$phylum_name_1),], aes(x=Sample, y=Abundance, fill = phylum_name_1,))+
  geom_bar(stat="identity")+color_scale + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
                                                axis.text.x = element_text(angle=0, hjust=0.5, size =11),
                                                legend.background = element_rect(fill="transparent", color= NA),
                                                legend.box.background = element_rect(fill = "transparent", color=NA),
                                                legend.key = element_rect(fill = "transparent",color=NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  guides(fill=guide_legend(ncol=1))+ ylab("Relative abundance (%)")+
  ggtitle("Environmental DNA - water (97% identity)")

## 5.1.3 Get fish taxa 
merged_phyloseq_97_fish = 
  subset_taxa(merged_phyloseq_97, class_name == "Actinopteri" & id != "Anguilla") # Anguilla is a freshwater eel that can live in brackish waters, but we remove here to be conservative 
merged_phyloseq_97_fish_glom = tax_glom(merged_phyloseq_97_fish, taxrank = "id",NArm=FALSE) # species level id's
# Get rid of 0 abundance taxonomic groups
merged_phyloseq_97_fish_glom = prune_taxa(taxa_sums(merged_phyloseq_97_fish_glom) > 0, merged_phyloseq_97_fish_glom)
merged_phyloseq_97_fish_export = as.data.frame(psmelt(merged_phyloseq_97_fish_glom)) 
merged_phyloseq_97_fish_export = aggregate(Abundance ~ id + region + sampleType, data = merged_phyloseq_97_fish_export, FUN = sum) # species level id's

## 5.1.2 Get shark taxa
merged_phyloseq_97_shark = subset_taxa(merged_phyloseq_97, class_name == "Chondrichthyes")
merged_phyloseq_97_shark_glom = tax_glom(merged_phyloseq_97_shark, taxrank = "id",NArm=FALSE) # Species level id's
# Get rid of 0 abundance taxonomic groups
merged_phyloseq_97_shark_glom = prune_taxa(taxa_sums(merged_phyloseq_97_shark_glom) > 0, merged_phyloseq_97_shark_glom)
merged_phyloseq_97_shark_export = as.data.frame(psmelt(merged_phyloseq_97_shark_glom)) 
merged_phyloseq_97_shark_export = aggregate(Abundance ~ id + region + sampleType, data = merged_phyloseq_97_shark_export, FUN = sum) # species level id's

## 5.1.4 Get benthic taxa
merged_phyloseq_97_corals = subset_taxa(merged_phyloseq_97, order_name == "Scleractinia")
merged_phyloseq_97_corals_glom = tax_glom(merged_phyloseq_97_corals, taxrank = "genus_name",NArm=FALSE)

merged_phyloseq_97_corals_glom = prune_taxa(taxa_sums(merged_phyloseq_97_corals_glom) > 0, merged_phyloseq_97_corals_glom)
merged_phyloseq_97_corals_export = as.data.frame(psmelt(merged_phyloseq_97_corals_glom)) 
merged_phyloseq_97_corals_export = aggregate(Abundance ~ family_name + region + sampleType + family_confidence, data = merged_phyloseq_97_corals_export, FUN = sum)

# Get phyla-level abundances
merged_phyloseq_phyla = tax_glom(merged_phyloseq_97, taxrank = "phylum_name", NArm=FALSE)
merged_phyloseq_phyla_glom = prune_taxa(taxa_sums(merged_phyloseq_phyla)>0, merged_phyloseq_phyla)
merged_phyloseq_phyla_export = as.data.frame(psmelt(merged_phyloseq_phyla_glom))
merged_phyloseq_phyla_export = aggregate(Abundance~phylum_name + sampleType, data = merged_phyloseq_phyla_export, FUN=sum)

#### 5.2 Underwater Visual Census data #### 
uvc_region = prune_taxa(taxa_sums(uvc_phylo) > 0, uvc_phylo)
uvc_region = merge_samples(uvc_region,"region")
uvc_region_relabund = transform_sample_counts(uvc_region, function(x) x / sum(x) )

uvc_region_plot = plot_bar(uvc_region_relabund, fill="family") + 
  geom_bar(stat="identity", position = "stack")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5, size =11),
        legend.background = element_rect(fill="transparent", color= NA),
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+ color_scale_fish +
  scale_x_discrete(labels=c("Lombok", "Misool", "Waigeo"))+
  guides(fill=guide_legend(ncol=2))+ ylab("Relative abundance (%)")+
  ggtitle("Family-level relative abundances of fishes (UVC)")

# Get fish taxa - uvc
# Get rid of 0 abundance taxonomic groups
uvc_region_fish_glom = prune_taxa(taxa_sums(uvc_phylo1) > 0, uvc_phylo1)
uvc_region_fish_export = as.data.frame(psmelt(uvc_region_fish_glom)) 
uvc_region_fish_export = aggregate(Abundance ~ OTU+region, data = uvc_region_fish_export, FUN = sum)

uvc_region_fish_export$sampleType = "UVC"
uvc_region_fish_export =uvc_region_fish_export[, c("OTU", "region", "sampleType", "Abundance")]
colnames(uvc_region_fish_export) = colnames(merged_phyloseq_97_fish_export)

# Compare detections between eDNA and UVC for a SI file
fish_detections = rbind(uvc_region_fish_export, merged_phyloseq_97_fish_export)
fish_detections$id = sub("_", " ", fish_detections$id)
fish_detections$region_sampleType = paste(fish_detections$region, fish_detections$sampleType)
fish_detections = subset(fish_detections, select = -c(region, sampleType))
fish_detections_wide = spread(fish_detections, region_sampleType, Abundance)

# Get giant clam data:
merged_phyloseq_97_clams = subset_taxa(merged_phyloseq_97, phylum_name == "Mollusca")
merged_phyloseq_97_clams_glom = tax_glom(merged_phyloseq_97_clams, taxrank = "id",NArm=FALSE)

merged_phyloseq_97_clams_glom = prune_taxa(taxa_sums(merged_phyloseq_97_clams_glom) > 0, merged_phyloseq_97_clams_glom)
merged_phyloseq_97_clams_export = as.data.frame(psmelt(merged_phyloseq_97_clams_glom)) 
merged_phyloseq_97_clams_export = aggregate(Abundance ~ id + region + sampleType, data = merged_phyloseq_97_clams_export, FUN = sum)

#### 5.3 Benthic data #### 
benthic_region = merge_samples(benthic_phylo,"region")
benthic_region_relabund = transform_sample_counts(benthic_region, function(x) x / sum(x) )

benthic_region_plot = plot_bar(benthic_region_relabund,fill="family") + 
  geom_bar(stat="identity", position = "stack")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),                       
                                                      panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"),
                                                      axis.title.x=element_blank(),
                                                      axis.text.x = element_text(angle=0, hjust=0.5, size =11),
                                                      legend.background = element_rect(fill="transparent", color= NA),
                                                      legend.box.background = element_rect(fill = "transparent", color=NA),
                                                      legend.key = element_rect(fill = "transparent",color=NA))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+ color_scale_benthic +
  scale_x_discrete(labels=c("Lombok", "Misool", "Waigeo"))+ ylab("Relative abundance (%)")+
  guides(fill=guide_legend(ncol=1))+
  ggtitle("Family-level relative abundance of stony corals (benthic)")

# Get benthic data
merged_benthic_corals_glom = tax_glom(benthic_phylo, taxrank = "family",NArm=FALSE)

merged_benthic_corals_glom = prune_taxa(taxa_sums(merged_benthic_corals_glom) > 0, merged_benthic_corals_glom)
merged_benthic_corals_export = as.data.frame(psmelt(merged_benthic_corals_glom)) 
merged_benthic_corals_export = aggregate(Abundance ~ family + region, data = merged_benthic_corals_export, FUN = sum)

#### 5.4 Arrange into one plot ####
combo_relabund_eDNA = ggarrange(water_97_plot,sed_97_plot,labels = c("(a)","(b)"), heights = c(1,1), common.legend = TRUE, legend="right")
ggsave("relabund_eDNA_v4.png",combo_relabund_eDNA, width=15, height=9, bg="white")

combo_relabund_UVC = ggarrange(uvc_region_plot, benthic_region_plot, widths = c(1.2,1), labels = c("(a)","(b)"), heights = c(1,1))
ggsave("relabund_uvc_v4.png",combo_relabund_UVC, width=15, height=9, bg="transparent")

####### 6. ORDINATIONS/PERMANOVAS  #######
#### 6.1 Environmental DNA data (full) #### 

##### 6.1.0 Set up distance matrix for OTU data, run PERMANOVA on region and sample type effects #####
otu_combo = merge_phyloseq(otu_water, otu_sediment)

otu_combo_perm = prune_samples(sample_sums(otu_combo)>0, otu_combo)
otu_combo_perm = subset_samples(otu_combo_perm, sampleType != "gut")

sample_data(otu_combo_perm)$name_sample = paste0(sample_data(otu_combo_perm)$region_code,"_", sample_data(otu_combo_perm)$sampleType)
sample_data(otu_combo_perm)$region_sample = paste0(sample_data(otu_combo_perm)$region,"_", sample_data(otu_combo_perm)$sampleType)
sample_data(otu_combo_perm)$sample_sums = sample_sums(otu_combo_perm)

sort((as.data.frame(sample_data((otu_combo_perm)))$region_code))
otu_perm_site_abund = as.matrix(sample_data(otu_combo_perm)[with(sample_data(otu_combo_perm), order(region,sample_sums))])

# Make sampling design balanced by taking out 3 Waigeo sites in both sediment and water, and 2 Misool sites only in water - now all regions should have 13 eDNA samples. Calculated by ensuring most sites have paired water and sediment samples, then removed sites based on their abundance
samples_to_remove = c("WG_YENW_sediment","WG_YENW_water", "WG_CKRI_sediment", "WG_CKRI_water")
otu_combo_perm = subset_samples(otu_combo_perm, !(name_sample %in% samples_to_remove))

# Get rid of ASVs that have zero sums across all samples
otu_combo_perm = prune_taxa(taxa_sums(otu_combo_perm) > 0, otu_combo_perm)

# Do a Hellinger transformation of the otu table 
otu_combo_perm = transform_sample_counts(otu_combo_perm, function(x) sqrt(x/sum(x)))

# add other environmental and social variables from sesync to eDNA metadata to use in permanova - here, use just the year
otu_combo_metadata = as(sample_data(otu_combo_perm),"data.frame")
site_sesync_no_t = subset(site_sesync, select = -c(site_code_t))
otu_combo_metadata$eDNA_name = row.names(otu_combo_metadata)
otu_combo_metadata1 = unique(merge(otu_combo_metadata, site_sesync_no_t))
rownames(otu_combo_metadata1) = otu_combo_metadata1$eDNA_name
sample_data(otu_combo_perm) = otu_combo_metadata1

# tried these variables as factors
otu_combo_metadata1$date = as.character(otu_combo_metadata1$date)

distance_method = distance(otu_combo_perm, "bray")

set.seed(617)
# This tests for regional effects, to see if there is spatial autocorrelation in community data
permanova_region  = adonis2(formula = distance_method ~ region + date+ sampleType, data = otu_combo_metadata1,perm = 999)

##### 6.1.0.1 PERMADISP to examine regional dispersions ##### 
permadisp_site = betadisper(distance_method,otu_combo_metadata$region)
permutest(permadisp_site, pairwise=TRUE)

edna_permadisp = gg_ordiplot(permadisp_site,groups = otu_combo_metadata$region, spiders = FALSE, ellipse = FALSE , hull = TRUE)
edna_permadisp_ord = edna_permadisp$df_ord
edna_permadisp_ord$sampleType = otu_combo_metadata$sampleType
edna_permadisp_ellipse = edna_permadisp$df_ellipse

edna_permadisp_ellipse_lbk = subset(edna_permadisp_ellipse, Group == "Lombok", droplevels = TRUE)
edna_permadisp_ellipse_ms = subset(edna_permadisp_ellipse, Group == "Misool", droplevels = TRUE)
edna_permadisp_ellipse_wg = subset(edna_permadisp_ellipse, Group == "Waigeo", droplevels = TRUE)

# used in MS - plotting just region on plot
edna_permadisp_plot = ggplot(edna_permadisp_ord, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = edna_permadisp_ellipse_lbk, color = "#CC3311")+geom_path(data=edna_permadisp_ellipse_ms, color ="#009988" )+geom_path(data=edna_permadisp_ellipse_wg, color = "#0077BB")+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" )) + scale_shape_manual(values = c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = edna_permadisp$plot$labels$x, y = edna_permadisp$plot$labels$y, fill="Region", shape = "Region",color="Region")+ coord_equal()

edna_boxplot_table = data.frame(Distance_to_centroid=permadisp_site$distances,Region=permadisp_site$group)
edna_boxplot_groups = permadisp_site$group
edna_boxplot = ggplot(data=edna_boxplot_table,aes(x=Region,y=Distance_to_centroid,colour=edna_boxplot_groups, shape=edna_boxplot_groups)) + geom_boxplot(alpha=0)+ scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" ))+
  scale_shape_manual(values=c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA), legend.position = "none")+
  labs(y = "Distance to centroid")+geom_jitter()

##### 6.1.0.2. PERMADISP to examine sample type dispersions (full dataset) #####
permadisp_method = betadisper(distance_method, otu_combo_metadata$sampleType)
permutest(permadisp_method, pairwise=TRUE)

edna_permadisp_method = gg_ordiplot(permadisp_method,groups = otu_combo_metadata$sampleType, spiders = FALSE, ellipse = FALSE , hull = TRUE)
edna_permadisp_method_ord = edna_permadisp_method$df_ord
edna_permadisp_method_ellipse = edna_permadisp_method$df_ellipse

edna_permadisp_water_ellipse = subset(edna_permadisp_method_ellipse, Group == "water", droplevels = TRUE)
edna_permadisp_sediment_ellipse = subset(edna_permadisp_method_ellipse, Group == "sediment", droplevels = TRUE)

edna_permadisp_method_plot = ggplot(edna_permadisp_method_ord, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = edna_permadisp_water_ellipse, color = "#439BA9")+geom_path(data=edna_permadisp_sediment_ellipse, color ="#DAB851" )+
  scale_color_manual(values = c("#DAB851","#439BA9")) + scale_shape_manual(values = c(0,3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = edna_permadisp_method$plot$labels$x, y = edna_permadisp_method$plot$labels$y, fill="Sample Type", shape = "Sample Type",color="Sample Type")+ coord_equal()

edna_method_boxplot_table = data.frame(Distance_to_centroid=permadisp_method$distances,SampleType=permadisp_method$group)
edna_method_boxplot_groups = permadisp_method$group
edna_method_boxplot = ggplot(data=edna_method_boxplot_table,aes(x=SampleType,y=Distance_to_centroid,colour=edna_method_boxplot_groups, shape = edna_method_boxplot_groups)) + 
  geom_boxplot(alpha = 0)+ scale_color_manual(values = c("#DAB851","#439BA9")) + scale_shape_manual(values = c(0,3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA), legend.position = "none")+
  labs(y = "Distance to centroid", x = "Sample Type")+geom_jitter()

##### 6.1.0.3. PERMADISP to examine sample type dispersions - Lombok #####
otu_combo_perm_lbk = subset_samples(otu_combo_perm, region == "Lombok")
distance_lombok = distance(otu_combo_perm_lbk, "bray")
permadisp_lombok = betadisper(distance_lombok, otu_combo_metadata[otu_combo_metadata$region=="Lombok",]$sampleType)
permutest(permadisp_lombok) # significant dispersion differences between Lombok water and sediment samples
ordiplot_lombok = gg_ordiplot(permadisp_lombok,groups = otu_combo_metadata[otu_combo_metadata$region=="Lombok",]$sampleType, spiders = FALSE, ellipse = FALSE , hull = FALSE)
lombok_ord_method = ordiplot_lombok$df_ord
ellipse_method_lombok_water = ordiplot_lombok$df_ellipse[ordiplot_lombok$df_ellipse$Group == "water",]
ellipse_method_lombok_sediment = ordiplot_lombok$df_ellipse[ordiplot_lombok$df_ellipse$Group == "sediment",]

lombok_permadisp_method_plot = ggplot(lombok_ord_method, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = ellipse_method_lombok_water, color = "#439BA9")+geom_path(data=ellipse_method_lombok_sediment, color ="#DAB851" )+
  scale_color_manual(values = c("#DAB851","#439BA9")) + scale_shape_manual(values = c(0,3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = ordiplot_lombok$plot$labels$x, y = ordiplot_lombok$plot$labels$y, fill="Sample Type", shape = "Sample Type",color="Sample Type")+ coord_equal()

lombok_method_boxplot_table = data.frame(Distance_to_centroid=permadisp_lombok$distances,SampleType=permadisp_lombok$group)
lombok_method_boxplot_groups = permadisp_lombok$group
lombok_method_boxplot = ggplot(data=lombok_method_boxplot_table,aes(x=SampleType,y=Distance_to_centroid,colour=lombok_method_boxplot_groups, shape = lombok_method_boxplot_groups)) + 
  geom_boxplot(alpha = 0)+ scale_color_manual(values = c("#DAB851","#439BA9")) + scale_shape_manual(values = c(0,3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA), legend.position = "none")+
  labs(y = "Distance to centroid", x = "Sample Type")+geom_jitter()

##### 6.1.0.4. PERMADISP to examine sample type dispersions - Misool #####
otu_combo_perm_ms = subset_samples(otu_combo_perm, region == "Misool")
distance_misool = distance(otu_combo_perm_ms, "bray")
permadisp_misool = betadisper(distance_misool, otu_combo_metadata[otu_combo_metadata$region=="Misool",]$sampleType)
permutest(permadisp_misool) # significant dispersion differences between Misool water and sediment samples
ordiplot_misool = gg_ordiplot(permadisp_misool, groups = otu_combo_metadata[otu_combo_metadata$region=="Misool",]$sampleType, spiders = FALSE, ellipse = TRUE , hull = FALSE)
misool_ord_method = ordiplot_misool$df_ord
ellipse_method_misool_water = ordiplot_misool$df_ellipse[ordiplot_misool$df_ellipse$Group=="water",]
ellipse_method_misool_sediment = ordiplot_misool$df_ellipse[ordiplot_misool$df_ellipse$Group=="sediment",]

misool_permadisp_method_plot =ggplot(misool_ord_method, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = ellipse_method_misool_water, color = "#439BA9")+geom_path(data=ellipse_method_misool_sediment, color ="#DAB851" )+
  scale_color_manual(values = c("#DAB851","#439BA9")) + scale_shape_manual(values = c(0,3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = ordiplot_misool$plot$labels$x, y = ordiplot_misool$plot$labels$y, fill="Sample Type", shape = "Sample Type",color="Sample Type")+ coord_equal()

misool_method_boxplot_table = data.frame(Distance_to_centroid=permadisp_misool$distances,SampleType=permadisp_misool$group)
misool_method_boxplot_groups = permadisp_misool$group
misool_method_boxplot = ggplot(data=misool_method_boxplot_table,aes(x=SampleType,y=Distance_to_centroid,colour=misool_method_boxplot_groups, shape = misool_method_boxplot_groups)) + 
  geom_boxplot(alpha = 0)+ scale_color_manual(values = c("#DAB851","#439BA9")) + scale_shape_manual(values = c(0,3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA), legend.position = "none")+
  labs(y = "Distance to centroid", x = "Sample Type")+geom_jitter()

##### 6.1.0.5. PERMADISP to examine sample type dispersions - Waigeo #####
otu_combo_perm_wg = subset_samples(otu_combo_perm, region == "Waigeo")
distance_waigeo = distance(otu_combo_perm_wg, "bray")
permadisp_waigeo = betadisper(distance_waigeo, otu_combo_metadata[otu_combo_metadata$region=="Waigeo",]$sampleType)
permutest(permadisp_waigeo) #  no significant dispersion
ordiplot_waigeo = gg_ordiplot(permadisp_waigeo, groups = otu_combo_metadata[otu_combo_metadata$region=="Waigeo",]$sampleType, spiders = FALSE, ellipse = TRUE , hull = FALSE)
waigeo_ord_method = ordiplot_waigeo$df_ord
ellipse_method_waigeo_water = ordiplot_waigeo$df_ellipse[ordiplot_waigeo$df_ellipse$Group=="water",]
ellipse_method_waigeo_sediment = ordiplot_waigeo$df_ellipse[ordiplot_waigeo$df_ellipse$Group=="sediment",]

waigeo_permadisp_method_plot =ggplot(waigeo_ord_method, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = ellipse_method_waigeo_water, color = "#439BA9")+geom_path(data=ellipse_method_waigeo_sediment, color ="#DAB851" )+
  scale_color_manual(values = c("#DAB851","#439BA9")) + scale_shape_manual(values = c(0,3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = ordiplot_waigeo$plot$labels$x, y = ordiplot_waigeo$plot$labels$y, fill="Sample Type", shape = "Sample Type",color="Sample Type")+ coord_equal()

waigeo_method_boxplot_table = data.frame(Distance_to_centroid=permadisp_waigeo$distances,SampleType=permadisp_waigeo$group)
waigeo_method_boxplot_groups = permadisp_waigeo$group
waigeo_method_boxplot = ggplot(data=waigeo_method_boxplot_table,aes(x=SampleType,y=Distance_to_centroid,colour=waigeo_method_boxplot_groups, shape = waigeo_method_boxplot_groups)) + 
  geom_boxplot(alpha = 0)+ scale_color_manual(values = c("#DAB851","#439BA9")) + scale_shape_manual(values = c(0,3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA), legend.position = "none")+
  labs(y = "Distance to centroid", x = "Sample Type")+geom_jitter()



##### 6.1.1. PERMADISP of regional effects on water samples #####
otu_water_filt = subset_samples(otu_combo_perm, sampleType == "water")

set.seed(617)

otu_water_filt_dist = distance(otu_water_filt, "bray")
otu_water_filt_meta = as.data.frame(as.matrix(sample_data(otu_water_filt)))
otu_water_filt_perm = adonis2(formula = otu_water_filt_dist ~ region + date, data = otu_water_filt_meta, perm = 999)

permadisp_site_otu_water = betadisper(otu_water_filt_dist,otu_water_filt_meta$region)

permutest(permadisp_site_otu_water, pairwise=TRUE)

plot(permadisp_site_otu_water, main ="Ordination centroids of Bray-Curtis distances and dispersion effects of region on water samples only", col = c("#C87A8A", "#6B9D59", "#5F96C2")) 
boxplot(permadisp_site_otu_water, main = "Dispersion effects by region in eDNA water samples only", xlab = "")

edna_permadisp_water = gg_ordiplot(permadisp_site_otu_water, groups = otu_water_filt_meta$region, spiders = FALSE, ellipse=FALSE, hull=TRUE)
edna_permadisp_water_ord = edna_permadisp_water$df_ord
edna_permadisp_water_ellipse = edna_permadisp_water$df_ellipse

edna_pd_water_ellipse_lbk = subset(edna_permadisp_water_ellipse, Group == "Lombok", droplevels = TRUE)
edna_pd_water_ellipse_ms = subset(edna_permadisp_water_ellipse, Group == "Misool", droplevels = TRUE)
edna_pd_water_ellipse_wg = subset(edna_permadisp_water_ellipse, Group == "Waigeo", droplevels = TRUE)

edna_permadisp_water_plot = ggplot(edna_permadisp_water_ord, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = edna_pd_water_ellipse_lbk, color = "#CC3311")+geom_path(data=edna_pd_water_ellipse_ms, color ="#009988" )+geom_path(data=edna_pd_water_ellipse_wg, color = "#0077BB")+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" )) + scale_shape_manual(values = c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = edna_permadisp_water$plot$labels$x, y = edna_permadisp_water$plot$labels$y, fill="Region", shape = "Region",color="Region")+ coord_equal()

edna_boxplot_water_table = data.frame(Distance_to_centroid=permadisp_site_otu_water$distances,Region=permadisp_site_otu_water$group)
edna_boxplot_water_groups = permadisp_site_otu_water$group
edna_boxplot_water = ggplot(data=edna_boxplot_water_table,aes(x=Region,y=Distance_to_centroid,colour=edna_boxplot_water_groups, shape=edna_boxplot_water_groups)) + geom_boxplot(alpha=0)+ scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" ))+
  scale_shape_manual(values=c(5,6,7))+ 
  color_scale_region +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA), legend.position="none")+
  labs(y = "Distance to centroid")+geom_jitter()

##### 6.1.2 PERMADISP of regional effect on sediment samples ##### 
otu_sediment_filt = subset_samples(otu_combo_perm, sampleType == "sediment")


set.seed(617)

otu_sediment_filt_dist = distance(otu_sediment_filt, "bray")
otu_sediment_filt_meta = as.data.frame(as.matrix(sample_data(otu_sediment_filt)))
otu_sediment_filt_perm = adonis2(formula = otu_sediment_filt_dist ~ region + date, data = otu_sediment_filt_meta, perm = 999)

permadisp_site_otu_sediment = betadisper(otu_sediment_filt_dist,otu_sediment_filt_meta$region)

permutest(permadisp_site_otu_sediment, pairwise=TRUE)

plot(permadisp_site_otu_sediment, main ="Ordination centroids of Bray-Curtis distances and dispersion effects of region on sediment samples only", col = c("#C87A8A", "#6B9D59", "#5F96C2")) 
boxplot(permadisp_site_otu_sediment, main = "Dispersion effects by region in eDNA sediment samples only", xlab = "")

edna_permadisp_sediment = gg_ordiplot(permadisp_site_otu_sediment, groups = otu_sediment_filt_meta$region, spiders = FALSE, ellipse=FALSE, hull=TRUE)
edna_permadisp_sediment_ord = edna_permadisp_sediment$df_ord
edna_permadisp_sediment_ellipse = edna_permadisp_sediment$df_ellipse

edna_pd_sediment_ellipse_lbk = subset(edna_permadisp_sediment_ellipse, Group == "Lombok", droplevels=TRUE)
edna_pd_sediment_ellipse_ms = subset(edna_permadisp_sediment_ellipse, Group == "Misool", droplevels=TRUE)
edna_pd_sediment_ellipse_wg = subset(edna_permadisp_sediment_ellipse, Group == "Waigeo", droplevels=TRUE)

edna_permadisp_sediment_plot = ggplot(edna_permadisp_sediment_ord, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = edna_pd_sediment_ellipse_lbk, color = "#CC3311")+geom_path(data=edna_pd_sediment_ellipse_ms, color ="#009988" )+geom_path(data=edna_pd_sediment_ellipse_wg, color = "#0077BB")+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" )) + scale_shape_manual(values = c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = edna_permadisp_sediment$plot$labels$x, y = edna_permadisp_sediment$plot$labels$y, fill="Region", shape = "Region",color="Region")+ coord_equal()

edna_boxplot_sediment_table = data.frame(Distance_to_centroid=permadisp_site_otu_sediment$distances,Region=permadisp_site_otu_sediment$group)
edna_boxplot_sediment_groups = permadisp_site_otu_sediment$group
edna_boxplot_sediment = ggplot(data=edna_boxplot_sediment_table,aes(x=Region,y=Distance_to_centroid,colour=edna_boxplot_sediment_groups, shape=edna_boxplot_sediment_groups)) + geom_boxplot(alpha=0)+ scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" ))+
  scale_shape_manual(values=c(5,6,7))+ 
  color_scale_region +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA), legend.position="none")+
  labs(y = "Distance to centroid")+geom_jitter()

##### 6.1.3.  PERMANOVA and PERMADISP region and method effects on non-metazoans ##### 
nonanimal_combo_perm = prune_samples(sample_sums(merged_phyloseq_97_nonanimal)>0, merged_phyloseq_97_nonanimal)

sample_data(nonanimal_combo_perm)$name_sample = paste0(sample_data(nonanimal_combo_perm)$region_code,"_", sample_data(nonanimal_combo_perm)$sampleType)
sample_data(nonanimal_combo_perm)$region_sample = paste0(sample_data(nonanimal_combo_perm)$region,"_", sample_data(nonanimal_combo_perm)$sampleType)
sample_data(nonanimal_combo_perm)$sample_sums = sample_sums(nonanimal_combo_perm)

sort((as.data.frame(sample_data((nonanimal_combo_perm)))$region_code))
nonanimal_perm_site_abund = as.matrix(sample_data(nonanimal_combo_perm)[with(sample_data(nonanimal_combo_perm), order(region,sample_sums))])

# Make sampling design balanced - same as OTUs above, except more was removed so it's 13 sites per region. Additional sites were removed based on abundance (most-low abundance sites removed)
samples_to_remove = c("WG_YENW_sediment","WG_YENW_water", "WG_CKRI_sediment", "WG_CKRI_water")
nonanimal_combo_perm = subset_samples(nonanimal_combo_perm, !(name_sample %in% samples_to_remove))

# Get rid of ASVs that have zero sums across all samples
nonanimal_combo_perm = prune_taxa(taxa_sums(nonanimal_combo_perm) > 0, nonanimal_combo_perm)

# Do a Hellinger transformation of the nonanimal table 
nonanimal_combo_perm = transform_sample_counts(nonanimal_combo_perm, function(x) sqrt(x/sum(x)))

# add other environmental and social variables from sesync to eDNA metadata to use in permanova - here, use just the year
nonanimal_combo_metadata = as(sample_data(nonanimal_combo_perm),"data.frame")
site_sesync_no_t = subset(site_sesync, select = -c(site_code_t))
nonanimal_combo_metadata$eDNA_name = row.names(nonanimal_combo_metadata)
nonanimal_combo_metadata1 = unique(merge(nonanimal_combo_metadata, site_sesync_no_t))
rownames(nonanimal_combo_metadata1) = nonanimal_combo_metadata1$eDNA_name
sample_data(nonanimal_combo_perm) = nonanimal_combo_metadata1

# tried these variables as factors
nonanimal_combo_metadata1$date = as.character(nonanimal_combo_metadata1$date)

distance_nonanimal_method = distance(nonanimal_combo_perm, "bray")

set.seed(617)
# This tests for regional effects, to see if there is spatial autocorrelation in community data
permanova_nonanimal_region  = adonis2(formula = distance_nonanimal_method ~ region + site_name + date, data = nonanimal_combo_metadata1,perm = 999)

# This blocks for region, accounts for spatial autocorrelation
permanova_nonanimal_block = how(within = Within(type="free", mirror = TRUE),
                                plots = Plots(strata = nonanimal_combo_metadata1$region, type = "free", mirror = TRUE))

permanova_nonanimal_method = adonis2(formula = distance_nonanimal_method ~ sampleType + pop2020_25km + npp_mean + management_type , permutations = permanova_nonanimal_block, data = nonanimal_combo_metadata1, by="margin", perm=999)

##### 6.1.3.1 PERMADISP to examine sample type dispersion on non-metazoans ##### 
permadisp_nonanimal_site = betadisper(distance_nonanimal_method,nonanimal_combo_metadata$region)

permutest(permadisp_nonanimal_site, pairwise=TRUE)
#plot(permadisp_nonanimal_site, main ="Ordination centroids of Bray-Curtis distances and dispersion effects of region for eDNA nonanimals", col = c("#C87A8A", "#6B9D59", "#5F96C2")) 
#boxplot(permadisp_nonanimal_site, main = "Dispersion effects by region in eDNA nonanimal data", xlab = "")
edna_permadisp_nonanimal = gg_ordiplot(permadisp_nonanimal_site,groups = nonanimal_combo_metadata1$region, spiders = FALSE, ellipse = FALSE , hull = TRUE)
edna_permadisp_nonanimal_ord = edna_permadisp_nonanimal$df_ord
edna_permadisp_nonanimal_ellipse = edna_permadisp_nonanimal$df_ellipse

edna_pd_nonanimal_ellipse_lbk = subset(edna_permadisp_nonanimal_ellipse, Group =="Lombok", droplevels=TRUE)
edna_pd_nonanimal_ellipse_ms = subset(edna_permadisp_nonanimal_ellipse, Group =="Misool", droplevels=TRUE)
edna_pd_nonanimal_ellipse_wg = subset(edna_permadisp_nonanimal_ellipse, Group =="Waigeo", droplevels=TRUE)

edna_permadisp_nonanimal_plot = ggplot(edna_permadisp_nonanimal_ord, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = edna_pd_nonanimal_ellipse_lbk, color = "#CC3311")+geom_path(data=edna_pd_nonanimal_ellipse_ms, color ="#009988" )+geom_path(data=edna_pd_nonanimal_ellipse_wg, color = "#0077BB")+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" )) + scale_shape_manual(values = c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = edna_permadisp_nonanimal$plot$labels$x, y = edna_permadisp_nonanimal$plot$labels$y, fill="Region", shape = "Region",color="Region")+ coord_equal()

edna_boxplot_nonanimal_table = data.frame(Distance_to_centroid=permadisp_nonanimal_site$distances,Region=permadisp_nonanimal_site$group)
edna_boxplot_nonanimal_groups = permadisp_nonanimal_site$group
edna_boxplot_nonanimal = ggplot(data=edna_boxplot_nonanimal_table,aes(x=Region,y=Distance_to_centroid,colour=edna_boxplot_nonanimal_groups, shape = edna_boxplot_nonanimal_groups)) + geom_boxplot(alpha = 0)+ color_scale_region + scale_shape_manual(values=c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA), legend.position="none")+labs(y = "Distance to centroid")+geom_jitter()

permadisp_nonanimal_method = betadisper(distance_nonanimal_method, nonanimal_combo_metadata$sampleType)

permutest(permadisp_nonanimal_method, pairwise=TRUE)
plot(permadisp_nonanimal_method, main ="Ordination centroids of Bray-Curtis distances and dispersion effects of eDNA sample type on nonanimals") 
boxplot(permadisp_nonanimal_method, main = "Dispersion effects by sample type in eDNA nonanimal data", xlab = "")

##### 6.1.4. Combine plots #####
permanova_arranged_nonanimal = ggarrange(edna_permadisp_nonanimal_plot + coord_fixed(ratio=1)+ggtitle("Environmental DNA - non-animals"), edna_boxplot_nonanimal, widths = c(1,0.5), align="hv", common.legend = TRUE, legend = "right")
ggsave("permanova_nonanimal.png", permanova_arranged_nonanimal,bg="white")

#### 6.2 PERMANOVA and PERMADISP of regional effects on Underwater Visual Census #### 
samples_to_keep = sort(unique((as.data.frame(sample_data((otu_combo_perm)))$region_code)))
uvc_phylo1_perm = subset_samples(uvc_phylo1, region_code %in% samples_to_keep)

uvc_phylo1_perm = prune_taxa(taxa_sums(uvc_phylo1_perm) > 0, uvc_phylo1_perm) # get rid of zero-abundance counts
meta_region_uvc = as(sample_data(uvc_phylo1_perm), "data.frame")

# Transform the data using a hellinger transformation 
uvc_phylo1_perm = transform_sample_counts(uvc_phylo1_perm, function(x) sqrt(x/sum(x)))

# Look at regional effect
set.seed(617)
distance_region_uvc = phyloseq::distance(uvc_phylo1_perm, method = "bray")
permanova_region_uvc = adonis2(formula = distance_region_uvc ~ region + date + region*date, data = meta_region_uvc, perm=999)

# Look at dispersion
permadisp_site_uvc = betadisper(distance_region_uvc,meta_region_uvc$region, add="cailliez")
permadisp_site_uvc$vectors[,2] = permadisp_site_uvc$vectors[,2]*-1
permutest(permadisp_site_uvc, pairwise=TRUE)

# plot(permadisp_site_uvc,main ="PCoA of hellinger-transformed UVC data using Bray-Curtis distances", col = c("#C87A8A", "#6B9D59", "#5F96C2"))
uvc_permadisp = gg_ordiplot(permadisp_site_uvc, groups = meta_region_uvc$region, spiders = FALSE, ellipse = FALSE, hull = TRUE)
uvc_permadisp_ord = uvc_permadisp$df_ord
uvc_permadisp_ellipse = uvc_permadisp$df_ellipse

uvc_pd_ellipse_lbk = subset(uvc_permadisp_ellipse, Group == "Lombok", droplevels=TRUE)
uvc_pd_ellipse_ms = subset(uvc_permadisp_ellipse, Group == "Misool", droplevels=TRUE)
uvc_pd_ellipse_wg = subset(uvc_permadisp_ellipse, Group == "Waigeo", droplevels=TRUE)

uvc_permadisp_plot = ggplot(uvc_permadisp_ord, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = uvc_pd_ellipse_lbk, color = "#CC3311")+geom_path(data=uvc_pd_ellipse_ms, color ="#009988" )+geom_path(data=uvc_pd_ellipse_wg, color = "#0077BB")+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" )) + scale_shape_manual(values = c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = uvc_permadisp$plot$labels$x, y = uvc_permadisp$plot$labels$y, fill="Region", shape = "Region",color="Region")+ coord_equal()

uvc_boxplot_table = data.frame(Distance_to_centroid = permadisp_site_uvc$distances, Region = permadisp_site_uvc$group)
uvc_boxplot_groups = permadisp_site_uvc$group
uvc_boxplot = ggplot(data = uvc_boxplot_table, aes(x=Region, y=Distance_to_centroid, colour = uvc_boxplot_groups, shape = uvc_boxplot_groups))+geom_boxplot(alpha=0)+ color_scale_region + scale_shape_manual(values=c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA),
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA),legend.position="none")+labs(y = "Distance to centroid")+geom_jitter()

#### 6.3 PERMANOVA and PERMADISP of regional effect on benthic data #### 
# If code doesn't work, re-load and merge benthic_phylo1 from step 1 (loading data)
# Balance sampling design
samples_to_keep = sort(unique((as.data.frame(sample_data((uvc_phylo1_perm)))$region_code)))

benthic_phylo1_perm = subset_samples(benthic_phylo1, sample_names(benthic_phylo1) %in% samples_to_keep)

# Hellinger transformation
benthic_phylo1_perm = transform_sample_counts(benthic_phylo1_perm , function(x) sqrt(x/sum(x)))

# Permanova - significant regional differences
meta_region_benthic = as(sample_data(benthic_phylo1_perm),"data.frame")
distance_region_benthic = phyloseq::distance(benthic_phylo1_perm, "bray")

set.seed(617)
permanova_region_benthic = adonis2(formula = distance_region_benthic ~ region + date + region*date , data = meta_region_benthic, perm=999)

# Look at dispersion
permadisp_site_benthic = betadisper(distance_region_benthic,meta_region_benthic$region, add = TRUE)
permadisp_site_benthic$vectors[,c(1,2)] = permadisp_site_benthic$vectors[,c(1,2)]*-1
permutest(permadisp_site_benthic, pairwise=TRUE)

benthic_permadisp = gg_ordiplot(permadisp_site_benthic, groups = meta_region_benthic$region, spiders = FALSE, ellipse = FALSE, hull = TRUE )
benthic_permadisp_ord = benthic_permadisp$df_ord
benthic_permadisp_ellipse = benthic_permadisp$df_ellipse

benthic_pd_ellipse_lbk = subset(benthic_permadisp_ellipse, Group =="Lombok", droplevels=TRUE)
benthic_pd_ellipse_ms = subset(benthic_permadisp_ellipse, Group =="Misool", droplevels=TRUE)
benthic_pd_ellipse_wg = subset(benthic_permadisp_ellipse, Group =="Waigeo", droplevels=TRUE)

benthic_permadisp_plot = ggplot(benthic_permadisp_ord, aes(x=x,y=y), color = Group) + geom_point(aes(color = Group, shape = Group))+
  geom_path(data = benthic_pd_ellipse_lbk, color = "#CC3311")+geom_path(data=benthic_pd_ellipse_ms, color ="#009988" )+geom_path(data=benthic_pd_ellipse_wg, color = "#0077BB")+
  scale_color_manual(values = c("#CC3311", "#009988", "#0077BB" )) + scale_shape_manual(values = c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))+scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = benthic_permadisp$plot$labels$x, y = benthic_permadisp$plot$labels$y, fill="Region", shape = "Region",color="Region")+ coord_equal()

benthic_boxplot_table = data.frame(Distance_to_centroid = permadisp_site_benthic$distances, Region = permadisp_site_benthic$group)
benthic_boxplot_groups = permadisp_site_benthic$group
benthic_boxplot = ggplot(data = benthic_boxplot_table, aes(x=Region, y=Distance_to_centroid, colour = benthic_boxplot_groups, shape = benthic_boxplot_groups))+geom_boxplot(alpha=0)+ color_scale_region + scale_shape_manual(values=c(5,6,7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA),
        axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA), legend.position="none")+
  labs(y = "Distance to centroid")+geom_jitter()

#### 6.4 Combine PERMANOVA and PERMADISP plots ####
permanova_arranged_edna = ggarrange(edna_permadisp_plot+ggtitle("(a) Environmental DNA (full)"), edna_boxplot, 
                                    edna_permadisp_water_plot+ggtitle("(b) Environmental DNA (water samples)"), edna_boxplot_water, 
                                    edna_permadisp_sediment_plot+ggtitle("(c) Environmental DNA (sediment samples)"), edna_boxplot_sediment,
                                    nrow=3, ncol=2, widths = c(1,0.7), align = "v", common.legend = TRUE, legend = "bottom")
ggsave("permanova_edna_v4.png", permanova_arranged_edna, bg="white", scale = 1.3, width=6)

permanova_arranged_uvc = ggarrange(uvc_permadisp_plot+ggtitle("(a) Visual surveys of fishes"), uvc_boxplot, benthic_permadisp_plot+ggtitle("(b) Visual surveys of benthos"), benthic_boxplot, nrow = 2, ncol = 2, heights = c(1,1), widths = c(1,0.7), align = "v", common.legend=TRUE, legend = "bottom")
ggsave("permanova_uvc.png", permanova_arranged_uvc, bg = "white", scale = 1.3, width = 6)

permanova_arranged_edna_method = ggarrange(edna_permadisp_method_plot + ggtitle("(a) Environmental DNA (full) - Sample Type"), edna_method_boxplot, 
                                           lombok_permadisp_method_plot + ggtitle("(b) Environmental DNA (Lombok)"), lombok_method_boxplot,
                                           misool_permadisp_method_plot + ggtitle("(c) Environmental DNA (Misool)"), misool_method_boxplot,
                                           waigeo_permadisp_method_plot + ggtitle("(d) Environmental DNA (Waigeo)"), waigeo_method_boxplot,
                                           nrow=4, ncol=2, widths = c(1,0.6), align="h", common.legend=TRUE, legend = "bottom" )
ggsave("permanova_edna_method_v4.png", permanova_arranged_edna_method,dpi=800, bg = "white", scale = 1.3, width = 5)

####### 7. ALPHA DIVERSITY  #######
#### 7.1 Environmental DNA #####
## 6.1.1. Water samples ##
otu_div = microbiome::alpha(otu_water,index = c("observed","diversity_shannon", "evenness_pielou"))
meta_otu_water = meta(otu_water)
meta_otu_water$sam_name = rownames(meta_otu_water)
otu_div$sam_name <- rownames(meta_otu_water)
otu_div = merge(otu_div, meta_otu_water, by ="sam_name")
otu_div = otu_div[,c(1:7)]

otu_div_melt <- reshape2::melt(otu_div)

# Plot
p <- ggboxplot(otu_div_melt, x = "region", y = "value",
               order = c("Lombok","Misool", "Waigeo"), fill = "region",
               legend= "right",
               facet.by = "variable", 
               scales = "free", palette = color_table_region) + rremove("x.text")

lev <- c("Lombok","Misool", "Waigeo")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif") + color_scale_region +  ggtitle("Alpha diversity - eDNA water samples")
p2

print(p2)

## 6.1.2 Sediment samples ## 
otu_div_sed = microbiome::alpha(otu_sediment,index = c("observed","diversity_shannon", "evenness_pielou"))
meta_otu_sed = meta(otu_sediment)
meta_otu_sed$sam_name = rownames(meta_otu_sed)
otu_div_sed$sam_name <- rownames(meta_otu_sed)
otu_div_sed = merge(otu_div_sed, meta_otu_sed, by ="sam_name")
otu_div_sed = otu_div_sed[,c(1:7)]

otu_div_melt <- reshape2::melt(otu_div_sed)

# Plot
p <- ggboxplot(otu_div_melt, x = "region", y = "value",
               fill = "region",
               order = c("Lombok","Misool", "Waigeo"),
               legend= "right",
               facet.by = "variable", 
               scales = "free",palette = color_table_region) + rremove("x.text")

lev <- c( "Lombok","Waigeo", "Misool")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif") + ggtitle("Alpha diversity - eDNA sediment samples")

print(p2)

#### 7.2 Underwater Visual Census #####
uvc_phylo_div = microbiome::alpha(uvc_phylo1,index = c("observed","diversity_shannon", "evenness_pielou"))
meta_uvc_phylo = meta(uvc_phylo1)
meta_uvc_phylo$region_code = rownames(meta_uvc_phylo)
uvc_phylo_div$region_code <- rownames(meta_uvc_phylo)
uvc_phylo_div = merge(uvc_phylo_div, meta_uvc_phylo, by ="region_code")
uvc_phylo_div = uvc_phylo_div[,c(1:4,8:9)]

uvc_phylo_div_melt <- reshape2::melt(uvc_phylo_div)

# Plot
p <- ggboxplot(uvc_phylo_div_melt, x = "region", y = "value",
               fill = "region",
               order = c("Lombok","Misool", "Waigeo"),
               legend= "right",
               facet.by = "variable", 
               scales = "free", palette = color_table_region) + rremove("x.text")

lev <- c("Lombok","Misool", "Waigeo")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif") + ggtitle("Alpha diversity - Underwater Visual Census")

print(p2)

#### 7.3 Benthic Surveys ####
benthic_div = microbiome::alpha(benthic_phylo1,index = c("observed","diversity_shannon", "evenness_pielou"))
meta_benthic_div = meta(benthic_phylo1)
meta_benthic_div$sam_name = rownames(meta_benthic_div)
benthic_div$sam_name <- rownames(meta_benthic_div)
benthic_div = merge(benthic_div, meta_benthic_div, by ="sam_name")
benthic_div = benthic_div[,c(1:6)]

benthic_div_melt <- reshape2::melt(benthic_div)
benthic_div_melt$region = str_to_title(benthic_div_melt$region)

# Plot
p <- ggboxplot(benthic_div_melt, x = "region", y = "value",
               order = c("Lombok","Misool", "Waigeo"), fill = "region",
               legend= "right",
               facet.by = "variable", 
               scales = "free", palette = color_table_region) + rremove("x.text")

lev <- c("Waigeo", "Lombok","Misool")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif")+ color_scale_region + ggtitle("Alpha diversity - Benthic Surveys")
p2

print(p2)

#### 7.4 Combine based on metric ####
otu_div$method = "Environmental DNA - water"
otu_div = subset(otu_div, select = -c(sampleID))
otu_div_sed$method = "Environmental DNA - sediment"
otu_div_sed = subset(otu_div_sed, select = -c(sampleID))
uvc_phylo_div$method = "Visual - fishes"
uvc_phylo_div$sam_name = uvc_phylo_div$site_code_1 
colnames(uvc_phylo_div)[1] = "sam_name"
benthic_div$method = "Visual - benthos"

alpha_div = rbind(otu_div, otu_div_sed)
alpha_div_vis = rbind(uvc_phylo_div,benthic_div)
alpha_div$region = str_to_title(alpha_div$region)
alpha_div_vis$region = str_to_title(alpha_div_vis$region)

alpha_div_melt = melt(alpha_div)
alpha_div_melt_vis = melt(alpha_div_vis)
alpha_div_melt$facet = factor(alpha_div_melt$method, levels = c("Environmental DNA - water", "Environmental DNA - sediment"))
alpha_div_melt_vis$facet = factor(alpha_div_melt_vis$method, levels = c("Visual - fishes", "Visual - benthos"))
levels(alpha_div_melt$variable) = c("Observed", "Diversity (Shannon)", "Evenness (Pielou)")
levels(alpha_div_melt_vis$variable) = c("Observed", "Diversity (Shannon)", "Evenness (Pielou)")

p = ggboxplot(alpha_div_melt, x = "region", y = "value", 
              order = c("Lombok","Misool", "Waigeo"), color = "region",fill="transparent",palette = c("#CC3311", "#009988", "#0077BB"),
              facet.by = c("variable", "facet"),
              scales = "free") + rremove("axis") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  geom_point(aes(col=region, shape=region), position="jitter")+scale_shape_manual(values=c(5,6,7))

lev <- c("Lombok","Misool", "Waigeo")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif")+scale_y_continuous(expand = expansion(mult=c(0.05,0.15)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))

p2 

ggsave("eDNA_alphadiv_v4.png",p2, bg="white")

# For visual surveys
p_vis = ggboxplot(alpha_div_melt_vis, x = "region", y = "value", 
                  order = c("Lombok","Misool", "Waigeo"), color = "region",fill="transparent", palette = c("#CC3311", "#009988", "#0077BB"),
                  facet.by = c("variable", "facet"),
                  scales = "free") + rremove("axis")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  geom_point(aes(col=region, shape=region), position="jitter")+scale_shape_manual(values=c(5,6,7))


lev <- c("Lombok","Misool", "Waigeo")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

p2_vis <- p_vis + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif") + scale_y_continuous(expand = expansion(mult=c(0.05,0.15)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent", color = NA), plot.background = element_rect(fill="transparent", color = NA), axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill="transparent", color= NA),
        legend.box.background = element_rect(fill = "transparent", color=NA),
        legend.key = element_rect(fill = "transparent",color=NA))

p2_vis

ggsave("uvc_alphadiv_v4.png", p2_vis, bg="white")

#### 7.5 Metazoans vs non-Metazoans ####
# Use ASV's for this 
animal_water_div = microbiome::alpha(merged_phyloseq_97_animal_w,index = c("observed","diversity_shannon", "evenness_pielou"))
meta_animal_water = meta(merged_phyloseq_97_animal_w)
meta_animal_water$sam_name = rownames(meta_animal_water)
animal_water_div$sam_name <- rownames(meta_animal_water)
animal_water_div = merge(animal_water_div, meta_animal_water, by ="sam_name")
animal_water_div = animal_water_div[,c(1:7)]
animal_water_div$method = "eDNA water - Metazoans"

animal_sed_div = microbiome::alpha(merged_phyloseq_97_animal_s,index = c("observed","diversity_shannon", "evenness_pielou"))
meta_animal_sed = meta(merged_phyloseq_97_animal_s)
meta_animal_sed$sam_name = rownames(meta_animal_sed)
animal_sed_div$sam_name <- rownames(meta_animal_sed)
animal_sed_div = merge(animal_sed_div, meta_animal_sed, by ="sam_name")
animal_sed_div = animal_sed_div[,c(1:7)]
animal_sed_div$method = "eDNA sediment - Metazoans"

nonanimal_water_div = microbiome::alpha(merged_phyloseq_97_nonanimal_w,index = c("observed","diversity_shannon", "evenness_pielou"))
meta_nonanimal_water = meta(merged_phyloseq_97_nonanimal_w)
meta_nonanimal_water$sam_name = rownames(meta_nonanimal_water)
nonanimal_water_div$sam_name <- rownames(meta_nonanimal_water)
nonanimal_water_div = merge(nonanimal_water_div, meta_nonanimal_water, by ="sam_name")
nonanimal_water_div = nonanimal_water_div[,c(1:7)]
nonanimal_water_div$method = "eDNA water - Non-metazoans"

nonanimal_sed_div = microbiome::alpha(merged_phyloseq_97_nonanimal_s,index = c("observed","diversity_shannon", "evenness_pielou"))
meta_nonanimal_sed = meta(merged_phyloseq_97_nonanimal_s)
meta_nonanimal_sed$sam_name = rownames(meta_nonanimal_sed)
nonanimal_sed_div$sam_name <- rownames(meta_nonanimal_sed)
nonanimal_sed_div = merge(nonanimal_sed_div, meta_nonanimal_sed, by ="sam_name")
nonanimal_sed_div = nonanimal_sed_div[,c(1:7)]
nonanimal_sed_div$method = "eDNA sediment - Non-metazoans"

alpha_div2 = rbind(animal_water_div, animal_sed_div,nonanimal_water_div,nonanimal_sed_div)
alpha_div2$region = str_to_title(alpha_div2$region)

alpha_div2_melt = melt(alpha_div2)
alpha_div2_melt$group = factor(alpha_div2_melt$method, levels = c("eDNA water - Metazoans", "eDNA water - Non-metazoans", "eDNA sediment - Metazoans","eDNA sediment - Non-metazoans"))

p = ggboxplot(alpha_div2_melt, x = "region", y = "value", 
              order = c("Lombok","Misool", "Waigeo"), color = "region",palette = c("#CC3311", "#009988", "#0077BB"),
              facet.by = c("variable", "group"),
              scales = "free") + rremove("axis") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

lev <- c("Lombok","Misool", "Waigeo")
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif") + scale_y_continuous(expand = expansion(mult=c(0.05,0.15)))
p2

ggsave("alpha_div_metazoans_edna.png",p2, width=10)

####### 8. SEQUENCING DEPTH  #######
#### 8.0 ggrare() function from ranacapa ####
ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}
#### 8.1 Environmental DNA - OTUs water ####
sampling_depth_otu = ggrare(otu_water, step = 100, label = "sample") 
sampling_depth_otu_facet = sampling_depth_otu + facet_grid(region ~ .) +theme_bw() + ylab("OTU Richness") + ggtitle("a) Environmental DNA - water")

#### 8.2 Environmental DNA - OTUs sediment ####
sampling_depth_otu_s = ggrare(otu_sediment, step = 100, label = "sample")
sampling_depth_otu_s_facet = sampling_depth_otu_s + facet_grid(region ~ .) +theme_bw() + ylab("OTU Richness") + ggtitle("b) Environmental DNA - sediment")

#### 8.3 Combine into one plot #### 
sampling_depth_combo = ggarrange(sampling_depth_otu_facet, sampling_depth_otu_s_facet, legend="none")
ggsave("sampling_depth_combo_v4.png", width=10)


####### 9. SITE MAP #######
# Libraries
library(sf)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(ggsn)
library(maps)
library(ggspatial)
library(cowplot)
library(ggsflabel)

site_meta  = meta(benthic_phylo1)
site_meta$eDNA_water = NA
site_meta$eDNA_sed = NA
site_meta$UVC = NA
site_meta$Benthic = NA

# Custom E-W and N-S lat/long scales
scale_x_longitude <- function(xmin=-180, xmax=180, step=1, ...) {
  ewbrks <- seq(xmin,xmax,step)
  ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(x, "°W"), ifelse(x > 0, paste(x, "°E"),x))))
  return(scale_x_continuous( breaks = ewbrks, labels = ewlbls, expand = c(0, 0), ...))
}
scale_y_latitude <- function(ymin=-90, ymax=90, step=0.5, ...) {
  nsbrks <- seq(ymin,ymax,step)
  nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "°S"), ifelse(x > 0, paste(x, "°N"),x))))
  return(scale_y_continuous( breaks = nsbrks, labels = nslbls, expand = c(0, 0), ...))
}

#### 9.1 Fill in site_meta table with which sites had which methods sampled #### 
## 9.1.1 eDNA water samples - use summary table variables for this ## 
site_meta[site_meta$region_code %in% unique(sample_data(merged_ps_water_lbk)$region_code), ]$eDNA_water = TRUE
site_meta[site_meta$region_code %in% unique(sample_data(merged_ps_water_mi)$region_code), ]$eDNA_water = TRUE
site_meta[site_meta$region_code %in% unique(sample_data(merged_ps_water_wg)$region_code), ]$eDNA_water = TRUE

## 9.1.2 eDNA sediment samples ## 
site_meta[site_meta$region_code %in% unique(sample_data(merged_ps_sed_lbk)$region_code), ]$eDNA_sed = TRUE
site_meta[site_meta$region_code %in% unique(sample_data(merged_ps_sed_mi)$region_code), ]$eDNA_sed = TRUE
site_meta[site_meta$region_code %in% unique(sample_data(merged_ps_sed_wg)$region_code), ]$eDNA_sed= TRUE

## 9.1.3 UVC sites ## 
site_meta[site_meta$region_code %in% unique(sample_data(uvc_phylo_lbk)$region_code), ]$UVC= TRUE
site_meta[site_meta$region_code %in% unique(sample_data(uvc_phylo_mis)$region_code), ]$UVC= TRUE
site_meta[site_meta$region_code %in% unique(sample_data(uvc_phylo_wg)$region_code), ]$UVC= TRUE

## 9.1.4 Benthic sites ## 
site_meta[site_meta$region_code %in% unique(sample_data(benthic_phylo_lbk)$region_code), ]$Benthic= TRUE
site_meta[site_meta$region_code %in% unique(sample_data(benthic_phylo_mis)$region_code), ]$Benthic= TRUE
site_meta[site_meta$region_code %in% unique(sample_data(benthic_phylo_wg)$region_code), ]$Benthic= TRUE

#### 9.2 Plot regions ####
indonesia = readRDS("maps/gadm36_IDN_3_sf.rds")
site_coords = st_as_sf(site_meta, coords = c("long","lat"), crs=4326, agr = "constant", remove = FALSE)

## 9.2.1 Lombok 
lombok_coords = subset(site_coords, region == "Lombok")
lombok_map = dplyr::filter(indonesia, NAME_1 %in% c("Nusa Tenggara Barat"))
lombok_map = lombok_map[lombok_map$NAME_2 == "Lombok Barat"| lombok_map$NAME_2 =="Lombok Tengah"| lombok_map$NAME_2 == "Lombok Timur"| lombok_map$NAME_2 =="Lombok Utara",]
lombok_map = lombok_map %>%
  st_union %>%
  st_sf()

options(ggrepel.max.overlaps = Inf)

lombok_plot = ggplot() + 
  geom_sf(data = lombok_map, fill="grey70", lwd=0)+
  geom_sf(data = lombok_coords, aes(color="#CC3311",shape = region, stroke=1))+ scale_shape_manual(values = c(5,6,7,8))+
  scale_color_manual(values = c("#CC3311"))+
  scale_fill_manual(values =c("#CC3311"))+
  scale_x_continuous(breaks= pretty(lombok_coords$long,n=2))+
  annotation_scale()+
  theme(panel.grid.major = element_line(colour = "white", linetype = "solid", size = 0.2), axis.text.y = element_text(angle=45, vjust=1,hjust=1),
        panel.background = element_rect(fill = "transparent"), 
        panel.border = element_rect(fill = NA),
        rect = element_rect(fill="transparent"),plot.background = element_rect(fill = "transparent", colour = NA),
        legend.title=element_blank(), legend.position = "None")+
  coord_sf(xlim = c(116.314129-0.47,116.314129+0.47), ylim = c(-8.547707 - 0.47,-8.547707+0.47)) # find a center point and add same value up/down left/right of it



## 9.2.2 Misool
misool_coords = subset(site_coords, region == "Misool")

misool_map = dplyr::filter(indonesia, NAME_2 %in% c("Raja Ampat"))
misool_map = misool_map[misool_map$NAME_3 == "Misool" | misool_map$NAME_3 == "Misool Selatan" | misool_map$NAME_3 == "Misool Barat" | misool_map$NAME_3 == "Misool Timur",]
misool_map = misool_map %>%
  st_union %>%
  st_sf()

misool_plot = ggplot() + 
  geom_sf(data = misool_map, fill= "grey70", lwd=0)+
  geom_sf(data = misool_coords, aes(color="#009988",shape = region, stroke=1))+ scale_shape_manual(values = c(6))+
  scale_color_manual(values = c("#009988"))+
  scale_fill_manual(values =c("#009988"))+
  scale_x_longitude(xmin=130.2, xmax=130.8, step=0.2) +
  scale_y_latitude(ymin=-2.30, ymax=-1.70, step=0.25) +
  coord_sf(xlim = c(130.538775-0.37,130.538775+0.37), ylim = c(-2.148162 - 0.37,-2.148162+0.37))+ # find a center point and add same value up/down left/right of it
  theme(panel.grid.major = element_line(colour = "white", linetype = "solid", size = 0.2), axis.text.y = element_text(angle=45, vjust=1,hjust=1), 
        panel.background = element_rect(fill = "transparent"), 
        panel.border = element_rect(fill = NA),
        rect = element_rect(fill="transparent"),plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
        legend.title=element_blank(), legend.position = "None")+
  annotation_scale()

## 9.2.3 Waigeo
waigeo_coords = subset(site_coords, region == "Waigeo")

waigeo_map = dplyr::filter(indonesia, NAME_2 %in% c("Raja Ampat"))
waigeo_map = waigeo_map[waigeo_map$NAME_3 == "Waigeo Barat"| waigeo_map$NAME_3 =="Waigeo Barat Kepulau"| waigeo_map$NAME_3 == "Waigeo Selatan"| waigeo_map$NAME_3 =="Waigeo Timur" | waigeo_map$NAME_3 =="Waigeo Utara" | waigeo_map$NAME_3 =="Warwabomi"| waigeo_map$NAME_3 =="Kebar"| waigeo_map$NAME_3 =="Teluk Mayalibit" | waigeo_map$NAME_3 =="Meosmansar",]
waigeo_map = waigeo_map %>%
  st_union %>%
  st_sf()

waigeo_plot = ggplot() + 
  geom_sf(data = waigeo_map, fill="grey70",lwd=0)+
  geom_sf(data = waigeo_coords, aes(color="#0077BB",shape = region, stroke=1))+ scale_shape_manual(values = c(7))+
  scale_color_manual(values = c("#0077BB"))+
  scale_fill_manual(values =c("#0077BB"))+
  coord_sf(xlim = c(130.530069-0.33,130.530069+0.33), ylim = c(-0.565138 - 0.33,-0.565138+0.33))+ # find a center point and add same value up/down left/right of it
  scale_x_longitude(xmin=130.21, xmax=130.84, step=0.3) +
  scale_y_latitude(ymin=-0.8, ymax=-0.3, step=0.24) +
  annotation_scale()+
  theme(panel.grid.major = element_line(colour = "white", linetype = "solid", size = 0.2),  axis.text.y = element_text(angle=45, vjust=1,hjust=1),
        panel.background = element_rect(fill = "transparent"), 
        panel.border = element_rect(fill = "transparent"),
        rect = element_rect(fill="transparent"),plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
        legend.title=element_blank(), legend.position = "None")

#### 9.3 Indonesia map - inset world map ####
worldmap = map_data("world")
indonesia_fill = worldmap[worldmap$region=="Indonesia",]


indonesia_map = ggplot(data=worldmap, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill="grey95", colour = "darkgray", aes(long, lat, group=group))+
  geom_polygon(data=indonesia_fill, aes(long, lat, group=group),fill="grey84", color="black",size=0.1)+
  coord_sf(xlim= c(111,135), ylim = c(-11.2610,1), expand=FALSE)+  
  scale_x_longitude(xmin=115, xmax=130, step=5) +
  scale_y_latitude(ymin=-10, ymax=1, step=2) +
  theme(panel.grid.major = element_line(colour = "white", linetype = "solid", size = 0.2), axis.text.y = element_text(angle=45, vjust=1,hjust=1),
        panel.background = element_rect(fill = "transparent"), 
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        rect = element_rect(fill="transparent"),plot.background = element_rect(fill = "transparent", colour = NA),
        legend.title=element_blank(), legend.position = "None")+
  annotation_north_arrow()+
  geom_rect(xmin = 115.773, ymax=-8.011, xmax = 116.81, ymin = -9.23, fill=NA, color="#CC3311")+ # lombok box
  geom_rect(xmin = 129.591, ymax=-1.592, xmax = 130.690, ymin = -2.232, fill=NA, color="#009988") + # misool box
  geom_rect(xmin = 130.162, ymax=0.0877, xmax = 131.492, ymin = -0.6050, fill=NA, color="#0077BB") # waigeo box

world_map = ggplot(data = worldmap, aes(x=long, y=lat, group=group)) + 
  geom_polygon(fill="gray90",colour= "gray50", size = 0.08, aes(long, lat, group=group))+
  geom_polygon(data=indonesia_fill, fill="darkgray", size=0.05)+
  coord_sf(xlim= c(30,180), ylim = c(-50,50), expand=FALSE)+ 
  annotation_north_arrow()+
  theme_map()+ 
  theme(panel.border = element_rect(color="black",fill = NA, size=1), panel.background = element_rect(fill="white"))+
  geom_rect(xmin = 111, ymax=1, xmax = 135, ymin = -10, fill=NA, color="black")

#### 9.4 Export as pdf or arrange ####

combined_map = ggarrange(indonesia_map+ggtitle("(a) Indonesia"), ggarrange(lombok_plot+ggtitle("(b) Lombok"), misool_plot+ggtitle("(c) Misool"), waigeo_plot+ggtitle("(d) Waigeo"),ncol=3, widths = c(1,1,1)), nrow=2, heights = c(1,1))

ggsave("combined_map_v3.png", combined_map, bg="white",dpi=1000, width = 12)
ggsave("world_map_v3.png", world_map,bg="transparent", dpi=700)


