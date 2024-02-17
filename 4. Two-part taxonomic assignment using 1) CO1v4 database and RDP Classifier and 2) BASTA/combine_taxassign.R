#######################################################################################################################################################################################
#
# AUTHORS: Elaine Shen
# DATE CREATED: 04/16/2021
# DATE MODIFIED: 02/17/2024 (modified for brevity)
# DESCRIPTION: Combine taxonomic assignment files and create phyloseq object for downstream analysis
#
# INPUTS:
#   CO1v4.output : output from the RDP classifier
#   all_blast_basta.txt: output from BASTA (iterative BLAST)
#   combined-sample-data.txt: eDNA sample metadata
######################################################################################################################################################################################

setwd("/Users/elaineshen/Desktop/cox1_raw/combined/databases/CO1v4_training/mydata_trained")

##### LIBRARIES ######
library(dplyr)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(microbiome) # accessing more of phyloseq
library(vegan) # rarefaction curves
library(metagMisc) # format for iNEXT
library(iNEXT) # rarefaction curves


##### IMPORT RAW DATA ######

## Taxonomy files
seqs_raw = read.delim("CO1v4.output",header=FALSE)
colnames(seqs_raw) = c("name","V2", "root", "root2","root_confidence","domain_name","domain_rank", "domain_confidence", 
                       "kingdom_name", "kingdom_rank", "kingdom_confidence",
                       "phylum_name", "phylum_rank", "phylum_confidence",
                       "class_name", "class_rank","class_confidence",
                       "order_name", "order_rank", "order_confidence",
                       "family_name", "family_rank","family_confidence",
                       "genus_name", "genus_rank", "genus_confidence",
                       "species_name", "species_rank","species_confidence")

# Iterative BLAST output with LCA 
ib_output = read.delim("/Users/elaineshen/Desktop/cox1_raw/combined/iterative_blast/ib_output/all_blast_basta.txt", header=F)
ib_output = ib_output[!duplicated(ib_output[c(1,2)]),] # Remove duplicate assignments for an ASV if they are the same
ib_output = ib_output %>% filter(!str_detect(V2, 'Unknown')) # Get rid of "Unknown" IB outputs
colnames(ib_output) = c("name", "tax")

ib_output = separate(data = ib_output, col = "tax", into = c("domain_name", "phylum_name","class_name","order_name","family_name","genus_name","species_name"), sep = "\\;")
ib_output[ ,c("domain_confidence","phylum_confidence","class_confidence","order_confidence", "family_confidence","genus_confidence","species_confidence")] <- 0.9701010101 # Minimum possible confidence for the iterative blast and the LCA was 97, so set all of these rank confidences as 0.9701010101 to separate these from RDP outputs
ib_output$domain_rank = "superkingdom"
ib_output$phylum_rank = "phylum"
ib_output$class_rank = "class"
ib_output$order_rank = "order"
ib_output$family_rank = "family"
ib_output$genus_rank = "genus"
ib_output$species_rank = "species"

# Add kingdom info as undef_domain
ib_output$kingdom_name = paste0("undef_",ib_output$domain_name)
ib_output$kingdom_rank = "kingdom"
ib_output$kingdom_confidence = 0.9701010101 # This output did not provide kingdom-level information, only domain. But will assume 97% to be conservative

# Add misc columns for ease of merging
ib_output[,c("root","root2")] = "cellularOrganisms"
ib_output$V2 = ""
ib_output$root_confidence = 1 # all should be 1

# Reorder columns
col_order = c("name","V2", "root", "root2","root_confidence","domain_name","domain_rank", "domain_confidence", 
              "kingdom_name", "kingdom_rank", "kingdom_confidence",
              "phylum_name", "phylum_rank", "phylum_confidence",
              "class_name", "class_rank","class_confidence",
              "order_name", "order_rank", "order_confidence",
              "family_name", "family_rank","family_confidence",
              "genus_name", "genus_rank", "genus_confidence",
              "species_name", "species_rank","species_confidence")
ib_output = ib_output[ ,col_order]


# Replace the CO1v4 taxonomy outputs with the iterative blast output
seqs_raw[match(ib_output$name,seqs_raw$name), ] <- ib_output

## ASV table
asv_table = as.data.frame(readRDS("/Users/elaineshen/Desktop/cox1_raw/combined/combined.rds"))

asv_tax = read.csv("/Users/elaineshen/Desktop/cox1_raw/combined/asv_table.csv", header=TRUE)
asv_tax$name<-gsub(">","",as.character(asv_tax$name))
asv_tax = asv_tax[,c(2,3)]

# Need to replace asv table's columns, which are the sequences, with the ASV ID. I could also add the sequence to the taxonomy table (but that's less clean)
names(asv_table)[match(asv_tax[,"sequence"], names(asv_table))] = asv_tax[,"name"]

##### CREATING GLOBAL CONFIDENCE/TAXONOMY COLUMN, REMOVING IDs WITH CONFIDENCE < 0.8 ######

##### SHIFT IDS ######

# Inspect to see if any of the ranks are in the wrong spot due to missing taxonomic annotations - none are missing
rows_to_change_p <- !seqs_raw$phylum_rank %in% "phylum"
length(rows_to_change_p[rows_to_change_p==TRUE]) 

rows_to_change_c <- !seqs_raw$class_rank %in% "class"
length(rows_to_change_c[rows_to_change_c==TRUE]) 

rows_to_change_o <- !seqs_raw$order_rank %in% "order"
length(rows_to_change_o[rows_to_change_o==TRUE]) 

rows_to_change_f <- !seqs_raw$family_rank %in% "family"
length(rows_to_change_f[rows_to_change_f==TRUE])

rows_to_change_g <- !seqs_raw$genus_rank %in% "genus"
length(rows_to_change_g[rows_to_change_g==TRUE]) 

##### FILTER BY CONFIDENCE ######

## Remove any identifications that have a confidence < 0.8. Go from lowest taxonomic rank (species) to highest (phylum)

CO1v4_80 = seqs_raw %>% 
  mutate_at(c("species_name","species_confidence", "species_rank"), ~ifelse(species_confidence < 0.8,NA,.)) %>%
  mutate_at(c("genus_name","genus_confidence", "genus_rank"), ~ifelse(genus_confidence < 0.8,NA,.)) %>%
  mutate_at(c("family_name","family_confidence","family_rank"), ~ifelse(family_confidence < 0.8,NA,.)) %>%
  mutate_at(c("order_name","order_confidence","order_rank"), ~ifelse(order_confidence < 0.8,NA,.)) %>%
  mutate_at(c("class_name","class_confidence","class_rank"), ~ifelse(class_confidence < 0.8,NA,.)) %>%
  mutate_at(c("phylum_name","phylum_confidence","phylum_rank"), ~ifelse(phylum_confidence < 0.8,NA,.))%>%
  mutate_at(c("kingdom_name","kingdom_confidence","kingdom_rank"), ~ifelse(kingdom_confidence < 0.8, NA,.))


## And then, if there are blank names (e.g., not NA), change to NA
CO1v4_80 = CO1v4_80%>%
  mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))

CO1v4_80_only_d = seqs_raw[seqs_raw$domain_confidence >= 0.8,]
CO1v4_80_only_k = seqs_raw[seqs_raw$kingdom_confidence >= 0.8,]
CO1v4_80_only_p = seqs_raw[seqs_raw$phylum_confidence >= 0.8,]
CO1v4_80_only_c = seqs_raw[seqs_raw$class_confidence >= 0.8,]
CO1v4_80_only_f = seqs_raw[seqs_raw$family_confidence >= 0.8,]
CO1v4_80_only_g = seqs_raw[seqs_raw$genus_confidence >= 0.8,]
CO1v4_80_only_s = seqs_raw[seqs_raw$species_confidence >= 0.8,]

##### MAKE GLOBAL ID AND CONFIDENCE COLUMN ######
# This will be the lowest taxonomic rank and confidence of an ASV ID in its own column 
# Makes it easier to filter later
# Also, distinguish lowest taxonomic rank 
CO1v4_80$confidence <- ifelse(!is.na(CO1v4_80$species_confidence), CO1v4_80$species_confidence,
                              ifelse(!is.na(CO1v4_80$genus_confidence), CO1v4_80$genus_confidence,
                                     ifelse(!is.na(CO1v4_80$family_confidence), CO1v4_80$family_confidence,
                                            ifelse(!is.na(CO1v4_80$order_confidence), CO1v4_80$order_confidence,
                                                   ifelse(!is.na(CO1v4_80$class_confidence), CO1v4_80$class_confidence,
                                                          ifelse(!is.na(CO1v4_80$phylum_confidence), CO1v4_80$phylum_confidence,
                                                                 ifelse(!is.na(CO1v4_80$kingdom_confidence), CO1v4_80$kingdom_confidence,
                                                                        ifelse(!is.na(CO1v4_80$domain_confidence), CO1v4_80$domain_confidence,  NA))))))))

CO1v4_80$id <- ifelse(!is.na(CO1v4_80$species_name), CO1v4_80$species_name,
                      ifelse(!is.na(CO1v4_80$genus_name), CO1v4_80$genus_name,
                             ifelse(!is.na(CO1v4_80$family_name), CO1v4_80$family_name,
                                    ifelse(!is.na(CO1v4_80$order_name), CO1v4_80$order_name,
                                           ifelse(!is.na(CO1v4_80$class_name), CO1v4_80$class_name,
                                                  ifelse(!is.na(CO1v4_80$phylum_name), CO1v4_80$phylum_name,
                                                         ifelse(!is.na(CO1v4_80$kingdom_name), CO1v4_80$kingdom_name,
                                                                ifelse(!is.na(CO1v4_80$domain_name), CO1v4_80$domain_name, NA))))))))

CO1v4_80$rank <- ifelse(!is.na(CO1v4_80$species_rank), CO1v4_80$species_rank,
                        ifelse(!is.na(CO1v4_80$genus_rank), CO1v4_80$genus_rank,
                               ifelse(!is.na(CO1v4_80$family_rank), CO1v4_80$family_rank,
                                      ifelse(!is.na(CO1v4_80$order_rank), CO1v4_80$order_rank,
                                             ifelse(!is.na(CO1v4_80$class_rank), CO1v4_80$class_rank,
                                                    ifelse(!is.na(CO1v4_80$phylum_rank), CO1v4_80$phylum_rank,
                                                           ifelse(!is.na(CO1v4_80$kingdom_rank), CO1v4_80$kingdom_rank,
                                                                  ifelse(!is.na(CO1v4_80$domain_rank), CO1v4_80$domain_rank,  NA))))))))

# Reorder these columns to the front
CO1v4_80 = CO1v4_80%>%
  dplyr::select(c("id","rank", "confidence"), everything())

# Set ASV as row name
rownames(CO1v4_80) <- CO1v4_80[,"name"]
CO1v4_80 = CO1v4_80[,-4]

##### IMPORT PHYLOSEQ ######
metadata = read.delim("/Users/elaineshen/Desktop/cox1_raw/combined/combined-sample-data.txt")
merged_phyloseq_raw = phyloseq(otu_table(asv_table, taxa_are_rows=FALSE), 
                               sample_data(metadata), 
                               tax_table(as.matrix(CO1v4_80)))

# exporting this phyloseq object to be accessed easier - 07/06/2021
# saveRDS(merged_phyloseq_raw, "CO1v4_rdp_raw.rds")