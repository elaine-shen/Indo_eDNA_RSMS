##### MERGE RUNS, REMOVE CHIMERAS ######
# Merge multiple runs (if necessary)
st_EB100 <- readRDS("/Users/elaineshen/Desktop/cox1_raw/EB101-194/EB101-194-seqtab.rds")
st_EB800 <- readRDS("/Users/elaineshen/Desktop/cox1_raw/EB801-939/EB801-939-seqtab.rds")
st.all <- mergeSequenceTables(st_EB100, st_EB800)

# Remove chimeras, final sanity check

seqtab_nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE) # seqtab would be st.all if multiple runs were merged
sum(seqtab_nochim)/sum(st.all) 
# sum(seqtab_nochim)/sum(seqtab) #see how many reads remain after chimera removal - 

# Track reads through pipeline (will be run-specific)
getN <- function(x) sum(getUniques(x))
track = cbind(out, rowSums(seqtab),rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtering", "merging", "chimera")
rownames(track) <- sample.names

# Create a fasta file for the unique sequences and save final combined sequence table
uniquesToFasta(getUniques(seqtab_nochim), fout="/Users/elaineshen/Desktop/cox1_raw/combined/uniqueSeqs.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab_nochim)))))
saveRDS(seqtab_nochim, "/Users/elaineshen/Desktop/cox1_raw/combined.rds") # CHANGE ME to where you want sequence table saved

#### MAKE ASV FASTA FILE FOR TAXONOMIC ASSIGNMENT ####
#seqtab_nochim = readRDS("/Users/elaineshen/Desktop/cox1_raw/combined/combined.rds")
sq <- getSequences(seqtab_nochim)
asv_headers <- vector(dim(seqtab_nochim)[2], mode="character")
for (i in 1:dim(seqtab_nochim)[2]) {
  asv_headers[i] <- paste("ASV", i, sep="_")
}
names(sq) <- asv_headers
library(ShortRead)
writeFasta(sq, file="combined-asvs.fasta")

# Make taxonomy table with ASVID and sequence 
library (devtools)
library (tidyverse)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")
FastaToTabular("combined-asvs.fasta")