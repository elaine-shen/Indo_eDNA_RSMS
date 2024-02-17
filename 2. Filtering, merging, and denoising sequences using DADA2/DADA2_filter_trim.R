#######################################################################################################################################################################################
#
# AUTHORS: Elaine Shen
# DATE CREATED: 01/30/2021
# DATE MODIFIED: 02/17/2024 (modified for brevity)
# DESCRIPTION: Take trimmed primers with adapters removed and filter, merge, and denoise for taxonomy assignment
#
# INPUTS:
#    Set directory to trimmomatic output (depending on workflow, will need folder of all files or separated folders with fastq.gz or fastq files)
#    Follows the DADA2 Tutorial: https://benjjneb.github.io/dada2/tutorial.html and for big data and paired reads: https://benjjneb.github.io/dada2/bigdata_paired.html
#
######################################################################################################################################################################################

setwd("/Users/elaineshen/Desktop/cox1_raw/EB801-939/cutadapt_output")

##### LIBRARY ######
library(dada2)
library(stringr)

##### VIEW QUALITY PROFILES ######
# This uses the normal tutorial, so you need all the trimmed_paired sequences unzipped in the same directory
# Unlike in qiime, you can only view quality profiles one sample at a time
path <- "/Users/elaineshen/Desktop/cox1_raw/EB801-939/cutadapt_output/" # CHANGE ME to the directory containing the fastq files after unzipping.
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# Plot quality profiles
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

# re-zip and move files into respective trimmed_paired folders to continue with the big data analysis

##### FILTER ######
# This is using the big data workflow, so you need to move the F and R files into their respective directories
# File parsing
pathF <- "/Users/elaineshen/Desktop/cox1_raw/EB801-939/cutadapt_output/trimmed_F" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "/Users/elaineshen/Desktop/cox1_raw/EB801-939/cutadapt_output/trimmed_R" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
out = filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                    rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                    truncLen=c(190,140), trimRight=c(0,5), rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE, multithread=TRUE)

##### INFER SEQUENCE VARIANTS ######
# File parsing
filtpathF <- "/Users/elaineshen/Desktop/cox1_raw/EB801-939/cutadapt_output/trimmed_F/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "/Users/elaineshen/Desktop/cox1_raw/EB801-939/cutadapt_output/trimmed_R/filtered" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

##### ERROR RATES ######

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE) 

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR,maxMismatch=1)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

##### HAMMING DISTANCE ######
# See: https://github.com/Talitrus/bocas_eDNA for explanation
mergers1 = mergers # make a duplicate of mergers to do this step.
refseq = "ATTATCTGGGATTCAGGCTCATTCCGGGGGTTCGGTAGATTTGGTTATTTTTAGTTTACATTTAGCGGGTATTTCTTCTATATTGGCGGCTATGAATTTTATAACCACTATTATTAATATGAGGGCACCAGGGATAACAATGGATAGAACGCCATTGTTTGTTTGGTCAATTTTAGTAACTGCGGTTTTATTATTATTATCTTTACCAGTATTAGCAGGCGCAATTACTATGTTATTAACGGATAGAAATTTTAATACTGCTTTTTTTGATCCAGCAGGTGGAGGAGACCCGATTTTATATCAACATTTATTT"
for(sam in names(mergers1)) {
  seq.hd <- nwhamming(substr(refseq,1,70), substr(mergers1[[sam]][,1],1,70), vec=TRUE, endsfree=F)
  seq.rc <- nwhamming(substr(dada2:::rc(refseq),1,70), substr(mergers1[[sam]][,1],1,70), vec=TRUE, endsfree=F)
  seq.df <- seq.rc - seq.hd #positive = original orientation, negative = reverse orientation
  hist(seq.df, breaks=100)
  
  seq.r <- which(seq.df < 0)
  mergers1[[sam]][seq.r,1] <- dada2:::rc(mergers1[[sam]][seq.r,1])   
}

# If the histogram does not have negative values, then you don't have sequences that need to be re-oriented in the forward direction
# If it does, use mergers1 in the bottom code

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers1)
seqtab_table = table(nchar(getSequences(seqtab))) # Examine length of sequences
barplot(seqtab_table, ylab = "frequency", xlab = "seq length")

# Save seqtab to merge with EB100 samples
saveRDS(seqtab, "/Users/elaineshen/Desktop/cox1_raw/EB801-939/EB801-939-seqtab.rds")