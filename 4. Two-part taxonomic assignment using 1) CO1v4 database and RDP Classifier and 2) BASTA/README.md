# 4. Two-part taxonomic assignment using 1) CO1v4 database and RDP Classifier and 2) BASTA 

## 4.1 CO1v4 database and RDP classifier

- Input: the `combined-asvs.fasta` file from the previous step
- Used this database: [https://github.com/terrimporter/CO1Classifier](https://github.com/terrimporter/CO1Classifier), which contains COI sequences mined from GenBank (April 2019) and the BOLD data releases
    - Citation: Porter, T.M., & Hajibabaei, M. (2018) Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226.
- Used this taxonomic assignment algorithm: Wang et al. (2007) Naïve Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Applied and Environmental Microbiology, 73: 5261.
    - Installed using conda
```bash
############ Install the RDP classifier if you need it
# The easiest way to install the RDP classifier v2.13 is using conda - I used this
conda install -c bioconda rdp_classifier
# Alternatively, you can install from SourceForge and run with java if you don't use conda
wget https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.13.zip
# decompress it
unzip rdp_classifier_2.13
# record path to classifier.jar ex. /path/to/rdp_classifier_2.13/dist/classifier.jar

############ Get the latest COI training set
wget https://github.com/terrimporter/CO1Classifier/releases/download/v4/CO1v4_trained.tar.gz

# decompress it
tar -xzf CO1v4_trained.tar.gz

# record the path to the rRNAClassifier.properties file ex. /path/to/mydata_trained/rRNAClassifier.properties
/Users/elaineshen/Desktop/cox1_raw/combined/databases/CO1v4_training/mydata_trained/rRNAClassifier.properties

############ Run the RDP Classifier 
# If it was installed using conda, run it like this:
rdp_classifier -Xmx8g classify -t /path/to/mydata_trained/rRNAClassifier.properties -o rdp.output query.fasta
rdp_classifier -Xmx8g classify -t /Users/elaineshen/Desktop/cox1_raw/combined/databases/CO1v4_training/mydata_trained/rRNAClassifier.properties -o CO1v4.output /Users/elaineshen/Desktop/cox1_raw/combined/combined-asvs.fasta
```

## 4.2 Additional taxonomic assignment with BASTA

- Followed tutorials from here: [https://github.com/timkahlke/BASTA/wiki](https://github.com/timkahlke/BASTA/wiki#basta---basic-sequence-taxonomy-annotation)

### 4.2.0 Installation

- Need to create a new conda environment to install this, because it requires python 2.7

```bash
conda create -n basta-env python=2.7
```
```bash
conda activate basta-env
```
```bash
conda install -c bioconda -c bnoon -c timkahlke basta
```

- Test to see if installation worked

```bash
python -m unittest discover basta
```

### 4.2.1 Download and create NCBI taxonomy

```bash
basta taxonomy
```

### 4.2.2 Download NCBI mapping files

- Downloading gb - genbank-to-taxonID mapping file (for most nucleotide databases)

```bash
basta download gb
```

### 4.2.3 LCA for each query sequence

- Sample code:

```bash
./bin/basta sequence INPUT_FILE OUTPUT_FILE MAPPING_FILE_TYPE
```

- My code:

```bash
basta sequence /Users/elaineshen/Desktop/cox1_raw/combined/iterative_blast/ib_output/test_seqs_0/blasted_20210510_0809_e1e-50.txt  blasted_20210510_0809_e1e-50_basta.out gb
```

### 4.2.4 Add all BLAST output into one file and run BASTA

- Concatenate all the output from `ib_output` folder

```bash
cat test_seqs_*/*.txt > all_output.txt
```

### 4.2.5 Run BASTA with changes to parameters

- Changed the following default parameters:
    - -p : the default is 100% (which means it will return taxonomy that is shared by 100% of hits). We changed this to 80% to be less conservative (this gets at reasonable - aka occurring in the Indo-Pacific - taxonomic assignments)
    - -i: identity, the default is 80 but I will change to 97% to be consistent with other taxonomic identity.

```bash
conda activate basta-env
```

```bash
basta sequence /Users/elaineshen/Desktop/cox1_raw/combined/iterative_blast/ib_output/all_output.txt  all_blast_basta.txt gb -i 97 -p 80 
```

## 4.3. Combine taxonomic assignment files

- The input files for `combine_taxassign.R` are the two files from the two-step taxonomic assignment steps: `CO1v4.output` and `all_blast_basta.txt`
- Because the `iterative_blast` outputs (from BASTA) do not have rank-level confidence values like the CO1v4 output does, they were added in as 97.01010101. This is the minimum confidence value that I used in both the iterative blast search and LCA, so at minimum, the confidence at each rank has to be 97%. This numerical designation also makes enough of a distinction to remove them down the line if desired.
- Small note: I made a global confidence column and did not include confidence values ≤0.8. This was mostly for exploratory analysis, but the final analysis is filtered for confidence values ≥0.97 (which you’ll see in the next analysis steps).
- The resulting new and combined taxonomy file is what is imported into the R package phyloseq and used for downstream analysis. I export this phyloseq object at the end for convenience (it gets re-uploaded in the next analysis step).
