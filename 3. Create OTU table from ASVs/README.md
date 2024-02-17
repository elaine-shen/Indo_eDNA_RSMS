## 3.0 Overview

- This OTU table will be used for taxonomy-independent analysis down the line.
- Following scripts from [https://github.com/Talitrus/BCS_dada2/blob/master/09_vsearch_cluster.sh](https://github.com/Talitrus/BCS_dada2/blob/master/09_vsearch_cluster.sh)
- This was run locally.
- vsearch clusters ASVs into OTUs and LULU gets rid of erroneous OTUs. LULU can also curate ASV’s.
    - Here, LULU was run on ASVs to get rid of any erroneous ASVs. In our case, it did not get rid of any ASV’s. Nevertheless, we used that LULU output of ASV’s as input into vsearch, then ran LULU again on the OTUs.
    - Justfication:"The LULU algorithm was then used to curate ASV results, as it reduces the number of erroneous ASVs by merging “daughter” ASVs with consistently co-occurring more abundant “parent” ASVs (Frøslev et al., 2017)." from [https://onlinelibrary.wiley.com/doi/full/10.1002/edn3.199](https://onlinelibrary.wiley.com/doi/full/10.1002/edn3.199)

## 3.0. Setup and install vsearch and LULU

- LULU can be run in R: [https://github.com/tobiasgf/lulu](https://github.com/tobiasgf/lulu)
- vsearch can be executed in conda: [https://github.com/torognes/vsearch](https://github.com/torognes/vsearch)
- Inputs are ASV table and taxonomy table

## 3.1. Preparing ASV table for LULU and vsearch

- The R code `LULU_ASVs.R` embeds sequence/sample information in the data so that the fasta file can be made back into the original ASV table (with abundances and sample info preserved)
- Then finish out the rest of the script in terminal in the `lulu_vsearch` folder. The final output will be `concat.fasta`, which will be used in downstream analysis (vsearch)

```r
cat DADA2_extracted_samples/*.fas > DADA2_extracted_samples/concat.fasta
```

## 3.2. Make ASVs into OTUs using vsearch

- Dereplicate

```r
vsearch --derep_fulllength DADA2_extracted_samples/concat.fasta --sizein --sizeout --fasta_width 0 --minuniquesize 1 --output combined_derep_vsearch.fasta
```

- Get rid of extra semicolon added by VSearch

```r
sed -i -E "s/;.+;/;/g" combined_derep_vsearch.fasta
```

- Sort by size

```r
# Sort by size
vsearch --sortbysize combined_derep_vsearch.fasta --sizein --sizeout --output combined_derep_sorted_vsearch.fasta --fasta_width 0
```

- VSEARCH chimera checking

```r
vsearch --uchime_denovo combined_derep_sorted_vsearch.fasta --sizein --sizeout --qmask none --nonchimeras combined_derep_sorted_chimera.fasta 
```

- VSEARCH clustering at 97%

```r
# VSEARCH clustering
vsearch --cluster_size combined_derep_sorted_chimera.fasta --id 0.97 --sizein --centroids vsearch_otus_97.fasta --qmask none
```

```r
sed -i -E "s/;size=.*//g" vsearch_otus_97.fasta
```

- Map against original fasta file

```r
# Map raw reads against OTU representatives concat.fasta is from LULU_prep.sh, which should be run before this.
vsearch --usearch_global DADA2_extracted_samples/concat.fasta --db vsearch_otus_97.fasta --id 0.97 --maxaccepts 0 --dbmask none --qmask none --uc FINAL_FOR_LULU.uc
```

- Convert from UC format to OTU table format

```r
# Then convert UC format to OTU table format with uc2otutab.py
python drive5_py/uc2otutab.py FINAL_FOR_LULU.uc > FINAL_FOR_LULU_otutab.txt
```

- Before running this code, which was built in Python 2, you need to change the syntax to Python 3 using the command `2to3`, which is already built in on the command line. Do this for the entire folder `drive5_py`

```r
2to3 -w drive5_py/*.py
```

- Also had to change `time.clock()` command in `progress.py` file to `time.process_time()`
- Delete intermediate files `combined_derep_sorted_chimera.fasta` , `combined_derep_sorted_vsearch.fasta`

## 3.3. Quality check OTUs using LULU

- Input file: `FINAL_FOR_LULU_otutab.txt`
- Make matchlist for OTUs

```r
vsearch --usearch_global vsearch_otus_97.fasta --db vsearch_otus_97.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
```

- Run LULU in R using `LULU_OTUs.R`
