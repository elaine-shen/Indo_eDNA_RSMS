# 1. Trimming adapters and primers using Trimmomatic and Cutadapt

## 1.0 Overview:

- There are two MiSeq runs of water samples, the ones from 2018 (Lombok, Waigeo) and the one from Misool in 2019. The previous trips are denoted as EB100-196 samples and the Misool 2019 collections are titled EB801-939
    - These samples contain different filter sizes (12 um and 0.4 um), sample types (water, sediment, gut), but for the COI Leray primers
    - Currently, they are in `/Dekstop/cox1_raw/` in the folders `EB100-194` and `EB801-939`
        - They were processed bioinformatically as separate runs until taxonomy assignment, when the ASV table was combined. This is because the quality of sequences may be run-specific. This is also the recommended treatment of data from the DADA2 tutorial.
- Sequences from this study can be found here in the sequence read archive (SRA): https://www.ncbi.nlm.nih.gov/sra/PRJNA1076664

## 1.1 EB100-194 samples

### 1.1.0 Initial folder organization

```bash
cd /Users/elaineshen/Desktop/cox1_raw/EB101-194/
```

- This folder should contain a folder `raw` with the raw sequences as `fastq.gz` files with the following naming convention:`EB109_R1.fastq.gz` and `EB109_R2.fastq.gz`

### 1.1.1 Trimmomatic for adapter removal

- Download Trimmomatic-0.39 binary here: [http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) into the `EB101-194` folder
    - I used `Trimmomatic-0.39`
    - The output will be in the folder `Trimmomatic-0.39`
- The general structure of the Trimmomatic code is as follows:

```bash
java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] [-trimlog <logFile>] <input 1> <input 2> <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <step 1>
```

- Execute a for-loop for Trimmomatic while in the `raw` folder (outputs will be in the `Trimmomatic-0.39`) folder
    - On command line, the intermediate outputs will not be legible

```bash
for i in *_R1.fastq.gz
do
  SAMPLE=$(echo ${i} | sed "s/_R1.fastq\.gz//") 
  echo ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz
java -jar /Users/elaineshen/Desktop/cox1_raw/EB101-194/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 ${SAMPLE}_R1.fastq.gz  ${SAMPLE}_R2.fastq.gz /Users/elaineshen/Desktop/cox1_raw/EB101-194/Trimmomatic-0.39/trimmomatic_output/${SAMPLE}_paired_R1.fastq.gz /Users/elaineshen/Desktop/cox1_raw/EB101-194/Trimmomatic-0.39/trimmomatic_output/${SAMPLE}_unpaired_R1.fastq.gz /Users/elaineshen/Desktop/cox1_raw/EB101-194/Trimmomatic-0.39/trimmomatic_output/${SAMPLE}_paired_R2.fastq.gz /Users/elaineshen/Desktop/cox1_raw/EB101-194/Trimmomatic-0.39/trimmomatic_output/${SAMPLE}_unpaired_R2.fastq.gz ILLUMINACLIP:/Users/elaineshen/Desktop/cox1_raw/EB101-194/Trimmomatic-0.39/trimmomatic-0.39.jar:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50
done
```

- Move unpaired reads into their own folder

```bash
# In trimmomatic_output
mv *unpaired* unpaired
```

### 1.1.2 Check adapter removal using MultiQC

- MultiQC allows you to check the quality of all of your samples at once, providing information about the number of reads per sample, average sequence length, and adapter contamination
    - It can be downloaded here: [https://multiqc.info/](https://multiqc.info/)
- Download MultiQC into the `EB101-194` folder to execute (here it is called FastQC)
- The first for-loop makes a MultiQC for individual samples
- The second makes a summary MultiQC file for all of the samples and allows you to scroll interactively in a .html file

```bash
cd FastQC

for file in /Users/elaineshen/Desktop/cox1_raw/EB101-194/Trimmomatic-0.39/trimmomatic_output/*.fastq.gz
do
   ./fastqc $file 
done
```

```bash
cd trimmomatic_output
```
```bash
multiqc . --interactive
```

### 1.1.3 Cutadapt for paired-end primer removal

- Download cutadapt here using conda: [https://cutadapt.readthedocs.io/en/stable/installation.html](https://cutadapt.readthedocs.io/en/stable/installation.html) and put `cutadapt` folder into `EB101-194` folder

```bash
conda activate cutadapt
```
```bash
cd trimmomatic_output
```

- Because we have paired-end sequences, we need to remove the reverse complement of the primers for the reverse sequences as well
- For more information, see: [https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads](https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads)

```bash
for i in *_R1.fastq.gz
do
  SAMPLE=$(echo ${i} | sed "s/_R1.fastq\.gz//") 
  echo ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz
cutadapt -q 20 --trim-n --report=minimal -g GGWACWGGWTGAACWGTWTAYCCYCC -a TGRTTYTTTGGTCACCCTGAAGTTTA -G TAAACTTCAGGGTGACCAAARAAYCA -A GGRGGRTAWACWGTTCAWCCWGTWCC  -o /Users/elaineshen/Desktop/cox1_raw/EB101-194/cutadapt_output/${SAMPLE}_trimmed_R1.fastq.gz -p /Users/elaineshen/Desktop/cox1_raw/EB101-194/cutadapt_output/${SAMPLE}_trimmed_R2.fastq.gz ${SAMPLE}_R1.fastq.gz  ${SAMPLE}_R2.fastq.gz 
done
```

- The output files will go into a folder `cutadapt/cutadapt_output` . I recommend running MultiQC again to make sure the primers were properly removed.

## 1.2 EB801-939 samples

### 1.2.0 Initial file organization

```bash
cd /Users/elaineshen/Desktop/cox1_raw/EB801-939
```

- This folder should contain a folder `raw` with the raw sequences as `fastq.gz` files with the following naming conventions: `EB801_S1_L001_R1_001.fastq.gz`
`EB801_S1_L001_R2_001.fastq.gz`
- I copied the `Trimmomatic-0.39` folder from EB100-194 folder here, so it was easier to execute
    - But I did not copy MultiQC, just used the original filpath in EB100-194 folder to run MultiQC on these samples

### 1.2.1 Trimmomatic for adapter removal

- Because the input file naming convention is slightly different, the above for-loop was modified:

```bash
for i in *_R1_001.fastq.gz
do
  SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq\.gz//") 
  echo ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz
java -jar /Users/elaineshen/Desktop/cox1_raw/EB801-939/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 ${SAMPLE}_R1_001.fastq.gz  ${SAMPLE}_R2_001.fastq.gz /Users/elaineshen/Desktop/cox1_raw/EB801-939/Trimmomatic-0.39/trimmomatic_output/${SAMPLE}_paired_R1.fastq.gz /Users/elaineshen/Desktop/cox1_raw/EB801-939/Trimmomatic-0.39/trimmomatic_output/${SAMPLE}_unpaired_R1.fastq.gz /Users/elaineshen/Desktop/cox1_raw/EB801-939/Trimmomatic-0.39/trimmomatic_output/${SAMPLE}_paired_R2.fastq.gz /Users/elaineshen/Desktop/cox1_raw/EB801-939/Trimmomatic-0.39/trimmomatic_output/${SAMPLE}_unpaired_R2.fastq.gz ILLUMINACLIP:/Users/elaineshen/Desktop/cox1_raw/EB801-939/Trimmomatic-0.39/trimmomatic-0.39.jar:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50
done
```

- After this step, you can use MultiQC (as demonstrated above) to check that adapters were removed.

### 1.2.3 Cutadapt for paired-end primer removal

- Make a directory in `EB801-939` called `cutadapt_output`where the output files are stored

```bash
conda activate cutadapt
```
```bash
cd trimmomatic_output
```

```bash
for i in *_R1.fastq.gz
do
  SAMPLE=$(echo ${i} | sed "s/_R1.fastq\.gz//") 
  echo ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz
cutadapt -q 20 --trim-n --report=minimal -g GGWACWGGWTGAACWGTWTAYCCYCC -a TGRTTYTTTGGTCACCCTGAAGTTTA -G TAAACTTCAGGGTGACCAAARAAYCA -A GGRGGRTAWACWGTTCAWCCWGTWCC  -o /Users/elaineshen/Desktop/cox1_raw/EB801-939/cutadapt_output/${SAMPLE}_trimmed_R1.fastq.gz -p /Users/elaineshen/Desktop/cox1_raw/EB801-939/cutadapt_output/${SAMPLE}_trimmed_R2.fastq.gz ${SAMPLE}_R1.fastq.gz  ${SAMPLE}_R2.fastq.gz 
done
```

- Once more, the output files will go into a folder `cutadapt/cutadapt_output` . I recommend running MultiQC again to make sure the primers were properly removed.
