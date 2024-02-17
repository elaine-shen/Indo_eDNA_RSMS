# 2. Filtering, merging, and denoising sequences using DADA2
## 2.0 Overview

- This was done using the DADA2 package in R with reference to the following tutorials:
    - [https://benjjneb.github.io/dada2/tutorial.html](https://benjjneb.github.io/dada2/tutorial.html)
    - [https://benjjneb.github.io/dada2/bigdata.html](https://benjjneb.github.io/dada2/bigdata.html)
- The DADA2 code for filtering and trimming sequences (`DADA2_filter_trim.R`) is the same for both runs, but run twice. Some specific parameters were changed based on the error profiles throughout. At the end, both the ASV sequence tables were combined into one for denoising.
    - Specifically, the truncation lengths and parameters for `filterAndTrim()` are different between EB101-194 and EB801-939 since the qualie runs (and subsequently, sequences) varied slightly
- For EB100-194 samples, the following params were used for `filterAndTrim()`:

```r
out = filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                    rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                    truncLen=c(190,135), trimRight=c(0,13), rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE, multithread=TRUE)
```

- For EB801-939 samples, the following params were used for `filterAndTrim()`:

```r
out = filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                    rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                    truncLen=c(190,140), trimRight=c(0,5), rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE, multithread=TRUE)
```

- Following what was done by Bryan Nguyen in a previous eDNA paper (https://github.com/Talitrus/bocas_eDNA), to correct for the tendency of the adapter-ligation library preparation method to cause random read directions, sequence variants were re-oriented into the forward direction by reverse-complementing if the reverse complement had a better Hamming distance than the original sequence after Needleman-Wunsch alignment to a reference sequence from a 313-bp section of the COI gene from a Tedania ignis sequence, GenBank accession number DQ133904.1. **
    - If the histogram does not have negative values, then you don't have sequences that need to be re-oriented in the forward direction
    - If it does, then use the variable `mergers1`in the code that follows instead of `mergers`. In this study, `mergers1` was used for EB101-194 and `mergers` was used for EB801-939.
- The code to merge and continue with DADA2 denoising steps is in `DADA2_merge.R` and includes the last formatting step needed to get the final DADA2 output (here the variable `seqtab_nochim`) into the proper format for taxonomic assignment.
