# 5. Ecological and statistical analyses
## Overview

- There are a few inputs for this data:
    - The taxonomy-dependent input of ASVs is `CO1v4_rdp_raw.rds`, which is uploaded here as a phyloseq file.
    - The taxonomy-independent input of OTUs is `LULU_curation_otu_97.rds`, also uploaded here as a phyloseq file.
- Please reach out to me directly for the raw visual survey data, benthic data, and the site metadata associated with them. Briefly, these data were converted into species x site count tables and combined with site metadata so they could be imported into R as phyloseq files. This made the eDNA and visual survey analyses parallel. These data are summarized as taxonomic checklists (with aggregated abundances) in the Supplementary Materials.
- The `main_analysis.R` file is organized by the type of analysis, not figure #, as some main-text and supplemental figures/tables were generated at the same time. These are numbered 1-9 with headers you can navigate to in RStudio and are as follows:
    1. IMPORT AND PREPROCESS DATA
    2. SUMMARY TABLE
    3. RAREFACTION CURVES
    4. SETTING COLORS FOR TAXONOMY BARPLOTS
    5. TAXONOMY BARPLOTS - RELATIVE ABUNDANCE
    6. ORDINATIONS/PERMANOVAS
    7. ALPHA DIVERSITY
    8. SEQUENCING DEPTH
    9. SITE MAP
