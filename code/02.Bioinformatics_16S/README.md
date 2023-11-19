# Bacterial Bioinformatics Pipeline for Illumina Paired-end Reads Targeting the 16S V1-V3 Region

## Overview:

- I have tailored this pipeline for data from the [Australian Microbiome Initiative](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) that targets the 16S V1-V4 region.
- Australian Microbiome 16S data are from Illumina paired-end sequencing, 2x300bp amplicons, and they've been demultiplexed.
- This pipeline compares denoising with DADA2 and VSEARCH, which uses the UNOISE3 algorithm. 
- I process each sequencing run individually before denoising to enhance error rate estimations because each sequencing run has a unique error profile.

## Dependencies:

- You can access all the necessary dependencies by creating and activating the `shell` and `R` environments. For more details, check out the [envs subdirectory](https://github.com/LukeLikesDirt/MycorrhizaAusForests/tree/main/envs).

## Data Acquisition:

- To access the AMI data used in this repository, you will need to create a profile at [Australian Microbiome Initiative](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative).
