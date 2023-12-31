# Fungal Bioinformatics Pipeline for Illumina Forward Reads Targeting the ITS1 Region

## Overview:

- I have tailored this pipeline for data from the [Australian Microbiome Initiative](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) that targets the ITS region.
- Australian Microbiome ITS data are from Illumina paired-end sequencing, 2x300bp amplicons, and they've been demultiplexed.
- The Australian Microbiome ITS protocol captures the entire ITS region, which is longer than the 2x300bp amplicon, so we can't merge reads. As a result, we cannot merge reads, so I am processing only the R1 reads here.
- Australian Microbiome ITS data are notoriously noisy, so I apply relatively strict quality filtering.
- The primers could not be found or removed using `Cutadapt`, most likely due to quality issues. Therefore, I haven't strictly performed primer trimming, yet reads are trimmed and truncated using `ITSxpress` and `Trimmomatic`.
- This pipeline compares denoising with DADA2 and VSEARCH, which uses the UNOISE3 algorithm. 
- I process each sequencing run individually before denoising to enhance error rate estimations because each sequencing run has a unique error profile.

## Dependencies:

- You can access all the necessary dependencies by creating and activating the `shell` and `R` environments. For more details, check out the [envs subdirectory](https://github.com/LukeLikesDirt/MycorrhizaAusForests/tree/main/envs).

## Data Acquisition:

- To access the AMI data used in this repository, you will need to create a profile at [Australian Microbiome Initiative](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative).
