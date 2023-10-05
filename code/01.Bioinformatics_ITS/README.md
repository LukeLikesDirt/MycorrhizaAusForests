# Fungal bioinformatics pipeline for Illumina forward reads targeting the ITS1 region:

## Overview:
 - This pipline has been optimised for [Australian Micropbiome Initiative (AMI)](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) dataset targeting the ITS (fungal) region.
 - AMI ITS data are Illumia paired-end, 2x300bp, and have been demultiplexed.
 - The AMI ITS protocol targets the entire ITS region, which is longer than the 2x300bp amplicon, therefore reads cannot be merged.
 - AMI data are netroiously noisy. Therefore quality filtering is relatively stringent.
 - Primer could not be found/removed using `Cutadapt`, probably due to quality quality issues. Therefore primer trimming has not been performed. Nonetheless reads are trimmed and truncated using `ITSxpress` and `Trimmomatic`.
 - Sequencing runs are processed induvidually to prior to denoising with `DADA2`
to improve error rate estimations, as each sequencing run has a unique error.

## Data aquisition:
 - To access the [AMI](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) data used in this repository, follow the instructions within the ['../data/AusMicrobime/ITS/BPA_*'](https://github.com/LukeLikesDirt/MycorrhizaAusForests/tree/main/data/AusMicrobiome/ITS) subdirectories.
