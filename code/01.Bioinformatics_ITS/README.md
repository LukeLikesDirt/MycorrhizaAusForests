# Fungal Bioinformatics Pipeline for Illumina Forward Reads Targeting the ITS1 Region

## Overview:

- I have tailored this pipeline for data from the [Australian Microbiome Initiative (AMI)](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) that targets the ITS region.
- AMI ITS data are from Illumina paired-end sequencing, 2x300bp amplicons, and they've been demultiplexed.
- The AMI ITS protocol captures the entire ITS region, which is longer than the 2x300bp amplicon, so we can't merge reads. As a result, we cannot merge reads, so I am processing only the R1 reads here.
- AMI ITS data are notoriously noisy, so I apply relatively strict quality filtering.
- The primers could not be found or remove primers using `Cutadapt`, most likely due to quality issues. Therefore, I haven't performed primer trimming. Nevertheless, reads are trimmed and truncated using `ITSxpress` and `Trimmomatic`.
- I process each sequencing run individually before denoising with `DADA2` to enhance error rate estimations because each sequencing run has a unique error profile.

## Dependencies:

- You can access all the necessary dependencies by creating and activating the `shell` and `R` environments. For more details, check out the [envs subdirectory](https://github.com/LukeLikesDirt/MycorrhizaAusForests/tree/main/envs).

## Data Acquisition:

- To access the AMI data used in this repository, follow the instructions within the ['../data/AusMicrobiome/ITS/BPA_*'](https://github.com/LukeLikesDirt/MycorrhizaAusForests/tree/main/data/AusMicrobiome/ITS) subdirectories. Specifically, execute the 'download.sh' scripts.
