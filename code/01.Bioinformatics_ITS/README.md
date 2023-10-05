# Fungal Bioinformatics Pipeline for Illumina Forward Reads Targeting the ITS1 Region

## Overview:

- This pipeline has been optimized for [Australian Microbiome Initiative (AMI)](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) data targeting the ITS region.
- AMI ITS data are Illumina paired-end, 2x300bp, and have been demultiplexed.
- The AMI ITS protocol targets the entire ITS region, which is longer than the 2x300bp amplicon, therefore, reads cannot be merged.
- AMI data are notoriously noisy. As such, quality filtering is relatively stringent.
- Primers could not be found/removed using `Cutadapt`, probably due to quality issues. Therefore, primer trimming has not been performed. Nonetheless, reads are trimmed and truncated using `ITSxpress` and `Trimmomatic`.
- Sequencing runs are processed individually to prior to denoising with `DADA2` to improve error rate estimations, as each sequencing run has a unique error.

## Dependencies:

- All dependencies can be called by creating and activating the `shell` and `R` environments. See the [envs subdirectory](https://github.com/LukeLikesDirt/MycorrhizaAusForests/tree/main/envs) for more details.

## Data Acquisition:

- To access the [AMI](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) data used in this repository, follow the instructions within the ['../data/AusMicrobiome/ITS/BPA_*'](https://github.com/LukeLikesDirt/MycorrhizaAusForests/tree/main/data/AusMicrobiome/ITS) subdirectories. Specifically, execute the 'download.sh' scripts.





# Fungal Bioinformatics Pipeline for Illumina Forward Reads Targeting the ITS1 Region

## Overview:

- I've tailored this pipeline for data from the [Australian Microbiome Initiative (AMI)](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) that targets the ITS region.
- AMI ITS data are from Illumina paired-end sequencing, 2x300bp reads, and they've been demultiplexed.
- The AMI ITS protocol captures the entire ITS region, which is longer than the 2x300bp amplicon, so we can't merge reads.
- AMI ITS data are notoriously noisy, so I apply relatively strict quality filtering.
- The primers could not be found or remove primers using `Cutadapt`, most likely due to quality issues. Therefore, I haven't performed primer trimming. Nevertheless, reads are trimmed and truncated using `ITSxpress` and `Trimmomatic`.
- I process each sequencing run individually before denoising with `DADA2` to enhance error rate estimations because each sequencing run has its unique errors.

## Dependencies:

- You can access all the necessary dependencies by creating and activating the `shell` and `R` environments. For more details, check out the [envs subdirectory](https://github.com/LukeLikesDirt/MycorrhizaAusForests/tree/main/envs).

## Data Acquisition:

- To access the data from the [Australian Microbiome Initiative (AMI)](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) used in this repository, follow the instructions within the ['../data/AusMicrobiome/ITS/BPA_*'](https://github.com/LukeLikesDirt/MycorrhizaAusForests/tree/main/data/AusMicrobiome/ITS) subdirectories. Specifically, execute the 'download.sh' scripts.
