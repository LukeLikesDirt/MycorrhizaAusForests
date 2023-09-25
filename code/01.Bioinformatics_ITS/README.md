# Fungal bioinformatics pipeline for Illumina forward reads targeting the ITS1 region:

## Data aquisition:
 - To access the [Australian Micropbiome Initiative (AMI)](https://data.bioplatforms.com/organization/about/australian-microbiome-initiative) data used here, follow the instructions within the ['../data/AusMicrobime/ITS/BPA_*'](https://github.com/LukeLikesDirt/MycorrhizasAustralianForests/tree/main/data/AusMicrobiome/ITS) folders.
 - At the time of processing the UNITE reference datasets could not be accessed using `wget` or `curl` and were alternitively acceesed from here for [chimera detection](https://dx.doi.org/10.15156/BIO/2483933) and here for [taxanomic assignement](https://dx.doi.org/10.15156/BIO/2938069)

 ## General comments:
 - This script has been optimised for the AMI ITS (fungal) dataset
 - AMI ITS data are Illumia paired-end, 2x300bp, and have been demultiplexed.
 - The AMI ITS protocol targets the entire ITS region, which is longer than the 2x300bp amplicon, therefore reads cannot be merged.
 - AMI data are netroiously noisy, and require stringent quality filtering.
 - Primer could not be found or removed using `Cutadapt`, probably due to quality quality issues. Therefore primer trimming has not been performed despite representing an important quality filtering step. Nonetheless reads are trimmed and truncated using ITSxpress.
 - To improve error rate approximations (that is, to improve the ability of the `DADA2` algorithum to identify and correct errors) these data are processed in 'seuqencing run' subsets until after denoising because: (1) the AMI dataset spans multiple seuqencing runs, and (2) different seuqencing runs have unique error profiles.
