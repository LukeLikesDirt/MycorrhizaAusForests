#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --mem-per-cpu=32G
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/%x.%j.out

## This chimera filtering script is adapted from: https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
##
## Pre-requisites:
## ---------------
## At this stage, there should be one '.fasta' file formatted for VSEARCH, 
## containing samples from all sequencing runs. Reads should be ITS extracted,
## quality-trimmed, quality-filtered, denoised and dereplicated within samples.
##
## Script Purpose:
## ---------------
## Here, I will dereplicate across samples and perform both de-novo and
## reference-based chimera filtering in VSEARCH. I have chosen to perform
## chimera filtering using VSEARCH rather than DADA2 because:
##  (1) The DADA2 de novo chimera detection has high rates of false positives
##      compared to the UCHIME3 de novo algorithm that are implemented in
##      VSEARCH (Edgar 2016, bioRxiv 074252; https://doi.org/10.1101/074252).
##  (2) DADA2 does not implement reference-based chimera detection. Although
##      de novo chimera detection outperforms reference-based chimera detection,
##      in terms of true positive chimera detection, referenced-based chimera
##      detection is recommended as most referenced-based chimeras are true
##      chimeras (Tedersoo et al. 2022, Molecular ecology 31, no. 10 (2022): 2769-2795).

## Constant variables and file paths
readonly THREADS=8                                                                      ## Set the number of threads
readonly path="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"            ## The path to the project directory
readonly map="$path/data/AusMicrobiome/ITS/map.pl"                                      ## Map file for fasta reconstruction
readonly refseqs="$path/data/AusMicrobiome/ITS/06.Reference_dataset/ITS1/ITS1.fasta"    ## Path to UNITE reference dataset
readonly denoised_dir="$path/data/AusMicrobiome/ITS/05.Denoised"                        ## Path to denoised fasta file
readonly chimeraFiltered_dir="$path/data/AusMicrobiome/ITS/07.Chimera_filtered"         ## Path for chimera filtered fasta file
mkdir -p "$chimeraFiltered_dir"                                                         ## Chimera filtered directory

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## Function for dereplication across samples, de novo and referenced-based chimera filtering, and mapping reads back to dereplicated non-chimeric sequences
chimera_filter() {

    log 'Dereplicating across samples at:'

    vsearch \
        --derep_fulllength "$denoised_dir/all.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$chimeraFiltered_dir/all.derep.uc" \
        --output "$chimeraFiltered_dir/all.derep.fasta"

    log 'De novo chimera detection at:'

    vsearch \
        --uchime3_denovo "$chimeraFiltered_dir/all.derep.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$chimeraFiltered_dir/all.denovo.nonchimeras.fasta"

    log 'Reference-based chimera detection at:'

    vsearch \
        --uchime_ref "$chimeraFiltered_dir/all.denovo.nonchimeras.fasta" \
        --threads "$THREADS" \
        --db "$refseqs" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$chimeraFiltered_dir/all.ref.nonchimeras.fasta"

    ## Map reads back to dereplicated non-chimeric sequences
    log 'Extract all dereplicated non-chimeric sequences at:'

    perl "$map" "$chimeraFiltered_dir/all.derep.fasta" "$chimeraFiltered_dir/all.derep.uc" "$chimeraFiltered_dir/all.ref.nonchimeras.fasta" > "$chimeraFiltered_dir/all.nonchimeras.derep.fasta"

    ## Track reads across the pipeline
    printf 'Track reads across the pipeline:\n'
    printf 'Unique sequence input: %s\n' "$(grep -c "^>" "$chimeraFiltered_dir/all.derep.fasta")"
    printf 'Unique sequences after de novo chimera detection: %s\n' "$(grep -c "^>" "$chimeraFiltered_dir/all.denovo.nonchimeras.fasta")"
    printf 'Unique sequences after reference-based chimera detection: %s\n' "$(grep -c "^>" "$chimeraFiltered_dir/all.ref.nonchimeras.fasta")"
    printf 'Unique non-chimeric sequence output: %s\n' "$(grep -c "^>" "$chimeraFiltered_dir/all.nonchimeras.derep.fasta")"

}

## Execute the chimera filter function

log 'Starting at:'

## Activate the conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell
## Dereplicate and remove chimeras
chimera_filter
## Deactivate the conda environment
conda deactivate

log 'Finishing at:'
