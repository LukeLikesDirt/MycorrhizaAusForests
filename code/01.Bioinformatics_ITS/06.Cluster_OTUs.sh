#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=day
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/%x.%j.out

## Script: Cluster OTUs.
## Purpose: To cluster OTUs prior to the selection of OTU representative sequences.
## Credit: This script is adapted from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline.
## Author: Luke Florence.
## Date: 5th October 2023.

## Constants and file paths
readonly THREADS=8                                                                       ## Set the number of threads
readonly IDENTITY=0.97                                                                   ## Set the identity threshold for clustering
readonly PROJECT_PATH="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"     ## Path to the project directory
readonly CHIMERA_FILTERED_DIR="$PROJECT_PATH/data/AusMicrobiome/ITS/07.Chimera_filtered" ## Path to chimera filtered fasta file
readonly CLUSTERED_DIR="$PROJECT_PATH/data/AusMicrobiome/ITS/08.Clustered"               ## Path to clustered fasta file and OTU table
mkdir -p "$CLUSTERED_DIR"                                                                ## Make the clustered subdirectory

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## Function for clustering OTUs and formatting the OTU table
cluster_OTUs() {

    ## Cluster at 97% and generate OTU table
    log 'Clustering at 97% at:'

    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta" \
        --threads "$THREADS" \
        --id "$IDENTITY" \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.clustered.uc" \
        --relabel_sha \
        --centroids "$CLUSTERED_DIR/OTUs.fasta" \
        --otutabout "$CLUSTERED_DIR/OTUs.txt"

    ## Format OTU table
    ## Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "$CLUSTERED_DIR/OTUs.txt"
    ## Convert to .csv and save to the output directory
    sed -e 's/\s\+/,/g' "$CLUSTERED_DIR/OTUs.txt" > "$OUTPUT_DIR/OTUs.csv"

    printf '\nNumber of unique sequences and OTUs\n'
    printf '    Unique non-chimeric sequence: %s\n' "$(grep -c "^>" "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta")"
    printf '    Clustered OTUs: %s\n' "$(grep -c "^>" "$CLUSTERED_DIR/OTUs.fasta")"

}

###############################################################################
## Main script ################################################################
###############################################################################

log 'Starting at:'

source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

cluster_OTUs

conda deactivate

log 'Finished at:'
