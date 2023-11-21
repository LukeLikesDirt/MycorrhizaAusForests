#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --partition=day
#SBATCH --mem=100G
#SBATCH --output=../01.Bioinformatics_ITS/slurm/%x.%j.out

## Script: Cluster OTUs.
## Author: Luke Florence.
## Date: 5th November 2023.

## Constants and file paths
readonly THREADS=16                                                     ## Set the number of threads
readonly IDENTITY=0.97                                                  ## Set the identity threshold for clustering
CHIMERA_METHODS_DIR=../../data/AusMicrobiome/ITS/07.Chimera_filtered    ## Path to chimera filtered fasta file
CLUSTERED_METHODS_DIR="../../data/AusMicrobiome/ITS/08.OTUs"            ## Path to clustered fasta file and OTU table

# Define extensions for different methods used for denoising
readonly method_UNOISE3="_UNOISE3"
readonly method_DADA2="_DADA2"

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## Function for clustering OTUs and formatting the OTU table
cluster_OTUs() {
    local method=$1
    local CHIMERA_FILTERED_DIR="$CHIMERA_METHODS_DIR$method"
    local CLUSTERED_DIR="$CLUSTERED_METHODS_DIR$method"

    # Create the directory if it doesn't exist
    mkdir -p "$CLUSTERED_DIR"

    log "Clustering at 97% for $method at:"

    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta" \
        --threads "$THREADS" \
        --id "$IDENTITY" \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --relabel_sha \
        --uc "$CHIMERA_FILTERED_DIR/OTUs.uc" \
        --centroids "$CLUSTERED_DIR/OTUs.fasta" \
        --biomout "$CLUSTERED_DIR/OTUs.biom" \
        --otutabout "$CLUSTERED_DIR/OTUs.txt"

    # Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "$CLUSTERED_DIR/OTUs.txt"

    printf '\nNumber of unique sequences and OTUs\n'
    printf '    Unique non-chimeric sequence: %s\n' "$(grep -c "^>" "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta")"
    printf '    Clustered OTUs: %s\n' "$(grep -c "^>" "$CLUSTERED_DIR/OTUs.fasta")"
}

# Function to build sequence matching table for post-clustering curation with LULU
match_list() {
    local method=$1
    local CHIMERA_FILTERED_DIR="$CHIMERA_METHODS_DIR$method"
    local CLUSTERED_DIR="$CLUSTERED_METHODS_DIR$method"

    log "Building sequence matching table for $method at:"

    vsearch \
        --usearch_global "$CLUSTERED_DIR/OTUs.fasta" \
        --db "$CLUSTERED_DIR/OTUs.fasta" \
        --self --id .84 \
        --iddef 1 --userout "$CLUSTERED_DIR/match_list.txt" \
        -userfields query+target+id \
        --maxaccepts 0 --query_cov 0.9 --maxhits 10
}

###############################################################################
## Main script ################################################################
###############################################################################

log 'Starting at:'

source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

# Run cluster_OTUs for both methods concurrently
(cluster_OTUs "$method_DADA2") &  # Run the first method in the background
(cluster_OTUs "$method_UNOISE3") &  # Run the second method in the background
wait  # Wait for background processes to finish

# Run match_list for both methods concurrently
(match_list "$method_DADA2") &  # Run the first method in the background
(match_list "$method_UNOISE3") &  # Run the second method in the background
wait  # Wait for background processes to finish


conda deactivate

log 'Finished at:'
