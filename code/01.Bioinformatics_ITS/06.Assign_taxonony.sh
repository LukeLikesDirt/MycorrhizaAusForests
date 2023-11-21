#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=long
#SBATCH --output="../01.Bioinformatics_ITS/slurm/%x.%j.out"

## Script: Taxonomic assignment with BLASTn
## Author: Luke Florence
## Date: 5th November 2023
## Software: BLAST v2.14.1 - https://blast.ncbi.nlm.nih.gov/Blast.cgi

# Constants and subdirectories
readonly THREADS=8
readonly REFERENCE_SEQUENCES="../../data/AusMicrobiome/ITS/06.Reference_dataset/ITS1/ITS1"
CLUSTERED_METHODS_DIR="../../data/AusMicrobiome/ITS/08.OTUs"
TAXA_METHODS_DIR="../../data/AusMicrobiome/ITS/09.Taxonomy"

# Define extensions for different methods used for denoising
readonly method_DADA2="_DADA2"
readonly method_UNOISE3="_UNOISE3"

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## BLAST function
blast_best_hit() {
    local method=$1
    local OTU_FASTA="$CLUSTERED_METHODS_DIR$method/OTUs.fasta"
    local TAXA_DIR="$TAXA_METHODS_DIR$method"

    # Create the directory if it doesn't exist
    mkdir -p "$TAXA_DIR"

    log 'BLAST best hit starting at'

    # Blast best hit
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES" \
        -max_target_seqs 1 \
        -max_hsps 1 \
        -out "$TAXA_DIR/BLAST_best_hit.txt" \
        -num_threads "$THREADS"

    # Reformat the taxa table
    sed -i '1s/^/OTU_ID;abundance\treference;kingdom;phylum;class;order;family;genus;species\tpident\tlength\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n/' "$TAXA_DIR/BLAST_best_hit.txt"
    sed 's/[[:space:]]\{1,\}/;/g' "$TAXA_DIR/BLAST_best_hit.txt" > "$TAXA_DIR/BLAST_best_hit.csv"
    rm "$TAXA_DIR/BLAST_best_hit.txt"

}

blast_best_ten_hits() {
    local method=$1
    local OTU_FASTA="$TAXA_METHODS_DIR$method/OTUs.fasta"
    local TAXA_DIR="$TAXA_METHODS_DIR$method"

    # Create the directory if it doesn't exist
    mkdir -p "$TAXA_DIR"

    # blast ten hits
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -strand both \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -out "$TAXA_DIR/BLAST_best_10.txt" \
        -num_threads "$THREADS"

    
    sed -i '1s/^/OTU_ID;abundance\treference;kingdom;phylum;class;order;family;genus;species\tpident\tlength\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n/' "$TAXA_DIR/BLAST_best_10.txt"
    sed 's/[[:space:]]\{1,\}/;/g' "$TAXA_DIR/BLAST_best_10.txt" > "$TAXA_DIR/BLAST_best_10.csv"
    rm "$TAXA_DIR/BLAST_best_10.txt"
        
}

###############################################################################
## Main script ################################################################
###############################################################################

log 'Starting at:'

source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

(blast_best_hit "$method_DADA2") & 
(blast_best_ten_hits "$method_DADA2") &
(blast_best_hit "$method_UNOISE3") &
(blast_best_ten_hits "$method_UNOISE3") &

# Wait for all background processes to finish
wait

conda deactivate

log 'Finished at:'
