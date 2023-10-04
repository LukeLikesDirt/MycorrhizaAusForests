#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_16S/slurm/%x.%j.out

# Script: code/02.Bioinformatics_16S/04.Prepare_SILVA.sh
# Purpose: Prepare the SILVA references dataset for chimera detection and taxonomic for V1-V3 amplicons using the 27F 519R primer pair.
# Author: Luke Florence
# Date: 3rd October 2023

# Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

# Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "[$timestamp] %s\n" "$1"
}

# Function to construct reverse-complement sequences: https://github.com/Mycology-Microbiology-Center/GSMc/blob/main/02.Extract_ITS.sh
RC () {
  echo "$1" | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBb]" | rev
}
export -f RC

# Constant variables
readonly THREADS=8                          # Number of threads
readonly PRIMER_F="AGAGTTTGATCMTGGCTCAG"    # Forwad primers for the V1-V3 region
readonly PRIMER_R="GWATTACCGCGGCKGCTG"      # Reverse primers for the V1-V3 region
readonly PRIMER_Rr=$( RC "$PRIMER_R")       # Reverse complement of the reverse primers
export OVERLAP=16                           # Minimum overlap for the forward and reverse primers

# File paths and names: 
readonly path="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"    # Main project directory
readonly SILVA_dir="$path/data/AusMicrobiome/16S/06.Reference_datasets"         # Full 16S reference dataset subdirectory
readonly V1V3_dir="$SILVA_dir/V1V3/"                                            # V1-V3 reference dataset subdirectory
readonly refseqs="$V1V3_dir/SILVA_V1V3.fasta"                                   # V1-V3 FASTA file

# URL to the latest release (v138; released 27.08.2020) of the SILVA SSU dataset with a 99% criterion applied to remove redundant sequences
readonly SILVA_138_NR99="https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"

# Make directories
mkdir -p "$SILVA_dir"
mkdir -p "$V1V3_dir"

log "Starting at: $(date)"

# Download SILVA database
log "Downloading SILVA database"
wget -P "$SILVA_dir" "$SILVA_138_NR99"

# Extract the V1-V3 region from the SILVA database using the 27F and 519R primer pair (the primers used in the AusMicrobiome project)
# The typical length of the V1-V3 is 490bp, but I allow for a range of 440-540bp
log "Extracting V1-V3 region from SILVA database"
cutadapt \
     -a "$PRIMER_F"..."$PRIMER_Rr" \
     -e 0.1 \
     -O "$OVERLAP" \
     --minimum-length 440 \
     --maximum-length 540 \
     --discard-untrimmed \
     --cores 2 \
     -o "$V1V3_dir/SILVA_V1V3.fasta.gz" "$SILVA_dir/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz" \
     >> "$SILVA_dir/log.txt"

# Decompress the V1-V3 FASTA file
gzip -d "$V1V3_dir/SILVA_V1V3.fasta.gz"

if [ ! -f "$refseqs" ]; then
    log "Error: Input FASTA file not found."
    exit 1
fi

# Extract sequence lengths and calculate stats
awk '/^>/ {if (seq!="") {print length(seq); seq="";} next} {seq = seq $0} END {if (seq!="") print length(seq)}' "$refseqs" > "$V1V3_dir/sequence_lengths.txt"

min_length=$(sort -n "$V1V3_dir/sequence_lengths.txt" | head -n 1)
max_length=$(sort -n "$V1V3_dir/sequence_lengths.txt" | tail -n 1)
total_length=$(awk '{sum += $efseqs} END {print sum}' "$V1V3_dir/sequence_lengths.txt")
sequence_count=$(wc -l "$V1V3_dir/sequence_lengths.txt" | awk '{print $refseqs}')
sequences_greater_than_520=$(awk '$refseqs > 520 {count++} END {print count}' "$V1V3_dir/sequence_lengths.txt")

# Calculate the average length
average_length=$(awk -v total_length="$total_length" -v sequence_count="$sequence_count" 'BEGIN { printf "%.2f\n", total_length / sequence_count }')

# Calculate the proportion of sequences with lengths less than and greater than 520bp
proportion_gt_520=$(awk -v sequences_greater_than_520="$sequences_greater_than_520" -v sequence_count="$sequence_count" 'BEGIN { printf "%.2f\n", sequences_greater_than_520 / sequence_count }')
proportion_lt_520=$(awk -v sequences_less_than_520="$sequences_less_than_520" -v sequence_count="$sequence_count" 'BEGIN { printf "%.2f\n", sequences_less_than_520 / sequence_count }')

# Log results
log "Minimum Length: $min_length"
log "Maximum Length: $max_length"
log "Average Length: $average_length"
log "Proportion of reads > 520bp: $proportion_gt_520"
log "Proportion of reads < 520bp: $proportion_gt_520"

# Clean up temporary file
rm "$V1V3_dir/sequence_lengths.txt"

log "Finished at: $(date)"
