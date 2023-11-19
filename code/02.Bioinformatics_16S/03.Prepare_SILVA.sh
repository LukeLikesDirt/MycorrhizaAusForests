#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --output=../02.Bioinformatics_16S/slurm/%x.%j.out

# Script: Prepare the SILVA references dataset for chimera detection and
#         taxonomic for V1-V3 amplicons using the 27F 519R primer pair.
# Author: Luke Florence
# Date:   3rd October 2023

# Constants and subdirectroies
readonly THREADS=8                          # Number of threads
readonly PRIMER_F="AGAGTTTGATCMTGGCTCAG"    # Forwad primers for the V1-V3 region
readonly PRIMER_R="GWATTACCGCGGCKGCTG"      # Reverse primers for the V1-V3 region
readonly SILVA_dir="../../data/AusMicrobiome/16S/05.Reference_dataset"  # Full 16S reference dataset subdirectory
readonly V1V3_dir="$SILVA_dir/V1V3/"                                    # V1-V3 reference dataset subdirectory
readonly refseqs="$V1V3_dir/SILVA_V1V3.fasta"                           # V1-V3 FASTA file
mkdir -p "$SILVA_dir" "$V1V3_dir"

# URL to the latest release (v138; released 27.08.2020) of the SILVA SSU dataset
# with a 99% criterion applied to remove redundant sequences
readonly SILVA_138_NR99="https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"

# Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

# Function to construct reverse-complement sequences, sourced from:
# https://github.com/Mycology-Microbiology-Center/GSMc/blob/main/02.Extract_ITS.sh
RC () {
  echo "$1" | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBb]" | rev
}
# Define the reverse complement of the reverse primers
readonly PRIMER_Rr=$( RC "$PRIMER_R")

#### Main script ###############################################################

log "Starting at:"

# Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

# Download SILVA database
log "Downloading SILVA database at"
wget -P "$SILVA_dir" "$SILVA_138_NR99"

# Extract the V1-V3 region from the SILVA database using the 27F and 519R primer pair (the primers used in the AusMicrobiome project)
# The typical length of the V1-V3 is 490bp
log "Extracting V1-V3 region from SILVA database"

cutadapt \
     -a "$PRIMER_F;min_overlap=19"..."$PRIMER_Rr;min_overlap=17" \
     -e 0.1 \
     --minimum-length 390 \
     --maximum-length 590 \
     --discard-untrimmed \
     --cores 2 \
     -o "$V1V3_dir/SILVA_V1V3.fasta.gz" "$SILVA_dir/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz" \
     >> "$SILVA_dir/log.txt"

# Decompress the V1-V3 FASTA file
gzip -d "$V1V3_dir/SILVA_V1V3.fasta.gz"

# Reformat the V1-V3 FASTA file:
#   (1) Replace the first occurrence of a spaces with a semicolon
#   (2) Eukaryota annotations are long and bizarre. Because I'm not interested
#       in Eukaryotes, I'll change all "Eukaryota" annotations to uninformative
#       annotatins that fit my taxonomy hierarchy after BLASTing the sequences
#   (3) Species names contain spaces so I'll replace spaces with underscores
log "Reformatting V1-V3 FASTA file at"
cat "$V1V3_dir/SILVA_V1V3.fasta" |
  sed -e 's/ /;/' -e 's/Eukaryota.*$/k__ Eukaryota;p__ Eukaryota;o__ Eukaryota;c__ Eukaryota;f__ Eukaryota;g__ Eukaryota;s__ Eukaryota/' |
  sed 's/ /_/g' > "$V1V3_dir/SILVA_V1V3_reformatted.fasta"

# Rename the reformatted file to the original name
mv "$V1V3_dir/SILVA_V1V3_reformatted.fasta" "$V1V3_dir/SILVA_V1V3.fasta"

if [ ! -f "$refseqs" ]; then
    log "Error: Input FASTA file not found."
    exit 1
fi

## Build the SILVA V1-V3 reference database
makeblastdb \
    -in "$V1V3_dir/SILVA_V1V3.fasta" \
    -out "$V1V3_dir/SILVA_V1V3" \
    -dbtype 'nucl' \
    -hash_index

# Extract sequence lengths and calculate stats
awk '/^>/ {if (seq!="") {print length(seq); seq="";} next} {seq = seq $0} END {if (seq!="") print length(seq)}' "$refseqs" > "$V1V3_dir/sequence_lengths.txt"

min_length=$(sort -n "$V1V3_dir/sequence_lengths.txt" | head -n 1)
max_length=$(sort -n "$V1V3_dir/sequence_lengths.txt" | tail -n 1)
total_length=$(awk '{sum += $1} END {print sum}' "$V1V3_dir/sequence_lengths.txt")
sequence_count=$(wc -l < "$V1V3_dir/sequence_lengths.txt")
sequences_greater_than_520=$(awk '$1 > 520 {count++} END {print count}' "$V1V3_dir/sequence_lengths.txt")
sequences_less_than_520=$(awk '$1 < 520 {count++} END {print count}' "$V1V3_dir/sequence_lengths.txt")

# Calculate the average length
average_length=$(awk -v total_length="$total_length" -v sequence_count="$sequence_count" 'BEGIN { printf "%.2f\n", total_length / sequence_count }')

# Calculate the proportion of sequences with lengths less than and greater than 520bp
proportion_gt_520=$(awk -v sequences_greater_than_520="$sequences_greater_than_520" -v sequence_count="$sequence_count" 'BEGIN { printf "%.2f\n", sequences_greater_than_520 / sequence_count }')
proportion_lt_520=$(awk -v sequences_less_than_520="$sequences_less_than_520" -v sequence_count="$sequence_count" 'BEGIN { printf "%.2f\n", sequences_less_than_520 / sequence_count }')

# Log results
printf "Minimum Length: %s\n" "$min_length"
printf "Maximum Length: %s\n" "$max_length"
printf "Average Length: %s\n" "$average_length"
printf "Proportion of reads > 520bp: %s\n" "$proportion_gt_520"
printf "Proportion of reads < 520bp: %s\n" "$proportion_lt_520"

# Clean up temporary files
rm "$V1V3_dir/sequence_lengths.txt"

# Deactivate conda environment
conda deactivate

log "Finished at:"
