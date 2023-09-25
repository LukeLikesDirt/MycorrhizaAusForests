#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1440:00:00
#SBATCH --partition=long
#SBATCH --mem-per-cpu=32G

echo "Starting at: $(date)"

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda
conda activate shell

## Take advantage of all the available threads
OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
echo
echo "Number of available threads: $OMP_NUM_THREADS"

## Organize subdirectories
path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/16S    ## Path to the main 16S data directory
mkdir $data/reference_datasets

## Download SILVA database
wget -P $data/reference_datasets https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

## Function to construct reverse-complement sequence from https://github.com/Mycology-Microbiology-Center/GSMc/blob/main/02.Extract_ITS.sh
RC () {
  echo "$1" | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBb]" | rev
}
export -f RC

export PRIMER_F="AGAGTTTGATCMTGGCTCAG"
export PRIMER_R="GWATTACCGCGGCKGCTG"
export PRIMER_Fr=$( RC "$PRIMER_F")
export PRIMER_Rr=$( RC "$PRIMER_R")
export OVERLAP_F=$(( ${#PRIMER_F} * 2 / 3 ))
export OVERLAP_R=$(( ${#PRIMER_R} * 2 / 3 ))
export OVERLAP=$(( ($OVERLAP_F + $OVERLAP_R) / 2 ))

## Extract the V1-V3 region from the SILVA database using the 27F and 519R primer pair (the primers using in the AusMicrobiome project)

cutadapt \
     -a $PRIMER_F...$PRIMER_Rr \
     -e 0.2 \
     --overlap $OVERLAP_R \
     --minimum-length 300 \
     --maximum-length 600 \
     --discard-untrimmed \
     --cores 2 \
     -o $path/reference_datasets/SILVA_V1V3.fasta.gz $path/reference_datasets/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz \
     >> $path/reference_datasets/log.txt


gzip -d $path/reference_datasets/SILVA_V1V3.fasta.gz

fasta=$path/reference_datasets/SILVA_V1V3.fasta

if [ ! -f "$fasta" ]; then
    echo "Error: Input FASTA file not found."
    exit 1
fi

# Extract sequence lengths and calculate stats
awk '/^>/ {if (seq!="") {print length(seq); seq="";} next} {seq = seq $0} END {if (seq!="") print length(seq)}' "$fasta" > sequence_lengths.txt

min_length=$(sort -n sequence_lengths.txt | head -n 1)
max_length=$(sort -n sequence_lengths.txt | tail -n 1)
total_length=$(awk '{sum += $fasta} END {print sum}' sequence_lengths.txt)
sequence_count=$(wc -l sequence_lengths.txt | awk '{print $fasta}')
sequences_greater_than_500=$(awk '$fasta > 520 {count++} END {print count}' sequence_lengths.txt)

# Calculate average length
average_length=$(awk -v total_length="$total_length" -v sequence_count="$sequence_count" 'BEGIN { printf "%.2f\n", total_length / sequence_count }')

# Calculate the proportion of sequences with lengths greater than 500
proportion_gt_500=$(awk -v sequences_greater_than_500="$sequences_greater_than_500" -v sequence_count="$sequence_count" 'BEGIN { printf "%.2f\n", sequences_greater_than_500 / sequence_count }')

# Output results
echo "Minimum Length: $min_length"
echo "Maximum Length: $max_length"
echo "Average Length: $average_length"
echo "Proportion > 500: $proportion_gt_500"

# Clean up temporary file
rm sequence_lengths.txt



