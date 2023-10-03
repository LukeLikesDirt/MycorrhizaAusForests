#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=720:00:00
#SBATCH --partition=month
#SBATCH --mem-per-cpu=32G
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/%x.%j.out

# Script: code/01.Bioinformatics_ITS/02.Extract_ITS.sh
# Purpose: Prepare Illumina forwards reads targeting the ITS1 for denoising
# with DADA2 by extracting the ITS1 subregion using ITSxpress and quality
# truncating reads using Trimmomatic.
# Author: Luke Florence
# Date: 28th September 2023

# Script Overview:
# ---------------
# The sequencing data processed here span 44 runs. These data are subset by
# Sequencing prior to denoising with DADA2 to improve error rate estimations,
# as each sequencing run has a unique error profile.
# This script performs the following tasks:
#   (1) Extract the ITS1 subregion with ITSxpress
#   (2) Quality truncate reads with Trimmomatic
#   (3) Quality check trimmed reads with fastQC
#   (4) Calculate statistics for tracking read survival across the pipeline

# Constant variables
readonly THREADS=8          # The number of threads to use for parallel processing
readonly first_run=37       # The first run to process
readonly last_run=40        # The last run to process
readonly WINDOW=4           # Trimmomatic: The sliding window size for averaging quality scores
readonly QUAL=13            # Trimmomatic: The quality threshold for sliding window trimming
readonly MINLEN=80          # Trimmomatic: The minimum read length after trimming

# Directory paths and file names
readonly path="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"    # The path to the project directory
readonly raw_data="$path/data/AusMicrobiome/ITS/01.Raw_data"               # Path to demultiplexed reads
readonly ITS_extracted="$path/data/AusMicrobiome/ITS/02.ITS_extracted"     # Path for the extracted ITS reads
readonly quality_trimmed="$path/data/AusMicrobiome/ITS/03.Quality_trimmed" # Path for the quality truncated reads
readonly track_reads="$path/data/AusMicrobiome/ITS/trim_summary.txt"       # Path to the file for tracking reads across the pipline

# Create subdirectories if they don't exist
mkdir -p "$path/02.ITS_extracted"
mkdir -p "$path/03.Quality_trimmed"
for dir in "$ITS_extracted" "$quality_trimmed"; do
  for run in $(seq "$first_run" "$last_run"); do
    mkdir -p "$dir/run$run"
  done
done

# Initialise read tracking file if it doesn't exist
if [ ! -f "$track_reads" ]; then
    printf "Run\tInput_Reads\tOutput_Reads\tRetained_Reads_%%\tDropped_Reads_%%\tMean_Reads_Per_Sample\tMin_Reads\tMax_Reads\tSamples_LT_10K_Reads\tSamples_10-20K_Reads\tSamples_MT_20K_Reads\n" > "$track_reads"
fi

# Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n[$timestamp] %s\n" "$1"
}

# Function to extract ITS1 with ITSxpress
extract_its1() {
    local run_number="$1"
    log "Run $run_number: Extracting ITS1 with ITSxpress"
    cd "$raw_data/run$run_number"

    for f in *_R1.fastq.gz; do
        itsxpress \
            --single_end \
            --fastq "${f}" \
            --cluster_id 1.0 \
            --region ITS1 \
            --taxa All \
            --log "$ITS_extracted/run${run_number}/logfile.txt" \
            --outfile "$ITS_extracted/run${run_number}/${f}" \
            --threads "$THREADS"
    done
}

# Function to quality trim reads with Trimmomatic
quality_trim() {
    local run_number="$1"
    log "Run $run_number: Quality trimming with Trimmomatic"
    cd "$ITS_extracted/run$run_number"

    for f in *.fastq.gz; do
        trimmomatic SE -threads "$THREADS" "$f" "$quality_trimmed/run$run_number/$f" \
            SLIDINGWINDOW:"$WINDOW":"$QUAL" MINLEN:"$MINLEN" 
    done
}

# Function to generate read quality report
generate_quality_report() {
    local run_number="$1"
    log "Run $run_number: Generating quality report"
    fastqc "$quality_trimmed/run$run_number"/*R1.fastq.gz -o "$quality_trimmed/run$run_number"
    multiqc "$quality_trimmed/run$run_number" -o "$quality_trimmed/run$run_number"
    
    # Remove intermediate files and 'multiqc_data' directory
    rm "$quality_trimmed/run$run_number"/*fastqc.zip "$quality_trimmed/run$run_number"/*fastqc.html
    rm "$quality_trimmed/run$run_number"/multiqc_data/*
    rmdir "$quality_trimmed/run$run_number"/multiqc_data
}

# Function to calculate statistics for tracking reads
calculate_statistics() {
    local run_number="$1"
    log "Run $run_number: Calculating trimming statistics"

    input_dir="$raw_data/run$run_number"
    output_dir="$quality_trimmed/run$run_number"

    # Calculate the number of input, output (retained), and dropped reads during trimming
    input_reads=$(($(zcat "$input_dir/s"*.fastq.gz| wc -l) / 4))
    output_reads=$(($(zcat "$output_dir/s"*.fastq.gz| wc -l) / 4))
    dropped_reads=$((input_reads - output_reads))

    # Calculate the percentage of reads retained and dropped during trimming
    percentage_retained=$(awk "BEGIN { printf \"%.2f\", ($output_reads / $input_reads) * 100 }")
    percentage_dropped=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads / $input_reads) * 100 }")
    
    # Calculate mean reads per sample per run after trimming
    num_samples=$(find "$output_dir/" -name "s*.fastq.gz" | wc -l)
    mean_reads_per_sample=$((output_reads / num_samples))

    # Find the minimum and maximum number of reads within a sample
    min_reads_in_sample=999999999
    max_reads_in_sample=0

    for sample_file in "$output_dir/s"*.fastq.gz; do
        total_reads=$(( $(zcat "$sample_file" | wc -l) / 4 ))
    
        if [ "$total_reads" -lt "$min_reads_in_sample" ]; then
            min_reads_in_sample="$total_reads"
        fi
    
        if [ "$total_reads" -gt "$max_reads_in_sample" ]; then
            max_reads_in_sample="$total_reads"
        fi
    done

    # Calculate the number of samples with <10,000 reads, 10,000-20,000 reads, and >20,000 reads within the run
    samples_less_than_10k=0
    samples_10K_to_20K=0
    samples_more_than_20k=0

    for sample_file in "$output_dir/s"*.fastq.gz; do
        total_reads=$(($(zcat "$sample_file" | wc -l) / 4))

        if [ "$total_reads" -lt 10000 ]; then
            ((samples_less_than_10k++))
        elif [ "$total_reads" -ge 10000 ] && [ "$total_reads" -lt 20000 ]; then
            ((samples_10K_to_20K++))
        else
            ((samples_more_than_20k++))
        fi
    done

    # Append the results to the result file
    printf "$run_number\t$input_reads\t$output_reads\t$percentage_retained\t$percentage_dropped\t$mean_reads_per_sample\t$min_reads_in_sample\t$max_reads_in_sample\t$samples_less_than_10k\t$samples_10K_to_20K\t$samples_more_than_20k\n" >> "$track_reads"
}

log "Starting at: $(date)"

# Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

# Execute the functions
for run_number in $(seq "$first_run" "$last_run"); do
    extract_its1 "$run_number"
    quality_trim "$run_number"
    generate_quality_report "$run_number"
    calculate_statistics "$run_number"
done

# Deactivate conda environment
conda deactivate

log "Finished at: $(date)"
