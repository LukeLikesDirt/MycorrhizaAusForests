#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00
#SBATCH --partition=week
#SBATCH --mem-per-cpu=32G
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/02.Bioinformatics_16S/slurm/%x.%j.out

# Script: code/02.Bioinformatics_16S/02.Quality_trim.sh
# Purpose: Prepare Illumina paired-end 16S amplicons (2 x 300bp) for denoising with DADA2 by quality truncating reads using Trimmomatic.
# Author: Luke Florence
# Date: 3rd October 2023

# Script Purpose:
# ----------------
# The Australian Microbiome 16S amplicon covers the V1-V3 region, which is 
# typically around 490 bases long (Allen, et al. 2016, BMC Res Notes 9, 380).
# In the current SILVA reference dataset, approximately 97% of V1-V3 amplicons
# are shorter than 520 bases. Therefore, the goal of this script is to retain
# good sequencing depth (>10,000 reads per sample) using a 532 base threshold,
# which equals 520 base length amplicons when accounting for the 12 base
# overlap required for merging with DADA2. Because the primers were not well
# detected in the raw reads, primers are trimmed using a fixed length in this 
# script.

# Script Overview:
# ---------------
# The '02.Quality_trim.sh' script performs the following tasks:
#   - Removes 20 bases from the 5' end of the forward reads to remove the
#     forward primers (27F).
#   - Removes 18 bases from the 5' end of the reverse reads to remove the
#     reverse primers (519R).
#   - Hard truncate reverse reads to 252 to require a maximum potential amplicon
#     length of 520 bases after merging with DADA2: 
#     280bp (R1) + 252bp (R2) - 12bp (merging overlap) = 520bp amplicons.
#   - Quality truncates using a threshold of Q10, which equals a probable error
#     rate of 1 in 10 calls. Considering that DADA2 will handle many incorrect
#     calls when denoising and because reads will be merged, this 'relaxed' 
#     threshold seems reasonable. Tighten it to at least Q13 if using forward
#     reads only.
#   - Removes reads shorter than 280 bases for forward reads and 252 bases for
#     reverse reads after quality truncation.

# Activate the conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

# Set constant variables:
readonly NUM_RUNS=43        # The number of sequencing runs to be processed
readonly THREADS=8          # The number of threads to use for parallel processing
readonly WINDOW=4           # Trimmomatic: The sliding window size for averaging quality scores
readonly QUAL=10            # Trimmomatic: The quality threshold for sliding window trimming
readonly HEADCROP_FWD=20    # Trimmomatic: The number of bases to remove from the start of the forward read
readonly MINLEN_FWD=280     # Trimmomatic: The minimum length read length to retain forward reads
readonly CROP_REV=252       # Trimmomatic: Fixed truncation length of reverse reads after HEADCROP and QUAL trimming (removes bases from the distal, low-quality end of the read)
readonly HEADCROP_REV=18    # Trimmomatic: The number of bases to remove from the start of the reverse read
readonly MINLEN_REV=252     # Trimmomatic: The minimum length read length to retain reverse reads

# Directory paths and file names:
path="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests" # The path to the project directory0
raw_data="$path/data/AusMicrobiome/16S/01.Raw_data"                 # The path to demultiplexed reads
trimmed="$path/02.Quality_trimmed"                                  # The path for the trimmed reads
results_file_fwd="$trimmed/trim_summary_fwd.txt"                    # The path to the forward read quality trimming summary file
results_file_rev="$trimmed/trim_summary_rev.txt"                    # The path to the reverse read quality trimming summary file

# Logging function
log() {
    local timestamp
    timestamp=$(date)
    printf "[$timestamp] %s\n" "$1"
}

# Function to run Trimmomatic and generate quality reports for each sample using fastQC and multiQC
run_trimmomatic_and_report() {
    local run_number="$1"

    log "Run $run_number: Quality truncate reads"

    cd "$raw_data/run$run_number"

    for f in *R1.fastq.gz; do
        trimmomatic SE "$f" "$trimmed/run$run_number/fwd/$f" HEADCROP:"$HEADCROP_FWD" SLIDINGWINDOW:4:"$QUAL" MINLEN:"$MINLEN_FWD" -threads "$THREADS"
    done

    for f in *R2.fastq.gz; do
        trimmomatic SE "$f" "$trimmed/run$run_number/rev/$f" CROP:"$CROP_REV" HEADCROP:"$HEADCROP_REV" SLIDINGWINDOW:4:"$QUAL" MINLEN:"$MINLEN_REV" -threads "$THREADS"
    done

    log "Run $run_number: Generating quality report"
    fastqc "$trimmed/run$run_number/fwd"/*R1.fastq.gz -o "$trimmed/run$run_number/fwd"
    fastqc "$trimmed/run$run_number/rev"/*R2.fastq.gz -o "$trimmed/run$run_number/rev"
    multiqc "$trimmed/run$run_number/fwd" -o "$trimmed/run$run_number/fwd"
    multiqc "$trimmed/run$run_number/rev" -o "$trimmed/run$run_number/rev"
    
    # Remove intermediate files and 'multiqc_data' directories (add error handling)
    rm -f "$trimmed/run$run_number"/fwd/*fastqc.zip "$trimmed/run$run_number"/fwd/*fastqc.html
    rm -f "$trimmed/run$run_number"/rev/*fastqc.zip "$trimmed/run$run_number"/rev/*fastqc.html
    rm -rf "$trimmed/run$run_number"/fwd/multiqc_data
    rm -rf "$trimmed/run$run_number"/rev/multiqc_data
}

# Function to calculate read tracking statistics for each sequencing run
calculate_and_log_statistics() {
    local run_number="$1"

    log "Run $run_number: Generating trim statistics report"

    # Calculate statistics (add error handling)
    input_reads_count_fwd=$(($(zcat "$raw_data/run$run_number/s"*R1.fastq.gz | wc -l) / 4))
    surviving_read_count_fwd=$(($(zcat "$trimmed/run$run_number/fwd/s"*R1.fastq.gz | wc -l) / 4))
    percentage_retained_fwd=$(awk "BEGIN { printf \"%.2f\", ($surviving_read_count_fwd / $input_reads_count_fwd) * 100 }")
    dropped_reads_count_fwd=$((input_reads_count_fwd - surviving_read_count_fwd))
    percentage_dropped_fwd=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads_count_fwd / $input_reads_count_fwd) * 100 }")

    input_reads_count_rev=$(($(zcat "$raw_data/run$run_number/s"*R2.fastq.gz | wc -l) / 4))
    surviving_read_count_rev=$(($(zcat "$trimmed/run$run_number/rev/s"*R2.fastq.gz | wc -l) / 4))
    percentage_retained_rev=$(awk "BEGIN { printf \"%.2f\", ($surviving_read_count_rev / $input_reads_count_rev) * 100 }")
    dropped_reads_count_rev=$((input_reads_count_rev - surviving_read_count_rev))
    percentage_dropped_rev=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads_count_rev / $input_reads_count_rev) * 100 }")

    # Calculate mean reads per sample
    mean_reads_per_sample_fwd=$(awk "BEGIN { printf \"%.2f\", $surviving_read_count_fwd / $NUM_RUNS }")
    mean_reads_per_sample_rev=$(awk "BEGIN { printf \"%.2f\", $surviving_read_count_rev / $NUM_RUNS }")

    # Calculate minimum reads in a sample
    min_reads_in_sample_fwd=$(awk "BEGIN { min = $surviving_read_count_fwd; } $surviving_read_count_fwd < min { min = $surviving_read_count_fwd; } END { print min; }")
    min_reads_in_sample_rev=$(awk "BEGIN { min = $surviving_read_count_rev; } $surviving_read_count_rev < min { min = $surviving_read_count_rev; } END { print min; }")

    # Calculate maximum reads in a sample
    max_reads_in_sample_fwd=$(awk "BEGIN { max = 0; } $surviving_read_count_fwd > max { max = $surviving_read_count_fwd; } END { print max; }")
    max_reads_in_sample_rev=$(awk "BEGIN { max = 0; } $surviving_read_count_rev > max { max = $surviving_read_count_rev; } END { print max; }")

    # Calculate the number of samples with less than 10K reads
    samples_less_than_10k_fwd=$(awk "BEGIN { count = 0; } $surviving_read_count_fwd < 10000 { count++; } END { print count; }")
    samples_less_than_10k_rev=$(awk "BEGIN { count = 0; } $surviving_read_count_rev < 10000 { count++; } END { print count; }")

    # Calculate the number of samples with 10K to 20K reads
    samples_10K_to_20K_fwd=$(awk "BEGIN { count = 0; } $surviving_read_count_fwd >= 10000 && $surviving_read_count_fwd <= 20000 { count++; } END { print count; }")
    samples_10K_to_20K_rev=$(awk "BEGIN { count = 0; } $surviving_read_count_rev >= 10000 && $surviving_read_count_rev <= 20000 { count++; } END { print count; }")

    # Calculate the number of samples with more than 20K reads
    samples_more_than_20k_fwd=$(awk "BEGIN { count = 0; } $surviving_read_count_fwd > 20000 { count++; } END { print count; }")
    samples_more_than_20k_rev=$(awk "BEGIN { count = 0; } $surviving_read_count_rev > 20000 { count++; } END { print count; }")

    # Append results to the respective files (add error handling)
    printf "$run_number\t$input_reads_count_fwd\t$surviving_read_count_fwd\t$percentage_retained_fwd\t$dropped_reads_count_fwd\t$percentage_dropped_fwd\t$mean_reads_per_sample_fwd\t$min_reads_in_sample_fwd\t$max_reads_in_sample_fwd\t$samples_less_than_10k_fwd\t$samples_10K_to_20K_fwd\t$samples_more_than_20k_fwd\n" >> "$results_file_fwd"
    
    printf "$run_number\t$input_reads_count_rev\t$surviving_read_count_rev\t$percentage_retained_rev\t$dropped_reads_count_rev\t$percentage_dropped_rev\t$mean_reads_per_sample_rev\t$min_reads_in_sample_rev\t$max_reads_in_sample_rev\t$samples_less_than_10k_rev\t$samples_10K_to_20K_rev\t$samples_more_than_20k_rev\n" >> "$results_file_rev"
}

# Main script

log "Starting at: $(date)"

# Create trimmed read subdirectories for each run
for dir in "02.Quality_trimmed"; do
    for run in $(seq 1 "$NUM_RUNS"); do
        mkdir -p "$path/$dir/run$run/fwd"
        mkdir -p "$path/$dir/run$run/rev"
    done
done

# Run Trimmomatic and generate quality reports for each sample and read tracking statistics for each run
for run_number in $(seq 1 "$NUM_RUNS"); do
    run_trimmomatic_and_report "$run_number"
    calculate_and_log_statistics "$run_number"
done

log "Ending at: $(date)"
