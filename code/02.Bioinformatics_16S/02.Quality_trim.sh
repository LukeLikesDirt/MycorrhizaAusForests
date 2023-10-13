#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --mem=100G
#SBATCH --output="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/02.Bioinformatics_16S/slurm/%x.%j.out"

# Script: Quality truncate and trim Illumina paired-end reads using Trimmomatic.
# Purpose: Prepare Illumina paired-end 16 amplicons for denoising with DADA2.
# Author: Luke Florence
# Date: 3rd October 2023

# Script Summary:
# ---------------
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
source "/data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh"
conda activate shell

# Constants:
readonly NUM_RUNS=43            # The number of sequencing runs to be processed
readonly FILE_EXT=".fastq.gz"   # The file extension of the demultiplexed (raw) reads
readonly THREADS=8              # The number of threads to use for parallel processing
readonly WINDOW=4               # Trimmomatic: The sliding window size for averaging quality scores
readonly QUAL=10                # Trimmomatic: The quality threshold for sliding window trimming
readonly HEADCROP_FWD=20        # Trimmomatic: The number of bases to remove from the start of the forward read
readonly MINLEN_FWD=280         # Trimmomatic: The minimum length read length to retain forward reads
readonly CROP_REV=270           # Trimmomatic: Fixed truncation length of reverse reads before HEADCROP and QUAL trimming (removes bases from the distal, low-quality end of the read)
readonly HEADCROP_REV=18        # Trimmomatic: The number of bases to remove from the start of the reverse read
readonly MINLEN_REV=252         # Trimmomatic: The minimum length read length to retain reverse reads

# Directory paths and file names:
PROJECT_PATH="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"     # The path to the project directory
RAW_DATA="$PROJECT_PATH/data/AusMicrobiome/16S/01.Raw_data"                     # The path to demultiplexed reads
TRIMMED_DIR="$PROJECT_PATH/data/AusMicrobiome/16S/02.Quality_trimmed"           # The path for the trimmed reads
TRACKING_FILE_FWD="$PROJECT_PATH/data/AusMicrobiome/16S/summary_trimmed_fwd.txt"   # The path to the forward read quality trimming summary file
TRACKING_FILE_REV="$PROJECT_PATH/data/AusMicrobiome/16S/summary_trimmed_rev.txt"   # The path to the reverse read quality trimming summary file

# Create trimmed read subdirectories for each run
for run in $(seq 1 "$NUM_RUNS"); do
    mkdir -p "$TRIMMED_DIR/run$run/fwd"
    mkdir -p "$TRIMMED_DIR/run$run/rev"
done

# Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

# Function to run Trimmomatic and generate quality reports for each sample
quality_trim() {

    cd "$RAW_DATA/run$run_number"

    log "Run $run_number: Quality truncating forward reads at"

    for fwd_file in *R1$FILE_EXT; do
        trimmomatic SE "$fwd_file" "$TRIMMED_DIR/run$run_number/fwd/$fwd_file" HEADCROP:"$HEADCROP_FWD" SLIDINGWINDOW:4:"$QUAL" MINLEN:"$MINLEN_FWD" -threads "$THREADS"
    done

    log "Run $run_number: Quality truncating reverse reads at"

    for rev_file in *R2$FILE_EXT; do
        trimmomatic SE "$rev_file" "$TRIMMED_DIR/run$run_number/rev/$rev_file" CROP:"$CROP_REV" HEADCROP:"$HEADCROP_REV" SLIDINGWINDOW:4:"$QUAL" MINLEN:"$MINLEN_REV" -threads "$THREADS"
    done

}

quality_report() {

    log "Run $run_number: Generating quality report at"

    # Run fastQC to generate quality reports for each sample
    fastqc "$TRIMMED_DIR/run$run_number/fwd"/*"$FILE_EXT" -o "$TRIMMED_DIR/run$run_number/fwd"
    fastqc "$TRIMMED_DIR/run$run_number/rev"/*"$FILE_EXT" -o "$TRIMMED_DIR/run$run_number/rev"
    # Use multiQC to merge reports from all samples into a single "run" report
    multiqc "$TRIMMED_DIR/run$run_number/fwd" -o "$TRIMMED_DIR/run$run_number/fwd"
    multiqc "$TRIMMED_DIR/run$run_number/rev" -o "$TRIMMED_DIR/run$run_number/rev"
    
    # Remove intermediate files and directories
    rm -f "$TRIMMED_DIR/run$run_number"/fwd/*_fastqc.zip "$TRIMMED_DIR/run$run_number"/fwd/*fastqc.html
    rm -f "$TRIMMED_DIR/run$run_number"/rev/*_fastqc.zip "$TRIMMED_DIR/run$run_number"/rev/*fastqc.html
    rm -rf "$TRIMMED_DIR/run$run_number"/fwd/multiqc_data
    rm -rf "$TRIMMED_DIR/run$run_number"/rev/multiqc_data
}

# Function to track reads across the pipeline
track_reads() {
    
    log "Run $run_number: Calculating read tarcking statistics at"

    input_dir_fwd="$RAW_DATA/run$run_number/"
    input_dir_rev="$RAW_DATA/run$run_number/"
    output_dir_fwd="$TRIMMED_DIR/run$run_number/fwd"
    output_dir_rev="$TRIMMED_DIR/run$run_number/rev"

    # Initialise read tracking files if they don't exist
    if [ ! -f "$TRACKING_FILE_FWD" ]; then
        printf "Run\tInput_Reads\tOutput_Reads\tRetained_Reads_%%\tDropped_Reads_%%\tMean_Reads_Per_Sample\tMin_Reads\tMax_Reads\tSamples_LT_10K_Reads\tSamples_10-20K_Reads\tSamples_MT_20K_Reads\n" > "$TRACKING_FILE_FWD"
    fi

    if [ ! -f "$TRACKING_FILE_REV" ]; then
        printf "Run\tInput_Reads\tOutput_Reads\tRetained_Reads_%%\tDropped_Reads_%%\tMean_Reads_Per_Sample\tMin_Reads\tMax_Reads\tSamples_LT_10K_Reads\tSamples_10-20K_Reads\tSamples_MT_20K_Reads\n" > "$TRACKING_FILE_REV"
    fi

    # Calculate the number of input, output (retained), and dropped reads during trimming for both forward and reverse reads
    # The number of lines are divided by 4 because each read occupies 4 lines in a fastq file
    input_reads_fwd=$(($(zcat "$input_dir_fwd/s"*"R1$FILE_EXT" | wc -l) / 4))
    output_reads_fwd=$(($(zcat "$output_dir_fwd/s"*"$FILE_EXT" | wc -l) / 4))
    dropped_reads_fwd=$((input_reads_fwd - output_reads_fwd))

    input_reads_rev=$(($(zcat "$input_dir_rev/s"*"R2$FILE_EXT" | wc -l) / 4))
    output_reads_rev=$(($(zcat "$output_dir_rev/s"*"$FILE_EXT" | wc -l) / 4))
    dropped_reads_rev=$((input_reads_rev - output_reads_rev))

    # Calculate the percentage of reads retained and dropped during trimming for both forward and reverse reads
    percentage_retained_fwd=$(awk "BEGIN { printf \"%.2f\", ($output_reads_fwd / $input_reads_fwd) * 100 }")
    percentage_dropped_fwd=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads_fwd / $input_reads_fwd) * 100 }")

    percentage_retained_rev=$(awk "BEGIN { printf \"%.2f\", ($output_reads_rev / $input_reads_rev) * 100 }")
    percentage_dropped_rev=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads_rev / $input_reads_rev) * 100 }")

    # Calculate mean reads per sample for both forward and reverse reads per run after trimming
    num_samples_fwd=$(find "$output_dir_fwd/" -name s*"$FILE_EXT" | wc -l)
    mean_reads_per_sample_fwd=$((output_reads_fwd / num_samples_fwd))

    num_samples_rev=$(find "$output_dir_rev/" -name s*"$FILE_EXT" | wc -l)
    mean_reads_per_sample_rev=$((output_reads_rev / num_samples_rev))

    # Find the minimum and maximum number of reads within a sample for both forward and reverse reads
    # Initialise the 'min_reads_in_sample' variable with a large number to ensure the first sample read count is less than this number
    min_reads_in_sample_fwd=999999999
    max_reads_in_sample_fwd=0

    min_reads_in_sample_rev=999999999
    max_reads_in_sample_rev=0

    for sample_file_fwd in "$output_dir_fwd/s"*"$FILE_EXT"; do
        # The number of lines are divided by 4 because each read occupies 4 lines in a fastq file
        total_reads_fwd=$(( $(zcat "$sample_file_fwd" | wc -l) / 4 ))
    
        if [ "$total_reads_fwd" -lt "$min_reads_in_sample_fwd" ]; then
            min_reads_in_sample_fwd="$total_reads_fwd"
        fi
    
        if [ "$total_reads_fwd" -gt "$max_reads_in_sample_fwd" ]; then
            max_reads_in_sample_fwd="$total_reads_fwd"
        fi
    done

    for sample_file_rev in "$output_dir_rev/s"*"$FILE_EXT"; do
        # The number of lines are divided by 4 because each read occupies 4 lines in a fastq file
        total_reads_rev=$(( $(zcat "$sample_file_rev" | wc -l) / 4 ))
    
        if [ "$total_reads_rev" -lt "$min_reads_in_sample_rev" ]; then
            min_reads_in_sample_rev="$total_reads_rev"
        fi
    
        if [ "$total_reads_rev" -gt "$max_reads_in_sample_rev" ]; then
            max_reads_in_sample_rev="$total_reads_rev"
        fi
    done

    # Calculate the number of samples with <10,000 reads, 10,000-20,000 reads, and >20,000 reads within the run for both forward and reverse reads
    samples_less_than_10k_fwd=0
    samples_10K_to_20K_fwd=0
    samples_more_than_20k_fwd=0

    samples_less_than_10k_rev=0
    samples_10K_to_20K_rev=0
    samples_more_than_20k_rev=0

    for sample_file_fwd in "$output_dir_fwd/s"*"$FILE_EXT"; do
        # The number of lines are divided by 4 because each read occupies 4 lines in a fastq file
        total_reads_fwd=$(($(zcat "$sample_file_fwd" | wc -l) / 4))

        if [ "$total_reads_fwd" -lt 10000 ]; then
            ((samples_less_than_10k_fwd++))
        elif [ "$total_reads_fwd" -ge 10000 ] && [ "$total_reads_fwd" -lt 20000 ]; then
            ((samples_10K_to_20K_fwd++))
        else
            ((samples_more_than_20k_fwd++))
        fi
    done

    for sample_file_rev in "$output_dir_rev/s"*"$FILE_EXT"; do
        # The number of lines are divided by 4 because each read occupies 4 lines in a fastq file
        total_reads_rev=$(($(zcat "$sample_file_rev" | wc -l) / 4))

        if [ "$total_reads_rev" -lt 10000 ]; then
            ((samples_less_than_10k_rev++))
        elif [ "$total_reads_rev" -ge 10000 ] && [ "$total_reads_rev" -lt 20000 ]; then
            ((samples_10K_to_20K_rev++))
        else
            ((samples_more_than_20k_rev++))
        fi
    done

    # Append the results to the result files for both forward and reverse reads
    printf "$run_number\t$input_reads_fwd\t$output_reads_fwd\t$percentage_retained_fwd\t$percentage_dropped_fwd\t$mean_reads_per_sample_fwd\t$min_reads_in_sample_fwd\t$max_reads_in_sample_fwd\t$samples_less_than_10k_fwd\t$samples_10K_to_20K_fwd\t$samples_more_than_20k_fwd\n" >> "$TRACKING_FILE_FWD"

    printf "$run_number\t$input_reads_rev\t$output_reads_rev\t$percentage_retained_rev\t$percentage_dropped_rev\t$mean_reads_per_sample_rev\t$min_reads_in_sample_rev\t$max_reads_in_sample_rev\t$samples_less_than_10k_rev\t$samples_10K_to_20K_rev\t$samples_more_than_20k_rev\n" >> "$TRACKING_FILE_REV"
}

###############################################################################
### Main script ###############################################################
###############################################################################

log "Starting at:"

# Run Trimmomatic and generate quality reports for each sample and read tracking statistics for each run
for run_number in $(seq 1 "$NUM_RUNS"); do
    quality_trim "$run_number"
    quality_report "$run_number"
done

for run_number in $(seq 1 "$NUM_RUNS"); do
    track_reads "$run_number"
done

log "Ending at:"
