#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --mem=100G
#SBATCH --output="../02.Bioinformatics_16S/slurm/%x.%j.out"

# Script:   Quality truncate Illumina paired-end reads using Trimmomatic.
# Purpose:  Prepare Illumina paired-end 16S amplicons for denoising with DADA2.
# Author:   Luke Florence
# Date:     3rd November 2023

# Summary:
# --------
# The Australian Microbiome 16S amplicon covers the V1-V3 region, which is 
# typically around 490 bases long (Allen, et al. 2016, BMC Res Notes 9, 380).
# The goal of this script is to retain good sequencing depth (>10,000 reads per
# sample) using a 532 base length threshold; up to 520 base length amplicons when
# accounting for the 12 base overlap required for merging with DADA2. Because
# the primers were not well detected in the raw reads, primers are trimmed 
# using fixed lengths.

# Script Overview:
# ---------------
#   1.  Removes 20 bases from the 5' end of the forward reads to trim the
#       forward primers (27F).
#   2.  Removes 18 bases from the 5' end of the reverse reads to trim the
#       reverse primers (519R).
#   3.  Hard truncate reverse reads to 252 to allow up to 520 amplicon after
#       merging.
#   4.  Quality truncates using a threshold of Q10 (1 in 10 error rate).

# Constants and subdirectories:
readonly NUM_RUNS=43            # The number of sequencing runs to be processed
readonly FILE_EXT=".fastq.gz"   # The file extension of the demultiplexed raw reads
readonly THREADS=8              # The number of threads
readonly WINDOW=4               # Trimmomatic: The sliding window size for averaging quality scores
readonly QUAL=10                # Trimmomatic: The quality threshold for sliding window trimming
readonly HEADCROP_FWD=20        # Trimmomatic: The number of bases to remove from the start of the forward read
readonly MINLEN_FWD=280         # Trimmomatic: The minimum length read length to retain forward reads
readonly CROP_REV=270           # Trimmomatic: Fixed truncation length of reverse reads before HEADCROP and QUAL trimming (removes bases from the distal, low-quality end of the read)
readonly HEADCROP_REV=18        # Trimmomatic: The number of bases to remove from the start of the reverse read
readonly MINLEN_REV=252         # Trimmomatic: The minimum length read length to retain reverse reads
DATA_DIR="../../data/AusMicrobiome/16S"                                  # 16S data directory
RAW_DATA="../../data/AusMicrobiome/16S/01.Raw_data"                      # The path to demultiplexed reads
QUALITY_TRIMMED_DIR="../../data/AusMicrobiome/16S/02.Quality_trimmed"    # The path for the trimmed reads

# Log function
log() {
    local timestamp
    timestamp=$(date)
    echo -e "\n$1 $timestamp\n"
}

# Function to check if a directory exists, and if not, create it
create_directory() {
    local directory="$1"
    if [ ! -d "$directory" ]; then
        mkdir -p "$directory"
    fi
}

# Function to perform quality trimming with Trimmomatic
quality_trim() {
    local run_number="$1"
    log "Starting quality trimming with Trimmomatic for run $run_number"

    local input_dir="$RAW_DATA/run$run_number"
    local output_dir_fwd="$QUALITY_TRIMMED_DIR/run$run_number/fwd"
    local output_dir_rev="$QUALITY_TRIMMED_DIR/run$run_number/rev"

    create_directory "$output_dir_fwd"
    create_directory "$output_dir_rev"

    for fwd_file in "$input_dir"/*"R1$FILE_EXT"; do
    
        local fwd_filename=$(basename "$fwd_file")
        log "Trimming file $fwd_filename"

        trimmomatic SE -threads "$THREADS" "${fwd_file}" "$output_dir_fwd/${fwd_filename}" \
        HEADCROP:"$HEADCROP_FWD" SLIDINGWINDOW:"$WINDOW":"$QUAL" MINLEN:"$MINLEN_FWD"

    done

    for rev_file in "$input_dir"/*"R2$FILE_EXT"; do
    
        local rev_filename=$(basename "$rev_file")
        log "Trimming file $rev_filename"

        trimmomatic SE -threads "$THREADS" "${rev_file}" "$output_dir_rev/${rev_filename}" \
        CROP:"$CROP_REV" HEADCROP:"$HEADCROP_REV" SLIDINGWINDOW:"$WINDOW":"$QUAL" MINLEN:"$MINLEN_REV"

    done

}

# Function to generate read quality report for each sample using FastQC and MultiQC
generate_quality_report() {

    log "Run $run_number: Generating quality report at"

    # Run fastQC to generate quality reports for each sample
    fastqc "$QUALITY_TRIMMED_DIR/run$run_number/fwd"/*"$FILE_EXT" -o "$QUALITY_TRIMMED_DIR/run$run_number/fwd"
    fastqc "$QUALITY_TRIMMED_DIR/run$run_number/rev"/*"$FILE_EXT" -o "$QUALITY_TRIMMED_DIR/run$run_number/rev"
    # Use multiQC to merge reports from all samples into a single "run" report
    multiqc "$QUALITY_TRIMMED_DIR/run$run_number/fwd" -o "$QUALITY_TRIMMED_DIR/run$run_number/fwd"
    multiqc "$QUALITY_TRIMMED_DIR/run$run_number/rev" -o "$QUALITY_TRIMMED_DIR/run$run_number/rev"
    
    # Remove intermediate files and directories
    rm -f "$QUALITY_TRIMMED_DIR/run$run_number"/fwd/*_fastqc.zip "$QUALITY_TRIMMED_DIR/run$run_number"/fwd/*fastqc.html
    rm -f "$QUALITY_TRIMMED_DIR/run$run_number"/rev/*_fastqc.zip "$QUALITY_TRIMMED_DIR/run$run_number"/rev/*fastqc.html
    rm -rf "$QUALITY_TRIMMED_DIR/run$run_number"/fwd/multiqc_data
    rm -rf "$QUALITY_TRIMMED_DIR/run$run_number"/rev/multiqc_data
}

# Function to calculate statistics for tracking reads after quality trimming
track_reads_quality_trimmed() {
    local run_number="$1"
    local read_direction="$2"
    log "Run $run_number: Calculating trimming statistics"

    input_dir="$RAW_DATA/run$run_number"
    output_dir="$QUALITY_TRIMMED_DIR/run$run_number/$read_direction"

    # Initialise read tracking file if it doesn't exist
    local track_reads_file="$DATA_DIR/summary_trimmed_${read_direction,,}.txt"
    if [ ! -f "$track_reads_file" ]; then
        printf "Run\tInput_Reads\tOutput_Reads\tRetained_Reads_%%\tDropped_Reads_%%\tMean_Reads_Per_Sample\tMin_Reads\tMax_Reads\tSamples_LT_10K_Reads\tSamples_10-20K_Reads\tSamples_MT_20K_Reads\n" > "$track_reads_file"
    fi

    # Calculate the number of input, output (retained), and dropped reads after quality trimming
    input_reads=$(($(zcat "$input_dir/s"*"R1$FILE_EXT"| wc -l) / 4))
    output_reads=$(($(zcat "$output_dir/s"*"$FILE_EXT"| wc -l) / 4))
    dropped_reads=$((input_reads - output_reads))

    # Calculate the percentage of reads retained and dropped after quality trimming
    percentage_retained=$(awk "BEGIN { printf \"%.2f\", ($output_reads / $input_reads) * 100 }")
    percentage_dropped=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads / $input_reads) * 100 }")

    # Calculate mean reads per sample per run after quality trimming
    num_samples=$(find "$output_dir/" -name "s"*"$FILE_EXT" | wc -l)
    mean_reads_per_sample=$((output_reads / num_samples))

    # Find the minimum and maximum number of reads within a sample
    min_reads_in_sample=999999999
    max_reads_in_sample=0

    for sample_file in "$output_dir/s"*"$FILE_EXT"; do
        total_reads=$(( $(zcat "$sample_file" | wc -l) / 4 ))
        if [ "$total_reads" -lt "$min_reads_in_sample" ]; then
            min_reads_in_sample="$total_reads"
        fi
        if [ "$total_reads" -gt "$max_reads_in_sample" ]; then
            max_reads_in_sample="$total_reads"
        fi
    done

    # Calculate the number of samples with <10,000 reads, 10,000-20,000 reads, and >20,000 reads within the run after quality trimming
    samples_less_than_10k=0
    samples_10K_to_20K=0
    samples_more_than_20k=0

    for sample_file in "$output_dir/s"*"$FILE_EXT"; do
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
    printf "$run_number\t$input_reads\t$output_reads\t$percentage_retained\t$percentage_dropped\t$mean_reads_per_sample\t$min_reads_in_sample\t$max_reads_in_sample\t$samples_less_than_10k\t$samples_10K_to_20K\t$samples_more_than_20k\n" >> "$track_reads_file"
}

### Main script ###############################################################

log "Starting at:"

# Activate the conda environment
conda activate shell

# Run Trimmomatic and generate quality reports for each sample and read tracking statistics for each run
for run_number in $(seq 1 "$NUM_RUNS"); do
    quality_trim "$run_number"
    quality_report "$run_number"
done

# Deactivate the conda environment
conda deactivate

# Track the reads
for run_number in $(seq 1 "$NUM_RUNS"); do
    track_reads_quality_trimmed "$run_number" "fwd"
    track_reads_quality_trimmed "$run_number" "rev"
done

log "Finished at:"