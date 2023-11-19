#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --mem=100G
#SBATCH --output=../01.Bioinformatics_ITS/slurm/%x.%j.out

# Script:   Quality truncate, extract the ITS and quality filter Illumina
#           single-end reads targeting the ITS region.
# Purpose:  Prepare Illumina single-end reads for denoising.
# Author:   Luke Florence
# Date:     28th October 2023

# Software:
# --------
# ITSxpress v2.0.0: https://github.com/USDA-ARS-GBRU/itsxpress
# VSEARCH v2.22.1: https://github.com/torognes/vsearch
# Trimmomatic v0.36: http://www.usadellab.org/cms/?page=trimmomatic
# FastQC v0.12.1: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# MultiQC v1.15: https://multiqc.info/

# Script Overview:
# ---------------
# This script performs the following tasks:
#   (1) Quality truncate reads with Trimmomatic
#   (2) Extract the ITS region with ITSxpress
#   (3) Quality filter reads with VSEARCH
#   (4) Quality check reads across sample with FastQC
#   (5) Track and log reads across the pipeline: Reads are grouped by sequencing
#       run.

# Pre-requisites:
# --------------
# The raw reads should be demultiplexed and separated into unique 
# subdirectories based on sequencing run. The sequencing run subdirectories
# should be named as 'run1', 'run2', 'run3', etc. The sequencing runs are
# processed separately prior to denoising to improve error rate estimations,
# as each sequencing run has a unique error profile.

# Notes:
# -----
# A full sequencing run will take about a day to be processed. Therefore, run
# the script in batches if required.
#
# Primers could not be detected in the AusMicrobiome dataset due to the overall
# low-quality of the reads. Add a primer trimming step prior to ITSxpress if
# required.
#
# The read tracking function expects the raw read file names for test samples
# to begin with the letter 's', for example, 's1_L001_R1_001.fastq.gz'. This
# is because I want to ignore the controls, particularly the negative controls,
# which begin with 'n' in my case, as negative controls are expected to have
# low read counts that will bias the read tracking calculations. Either rename
# test sample files or amend the 'track_reads' function if required. Otherwise,
# remove the 'track_reads' function from the main script.

# Constants and subdirectories
readonly NUM_THREADS=8            # The number of threads to use for parallel processing
readonly FIRST_RUN=1              # The first run to process
readonly LAST_RUN=4               # The last run to process
readonly FILE_EXT=".fastq.gz"     # The file extension of the raw reads
readonly ITS_REGION="ITS1"        # ITSxpress: The ITS region to extract. Options: ITS1, ITS2 or All
readonly TAXA="All"               # ITSxpress: The target taxonomic group: Alveolata, Bryophyta, Bacillariophyta, Amoebozoa, Euglenozoa, Fungi, Chlorophyta, Rhodophyta, Phaeophyceae, Marchantiophyta, Metazoa, Oomycota, Haptophyceae, Raphidophyceae, Rhizaria, Synurophyceae, Tracheophyta, Eustigmatophyceae, All. Default Fungi.
readonly CLUSTER=1.0              # ITSxpress: The identity threshold for clustering reads. Range 0.99-1.0. Set to 1 for exact de-replication. Default 1.0.
readonly WINDOW=4                 # Trimmomatic: The sliding window size for averaging quality scores
readonly QUAL=13                  # Trimmomatic: The quality threshold for sliding window trimming
readonly MAXEE=1                  # VSEARCH: Maximum expected error rate
readonly MAXN=0                   # VSEARCH: Maximum number of mismatches
readonly QMAX=41                  # VSEARCH: Maximum quality score
readonly MINLEN=80                # Various functions: Minimum length of reads
RAW_DATA="../../data/AusMicrobiome/ITS/01.Raw_data"
QUALITY_TRIMMED_DIR="../../data/AusMicrobiome/ITS/02.Quality_trimmed_test"
ITS_EXTRACTED_DIR="../../data/AusMicrobiome/ITS/03.ITS_extracted_test"
QUALITY_FILTERED_DIR="../../data/AusMicrobiome/ITS/04.Quality_filtered_test"

# Define summary files for tracking reads
TRACK_READS_TRIMMED="../../data/AusMicrobiome/ITS/summary_trimmed_test.txt"
TRACK_READS_ITS="../../data/AusMicrobiome/ITS/summary_its_test.txt"
TRACK_READS_QUAL_FILTERED="../../data/AusMicrobiome/ITS/summary_quality_filtered_test.txt"


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
    local output_dir="$QUALITY_TRIMMED_DIR/run$run_number"

    create_directory "$output_dir"

    for file in "$input_dir"/*"$FILE_EXT"; do
        local filename=$(basename "$file")
        log "Trimming file $filename"
        
        trimmomatic SE -threads "$NUM_THREADS" "${file}" "$output_dir/${filename}" \
            SLIDINGWINDOW:"$WINDOW":"$QUAL" MINLEN:"$MINLEN"
        
        log "Completed trimming file $filename"
    done

    log "Completed quality trimming for run $run_number"
}

# Function to extract ITS region with ITSxpress
extract_its() {
    local run_number="$1"
    log "Starting ITS extraction with ITSxpress for run $run_number"

    local input_dir="$QUALITY_TRIMMED_DIR/run$run_number"
    local output_dir="$ITS_EXTRACTED_DIR/run$run_number"

    create_directory "$output_dir"

    for file in "$input_dir"/*"$FILE_EXT"; do
        local filename=$(basename "$file")
        log "Extracting ITS region from file $filename"

        itsxpress \
            --single_end \
            --fastq "${file}" \
            --cluster_id "$CLUSTER" \
            --region "$ITS_REGION" \
            --taxa "$TAXA" \
            --log "$output_dir/logfile.txt" \
            --outfile "$output_dir/${filename}" \
            --threads "$NUM_THREADS"

        log "Completed ITS extraction from file $filename"
    done

    log "Completed ITS extraction with ITSxpress for run $run_number"
}

# Function to perform quality filtering with VSEARCH
quality_filter() {
    local run_number="$1"
    log "Starting quality filtering with VSEARCH for run $run_number"

    local input_dir="$ITS_EXTRACTED_DIR/run$run_number"
    local output_dir="$QUALITY_FILTERED_DIR/run$run_number"

    create_directory "$output_dir"

    for file in "$input_dir"/*.fastq.gz; do
        local filename=$(basename "$file")

        vsearch \
            --fastq_filter "${file}" \
            --fastq_maxee $MAXEE \
            --fastq_minlen $MINLEN \
            --fastq_maxns $MAXN \
            --fastq_qmax $QMAX \
            --fastqout "$output_dir/${filename%$FILE_EXT}.fastq"

    done

    log "Completed quality filtering for run $run_number"

    log "Compressing the quality-filtered reads"

    for file in "$output_dir"/*"fastq"; do

        gzip -c "$file" > "${file}.gz"
        rm "$file"

    done

    log "Finised compressing the quality-filtered reads"

}

# Function to generate read quality report for each sample using FastQC and MultiQC
generate_quality_report() {
    local run_number="$1"

    log "Starting to generate the quality report for run $run_number"
    
    fastqc "$ITS_EXTRACTED_DIR/run$run_number"/*FILE_EXT -o "$ITS_EXTRACTED_DIR/run$run_number"
    multiqc "$ITS_EXTRACTED_DIR/run$run_number" -o "$ITS_EXTRACTED_DIR/run$run_number"
    
    # Remove intermediate files and 'multiqc_data' directory
    rm "$ITS_EXTRACTED_DIR/run$run_number"/*fastqc.zip "$ITS_EXTRACTED_DIR/run$run_number"/*fastqc.html
    rm -r "$ITS_EXTRACTED_DIR/run$run_number"/multiqc_data

    log "Completed generating a quality report for run $run_number"
}

# Function to calculate statistics for tracking reads after quality trimming
track_reads_quality_trim() {
    local run_number="$1"
    log "Run $run_number: Calculating trimming statistics"

    input_dir="$RAW_DATA/run$run_number"
    output_dir="$QUALITY_TRIMMED_DIR/run$run_number"

    # Initialise read tracking file if it doesn't exist
    if [ ! -f "$TRACK_READS_TRIMMED" ]; then
        printf "Run\tInput_Reads_QT\tOutput_Reads_QT\tRetained_Reads_%%_QT\tDropped_Reads_%%_QT\tMean_Reads_Per_Sample_QT\tMin_Reads_QT\tMax_Reads_QT\tSamples_LT_10K_Reads_QT\tSamples_10-20K_Reads_QT\tSamples_MT_20K_Reads_QT\n" > "$TRACK_READS_TRIMMED"
    fi

    # Calculate the number of input, output (retained), and dropped reads after quality trimming
    input_reads=$(($(zcat "$input_dir/s"*"$FILE_EXT"| wc -l) / 4))
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
    printf "$run_number\t$input_reads\t$output_reads\t$percentage_retained\t$percentage_dropped\t$mean_reads_per_sample\t$min_reads_in_sample\t$max_reads_in_sample\t$samples_less_than_10k\t$samples_10K_to_20K\t$samples_more_than_20k\t-\t-\t-\t-\n" >>  "$TRACK_READS_TRIMMED"
}

# Function to calculate statistics for tracking reads after ITS extraction
track_reads_extract_its() {
    local run_number="$1"
    log "Run $run_number: Calculating ITS extraction statistics"

    input_dir="$QUALITY_TRIMMED_DIR/run$run_number"
    output_dir="$ITS_EXTRACTED_DIR/run$run_number"

    # Initialise read tracking file if it doesn't exist
    if [ ! -f "$TRACK_READS_ITS" ]; then
        printf "Run\tInput_Reads_ITS\tOutput_Reads_ITS\tRetained_Reads_%%_ITS\tDropped_Reads_%%_ITS\tMean_Reads_Per_Sample_ITS\tMin_Reads_ITS\tMax_Reads_ITS\tSamples_LT_10K_Reads_ITS\tSamples_10-20K_Reads_ITS\tSamples_MT_20K_Reads_ITS\n" > "$TRACK_READS_ITS"
    fi

    # Calculate the number of input, output (retained), and dropped reads after ITS extraction
    input_reads=$(($(zcat "$input_dir/s"*"$FILE_EXT"| wc -l) / 4))
    output_reads=$(($(zcat "$output_dir/s"*"$FILE_EXT"| wc -l) / 4))
    dropped_reads=$((input_reads - output_reads))

    # Calculate the percentage of reads retained and dropped after ITS extraction
    percentage_retained=$(awk "BEGIN { printf \"%.2f\", ($output_reads / $input_reads) * 100 }")
    percentage_dropped=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads / $input_reads) * 100 }")

    # Calculate mean reads per sample per run after ITS extraction
    num_samples=$(find "$output_dir/" -name "s"*"$FILE_EXT" | wc -l)
    mean_reads_per_sample=$((output_reads / num_samples))

    # Find the minimum and maximum number of reads within a sample after ITS extraction
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

    # Calculate the number of samples with <10,000 reads, 10,000-20,000 reads, and >20,000 reads within the run after ITS extraction
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
    printf "$run_number\t$input_reads\t$output_reads\t$percentage_retained\t$percentage_dropped\t$mean_reads_per_sample\t$min_reads_in_sample\t$max_reads_in_sample\t$samples_less_than_10k\t$samples_10K_to_20K\t$samples_more_than_20k\t-\t-\t-\t-\n" >> "$TRACK_READS_ITS"
}

# Function to calculate statistics for tracking reads after quality filtering
track_reads_quality_filter() {
    local run_number="$1"
    log "Run $run_number: Calculating quality filtering statistics"

    input_dir="$ITS_EXTRACTED_DIR/run$run_number"
    output_dir="$QUALITY_FILTERED_DIR/run$run_number"

    # Initialize read tracking file if it doesn't exist
    if [ ! -f "$TRACK_READS_QUAL_FILTERED" ]; then
        printf "Run\tInput_Reads_QF\tOutput_Reads_QF\tRetained_Reads_%%_QF\tDropped_Reads_%%_QF\tMean_Reads_Per_Sample_QF\tMin_Reads_QF\tMax_Reads_QF\tSamples_LT_10K_Reads_QF\tSamples_10-20K_Reads_QF\tSamples_MT_20K_Reads_QF\n" > "$TRACK_READS_QUAL_FILTERED"
    fi

    # Calculate the number of input, output (retained), and dropped reads after quality filtering
    input_reads=$(($(zcat "$input_dir/s"*"$FILE_EXT"| wc -l) / 4))
    output_reads=$(($(zcat "$output_dir/s"*"$FILE_EXT"| wc -l) / 4))
    dropped_reads=$((input_reads - output_reads))

    # Calculate the percentage of reads retained and dropped after quality filtering
    percentage_retained=$(awk "BEGIN { printf \"%.2f\", ($output_reads / $input_reads) * 100 }")
    percentage_dropped=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads / $input_reads) * 100 }")

    # Calculate mean reads per sample per run after quality filtering
    num_samples=$(find "$output_dir/" -name "s"*"$FILE_EXT" | wc -l)
    mean_reads_per_sample=$((output_reads / num_samples))

    # Find the minimum and maximum number of reads within a sample after quality filtering
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

    # Calculate the number of samples with <10,000 reads, 10,000-20,000 reads, and >20,000 reads within the run after quality filtering
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
    printf "$run_number\t$input_reads\t$output_reads\t$percentage_retained\t$percentage_dropped\t$mean_reads_per_sample\t$min_reads_in_sample\t$max_reads_in_sample\t$samples_less_than_10k\t$samples_10K_to_20K\t$samples_more_than_20k\t-\t-\t-\t-\n" >> "$TRACK_READS_QUAL_FILTERED"
}

### Main script ###############################################################

log "Starting at:"

# Activate Conda environment
conda activate shell

# Process the fastq files
for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    quality_trim "$run_number"
    extract_its "$run_number"
    quality_filter "$run_number"
    generate_quality_report "$run_number"
done

# Deactivate Conda environment
conda deactivate

# Track the reads
for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    track_reads_quality_trim "$run_number"
    track_reads_extract_its "$run_number"
    track_reads_quality_filter "$run_number"
done

log "Finished at:"