#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=720:00:00
#SBATCH --partition=month
#SBATCH --mem-per-cpu=32G
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/%x.%j.out

# Script: Extract the ITS, quality trim and quality check Illumina single-end reads.
# Purpose: Prepare Illumina forwards reads targeting the ITS1 for denoising with DADA2.
# Author: Luke Florence
# Date: 28th September 2023

# Software:
# --------
# ITSxpress v2.0.0: https://github.com/USDA-ARS-GBRU/itsxpress
# Trimmomatic v0.36: http://www.usadellab.org/cms/?page=trimmomatic
# FastQC v0.12.1: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# MultiQC v1.15: https://multiqc.info/

# Script Overview:
# ---------------
# This script performs the following tasks:
#   (1) Extract the ITS region with ITSxpress
#   (2) Quality truncate reads with Trimmomatic
#   (3) Quality check trimmed reads for each sample with FastQC
#   (4) Track and log reads across the pipeline. Reads are grouped by sequencing run.
#
# The sequencing data processed here span 44 runs. Sequencing runs are
# processed separately prior to denoising with DADA2 to improve error rate
# estimations, as each sequencing run has a unique error profile.

# Pre-requisites:
# --------------
# The raw reads should be demultiplexed and separated into separate
# directories based on sequencing run. The sequencing run subdirectories
# should be named as 'run1', 'run2', 'run3', etc.
#
# The read tracking function expects the raw read file names for test samples
# to begin with the letter 's', for example, 's1_L001_R1_001.fastq.gz'. This
# is because I want to ignore the controls, particulaly the negative controls,
# which begin with 'n' in my case, as negative controls are expected to have
# very low read counts which will bias the read tracking calculations. Either 
# rename test sample files, amend the the 'track_reads' function by removing
# the letter 's', or remove the 'track_reads' function from the 'Main script'
# section if required.

# Constants
readonly THREADS=8                # The number of threads to use for parallel processing
readonly FIRST_RUN=29             # The first run to process
readonly LAST_RUN=32              # The last run to process
readonly FILE_EXT=".fastq.gz"     # The file extension of the raw reads
readonly ITSX="ITS1"              # ITSxpress: The ITS region to extract. Options: ITS1, ITS2 or All
readonly TAXA="All"               # ITSxpress: The target taxonomic group: Alveolata, Bryophyta, Bacillariophyta, Amoebozoa, Euglenozoa, Fungi, Chlorophyta, Rhodophyta, Phaeophyceae, Marchantiophyta, Metazoa, Oomycota, Haptophyceae, Raphidophyceae, Rhizaria, Synurophyceae, Tracheophyta, Eustigmatophyceae, All. Default Fungi.
readonly CLUSTER=1.0              # ITSxpress: The identity threshold for clustering reads. Range 0.99-1.0. Set to 1 for exact de-replication. Default 1.0.
readonly WINDOW=4                 # Trimmomatic: The sliding window size for averaging quality scores
readonly QUAL=13                  # Trimmomatic: The quality threshold for sliding window trimming
readonly MINLEN=80                # Trimmomatic: The minimum read length after trimming

# Directory paths and file names
PROJECT_PATH="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"    # The path to the project directory
RAW_DATA="$PROJECT_PATH/data/AusMicrobiome/ITS/01.Raw_data"                    # Path to demultiplexed reads
ITS_EXTRACTED="$PROJECT_PATH/data/AusMicrobiome/ITS/02.ITS_extracted"          # Path for the extracted ITS reads
QUALITY_TRIMMED="$PROJECT_PATH/data/AusMicrobiome/ITS/03.Quality_trimmed"      # Path for the quality truncated reads
TRACK_READS="$PROJECT_PATH/data/AusMicrobiome/ITS/summary_trimmed.txt"                    # Path to the file for tracking reads across the pipeline

# Create subdirectories if they don't exist
mkdir -p "$ITS_EXTRACTED"
mkdir -p "$QUALITY_TRIMMED"
for dir in "$ITS_EXTRACTED" "$QUALITY_TRIMMED"; do
  for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    mkdir -p "$dir/run$run_number"
  done
done

# Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

# Function to extract ITS region of interest with ITSxpress
extract_its() {

    log "Run $run_number: Extracting ITS1 with ITSxpress at"
    cd "$RAW_DATA/run$run_number"

    for f in *"$FILE_EXT"; do
        itsxpress \
            --single_end \
            --fastq "${f}" \
            --cluster_id "$CLUSTER" \
            --region "$ITSX" \
            --taxa "$TAXA" \
            --log "$ITS_EXTRACTED/run${run_number}/logfile.txt" \
            --outfile "$ITS_EXTRACTED/run${run_number}/${f}" \
            --threads "$THREADS"
    done
}

# Function to quality trim reads with Trimmomatic
quality_trim() {
    
    log "Run $run_number: Quality trimming with Trimmomatic at"
    cd "$ITS_EXTRACTED/run$run_number"

    for f in *"$FILE_EXT"; do
        trimmomatic SE -threads "$THREADS" "$f" "$QUALITY_TRIMMED/run$run_number/$f" \
            SLIDINGWINDOW:"$WINDOW":"$QUAL" MINLEN:"$MINLEN" 
    done
}

# Function to generate read quality report for each sample using fastQC and multiQC
generate_quality_report() {
    
    log "Run $run_number: Generating quality report at"
    fastqc "$QUALITY_TRIMMED/run$run_number"/*R1.fastq.gz -o "$QUALITY_TRIMMED/run$run_number"
    multiqc "$QUALITY_TRIMMED/run$run_number" -o "$QUALITY_TRIMMED/run$run_number"
    
    # Remove intermediate files and 'multiqc_data' directory
    rm "$QUALITY_TRIMMED/run$run_number"/*fastqc.zip "$QUALITY_TRIMMED/run$run_number"/*fastqc.html
    rm "$QUALITY_TRIMMED/run$run_number"/multiqc_data/*
    rmdir "$QUALITY_TRIMMED/run$run_number"/multiqc_data
}

# Function to calculate statistics for tracking reads
track_reads() {
    
    log "Run $run_number: Calculating trimming statistics at"

    input_dir="$RAW_DATA/run$run_number"
    output_dir="$QUALITY_TRIMMED/run$run_number"

    # Initialise read tracking file if it doesn't exist
    if [ ! -f "$TRACK_READS" ]; then
        printf "Run\tInput_Reads\tOutput_Reads\tRetained_Reads_%%\tDropped_Reads_%%\tMean_Reads_Per_Sample\tMin_Reads\tMax_Reads\tSamples_LT_10K_Reads\tSamples_10-20K_Reads\tSamples_MT_20K_Reads\n" > "$TRACK_READS"
    fi

    # Calculate the number of input, output (retained), and dropped reads during trimming
    # The number of lines are divided by 4 because each read occupies 4 lines in a fastq file
    input_reads=$(($(zcat "$input_dir/s"*"$FILE_EXT"| wc -l) / 4))
    output_reads=$(($(zcat "$output_dir/s"*"$FILE_EXT"| wc -l) / 4))
    dropped_reads=$((input_reads - output_reads))

    # Calculate the percentage of reads retained and dropped during trimming
    percentage_retained=$(awk "BEGIN { printf \"%.2f\", ($output_reads / $input_reads) * 100 }")
    percentage_dropped=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads / $input_reads) * 100 }")
    
    # Calculate mean reads per sample per run after trimming
    num_samples=$(find "$output_dir/" -name "s"*"$FILE_EXT" | wc -l)
    mean_reads_per_sample=$((output_reads / num_samples))

    # Find the minimum and maximum number of reads within a sample
    # Initialise the 'min_reads_in_sample' variable with a large number to ensure the first sample read count is less than this number
    min_reads_in_sample=999999999
    max_reads_in_sample=0

    for sample_file in "$output_dir/s"*"$FILE_EXT"; do
        # The number of lines are divided by 4 because each read occupies 4 lines in a fastq file
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

    for sample_file in "$output_dir/s"*"$FILE_EXT"; do
        # The number of lines are divided by 4 because each read occupies 4 lines in a fastq file
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
    printf "$run_number\t$input_reads\t$output_reads\t$percentage_retained\t$percentage_dropped\t$mean_reads_per_sample\t$min_reads_in_sample\t$max_reads_in_sample\t$samples_less_than_10k\t$samples_10K_to_20K\t$samples_more_than_20k\n" >> "$TRACK_READS"
}

###############################################################################
### Main script ###############################################################
###############################################################################

log "Starting at:"

# Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

# Execute the functions
for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    extract_its "$run_number"
    quality_trim "$run_number"
    generate_quality_report "$run_number"
    track_reads "$run_number"
done

# Deactivate conda environment
conda deactivate

log "Finished at:"
