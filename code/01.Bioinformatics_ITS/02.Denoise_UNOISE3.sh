#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --output="../01.Bioinformatics_ITS/slurm/%x.%j.out"

## Script:  Remove chimeras using de novo and reference-based chimera detection in VSEARCH.
## Purpose: To remove chimeras prior to clustering and generating the OTU table.
## Credit:  This script is adapted from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
## Author:  Luke Florence.
## Date:    4th November 2023.
# Software: VSEARCH v2.22.1: https://github.com/torognes/vsearch

# Constants and subdirectories
readonly FIRST_RUN=1                # The first run to process
readonly LAST_RUN=44                # The last run to process
readonly FILE_EXT=.fastq.gz         # File extension of the input files
readonly THREADS=8                  # The number of threads to use
readonly QUALITY_FILTERED_DIR="../../data/AusMicrobiome/ITS/04.Quality_filtered"
readonly DENOISED_DIR="../../data/AusMicrobiome/ITS/05.Denoised_DADA2"
readonly CHIMERA_FILTERED_DIR="../../data/AusMicrobiome/ITS/07.Chimera_filtered_UNOISE3"

# Define summary files for tracking reads
TRACK_READS_DENOISED="../../data/AusMicrobiome/ITS/summary_denoised_UNOISE3.txt"

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

# Denoise function
denoise() {
    local run_number="$1"
    log "Starting denoising for run $run_number"

    local input_dir="$QUALITY_FILTERED_DIR/run$run_number/"
    local output_dir="$DENOISED_DIR/run$run_number/"

    create_directory "$output_dir"

    for file in "$input_dir"/*"$FILE_EXT"; do

        local filename=$(basename "$file")
        local renamed_file=$(cut -d_ -f1 <<< "$filename") # Expects sequence names to be in the format: <sequence_name>_<sequening_information>.<file_extension>

        vsearch \
            --cluster_unoise "$file" \
            --threads $THREADS \
            --sizeout \
            --minsize 1 \
            --fasta_width 0 \
            --relabel "$renamed_file". \
            --centroids "$output_dir/${renamed_file%$FILE_EXT}.denoised.fasta" \
            --uc "$output_dir/${renamed_file%$FILE_EXT}.denoised.uc"

    done

    log "Finished denoising for run $run_number"

}

# Function to merge all fasta files in to one "all.fasta" file
merge_fasta_files() {
    local run_number="$1"
    log "Merging all denoised files"

    local input_dir="$DENOISED_DIR/run$run_number"
    local output_dir="$CHIMERA_FILTERED_DIR"

    create_directory "$CHIMERA_FILTERED_DIR"

    # Loop through each run and append denoised files to "all.fasta"
    for file in "$input_dir"/*.denoised.fasta; do

        cat "$file" >> "$CHIMERA_FILTERED_DIR/all.fasta"

    done

log "Finished merging denoised files"

}

# Function to calculate reads in a dereplicated fasta file
get_reads_count() {
    file="$1"
    count=$(grep -o 'size=[0-9]\+' "$file" | awk -F'=' '{ sum += $2 } END { print sum }')
    echo "$count"
}

# Function to calculate read tracking statistics across the denoising pipeline
track_reads_denoised() {
    local run_number="$1"
    log "Run $run_number: Calculating trimming statistics"

    input_dir="$QUALITY_FILTERED_DIR/run$run_number/"
    output_dir="$DENOISED_DIR/run$run_number/"

    # Initialise read tracking file if it doesn't exist
    if [ ! -f "$TRACK_READS_DENOISED" ]; then
        printf "Run\tInput_Reads\tOutput_Reads\tRetained_Reads_%%\tDropped_Reads_%%\tMean_Reads_Per_Sample\tMin_Reads\tMax_Reads\tSamples_LT_10K_Reads\tSamples_10-20K_Reads\tSamples_MT_20K_Reads\n" > "$TRACK_READS_DENOISED"
    fi

    # Calculate the number of input, output (retained), and dropped reads after quality trimming
    input_reads=$(($(zcat "$input_dir/s"*"$FILE_EXT"| wc -l) / 4))
    
    # Initialise to zero
    output_reads=0
    # Loop through each sample file and accumulate the output_reads count
    for sample_file in "$output_dir/s"*".fasta"; do
        output_reads=$((output_reads + $(get_reads_count "$sample_file")))
    done
    dropped_reads=$((input_reads - output_reads))

    # Calculate the percentage of reads retained and dropped after quality trimming
    percentage_retained=$(awk "BEGIN { printf \"%.2f\", ($output_reads / $input_reads) * 100 }")
    percentage_dropped=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads / $input_reads) * 100 }")

    # Calculate mean reads per sample per run after quality trimming
    num_samples=$(find "$output_dir/" -name "s"*.fasta | wc -l)
    mean_reads_per_sample=$((output_reads / num_samples))

    # Find the minimum and maximum number of reads within a sample
    min_reads_in_sample=999999999
    max_reads_in_sample=0

    for sample_file in "$output_dir/s"*".fasta"; do
        total_reads="$(get_reads_count "$sample_file")"
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

    for sample_file in "$output_dir/s"*".fasta"; do
         total_reads="$(get_reads_count "$sample_file")"
        if [ "$total_reads" -lt 10000 ]; then
            ((samples_less_than_10k++))
        elif [ "$total_reads" -ge 10000 ] && [ "$total_reads" -lt 20000 ]; then
            ((samples_10K_to_20K++))
        else
            ((samples_more_than_20k++))
        fi
    done

    # Append the results to the result file
    printf "$run_number\t$input_reads\t$output_reads\t$percentage_retained\t$percentage_dropped\t$mean_reads_per_sample\t$min_reads_in_sample\t$max_reads_in_sample\t$samples_less_than_10k\t$samples_10K_to_20K\t$samples_more_than_20k\n" >>  "$TRACK_READS_DENOISED"
}


### Main script ###############################################################

log "Starting at:"

# Activate conda environment
conda activate shell

# Run functions
for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    denoise "$run_number"
done

# Deactivate Conda environment
conda deactivate

# Merge all denoised fasta files into one file
# Inisiate the all.fasta file or overwrite it in case you are rerunning the analysis
cat /dev/null > "$CHIMERA_FILTERED_DIR/all.fasta"
# Run the function
for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    merge_fasta_files "$run_number"
done

# Track the reads
for run_number in $(seq "$FIRST_RUN" "$LAST_RUN"); do
    track_reads_denoised "$run_number"
done

log "Finished at:"
