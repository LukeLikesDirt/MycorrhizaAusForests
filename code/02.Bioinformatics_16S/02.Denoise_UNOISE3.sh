#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --mem=100G
#SBATCH --output="../02.Bioinformatics_16S/slurm/%x.%j.out"

# Script:   Denoise Illumina paired-end reads using the UNOISE3 algorithm in
#           VSEARCH
# Author:   Luke Florence
# Date:     3rd November 2023
# Software: VSEARCH v2.22.1: https://github.com/torognes/vsearch

# Constants and subdirectories:
readonly NUM_RUNS=43            # The number of sequencing runs to be processed
readonly FILE_EXT=".fastq.gz"   # The file extension of the demultiplexed (raw) reads
readonly THREADS=8              # The number of threads to use
QUALITY_FILTERED_DIR="../../data/AusMicrobiome/16S/03.Quality_filtered"
DENOISED_DIR="../../data/AusMicrobiome/16S/04.Denoised_UNOISE3"
CHIMERA_FILTERED_DIR="../../data/AusMicrobiome/16S/06.Chimera_filtered_UNOISE3"

# Define the file for tracking reads
TRACK_READS_DENOISED="../../data/AusMicrobiome/16S/summary_denoised_UNOISE3.txt"

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

# Function to merge forward and reverse reads
merge_reads() {
    local run_number="$1"
    log "Starting read merging for run $run_number"

    local input_dir_fwd="$QUALITY_FILTERED_DIR/run$run_number/fwd"
    local input_dir_rev="$QUALITY_FILTERED_DIR/run$run_number/rev"
    local output_dir_merged="$QUALITY_FILTERED_DIR/run$run_number/merged"

    create_directory "$output_dir_merged"

    for fwd_file in "$input_dir_fwd"/*R1$FILE_EXT; do

        # Extract the sample name from the forward file by cutting at the last "_" (remove R1 and R2 extensions)
        sample_name=$(basename "$fwd_file" | rev | cut -d '_' -f 2- | rev)

        # Construct the path to the corresponding reverse file in the input_dir_rev
        rev_file=$(find "$input_dir_rev" -type f -name "${sample_name}*R2$FILE_EXT" | head -n 1)

        vsearch \
            --fastq_mergepairs "${fwd_file}" \
            --reverse "${rev_file}" \
            --threads $THREADS \
            --fastq_minovlen 12 \
            --fastq_maxdiffpct 10 \
            --fastq_minmergelen 390 \
            --fastq_maxmergelen 530 \
            --fastqout "$output_dir_merged/${sample_name%$FILE_EXT}.fastq"
    done

    log "Completed read merging for run $run_number"

    log "Compressing the merged reads"
    
    for file in "$output_dir_merged"/*.fastq; do
        gzip -c "${file}" > "${file}.gz"
        rm "${file}"
    done

    log "Completed compressing the merged reads"
}

# Denoise function
denoise() {
    local run_number="$1"
    log "Starting denoising for run $run_number"

    local input_dir="$QUALITY_FILTERED_DIR/run$run_number/merged"
    local output_dir="$DENOISED_DIR/run$run_number"

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

    input_dir="$QUALITY_FILTERED_DIR/run$run_number/rev/"
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

# Activate Conda environment
conda activate shell

# Process the fastq files
for run_number in $(seq 1 "$NUM_RUNS"); do
    merge_reads "$run_number"
    denoise "$run_number"
done

# Deactivate Conda environment
conda deactivate

# Merge all denoised fasta files into one file
# Inisiate the all.fasta file or overwrite it in case you are rerunning the analysis
cat /dev/null > "$CHIMERA_FILTERED_DIR/all.fasta"
# Run the merge fatsa file function 
for run_number in $(seq 1 "$NUM_RUNS"); do
    merge_fasta_files "$run_number"
done

# Track the reads
for run_number in $(seq 1 "$NUM_RUNS"); do
    track_reads_denoised "$run_number"
done

log "Finished at:"