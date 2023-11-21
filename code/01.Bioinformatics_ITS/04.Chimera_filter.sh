#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --partition=day
#SBATCH --output=../01.Bioinformatics_ITS/slurm/%x.%j.out

# Constants and subdirectories
readonly THREADS=8
readonly IDENTITY=0.97
readonly MAP_SCRIPT="../01.Bioinformatics_ITS/map.pl"
readonly REFERENCE_SEQS="../../data/AusMicrobiome/ITS/06.Reference_dataset/ITS1/ITS1.fasta"
CHIMERA_METHODS_DIR=../../data/AusMicrobiome/ITS/07.Chimera_filtered
TRACK_REPSEQS_READS_DIR="../../data/AusMicrobiome/ITS/"

# Define extensions for different methods used for denoising
readonly method_DADA2="_DADA2"
readonly method_UNOISE3="_UNOISE3"

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp" | tee -a "$LOG_FILE"
}
LOG_FILE="slurm/%x.%j.out"

## Function for dereplication across samples, de novo and referenced-based chimera filtering
chimera_filter() {
    local method=$1
    local CHIMERA_FILTERED_DIR="$CHIMERA_METHODS_DIR$method"

    # Remove downstream files if re-running the pipeline
    local files_to_remove=("all.denovo.nonchimeras.fasta" "all.derep.fasta" "all.derep.uc" "all.nonchimeras.derep.fasta" "all.nonchimeras.fasta" "all.preclustered.fasta" "all.preclustered.uc" "all.ref.nonchimeras.fasta")

    for file in "${files_to_remove[@]}"; do
        rm -f "$CHIMERA_FILTERED_DIR/$file"
    done

    log 'Dereplicating across samples at:'
    vsearch \
        --derep_fulllength "$CHIMERA_FILTERED_DIR/all.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.derep.uc" \
        --output "$CHIMERA_FILTERED_DIR/all.derep.fasta"

    log 'Preclustering reads at:'
    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.derep.fasta" \
        --threads $THREADS \
        --id $IDENTITY \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.preclustered.uc" \
        --centroids "$CHIMERA_FILTERED_DIR/all.preclustered.fasta"

    log 'Starting de novo chimera detection at:'
    vsearch \
        --uchime3_denovo "$CHIMERA_FILTERED_DIR/all.preclustered.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta"

    log 'Starting reference-based chimera detection at:'
    vsearch \
        --uchime_ref "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta" \
        --threads "$THREADS" \
        --db "$REFERENCE_SEQS" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta"

    log 'Extracting all non-chimeric sequences at:'
    perl "$MAP_SCRIPT" "$CHIMERA_FILTERED_DIR/all.derep.fasta" "$CHIMERA_FILTERED_DIR/all.preclustered.uc" "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta"
    perl "$MAP_SCRIPT" "$CHIMERA_FILTERED_DIR/all.fasta" "$CHIMERA_FILTERED_DIR/all.derep.uc" "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta"
}

## Functions to track the number of representative sequences and reads across the pipeline

get_rep_seq_count() {
    file="$1"
    count=$(grep -c "^>" "$file")
    echo "$count"
}

track_representative_sequences() {
    local method=$1
    local CHIMERA_FILTERED_DIR="$CHIMERA_METHODS_DIR$method"
    local TRACK_REPSEQS_READS_FILE="$TRACK_REPSEQS_READS_DIR/summary_denoised$method.txt"

    log 'Track the number of representative sequences across the pipeline:' | tee -a "$TRACK_REPSEQS_READS_FILE"
    unique_sequences_denovo="$(get_rep_seq_count "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta")"
    unique_sequences_reference="$(get_rep_seq_count "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta")"
    total_cluatered_seqs="$(get_rep_seq_count "$CHIMERA_FILTERED_DIR/all.preclustered.fasta")"

    percentage_denovo_chimeras=$(bc <<< "scale=2; 100 - (100 * $unique_sequences_denovo / $total_cluatered_seqs)")
    percentage_reference_chimeras=$(bc <<< "scale=2; (100 - (100 * $unique_sequences_reference / $total_cluatered_seqs)) - $percentage_denovo_chimeras")
    total_chimeras_percentage=$(bc <<< "scale=2; $percentage_denovo_chimeras + $percentage_reference_chimeras")

    printf '    Number of unique sequences input: %s\n' "$(get_rep_seq_count "$CHIMERA_FILTERED_DIR/all.fasta")" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Unique sequences after dereplication and clustering: %s\n' "$total_cluatered_seqs" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Unique sequences after de novo chimera detection: %s (%.2f%% chimeras)\n' "$unique_sequences_denovo" "$percentage_denovo_chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Unique sequences after reference-based chimera detection: %s (%.2f%% chimeras)\n' "$unique_sequences_reference" "$percentage_reference_chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Percentage of total chimeric sequences relative to unique (clustered) sequences: %.2f%%\n' "$total_chimeras_percentage chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
}

get_reads_count() {
    file="$1"
    count=$(grep -o 'size=[0-9]\+' "$file" | awk -F'=' '{ sum += $2 } END { print sum }')
    echo "$count"
}

track_reads() {
    local method=$1
    local CHIMERA_FILTERED_DIR="$CHIMERA_METHODS_DIR$method"
    local TRACK_REPSEQS_READS_FILE="$TRACK_REPSEQS_READS_DIR/summary_denoised$method.txt"

    log 'Track the number of reads across the pipeline:' | tee -a "$TRACK_REPSEQS_READS_FILE"
    total_reads="$(get_reads_count "$CHIMERA_FILTERED_DIR/all.preclustered.fasta")"
    reads_denovo="$(get_reads_count "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta")"
    reads_reference="$(get_reads_count "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta")"
    
    percentage_denovo_chimeras=$(bc <<< "scale=2; 100 - (100 * $reads_denovo / $total_reads)")
    percentage_reference_chimeras=$(bc <<< "scale=2; (100 - (100 * $reads_reference / $total_reads)) - $percentage_denovo_chimeras")
    total_chimeras_percentage=$(bc <<< "scale=2; $percentage_denovo_chimeras + $percentage_reference_chimeras")
    
    printf '    Number of reads input: %s\n' "$total_reads" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Reads after de novo chimera detection: %s (%.2f%% chimeras)\n' "$reads_denovo" "$percentage_denovo_chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Reads after reference-based chimera detection: %s (%.2f%% chimeras)\n' "$reads_reference" "$percentage_reference_chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Percentage of total chimeric reads relative to total reads: %.2f%%\n' "$total_chimeras_percentage chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
}

###############################################################################
### Main script ###############################################################
###############################################################################

log 'Starting at:'

# Activate the conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

# Chimera detection for both methods in parallel
(chimera_filter "$method_DADA2") &
(chimera_filter "$method_UNOISE3") &

# Wait for background processes to finish
wait

# Deactivate the conda environment
conda deactivate

# Track the number of representative sequences across the pipeline
(track_representative_sequences "$method_DADA2") &
(track_representative_sequences "$method_UNOISE3") &

# Wait for background processes to finish
wait

# Track the number of reads across the pipeline
(track_reads "$method_DADA2") &
(track_reads "$method_UNOISE3") &

# Wait for background processes to finish
wait

for method in "$method_DADA2" "$method_UNOISE3"; do
    printf "Current method: $method\n"

    CHIMERA_FILTERED_DIR="$CHIMERA_METHODS_DIR$method"
    TRACK_REPSEQS_READS_FILE="$TRACK_REPSEQS_READS_DIR/summary_denoised$method.txt"

    # Activate the conda environment
    conda activate shell

    # Dereplicate across samples and remove chimeras in the background
    chimera_filter &

    # Deactivate the conda environment
    conda deactivate
done

# Wait for all background processes to finish
wait

# Track the number of representative sequences and reads across the pipeline
for method in "$method_DADA2" "$method_UNOISE3"; do
    printf "Current method: $method\n"

    CHIMERA_FILTERED_DIR="$CHIMERA_METHODS_DIR$method"
    TRACK_REPSEQS_READS_FILE="$TRACK_REPSEQS_READS_DIR/summary_chimeras$method.txt"

    track_representative_sequences
    track_reads
    
done

log 'Finishing at:'