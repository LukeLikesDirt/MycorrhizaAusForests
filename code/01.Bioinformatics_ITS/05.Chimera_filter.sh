#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=day
#SBATCH --mem-per-cpu=32G
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/%x.%j.out

## Script: Remove chimeras using de novo and reference-based chimera detection in VSEARCH.
## Purpose: To remove chimeras prior to clustering and generating the OTU table.
## Credit: This script is adapted from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
## Author: Luke Florence.
## Date: 4th October 2023.
##
## Script Overview:
## ---------------
##   (1) Dereplicate across samples.
##   (2) De novo chimera detection.
##   (3) Reference-based chimera detection.
##   (4) Extract all dereplicated non-chimeric sequences.
##
## Pre-requisites:
## ---------------
## By this point reads should be quality processed and dereplicated within
## samples. If processing ITS sequences, reads should be ITS extracted.
##
## There should be one file named 'all.fasta' formatted for VSEARCH. That is,
## each unique sequence should occupy two lines, formatted as follows:
##  (1) <sample_name>.fasta.<unique_ID>;size=<number_of_reads>;
##  (2) <the actual sequence>
##
## There should also be a mapping file named 'map.pl'. The mapping file can be
## found in this repository or on the VSEARCH wiki linked at the top of the 
## page. Otherwise here is a link to wiki with an example VSEARCH pipeline that
## doesn't use a map.pl file: https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline
##
## Additional notes:
## ----------------
## Here, I will dereplicate across samples and perform both de-novo and
## reference-based chimera filtering in VSEARCH. The previous quality and
## denosing steps in this pipeline were done with DADA2. I have chosen to
## switch from DADA2 to VSEARCH to perform chimera filtering because:
##  (1) De novo chimera detection with DADA2 results in has high rates
##      of false positives compared to the UCHIME3 algorithm that is used
##      by VSEARCH (Edgar 2016, bioRxiv 074252; https://doi.org/10.1101/074252).
##  (2) VSEARCH implements reference-based chimera detection whereas VSEARCH
##      does not. Despite de novo chimera detection outperforming the 
##      reference-based approach in terms of true positive detection, the
##      referenced-based approach is a recommended supplement as most 
##      referenced-based chimeras are true chimeras (Tedersoo et al. 2022, Molecular ecology 31, no. 10 (2022): 2769-2795).

readonly THREADS=8                                                                                  ## Set the number of threads
readonly PROJECT_PATH="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"                ## The path to the project directory
readonly MAP_SCRIPT="$PROJECT_PATH/data/AusMicrobiome/ITS/map.pl"                                   ## Map file for fasta reconstruction
readonly REFERENCE_SEQS="$PROJECT_PATH/data/AusMicrobiome/ITS/06.Reference_dataset/ITS1/ITS1.fasta" ## Path to UNITE reference dataset
readonly DENOISED_DIR="$PROJECT_PATH/data/AusMicrobiome/ITS/05.Denoised"                            ## Path to denoised fasta file
readonly CHIMERA_FILTERED_DIR="$PROJECT_PATH/data/AusMicrobiome/ITS/07.Chimera_filtered"            ## Path for chimera filtered fasta file
mkdir -p "$CHIMERA_FILTERED_DIR"                                                                    ## Chimera filtered directory
# File to track representative sequences and reads across the pipeline
LOG_FILE="$PROJECT_PATH/code/02.Bioinformatics_ITSS/slurm/%x.%j.out"
TRACK_REPSEQS_READS_FILE="$PROJECT_PATH/data/AusMicrobiome/ITS/summary_chimeraFiltered.txt"


## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp" | tee -a "$LOG_FILE"
}

## Function for dereplication across samples, de novo and referenced-based chimera filtering, and mapping reads back to dereplicated non-chimeric sequences
chimera_filter() {

    log 'Dereplicating across samples at:'

    vsearch \
        --derep_fulllength "$DENOISED_DIR/all.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.derep.uc" \
        --output "$CHIMERA_FILTERED_DIR/all.derep.fasta"

    log 'De novo chimera detection at:'

    vsearch \
        --uchime3_denovo "$CHIMERA_FILTERED_DIR/all.derep.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta"

    log 'Reference-based chimera detection at:'

    vsearch \
        --uchime_ref "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta" \
        --threads "$THREADS" \
        --db "$REFERENCE_SEQS" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta"

    ## Map reads back to dereplicated non-chimeric sequences
    log 'Extract all dereplicated non-chimeric sequences at:'

    perl "$MAP_SCRIPT" "$CHIMERA_FILTERED_DIR/all.derep.fasta" "$CHIMERA_FILTERED_DIR/all.derep.uc" "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta"

}

## Functions to track the number of representative sequences and reads across the pipeline

get_rep_seq_count() {
    file="$1"
    count=$(grep -c "^>" "$file")
    echo "$count"
}

track_representative_sequences() {
    log 'Track the number of representative sequences across the pipeline:' | tee -a "$TRACK_REPSEQS_READS_FILE"
    unique_sequences_deno="$(get_rep_seq_count "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta")"
    unique_sequences_ref="$(get_rep_seq_count "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta")"
    total_dereplicated_seqs="$(get_rep_seq_count "$CHIMERA_FILTERED_DIR/all.derep.fasta")"
    
     # Calculate the percentages of chimeras
    percentage_denovo_chimeras=$(bc <<< "scale=2; 100 * $unique_sequences_deno / $total_dereplicated_seqs")
    percentage_reference_chimeras=$(bc <<< "scale=2; 100 * $unique_sequences_ref / $total_dereplicated_seqs")
    total_chimeras_percentage=$(bc <<< "scale=2; 100 * ($unique_sequences_deno + $unique_sequences_ref) / $total_dereplicated_seqs")
    
    printf '    Unique sequence input: %s\n' "$(get_rep_seq_count "$DENOISED_DIR/all.fasta")" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Unique sequences after dereplication: %s\n' "$total_dereplicated_seqs" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Unique sequences after de novo chimera detection: %s (%.2f%% chimeras)\n' "$unique_sequences_deno" "$percentage_denovo_chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Unique sequences after reference-based chimera detection: %s (%.2f%% chimeras)\n' "$unique_sequences_ref" "$percentage_reference_chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '    Percentage of total chimeras relative to dereplicated unique sequences: %.2f%%\n' "$total_chimeras_percentage" | tee -a "$TRACK_REPSEQS_READS_FILE"
}

get_reads_count() {
    file="$1"
    count=$(grep -o 'size=[0-9]\+' "$file" | awk -F'=' '{ sum += $2 } END { print sum }')
    echo "$count"
}

track_reads() {
    log 'Track the number of reads across the pipeline:' | tee -a "$TRACK_REPSEQS_READS_FILE"
    total_dereplicated_reads="$(get_reads_count "$CHIMERA_FILTERED_DIR/all.derep.fasta")"
    reads_denovo="$(get_reads_count "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta")"
    reads_reference="$(get_reads_count "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta")"
    
    # Calculate the percentages of chimeras
    percentage_denovo_chimeras=$(bc <<< "scale=2; 100 * $reads_denovo / $total_dereplicated_reads")
    percentage_reference_chimeras=$(bc <<< "scale=2; 100 * $reads_reference / $total_dereplicated_reads")
    total_chimeras_percentage=$(bc <<< "scale=2; 100 * ($reads_denovo + $reads_reference) / $total_dereplicated_reads")
    
    printf '   Reads input: %s\n' "$(get_reads_count "$DENOISED_DIR/all.fasta")" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '   Reads after dereplication: %s\n' "$total_dereplicated_reads" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '   Reads after de novo chimera detection: %s (%.2f%% chimeras)\n' "$reads_denovo" "$percentage_denovo_chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '   Reads after reference-based chimera detection: %s (%.2f%% chimeras)\n' "$reads_reference" "$percentage_reference_chimeras" | tee -a "$TRACK_REPSEQS_READS_FILE"
    printf '   Percentage of total chimeras relative to dereplicated reads: %.2f%%\n' "$total_chimeras_percentage" | tee -a "$TRACK_REPSEQS_READS_FILE"
}

###############################################################################
### Main script ###############################################################
###############################################################################

log 'Starting at:'

## Activate the conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

## Dereplicate across samples and remove chimeras
chimera_filter
## Track the number of representative sequences and reads across the pipeline
track_representative_sequences
track_reads

## Deactivate the conda environment
conda deactivate

log 'Finishing at:'
