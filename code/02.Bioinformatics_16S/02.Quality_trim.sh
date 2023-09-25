#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week

echo "Starting at: $(date)"

## Prepare Illumina forward and reverse reads for denoising with DADA2 by
## quality truncating reads using Trimmomatic. The Australian Microbiome
## 16S amplicon covers the V1-V3 region, which is approximately 490 bases long,
## but varies from 450 to 600 bases. In the current SILVA reference dataset 
## used in this worklflow, ~95% of V1-V3 amplicons are shorter than 510 bases.
## Therefore, I aim to retain good sequencing depth (>10,000 reads per sample) 
## using a 512 base threshold, which equals 500 bases when accounting the 12 
## base overlap required by DADA2. General overview the '02.Quality_trim.sh'
## script:
##  - Remove 30 bases from the 5' end of the forward reads to remove the fwd
##    primers (27F), primer linker and primer pad.
##  - Remove 28 bases from the 5' end of the reverse reads to remove the rev
##    primers (519R), primer linker and primer pad.
##  - Remove 30 bases from the 3' end of the reverse reads to keep a maximum 
##    potential amplicon length of 510 bases after merging.
##  - Quality truncate using a threshold of Q10, which equals a probable error
##    rate of 1 in 10 calls. Considering that DADA2 will handel many incorrect
##    calls when denoising and because reads will be merged, this 'relaxed' 
##    threshold seems reasonable. Tighten to at least Q13 if using fwd reads 
##    only.
##  - Remove reads shorter than 270 (fwd) and 242 (rev) bases after quality 
##    truncation.

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

## Organize directories:
## Path to the main data directory
path="/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/16S"
## Create subdirectories for each run and for both fwd (R1) and rev (R2) reads
for dir in "02.Quality_trimmed"; do
  for run in {1..4}; do
    mkdir -p "$path/$dir/run$run/fwd"
    mkdir -p "$path/$dir/run$run/rev"
  done
done

## Define paths to subdirectory:
raw_data="$path/01.Raw_data"
trimmed="$path/02.Quality_trimmed"

## Define trimmomatic parameters:
num_runs=43                 # The number of runs to be processed
qual=10                     # The Q-score cut off
headcrop_fwd=30             # The number of bases to be trimmed from the 5' of the 'fwd' (R1) reads.
minlen_fwd=270              # Minimum length cut off
crop_rev=270                # Remove the available bases to remove from the distal end of the rev reads
headcrop_rev=28             # The number of bases to be trimmed from the 5' of the 'rev' (R2) read: This is being used because the primers were not found for most reads
minlen_rev=242              # Minimum length cut off
THREADS=$SLURM_CPUS_ON_NODE # Take advantage of all the available threads

###############################################################################
## Run trimmomatic and calculate run performance statistics ###################
###############################################################################

for run_number in $(seq 1 $num_runs); do
    echo
    echo "Run $run_number: Quality truncate reads"

    cd "$raw_data/run$run_number"

    for f in *R1.fastq.gz; do
        trimmomatic SE "$f" "$trimmed/run$run_number/fwd/$f" \
            HEADCROP:"$headcrop_fwd" SLIDINGWINDOW:4:"$qual" MINLEN:"$minlen_fwd" \
            -threads "$THREADS"
    done

    for f in *R2.fastq.gz; do
        trimmomatic SE "$f" "$trimmed/run$run_number/rev/$f" \
            CROP:"$crop_rev" HEADCROP:"$headcrop_rev" SLIDINGWINDOW:4:"$qual" MINLEN:"$minlen_rev" \
            -threads "$THREADS"
    done

    echo
    echo "Run $run_number: Generating quality report"

    ## Generate quality report for R1 (fwd) and R2 (rev) separately
    fastqc "$trimmed/run$run_number/fwd"/*R1.fastq.gz -o "$trimmed/run$run_number/fwd"
    fastqc "$trimmed/run$run_number/rev"/*R2.fastq.gz -o "$trimmed/run$run_number/rev"
    multiqc "$trimmed/run$run_number/fwd" -o "$trimmed/run$run_number/fwd"
    multiqc "$trimmed/run$run_number/rev" -o "$trimmed/run$run_number/rev"
    
    # Remove intermediate files and 'multiqc_data' directories
    rm "$trimmed/run$run_number"/fwd/*fastqc.zip "$trimmed/run$run_number"/fwd/*fastqc.html
    rm "$trimmed/run$run_number"/rev/*fastqc.zip "$trimmed/run$run_number"/rev/*fastqc.html
    rm "$trimmed/run$run_number"/fwd/multiqc_data/* "$trimmed/run$run_number"/rev/multiqc_data/*
    rmdir "$trimmed/run$run_number"/fwd/multiqc_data
    rmdir "$trimmed/run$run_number"/rev/multiqc_data
    
    # Calculate number and percentage of retained reads after trimming, along with additional statistics
    echo
    echo "Run $run_number: Generating trim statitics report"
    
    results_file_fwd="$trimmed/trim_summary_fwd.txt"
    results_file_rev="$trimmed/trim_summary_rev.txt"
    
    # Initialize the result files if they don't exist yet
    if [ ! -f "$results_file_fwd" ]; then
        echo -e "Run\tInput_Reads\tSurviving_Reads\tSurviving_Percent\tDropped_Reads\tDropped_Percent\tMean_Reads_Per_Sample\tMin_Reads_In_Sample\tMax_Reads_In_Sample\tSamples_Less_Than_10K\tSamples_10K_to_20K\tSamples_More_Than_20K" > "$results_file_fwd"
    fi
    
    if [ ! -f "$results_file_rev" ]; then
        echo -e "Run\tInput_Reads\tSurviving_Reads\tSurviving_Percent\tDropped_Reads\tDropped_Percent\tMean_Reads_Per_Sample\tMin_Reads_In_Sample\tMax_Reads_In_Sample\tSamples_Less_Than_10K\tSamples_10K_to_20K\tSamples_More_Than_20K" > "$results_file_rev"
    fi
    
    ### Calculate statistics for both R1 (fwd) and R2 (rev) files separately
    ### Note that I have included only test samples, i.e., those with a file name beginning with 's'
    
    # Calculate statistics for R1 (fwd) files
    input_reads_count_fwd=$(($(zcat "$raw_data/run$run_number/s"*R1.fastq.gz | wc -l) / 4)) # Divide by 4 because each read occupies four lines in a fastq file
    surviving_read_count_fwd=$(($(zcat "$trimmed/run$run_number/fwd/s"*R1.fastq.gz | wc -l) / 4))
    percentage_retained_fwd=$(awk "BEGIN { printf \"%.2f\", ($surviving_read_count_fwd / $input_reads_count_fwd) * 100 }")
    dropped_reads_count_fwd=$((input_reads_count_fwd - surviving_read_count_fwd))
    percentage_dropped_fwd=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads_count_fwd / $input_reads_count_fwd) * 100 }")
    
    # Calculate statistics for R2 (rev) files
    input_reads_count_rev=$(($(zcat "$raw_data/run$run_number/s"*R2.fastq.gz | wc -l) / 4)) # Divide by 4 because each read occupies four lines in a fastq file
    surviving_read_count_rev=$(($(zcat "$trimmed/run$run_number/rev/s"*R2.fastq.gz | wc -l) / 4))
    percentage_retained_rev=$(awk "BEGIN { printf \"%.2f\", ($surviving_read_count_rev / $input_reads_count_rev) * 100 }")
    dropped_reads_count_rev=$((input_reads_count_rev - surviving_read_count_rev))
    percentage_dropped_rev=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads_count_rev / $input_reads_count_rev) * 100 }")
    
    # Calculate mean reads per sample per run after trimming for R1 (fwd)
    num_samples_fwd=$(find "$trimmed/run$run_number/fwd" -name "s*R1.fastq.gz" | wc -l)
    mean_reads_per_sample_fwd=$((surviving_read_count_fwd / num_samples_fwd))

    # Calculate mean reads per sample per run after trimming for R2 (rev) files
    num_samples_rev=$(find "$trimmed/run$run_number/rev" -name "s*R2.fastq.gz" | wc -l)
    mean_reads_per_sample_rev=$((surviving_read_count_rev / num_samples_rev))

    # Find the minimum and maximum number of reads within a sample for R1 (fwd) and R2 (rev) 
    # Initialise variables to store the minimum and maximum read counts
    min_reads_in_sample_fwd=9999999999 # Start with a very large number so that the first sample will always be less than this number
    max_reads_in_sample_fwd=0
    min_reads_in_sample_rev=9999999999
    max_reads_in_sample_rev=0

    # Loop through each R1 (fwd) sample within individual runs
    for sample_file in "$trimmed/run$run_number/fwd/s"*R1.fastq.gz; do
        total_reads=$(( $(zcat "$sample_file" | wc -l) / 4 )) # Divide by 4 because each read occupies four lines in a fastq file
    
        if [ "$total_reads" -lt "$min_reads_in_sample_fwd" ]; then
            min_reads_in_sample_fwd="$total_reads"
        fi
    
        if [ "$total_reads" -gt "$max_reads_in_sample_fwd" ]; then
            max_reads_in_sample_fwd="$total_reads"
        fi
    done

    # Loop through each R2 (rev) sample within individual runs
    for sample_file in "$trimmed/run$run_number/rev/s"*R2.fastq.gz; do
        total_reads=$(( $(zcat "$sample_file" | wc -l) / 4 )) # Divide by 4 because each read occupies four lines in a fastq file
    
        if [ "$total_reads" -lt "$min_reads_in_sample_rev" ]; then
            min_reads_in_sample_rev="$total_reads"
        fi
    
        if [ "$total_reads" -gt "$max_reads_in_sample_rev" ]; then
            max_reads_in_sample_rev="$total_reads"
        fi
    done

    # Calculate the number of samples with <10,000 reads, 10,000-20,000 reads, and >20,000 reads within a run for R1 (fwd) and R2 (rev)
    # Initialize variables to count the samples in different read count ranges
    samples_less_than_10k_fwd=0
    samples_10K_to_20K_fwd=0
    samples_more_than_20k_fwd=0
    samples_less_than_10k_rev=0
    samples_10K_to_20K_rev=0
    samples_more_than_20k_rev=0

    # Loop through each R1 (fwd) sample within individual runs
    for sample_file in "$trimmed/run$run_number/fwd/s"*R1.fastq.gz; do
        total_reads=$(($(zcat "$sample_file" | wc -l) / 4)) # Divide by 4 because each read occupies four lines in a fastq file

        if [ "$total_reads" -lt 10000 ]; then
            ((samples_less_than_10k_fwd++))
        elif [ "$total_reads" -ge 10000 ] && [ "$total_reads" -lt 20000 ]; then
            ((samples_10K_to_20K_fwd++))
        else
            ((samples_more_than_20k_fwd++))
        fi
    done

    # Loop through each R2 (rev) sample within individual runs
    for sample_file in "$trimmed/run$run_number/rev/s"*R2.fastq.gz; do
        total_reads=$(($(zcat "$sample_file" | wc -l) / 4)) # Divide by 4 because each read occupies four lines in a fastq file

        if [ "$total_reads" -lt 10000 ]; then
            ((samples_less_than_10k_rev++))
        elif [ "$total_reads" -ge 10000 ] && [ "$total_reads" -lt 20000 ]; then
            ((samples_10K_to_20K_rev++))
        else
            ((samples_more_than_20k_rev++))
        fi
    done
    
    # Append the results to the respective files
    echo -e "$run_number\t$input_reads_count_fwd\t$surviving_read_count_fwd\t$percentage_retained_fwd\t$dropped_reads_count_fwd\t$percentage_dropped_fwd\t$mean_reads_per_sample_fwd\t$min_reads_in_sample_fwd\t$max_reads_in_sample_fwd\t$samples_less_than_10k_fwd\t$samples_10K_to_20K_fwd\t$samples_more_than_20k_fwd" >> "$results_file_fwd"
    
    echo -e "$run_number\t$input_reads_count_rev\t$surviving_read_count_rev\t$percentage_retained_rev\t$dropped_reads_count_rev\t$percentage_dropped_rev\t$mean_reads_per_sample_rev\t$min_reads_in_sample_rev\t$max_reads_in_sample_rev\t$samples_less_than_10k_rev\t$samples_10K_to_20K_rev\t$samples_more_than_20k_rev" >> "$results_file_rev"
done

echo "Ending at: $(date)"
