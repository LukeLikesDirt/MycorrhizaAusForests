#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=720:00:00
#SBATCH --partition=month
#SBATCH --mem-per-cpu=32G
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/%x.%j.out

printf '\nStarting at: %s\n' "$(date)"

## FILEPATH: /data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/02.Extract_ITS_b11.sh

## Description: Prepare Illumina forward reads for denoising with DADA2 by
## extracting the ITS1 subregion using ITSxpress and quality truncating reads
## using Trimmomatic.
##
## The sequenceing data processed here span 44 runs. Because denoising with
## DADA2 is improved when estamating error rates on induvidual runs these data
## are subset by sequencing run. This script will loop through each sequencing
## run and perform the following tasks: (1) extract the ITS1 subregion with
## ITSxpress, (2) quality truncate reads with Trimmomatic, (3) quality check
## trimmed reads with fastQC, and (4) calculate statistics for tracking read
## survival across the pipeline.

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

## First and last sequencing runs that will be looped through in this batch
readonly first_run=41
readonly last_run=44

## Number of threads to use
readonly THREADS=8

## Organise subdirectories
readonly data_dir="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/data/AusMicrobiome/ITS"  ## Path to the main 'ITS' data directory
mkdir -p "$data_dir/02.ITS_extracted"
mkdir -p "$data_dir/03.Quality_trimmed"
readonly raw_data="$data_dir/01.Raw_data"               ## Path to raw data
readonly ITS_extracted="$data_dir/02.ITS_extracted"     ## Path to the folder for the extracted ITS reads
readonly quality_trimmed="$data_dir/03.Quality_trimmed" ## Path to the folder for the quality truncated reads

## Create subdirectories for each run
for dir in "$ITS_extracted" "$quality_trimmed"; do
  for run in $(seq "$first_run" "$last_run"); do
    mkdir -p "$dir/run$run"
  done
done

## Define the summary file for tracking reads from raw to trimmed. This file
## will be grouped and by run and appended to after each loop
readonly results_file="$data_dir/trim_summary.txt"

## Initialise the trimmed reads result file if it doesn't exist yet
if [ ! -f "$results_file" ]; then
    printf "Run\tInput_Reads\tOutput_Reads\tRetained_Percent\tDropped_Percent\tMean_Reads_Per_Sample\tMin_Reads_In_Sample\tMax_Reads_In_Sample\tSamples_Less_Than_10K\tSamples_10K_to_20K\tSamples_More_Than_20K\n" > "$results_file"
fi

## Loop through each run: (1) extract the ITS1 with ITSxpress, (2) quality
## truncate reads with trimmomatic, (3) qulatity check trimmed reads with 
## fastQC, and (4) calculate read tracking statistics for each run, as follows:
## input reads, reads retained, percentage retained, percentage dropped, mean 
## reads per sample, minimum reads in a sample, maximum reads in a sample,
## number of samples with <10,000 reads, number of samples with 10,000-20,000 
## and number of samples with >20,000 reads.

for run_number in $(seq "$first_run" "$last_run"); do

    ## (1) Extract ITS1
    printf '\nRun %d: Extracting ITS1 with ITSxpress\n' "$run_number"

    cd "$raw_data/run$run_number"

    for f in *_R1.fastq.gz; do
        itsxpress \
            --single_end \
            --fastq "${f}" \
            --cluster_id 1.0 \
            --region ITS1 \
            --taxa All \
            --log "$ITS_extracted/run${run_number}/logfile.txt" \
            --outfile "$ITS_extracted/run${run_number}/${f}" \
            --threads "$THREADS"
    done

    ## (2) Quality trim reads
    printf '\nRun %d: Quality trimming with Trimmomatic\n' "$run_number"

    cd "$ITS_extracted/run$run_number"

    for f in *.fastq.gz; do
        trimmomatic SE -threads "$THREADS" "$f" "$quality_trimmed/run$run_number/$f" \
            SLIDINGWINDOW:4:13 MINLEN:80 
    done

    printf '\nRun %d: Generating quality report\n' "$run_number"

    ## (2) Generate read quality report
    fastqc "$quality_trimmed/run$run_number"/*R1.fastq.gz -o "$quality_trimmed/run$run_number"
    multiqc "$quality_trimmed/run$run_number" -o "$quality_trimmed/run$run_number"
    
    ## Remove intermediate files and 'multiqc_data' directory
    rm "$quality_trimmed/run$run_number"/*fastqc.zip "$quality_trimmed/run$run_number"/*fastqc.html
    rm "$quality_trimmed/run$run_number"/multiqc_data/*
    rmdir "$quality_trimmed/run$run_number"/multiqc_data

    ## (4) Calculate statistics for tracking reads:
    ##  - Divide reads by 4 because each read occupies four lines in a fastq file
    ##    When counting reads across a fasrtq file, I need to divide by 4 because
    ##    each read occupies four lines in a fastq file. 
    ##  - I have appended "s" to the start of each file seach to only include
    ##    test samples and exclude the negative controls, which have very few
    ##    reads and will skew the results.
    printf '\nRun %d: Calculating trimming statistics\n' "$run_number"

    input_dir="$raw_data/run$run_number"
    output_dir="$quality_trimmed/run$run_number"

    ## Calculate the number of input, output (retained) and dropped during trimming
    input_reads=$(($(zcat "$input_dir/s"*.fastq.gz| wc -l) / 4))
    output_reads=$(($(zcat "$output_dir/s"*.fastq.gz| wc -l) / 4))
    dropped_reads=$((input_reads - output_reads))

    ## Calculate the percentage of reads retained and dropped during trimming
    percentage_retained=$(awk "BEGIN { printf \"%.2f\", ($output_reads / $input_reads) * 100 }")
    percentage_dropped=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads / $input_reads) * 100 }")
    
    ## Calculate mean reads per sample per run after trimming
    num_samples=$(find "$output_dir/" -name "s*.fastq.gz" | wc -l)
    mean_reads_per_sample=$((output_reads / num_samples))

    ## Find the minimum and maximum number of reads within a sample
    ## For "min_reads_in_sample", start with a very large number to ensure that the first sample has fewer reads
    min_reads_in_sample=999999999
    max_reads_in_sample=0

    ## Divide reads count in "$sample_file" by 4 because each read occupies four lines in a fastq file
    for sample_file in "$output_dir/s"*.fastq.gz; do
        # Remember to divide by 4 because each read occupies four lines in a fastq file
        total_reads=$(( $(zcat "$sample_file" | wc -l) / 4 ))
    
        if [ "$total_reads" -lt "$min_reads_in_sample" ]; then
            min_reads_in_sample="$total_reads"
        fi
    
        if [ "$total_reads" -gt "$max_reads_in_sample" ]; then
            max_reads_in_sample="$total_reads"
        fi
    done

    ## Calculate the number of samples with <10,000 reads, 10,000-20,000 reads, and >20,000 reads within the run
    ## Initialize variables to count the samples in different read count ranges
    samples_less_than_10k=0
    samples_10K_to_20K=0
    samples_more_than_20k=0

    ## Loop through each sample within the run
    for sample_file in "$output_dir/s"*.fastq.gz; do
        # Remember to divide by 4 because each read occupies four lines in a fastq file
        total_reads=$(($(zcat "$sample_file" | wc -l) / 4))

        if [ "$total_reads" -lt 10000 ]; then
            ((samples_less_than_10k++))
        elif [ "$total_reads" -ge 10000 ] && [ "$total_reads" -lt 20000 ]; then
            ((samples_10K_to_20K++))
        else
            ((samples_more_than_20k++))
        fi
    done

    ## Append the results to the result file
    printf "$run_number\t$input_reads\t$output_reads\t$percentage_retained\t$percentage_dropped\t$mean_reads_per_sample\t$min_reads_in_sample\t$max_reads_in_sample\t$samples_less_than_10k\t$samples_10K_to_20K\t$samples_more_than_20k\n" >> "$results_file"

done   

conda deactivate
printf '\nFinished at: %s\n' "$(date)"