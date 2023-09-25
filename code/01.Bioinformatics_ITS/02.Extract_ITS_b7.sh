#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1440:00:00
#SBATCH --partition=long
#SBATCH --mem-per-cpu=32G

echo "Starting at: $(date)"

## Prepare Illumina forward reads for denoising with DADA2 by extracting the ITS1 subregion using ITSxpress and quality truncating reads using Trimmomatic
## Reads cannot be merged due to the high error rates: error rates are routine for Australian Microbiome data, hence only the forward reads are used
## Primers have not been removed, which is also due to the high error rates

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda
conda activate shell

## Take advantage of all the available threads
OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

## Organize subdirectories
path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/ITS  # Path to the main 'ITS' data directory
mkdir $path/02.ITS_extracted
mkdir $path/03.Quality_trimmed
raw_data=$path/01.Raw_data          # Path to raw data
ITSx=$path/02.ITS_extracted         # Path to the folder for the extracted ITS reads
trimmed=$path/03.Quality_trimmed    # Path to the folder for the quality truncated reads
## Create subdirectories for each run
for dir in {02.ITS_extracted,03.Quality_trimmed}; do
  for run in {1..4}; do
    mkdir -p "$path/$dir/run$run"
  done
done

## Extract ITS1 and quality truncate reads
## Define the number of runs
num_runs=7
## Loop through each run
for run_number in $(seq 1 4); do
    echo
    echo "Run $run_number: Extract ITS1 and quality truncate reads"

    ## Extract ITS1
    cd $raw_data/run$run_number

    for f in *_R1.fastq.gz; do
        itsxpress \
            --single_end \
            --fastq ${f} \
            --cluster_id 1.0 \
            --region ITS1 \
            --taxa All \
            --log $ITSx/run${run_number}/logfile.txt \
            --outfile $ITSx/run${run_number}/${f} \
            --threads $OMP_NUM_THREADS
    done

    ## Quality trim reads
    cd $ITSx/run$run_number

    for f in *.fastq.gz; do
        trimmomatic SE $f $trimmed/run$run_number/$f \
            SLIDINGWINDOW:4:13 MINLEN:80 \
            -threads OMP_NUM_THREADS
    done

    ## Calculate statistics for retained reads after trimming
    echo
    echo "Run $run_number: Calculating trimming statistics"

    input_dir="$raw_data/run$run_number"
    trimmed_dir="$trimmed/run$run_number"

    results_file="$trimmed_dir/trim_summary.txt"

    # Initialize the result file if it doesn't exist yet
    if [ ! -f "$results_file" ]; then
        echo -e "Run\tInput_Reads\tSurviving_Reads\tSurviving_Percent\tDropped_Reads\tDropped_Percent\tMean_Reads_Per_Sample\tMin_Reads_In_Sample\tMax_Reads_In_Sample\tSamples_Less_Than_10K\tSamples_10K_to_20K\tSamples_More_Than_20K" > "$results_file"
    fi

    # Calculate statistics for the trimmed reads
    input_reads_count=$(($(zcat "$input_dir"/*_R1.fastq.gz | wc -l) / 4)) # Divide by 4 because each read occupies four lines in a fastq file
    surviving_read_count=$(($(zcat "$trimmed_dir"/*.fastq.gz | wc -l) / 4))
    percentage_retained=$(awk "BEGIN { printf \"%.2f\", ($surviving_read_count / $input_reads_count) * 100 }")
    dropped_reads_count=$((input_reads_count - surviving_read_count))
    percentage_dropped=$(awk "BEGIN { printf \"%.2f\", ($dropped_reads_count / $input_reads_count) * 100 }")

    # Calculate mean reads per sample
    num_samples=$(find "$trimmed_dir" -name "*.fastq.gz" | wc -l)
    mean_reads_per_sample=$((surviving_read_count / num_samples))

    # Find the minimum and maximum number of reads within a sample
    min_reads_in_sample=9999999999 # Start with a very large number so that the first sample will always be less than this number
    max_reads_in_sample=0

    # Loop through each sample within the run
    for sample_file in "$trimmed_dir"/*.fastq.gz; do
        total_reads=$(( $(zcat "$sample_file" | wc -l) / 4 )) # Divide by 4 because each read occupies four lines in a fastq file

        if [ "$total_reads" -lt "$min_reads_in_sample" ]; then
            min_reads_in_sample="$total_reads"
        fi

        if [ "$total_reads" -gt "$max_reads_in_sample" ]; then
            max_reads_in_sample="$total_reads"
        fi
    done

    # Calculate the number of samples with <10,000 reads, 10,000-20,000 reads, and >20,000 reads within the run
    # Initialize variables to count the samples in different read count ranges
    samples_less_than_10k=0
    samples_10K_to_20K=0
    samples_more_than_20k=0

    # Loop through each sample within the run
    for sample_file in "$trimmed_dir"/*.fastq.gz; do
        total_reads=$(($(zcat "$sample_file" | wc -l) / 4)) # Divide by 4 because each read occupies four lines in a fastq file

        if [ "$total_reads" -lt 10000 ]; then
            ((samples_less_than_10k++))
        elif [ "$total_reads" -ge 10000 ] && [ "$total_reads" -lt 20000 ]; then
            ((samples_10K_to_20K++))
        else
            ((samples_more_than_20k++))
        fi
    done

    # Append the results to the result file
    echo -e "$run_number\t$input_reads_count\t$surviving_read_count\t$percentage_retained\t$dropped_reads_count\t$percentage_dropped\t$mean_reads_per_sample\t$min_reads_in_sample\t$max_reads_in_sample\t$samples_less_than_10k\t$samples_10K_to_20K\t$samples_more_than_20k" >> "$results_file"

done

conda deactivate
echo
echo "Finished at: $(date)"
