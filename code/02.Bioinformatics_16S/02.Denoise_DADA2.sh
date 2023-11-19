#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=720:00:00
#SBATCH --partition=month
#SBATCH --mem=100G
#SBATCH --output=../02.Bioinformatics_16S/slurm/%x.%j.out

## The R script executed here is adapted from: https://benjjneb.github.io/dada2/bigdata_paired.html
## See the R script for more details.

printf "\nStarting at: %s\n" "$(date)"

## Activate the R environment
conda activate R

# Run R script in a clean R instance, output a logfile
R_SCRIPT_PATH="../02.Bioinformatics_16S/02.Denoise_DADA2.R"
LOG_FILE="../02.Bioinformatics_16S/slurm/02.Denoise_DADA2.sh.${SLURM_JOBID}.Rout"
OUTPUT_LOG="../02.Bioinformatics_16S/slurm/02.Denoise_DADA2.sh.${SLURM_JOBID}.out"

Rscript --vanilla --verbose "$R_SCRIPT_PATH" > "$LOG_FILE" 2>&1 

# Append Rout log to the slurm out log
cat "$LOG_FILE" >> "$OUTPUT_LOG"

# Remove Rout log
rm "$LOG_FILE"

## Deactivate the R environment
conda deactivate

printf "\nFinished at: %s\n" "$(date)"
