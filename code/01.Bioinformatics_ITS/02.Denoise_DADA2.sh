#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --output=../01.Bioinformatics_ITS/slurm/%x.%j.out

## The R script executed here is adapted from: https://benjjneb.github.io/dada2/bigdata.html
## See the R script for more details

printf "\nStarting at: %s\n" "$(date)"

## Activate the R environment
conda activate R

# Run R script in a clean R instance, output a logfile
R_SCRIPT_PATH="../01.Bioinformatics_ITS/02.Denoise_DADA2.R"
LOG_FILE="../01.Bioinformatics_ITS/slurm/02.Denoise_DADA2.sh.${SLURM_JOBID}.Rout"
OUTPUT_LOG="../01.Bioinformatics_ITS/slurm/02.Denoise_DADA2.sh.${SLURM_JOBID}.out"

Rscript --vanilla --verbose "$R_SCRIPT_PATH" > "$LOG_FILE" 2>&1 

# Append Rout log to the slurm out log
cat "$LOG_FILE" >> "$OUTPUT_LOG"

# Remove Rout log
rm "$LOG_FILE"

## Deactivate the R environment
conda deactivate

printf "\nFinished at: %s\n" "$(date)"
