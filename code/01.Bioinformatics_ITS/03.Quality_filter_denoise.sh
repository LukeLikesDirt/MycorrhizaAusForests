#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --mem=100G
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/%x.%j.out

## The associated R script used here is adapted from: https://benjjneb.github.io/dada2/bigdata.html
## See the R script for details on script purpose.

printf "\nStarting at: %s\n" "$(date)"

## Activate the R environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate R

# Run Rscript in a clean R instance, output a logfile
R_SCRIPT_PATH="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/03.Quality_filter_denoise.R"
LOG_FILE="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/03.Quality_filter_denoise.sh.${SLURM_JOBID}.Rout"
OUTPUT_LOG="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/03.Quality_filter_denoise.sh.${SLURM_JOBID}.out"

Rscript --vanilla --verbose "$R_SCRIPT_PATH" > "$LOG_FILE" 2>&1 

# Append a logfile to this script's logfile
cat "$LOG_FILE" >> "$OUTPUT_LOG"

# Remove Rout log
rm "$LOG_FILE"

## Deactivate the R environment
conda deactivate

printf "\nFinished at: %s\n" "$(date)"