#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --mem=150G

printf "Starting at: $(date)"

## Activate the R environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda environment
conda activate R

## Run Rscript in a clean R instance, output a logfile
R_SCRIPT_PATH="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/02.Bioinformatics_16S/03.Quality_filter_denoise.R"
LOG_FILE="slurm-${SLURM_JOBID}.Rout"
OUTPUT_LOG="slurm-${SLURM_JOBID}.out"

Rscript --vanilla --verbose "$R_SCRIPT_PATH" > "$LOG_FILE" 2>&1 

## Append a logfile to this script's logfile
cat "$LOG_FILE" >> "$OUTPUT_LOG"

## Remove Rout log
rm "$LOG_FILE"

## Deactivate the R environment
conda deactivate

printf  "Finished at: $(date)"
