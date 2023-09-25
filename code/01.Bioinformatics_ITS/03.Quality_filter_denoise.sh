#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --mem-per-cpu=32G

echo "Starting at: $(date)"
echo

## Activate the R environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda environment
conda activate R

# Run Rscript in a clean R instance, output a logfile
path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/code/01.Bioinformatics_ITS/03.Quality_filter_denoise    ## Path to R script
Rscript --vanilla --verbose $path/Quality_filter_denoise.R > slurm-${SLURM_JOBID}.Rout 2>&1 

# Append a logfile to this scripts logfile
cat slurm-${SLURM_JOBID}.Rout >> slurm-${SLURM_JOBID}.out

# Remove Rout log
rm slurm-${SLURM_JOBID}.Rout

## Deactivate the R environment
conda deactivate

echo
echo "Finished at: $(date)"