#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=720:00:00
#SBATCH --partition=month
#SBATCH --mem-per-cpu=32G

echo "Starting at: $(date)"

## Organise directories
path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/ITS/reference_datasets/UNITE_general
UNITE=$path/sh_general_release_dynamic_all_eukaryotes_25.07.2023.fasta

## Set working directory
cd $path

## Activate the conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

## Extract the ITS regions from the UNITE dataset
ITSx \
  -i $UNITE \
  --complement T \
  --save_regions all \
  --graphical F \
  --positions T \
  -E 1e-1 \
  -t All \
  --cpu 8 \
  --preserve T \
  -o $path/

## Build UNITE ITS1 and ITS2 datasets for input into BLAST
makeblastdb \
   -in $path/ITS1.fasta \
   -out $path/ \
   -dbtype 'nucl' \
   -hash_index

makeblastdb \
   -in $path/ITS2.fasta \
   -out $path/ \
   -dbtype 'nucl' \
   -hash_index

## Deactivate the conda environmentput 
conda deactivate

echo "Finished at: $(date)"
