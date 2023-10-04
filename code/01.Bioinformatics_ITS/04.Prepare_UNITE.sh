#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --mem-per-cpu=32G
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/%x.%j.out

# Script: code/01.Bioinformatics_ITS/04.Prepare_UNITE.sh
# Purpose: Prepare the UNITE references dataset for chimera detection and taxonomic assignment for ITS1 and ITS2 amplicons
# Author: Luke Florence
# Date: 3rd October 2023

## Constants and file paths
readonly path="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"  ## The path to the project directory
readonly UNITE_dir="$path/data/AusMicrobiome/ITS/06.Reference_dataset"        ## Path to the UNITE fasta file
readonly ITS1_dir="$UNITE_dir/ITS1"                                           ## Path to the ITS1 fasta file
readonly ITS2_dir="$UNITE_dir/ITS2"                                           ## Path to the ITS2 fasta file
mkdir -p "$UNITE_dir"                                                         ## UNITE subdirectory
mkdir -p "$ITS1_dir"                                                          ## ITS1 subdirectory
mkdir -p "$ITS2_dir"                                                          ## ITS2 subdirectory

## UNITE dataset information: 
## CHANGE ME to the URL and file name for the desired UNITE version: https://unite.ut.ee/repository.php
## The URL to the latest release (25.07.2023) of UNITE general release in dynamic files for all eukaryotes with singletons set as RefS
UNITE_URL="https://files.plutof.ut.ee/public/orig/3D/FF/3DFF70D8FC94F5C331FA2AE7FCEDCC15181C70563641E9F4E8C3DC7869987C4D.tgz"
## UNITE dataset file names: 
##    - The stable version ends in '.fasta' at the development version ends in '_dev.fasta'.
##    - The 'prepare_UNITE_dataset' function uses the stable version by default.
##    - Amend "$UNITE_version" from '.fatsa' to '_dev.fasta' if you want to use the development version.
UNITE_dataset="$UNITE_dir/sh_general_release_dynamic_all_25.07.2023"
UNITE_version=".fasta"
UNITE_reformatted_dataset="$UNITE_dir/UNITE_reformatted.fasta"

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## Function to prepare the UNITE dataset for ITS extraction
prepare_UNITE_dataset() {

    ## Retrieve UNITE
    log 'Downloading UNITE at:'
    wget -O "$UNITE_dir/UNITE_compressed.tgz" "$UNITE_URL"

    ## Decompress UNITE
    log 'Decompressing UNITE at:'
    tar -zxvf "$UNITE_dir/UNITE_compressed.tgz" -C "$UNITE_dir/"

    ## Reformat to remove lowercase and spaces
    log 'Reformatting UNITE at:'
    awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' \
       "$UNITE_dataset$UNITE_version" | tr -d ' ' > \
       "$UNITE_reformatted_dataset"
    # Remove the intermediate files
    rm "$UNITE_dir/UNITE_compressed.tgz"
    rm ""$UNITE_dataset"_dev.fasta"
    rm ""$UNITE_dataset".fasta"
    # Rename the reformatted UNITE file
    mv "$UNITE_reformatted_dataset" "$UNITE_dataset$UNITE_version"

    ## Extract the ITS1 and ITS2 subregions from the UNITE dataset
    log 'Extracting the ITS region from UNITE at:'

    ITSx \
      -i "$UNITE_dataset$UNITE_version" \
      --complement T \
      --save_regions all \
      --graphical F \
      --positions T \
      -E 1e-1 \
      -t All \
      --cpu 8 \
      --preserve T \
      -o "$UNITE_dir/ITSx"

    ## Build UNITE ITS1 and ITS2 datasets for input into BLAST
    log 'Building blast databases for the ITS1 and ITS2 subregions at:'

    makeblastdb \
       -in "$UNITE_dir/ITSx.ITS1.fasta" \
       -out "$ITS1_dir/" \
       -dbtype 'nucl' \
       -hash_index

    makeblastdb \
       -in "$ITS2_dir/ITSx.ITS2.fasta" \
       -out "$ITS2_dir/" \
       -dbtype 'nucl' \
       -hash_index

}

## Execute the prepare_UNITE_dataset function

log 'Starting at:'

## Activate the conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell
## Execute the function
prepare_UNITE_dataset
## Deactivate the conda environment
conda deactivate

log 'Finished at:'
