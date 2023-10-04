#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --partition=week
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/01.Bioinformatics_ITS/slurm/%x.%j.out

## Constants and file paths
readonly THREADS=8                                                                      ## Set the number of threads
readonly IDENTITY=0.97                                                                  ## Set the identity threshold for clustering
readonly path="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"            ## Path to the project directory
readonly chimeraFiltered_dir="$path/data/AusMicrobiome/ITS/test/07.Chimera_filtered"    ## Path to chimera filtered fasta file
readonly clustered_dir="$path/data/AusMicrobiome/ITS/test/08.Clustered"                 ## Path to clustered fasta file and OTU table
readonly output_dir="$path/output/"                                                     ## Path to formatted OTU table
mkdir -p "$clustered_dir"                                                               ## Make the clustered subdirectory
mkdir -p "$output_dir"                                                                  ## Make the clustered subdirectory

## Log function
log() {
    local timestamp
    timestamp=$(date)
    printf "\n%s %s\n\n" "$1" "$timestamp"
}

## Function for clustering at 97% and formatting the OTU table
cluster_and_format_otu_table() {

    ## Cluster at 97% and generate OTU table
    log 'Clustering at 97% at:'

    vsearch \
        --cluster_size "$chimeraFiltered_dir/all.nonchimeras.derep.fasta" \
        --threads "$THREADS" \
        --id "$IDENTITY" \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$chimeraFiltered_dir/all.clustered.uc" \
        --relabel_sha \
        --centroids "$clustered_dir/OTUs.fasta" \
        --otutabout "$clustered_dir/OTUs.txt"

    ## Format OTU table
    ## Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "$clustered_dir/OTUs.txt"
    ## Convert to .csv and save to the output directory
    sed -e 's/\s\+/,/g' "$clustered_dir/OTUs.txt" > "$output_dir/OTUs.csv"

    printf '\nNumber of OTUs and reads:\n'
    printf '    Number of OTUs: %s\n' "$(grep -c "^>" "$clustered_dir/OTUs.fasta")"
    printf '    Number of reads: %s\n' "$(grep -c "^>" "$chimeraFiltered_dir/all.nonchimeras.derep.fasta")"

}

## Execute the cluster_and_format_otu_table function

log 'Starting at:'

source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

cluster_and_format_otu_table

conda deactivate

log 'Finished at:'