#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1440:00:00
#SBATCH --partition=long
#SBATCH --mem=150G
#SBATCH --output=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/code/02.Bioinformatics_16S/slurm/%x.%j.out

## Taxonomic assignment with BLAST

readonly THREADS=8
readonly PROJECT_PATH="/data/group/frankslab/project/LFlorence/MycorrhizaAusForests"
readonly REFERENCE_SEQUENCES="$PROJECT_PATH/data/AusMicrobiome/16S/05.Reference_dataset/V1V3/SILVA_V1V3"
readonly OTU_FASTA="$PROJECT_PATH/data/AusMicrobiome/16S/06.Chimera_filtered/OTUs.fasta"
readonly TAXA_DIR="$PROJECT_PATH/data/AusMicrobiome/16S/07.Taxonomy"
mkdir -p "$TAXA_DIR"

## BLAST function

my_blast() {

    log 'BLAST best hit starting at'

    # Best hit
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES" \
        -strand both \
        -evalue 0.001 \
        -word_size 7 \
        -reward 1 \
        -penalty -1 \
        -gapopen 1 \
        -gapextend 2 \
        -max_target_seqs 1 \
        -max_hsps 1 \
        -out "$TAXA_DIR/BLAST_best_hit.txt" \
        -num_threads "$THREADS"

    log 'BLAST best ten hits starting at'

    # Best 10
    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        -query "$OTU_FASTA" \
        -db "$REFERENCE_SEQUENCES" \
        -strand both \
        -evalue 0.001 \
        -word_size 7 \
        -reward 1 \
        -penalty -1 \
        -gapopen 1 \
        -gapextend 2 \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -out "$TAXA_DIR/BLAST_best_10.txt" \
        -num_threads "$THREADS"

}

## Reformat taxa table function

reformat_taxa_tables() {

    sed -i '1s/^/OTU_ID;abundance\treference;kingdom;phylum;class;order;family;genus;species\tpident\tlength\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n/' "$TAXA_DIR/BLAST_best_hit.txt"
    sed 's/[[:space:]]\{1,\}/;/g' "$TAXA_DIR/BLAST_best_hit.txt" > "$TAXA_DIR/BLAST_best_hit.csv"

    sed -i '1s/^/OTU_ID;abundance\treference;kingdom;phylum;class;order;family;genus;species\tpident\tlength\tslen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n/' "$TAXA_DIR/BLAST_best_10.txt"
    sed 's/[[:space:]]\{1,\}/;/g' "$TAXA_DIR/BLAST_best_10.txt" > "$TAXA_DIR/BLAST_best_10.csv"

}

## Main script

printf '\nStarting at: %s\n' "$(date)"

log 'Starting at:'

source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh
conda activate shell

my_blast
reformat_taxa_tables

conda deactivate

printf '\nFinished at: %s\n' "$(date)"