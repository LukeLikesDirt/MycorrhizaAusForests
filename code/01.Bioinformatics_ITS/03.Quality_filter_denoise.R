
# Script: Quality filter and denoise Illumina single-end reads, and merge
# sequences from multiple sequencing runs.
# Purpose: Prepare Illumina forwards reads for chimera detection in VSEARCH.
# Credit: This script is adapted from https://benjjneb.github.io/dada2/bigdata.html
# Author: Luke Florence.
# Date: 28th September 2023.
#
# Script Purpose:
# ---------------
# This R script is designed for the quality filtering and denoising of
# Illumina forward reads using DADA2. The script is designed for "big data"
# projects, and therefore, performs quality filtering and denoising on
# sequencing runs induvidually. Sequencing runs are processed induvidually to
# improve error rate estimations, as each sequencing run has a unique error
# profile. The script also converts the DADA2 'rds' sequence table to a 'fasta'
# file formatted for VESEARCH for the subsequent chimera detection and removal
# in VSEARCH.
#
# Script Overview:
# ---------------
#   (1) Quality filter reads from each sequencing run induvidually.
#   (2) Denoise the reads from each sequencing run induvidually.
#   (3) Merges the denoised sequence tables.
#   (4) Convert the merged sequence table to a fasta file for chimera detection
#       and removal in VSEARCH.
#
# Pre-processing:
# ---------------
# If processing ITS sequences, ITS extraction with ITSxpress is recomended
# prior to quality filtering and denoising.
#
# I have quality truncated reads with Trimmomatic prior to quality filtering
# to improve overall read quality (maxEE). This is because removing the low
# quality reads on the distal (3') end removes low-quality bases, which can
# have a disproportionately high impact on the overall quality score. Therefore,
# by removing these low quality bases more reads will have the potential to
# pass the initial quality quality filtering step in DADA2.
#
# For SSU or LSU regions, reads can be trimmed to the same length in Trimmomatic
# or DADA2.

# Load required libraries
library(dada2, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/envs/R-packages')
library(seqinr, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/envs/R-packages')
library(tidyverse, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/envs/R-packages')

# Constant variables
num_runs <- 44   # The number of sequencing runs that will be processed
mEE <- 1        # maxEE: maximum number of expected errors per read
tL <- 0         # truncLen: fixed length for truncating reads
tQ <- 0         # truncQ: quality threshold for truncating reads
mNs <- 0        # maxN: maximum number of Ns allowed in a read

# Set the working directory to the project directory
setwd("/data/group/frankslab/project/LFlorence/MycorrhizaAusForests")

# Define the path to the ITS directory
path <- "data/AusMicrobiome/ITS"

# Create subdirectories for quality-filtered and denoised files
for (run in 1:num_runs) {
  run_dir <- file.path(path, '04.Quality_filtered', paste0('run', run))
  dir.create(run_dir, recursive = TRUE) # Create parent directories if they don't exist
}
dir.create(file.path(path, '05.Denoised'))

# Create an empty data frame to store the read tracking summary
summary_track <- data.frame()

# Quality filter each run individually
for (run in 1:num_runs) {
    
    ## Input and output directories
    trimmed_dir <- file.path(path, paste0('03.Quality_trimmed/run', run, '/'))
    qualFilt_dir <- file.path(path, paste0('04.Quality_filtered/run', run))
    
    ## List of files to be filtered
    fns <- list.files(trimmed_dir, pattern = 'fastq.gz')
    sample.names <- sapply(strsplit(basename(fns), "_"), `[`, 1)  ## Assumes filename = samplename_XXX.fastq.gz
    names(fns) <- sample.names
    
    ## Quality filtering
    qualFilt <- filterAndTrim(file.path(trimmed_dir, fns), file.path(qualFilt_dir, basename(fns)),
                              maxEE = mEE, truncLen = tL, truncQ = tQ, maxN = mNs,
                              rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE,
                              matchIDs = TRUE)
    
    ## Track reads through the quality filtering pipeline
    colnames(qualFilt) <- c("input_qualFilter", "output_qualFilter")
    rownames(qualFilt) <- sample.names
    
    ## Save track file to denoised directory
    write.csv(qualFilt, file = file.path(path, "04.Quality_filtered", paste0("run", run), "stats.csv"))
}

# Denoise each run individually
for (run in 1:num_runs) {

    ## Input and output directories
    qualFilt_dir <- file.path(path, "04.Quality_filtered", paste0("run", run))
    denoised_dir <- file.path(path, "05.Denoised")

    ## Denoising
    filts <- list.files(qualFilt_dir, pattern = 'fastq.gz', full.names = TRUE)
    sample.names <- sapply(strsplit(basename(filts), '_'), `[`, 1)
    names(filts) <- sample.names
    
    ## Learn error rates
    set.seed(1986)
    err <- learnErrors(filts, nbases = 1e8, multithread = TRUE, randomize = TRUE)
    
    ## Infer sequence variants
    dds <- vector('list', length(sample.names))
    names(dds) <- sample.names
    for (sam in sample.names) {
        cat('Processing:', sam, '\n')
        derep <- derepFastq(filts[[sam]])
        dds[[sam]] <- dada(derep, err = err, multithread = TRUE)
    }
    
    ## Build sequence table
    seqtab <- makeSequenceTable(dds)
    saveRDS(seqtab, file.path(denoised_dir, paste0(run, '_seqtab.rds')))

    ## Track reads through the DADA2 pipeline for the current run
    quality <- read.csv(file.path(qualFilt_dir, "stats.csv"))
    getN <- function(x) sum(getUniques(x))
    
    ## Check if there are rows in the "stats.csv" file
    if (nrow(quality) > 0) {
        track <- data.frame(
            run_number = run,
            input_qualFilter = quality$input_qualFilter,
            output_qualFilter = quality$output_qualFilter,
            output_denoised = sapply(dds, getN)
        )
        
        ## Append the track information for the current run to the summary_track data frame
        summary_track <- rbind(summary_track, track)
    }
}

# Track reads across the pipeline, summarised by sequencing run
summary_track %>%
    rename("Run" = run_number) %>%
    group_by(Run) %>%
    summarise(
        Input_Reads = sum(input_qualFilter),
        Output_Reads = sum(output_qualFilter),
        "Retained_Reads_%" = round((sum(output_denoised) / sum(input_qualFilter)) * 100, digits = 2),
        "Dropped_Reads_%" = round((sum(input_qualFilter) - sum(output_qualFilter)) / sum(input_qualFilter), digits = 2),
        Mean_Reads_Per_Sample = round(sum(output_qualFilter) / n(), digits = 0),
        Min_Reads = min(output_qualFilter),
        Max_Reads = max(output_qualFilter),
        Samples_LT_10K_Reads = sum(output_qualFilter < 10000),
        "Samples_10-20K_Reads" = sum(output_qualFilter >= 10000 & output_qualFilter < 20000),
        Samples_MT_20K_Reads = sum(output_qualFilter >= 20000)) %>%
        ungroup() %>%
        write.csv(file.path(path, "summary_denoised.csv"), row.names = FALSE)

# Merge all sequence tables and convert rds to fasta format for chimera detection in VSEARCH

# NOTE: Some samples have been re-sequenced across multiple sequencing runs.
# Therefore, I sum read counts and merge re-sequenced samples when executing the mergeSequenceTables()function.

# Initialise an empty list to store the sequence tables
seqtab_list <- list()

# Loop through the sequencing runs to read and merge the sequence tables
for (run in 1:num_runs) {
  seqtab_file <- file.path(path, paste0('05.Denoised/', run, '_seqtab.rds'))
  seqtab <- readRDS(seqtab_file)
  seqtab_list[[run]] <- seqtab
}

# Merge all sequence tables with repeats = "sum" to combine abundance across samples that have been sequenced in multiple runs
all.seqtab <- mergeSequenceTables(tables = seqtab_list, repeats = "sum")

## Convert rds to fasta for chimera detection in VSEARCH

## Format data frame
fasta.tab <- as.data.frame(all.seqtab) %>%
    rownames_to_column('sample.names') %>%
    as_tibble() %>%
    pivot_longer(-sample.names, names_to = 'seq', values_to = 'size') %>%
    filter(size > 0) %>%
    ungroup() %>%
    mutate(seq.name = paste0(sample.names, '.fasta', '.', row_number(), ';size=', size)) %>%
    select(seq.name, seq)

## Save the rds file for the merged sequence table
saveRDS(all.seqtab, file.path(path, '05.Denoised/all_seqtab.rds'))

## Write and save the fasta file for the merged sequence table
write.fasta(as.list(fasta.tab$seq), fasta.tab$seq.name, file.path(path, '05.Denoised/all.fasta'), open = 'w', nbchar = 60, as.string = FALSE)
