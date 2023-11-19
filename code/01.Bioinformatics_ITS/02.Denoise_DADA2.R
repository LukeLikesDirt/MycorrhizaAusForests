
# Script:   Denoise Illumina single-end reads and merge sequence tables from
#           multiple sequencing runs.
# Credit:   https://benjjneb.github.io/dada2/bigdata.html
# Author:   Luke Florence.
# Date:     28th October 2023.
#
# Script Purpose:
# ---------------
# This R script is designed for the quality filtering and denoising of
# Illumina forward reads using DADA2. The script is designed for "big data"
# projects, and therefore, performs quality filtering and denoising on
# sequencing run individually. Sequencing runs are processed individually to
# improve error rate estimations, as each sequencing run has a unique error
# profile. The script also converts the DADA2 'rds' sequence table to a 'fasta'
# file formatted for VESEARCH for the subsequent chimera detection and removal
# in VSEARCH.
#
# Script Overview:
# ---------------
#   (1) Quality filter reads from each sequencing run individually.
#   (2) Denoise the reads from each sequencing run individually.
#   (3) Merges the denoised sequence tables.
#   (4) Convert the merged sequence table to a fasta file for chimera detection
#       and removal in VSEARCH.
#
# Note:
# -----
# Quality trimming and filtering have already been performed using
# Trimmomatic and VSEARCH.

# Load required libraries
library(dada2)
library(seqinr)
library(here)
library(tidyverse)

# The number of sequencing runs that will be processed
num_runs <- 44

# Define the path to the ITS data directory
path <- here("data/AusMicrobiome/ITS")

# Create subdirectories
dir.create(file.path(path, '05.Denoised_DADA2'))
dir.create(file.path(path, '07.Chimera_filtered_DADA2'))

# Create an empty data frame to store the read tracking summary
summary_track <- data.frame()

# Denoise each run individually
for (run in 1:num_runs) {

  ## Input and output directories
  qualFilt_dir <- file.path(path, "04.Quality_filtered",
                            paste0("run", run))
  denoised_dir <- file.path(path, "05.Denoised_DADA2")

  ## Denoising
  filts <- list.files(qualFilt_dir, pattern = 'fastq.gz', full.names = TRUE)
  sample.names <- sapply(strsplit(basename(filts), '_'), `[`, 1)
  names(filts) <- sample.names

  ## Learn error rates
  set.seed(1986)
  err <- learnErrors(filts, nbases = 1e8, multithread = TRUE,
                     randomize = TRUE)

  ## Infer sequence variants
  dds <- vector('list', length(sample.names))
  names(dds) <- sample.names
  for (sam in sample.names) {
    cat('Processing:', sam, '\n')
    derep <- derepFastq(filts[[sam]])
    dds[[sam]] <- dada(derep, err = err, multithread = TRUE,
                       DETECT_SINGLETONS = TRUE)
  }

  ## Build sequence table
  seqtab <- makeSequenceTable(dds)
  saveRDS(seqtab, file.path(denoised_dir, paste0(run, '_seqtab.rds')))

  ## Track reads through the DADA2 pipeline for the current run
  getN <- function(x) sum(getUniques(x))

  ## Check if there are rows in the "stats.csv" file
  track <- data.frame(Run = run,
                      sample = sample.names,
                      output = sapply(dds, getN)) %>%
    filter(substr(sample, 1, 1) == "s") %>%
    select(-sample) %>%
    as.data.frame()

  ## Append the track information for the current run to the summary_track file
  summary_track <- rbind(summary_track, track)

}

# Read in the summary of the quality filtered reads
quality <- read.table(file.path(path, "summary_quality_filtered.txt"),
                      header = TRUE, sep = "\t") %>%
  select(Run, input = Output_Reads) %>%
  as.data.frame()

inner_join(summary_track, quality, by = "Run") %>%
  group_by(Run) %>%
  summarise(
    Input_Reads = unique(input),
    Output_Reads = sum(output),
    "Retained_Reads_%" = round((sum(output) /
                                  sum(unique(input))) * 100, digits = 2),
    "Dropped_Reads_%" = round(((sum(unique(input)) - sum(output)) /
                                 sum(unique(input))) * 100, digits = 2),
    Mean_Reads_Per_Sample = round(sum(output) / n(), digits = 0),
    Min_Reads = min(output),
    Max_Reads = max(output),
    Samples_LT_10K_Reads = sum(output < 10000),
    "Samples_10-20K_Reads" = sum(output > 10000 & output < 20000),
    Samples_MT_20K_Reads = sum(output > 20000)
  ) %>%
  ungroup() %>%
  as.data.frame() %>%
  write.table(file.path(path, "summary_denoised_DADA2.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

# Merge all sequence tables and convert rds to fasta format for chimera
# detection in VSEARCH

# NOTE: Some samples have been re-sequenced across multiple sequencing runs.
# Therefore, I sum read counts and merge re-sequenced samples when executing the
# mergeSequenceTables()function.

# Initialise an empty list to store the sequence tables
seqtab_list <- list()

# Loop through the sequencing runs to read and merge the sequence tables
for (run in 1:num_runs) {
  seqtab_file <- file.path(path, paste0('05.Denoised_DADA2/',
                                        run, '_seqtab.rds'))
  seqtab <- readRDS(seqtab_file)
  seqtab_list[[run]] <- seqtab
}

# Merge all sequence tables with repeats = "sum" to combine abundance across
# samples that have been sequenced in multiple runs
all.seqtab <- mergeSequenceTables(tables = seqtab_list, repeats = "sum")

## Convert rds to fasta for chimera detection in VSEARCH

## Format data frame
fasta.tab <- as.data.frame(all.seqtab) %>%
  rownames_to_column('sample.names') %>%
  as_tibble() %>%
  pivot_longer(-sample.names, names_to = 'seq', values_to = 'size') %>%
  filter(size > 0) %>%
  ungroup() %>%
  mutate(seq.name = paste0(sample.names, '.fasta', '.', row_number(),
                           ';size=', size)) %>%
  select(seq.name, seq)

## Save the rds file for the merged sequence table
saveRDS(all.seqtab, file.path(path, '05.Denoised_DADA2/all_seqtab.rds'))

## Write and save the fasta file for the merged sequence table
write.fasta(as.list(fasta.tab$seq), fasta.tab$seq.name,
            file.path(path, '07.Chimera_filtered_DADA2/all.fasta'),
                      open = 'w', nbchar = 60, as.string = FALSE)
