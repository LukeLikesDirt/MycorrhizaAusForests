#
# Creadit: This script is adapted from https://benjjneb.github.io/dada2/bigdata_paired.html
#
# Script Purpose:
# ---------------
# This R script is designed for the quality filtering and denoising of
# illumina paired-end amplicons using DADA2. The script is designed
# for "big data" projects, and therefore, performs quality filtering and
# denoising on sequencing runs individually. The sequencing runs are
# individually denoised to improve error rate estimations, as each sequencing
# run has a unique error profile. The script also converts the DADA2 'rds'
# sequence table to a 'fasta' file formatted for VESEARCH, to run the
# subsequent chimera detection and removal process in VSEARCH.
#
# Script Overview:
# ---------------
#   (1) Quality filter reads from each sequencing run individually.
#   (2) Denoise each sequencing run individually and merge the R1 and R2 reads.
#   (3) Merges the denoised sequence tables.
#   (4) Convert the merged sequence table to a fasta file for chimera detection
#       and removal in VSEARCH.

# Load required libraries
library(dada2)
library(seqinr)
library(here)
library(tidyverse)

# Define the path to the data subdirectory containing the 16S fastq files
path <- here("data/AusMicrobiome/16S")

# Define the number of sequencing runs that will be processed
num_runs <- 43

# Define the parameters for quality filtering
mEE <- c(2,3) # maxEE: maximum number of expected errors per read
tL <- 0       # truncLen: fixed length for truncating reads
tQ <- 0       # truncQ: quality threshold for truncating reads
mNs <- 0      # maxN: maximum number of Ns allowed in a read
# NOTE: tL and tQ are set to 0 to because the reads are already quality
# truncated and trimmed to equal lengths using Trimmomatic.

# Create subdirectories for quality filtered .fastq and denoised .rda files
for (run in 1:num_runs) {
  run_dir = file.path(path, '03.Quality_filtered', paste0('run', run))
  dir.create(run_dir, recursive = TRUE)
  dir.create(file.path(run_dir, 'fwd'))
  dir.create(file.path(run_dir, 'rev'))
}
dir.create(file.path(path, '04.Denoised_DADA2'))
dir.create(file.path(path, '06.Chimera_filtered_DADA2'))

# Create an empty data frame to store the read tracking summary
summary_track <- data.frame()

# Quality filter each run individually
for (run in 1:num_runs) {
  
  # Define the input and output directories
  trim.fwd <- file.path(path, "02.Quality_trimmed", paste0("run", run), "fwd")
  trim.rev <- file.path(path, "02.Quality_trimmed", paste0("run", run), "rev")
  qualFilt.fwd <- file.path(path, "03.Quality_filtered",
                            paste0("run", run), "fwd")
  qualFilt.rev <- file.path(path, "03.Quality_filtered",
                            paste0("run", run), "rev")

  # List file names for 'fwd' and 'rev' reads
  fns.fwd <- sort(list.files(trim.fwd, pattern = "R1.fastq.gz",
                             full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fns.fwd),
                         "_"), `[`, 1) ## Assumes samplename_XXX.fastq.gz
  names(fns.fwd) <- sample.names
  fns.rev <- sort(list.files(trim.rev, pattern = "R2.fastq.gz",
                             full.names = TRUE))
  sample.names.rev <- sapply(strsplit(basename(fns.rev), "_"), `[`, 1)
  names(fns.rev) <- sample.names.rev

  # Check that forward and reverse files match for the current run
  stopifnot(identical(sample.names, sample.names.rev))

  # Quality filtering
  qualFilt <- filterAndTrim(
    fwd = fns.fwd, filt = file.path(qualFilt.fwd, basename(fns.fwd)),
    rev = fns.rev, filt.rev = file.path(qualFilt.rev, basename(fns.rev)),
    maxEE = mEE, truncLen = tL, truncQ = tQ, maxN = mNs,
    rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE,
    matchIDs = TRUE
  )

  # Track reads through the quality filtering pipeline
  colnames(qualFilt) <- c("input_quality_filter", "output_quality_filter")
  rownames(qualFilt) <- sample.names
  # Save track file to denoised directory
  write.csv(
    qualFilt, file = file.path(path, "03.Quality_filtered",
                               paste0("run", run), "stats.csv")
  )

}

# Denoise each run individually
for (run in 1:num_runs) {

  # Define the input and output directories
  qualFilt_dir <- file.path(path, "03.Quality_filtered", paste0("run", run))
  qualFilt.fwd <- file.path(qualFilt_dir, "fwd")
  qualFilt.rev <- file.path(qualFilt_dir, "rev")
  denoised_dir <- file.path(path, "04.Denoised_DADA2")

  # Rename sample files to samplename.fastq.gz and list names for denoising
  filts.fwd <- list.files(qualFilt.fwd, pattern = "R1.fastq.gz",
                          full.names = TRUE)
  sample.names <- sapply(strsplit(basename(filts.fwd), "_"), `[`, 1)
  names(filts.fwd) <- sample.names
  filts.rev <- list.files(qualFilt.rev, pattern = "R2.fastq.gz",
                          full.names = TRUE)
  sample.names.rev <- sapply(strsplit(basename(filts.rev), "_"), `[`, 1)
  names(filts.rev) <- sample.names.rev

  # Check that forward and reverse files match for the current run
  stopifnot(identical(sample.names, sample.names.rev))

  # Learn error rates for the current run
  set.seed(1986)
  err.fwd <- learnErrors(filts.fwd, nbases = 1e8,
                         multithread = TRUE, randomize = TRUE)
  err.rev <- learnErrors(filts.rev, nbases = 1e8,
                         multithread = TRUE, randomize = TRUE)

  # Denoise and merge paired-end reads for the current run
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  for(sam in sample.names) {
    message("Processing:", sam)
    derep.fwd <- derepFastq(filts.fwd[[sam]])
    dd.fwd <- dada(derep.fwd, err = err.fwd, multithread = TRUE,
                   DETECT_SINGLETONS = TRUE)
    derep.rev <- derepFastq(filts.rev[[sam]])
    dd.rev <- dada(derep.rev, err = err.rev, multithread = TRUE,
                   DETECT_SINGLETONS = TRUE)
    merger <- mergePairs(dd.fwd, derep.fwd, dd.rev, derep.rev,
                         maxMismatch = 2, minOverlap = 12)
    mergers[[sam]] <- merger
  }
  rm(derep.fwd, derep.rev)

  # Build sequence tables

  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, file.path(denoised_dir, paste0(run, '_seqtab.rds')))

  # Track reads through the DADA2 pipeline
  # Read in quality filtered stats file
  quality <- read.csv(file.path(qualFilt_dir, "stats.csv")) %>%
    mutate("Run" = run) %>%
    rename("sample" = X) %>%
    # Filter samples starting with "s" to retain only test samples and ignore
    # control samples
    filter(substr(sample, 1, 1) == "s") %>%
    select(-sample) %>%
    as.data.frame()

  # Denoised stats
  getN <- function(x) sum(getUniques(x))
  if (nrow(quality) > 0) {
    track <- data.frame(
      sample = sample.names,
      output_denoised = sapply(mergers, getN)
    ) %>%
      filter(substr(sample, 1, 1) == "s") %>%
      select(-sample) %>%
      mutate(Run = run,
             input_quality_filter = quality$input_quality_filter,
             output_quality_filter = quality$output_quality_filter) %>%
      as.data.frame()

  # Append the merged information for the current run to the summary_track file
  summary_track <- rbind(summary_track, track)
  }
}

# Track reads across the pipeline, summarised by sequencing run
summary_track %>%
  group_by(Run) %>%
  summarise(
            Input_Reads = sum(input_quality_filter),
            Output_Reads = sum(output_quality_filter),
            "Retained_Reads_%" = round((sum(output_quality_filter) /
                                          sum(input_quality_filter)) * 100,
                                       digits = 2),
            "Dropped_Reads_%" = round(((sum(input_quality_filter) -
                                          sum(output_quality_filter)) /
                                          sum(input_quality_filter)) * 100,
                                      digits = 2),
            Mean_Reads_Per_Sample = round(sum(output_quality_filter) / n(),
                                          digits = 0),
            Min_Reads = min(output_quality_filter),
            Max_Reads = max(output_quality_filter),
            Samples_LT_10K_Reads = sum(output_quality_filter < 10000),
            "Samples_10-20K_Reads" = sum(output_quality_filter > 10000 &
                                           output_quality_filter < 20000),
            Samples_MT_20K_Reads = sum(output_quality_filter > 20000)) %>%
  ungroup() %>%
  write.table(file.path(path, "summary_quality_filtered.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)


summary_track %>%
  group_by(Run) %>%
  summarise(
            Input_Reads = sum(output_quality_filter),
            Output_Reads = sum(output_denoised),
            "Retained_Reads_%" = round((sum(output_denoised) /
                                         sum(output_quality_filter)) * 100,
                                       digits = 2),
            "Dropped_Reads_%" = round(((sum(output_quality_filter) -
                                        sum(output_denoised)) /
                                        sum(output_quality_filter)) * 100,
                                      digits = 2),
        Mean_Reads_Per_Sample = round(sum(output_denoised) / n(), digits = 0),
        Min_Reads = min(output_denoised),
        Max_Reads = max(output_denoised),
        Samples_LT_10K_Reads = sum(output_denoised < 10000),
        "Samples_10-20K_Reads" = sum(output_denoised > 10000 &
                                     output_denoised < 20000),
        Samples_MT_20K_Reads = sum(output_denoised > 20000)) %>%
  ungroup() %>%
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
  seqtab_file <- file.path(path, paste0('04.Denoised_DADA2/',
                                        run, '_seqtab.rds'))
  seqtab <- readRDS(seqtab_file)
  seqtab_list[[run]] <- seqtab
}

# The vector length when using 'pivot_longer' exceeds the maximum vector length
# that is allowed in R, so I will process the sequence tables in chunks.
# Divide the sequencing runs into five chunks
chunk1 <- seqtab_list[1:10]
chunk2 <- seqtab_list[11:20]
chunk3 <- seqtab_list[21:30]
chunk4 <- seqtab_list[31:40]
chunk5 <- seqtab_list[41:43]

# Merge each chunk of sequence tables separately
# Use repeats = "sum" to combine abundance across samples that have been
# sequenced in multiple sequencing runs. This should only be used when you are
# certain that the samples have been re-sequenced.
all.seqtab_chunk1 <- mergeSequenceTables(tables = chunk1, repeats = "sum")
all.seqtab_chunk2 <- mergeSequenceTables(tables = chunk2, repeats = "sum")
all.seqtab_chunk3 <- mergeSequenceTables(tables = chunk3, repeats = "sum")
all.seqtab_chunk4 <- mergeSequenceTables(tables = chunk4, repeats = "sum")
all.seqtab_chunk5 <- mergeSequenceTables(tables = chunk5, repeats = "sum")

# Merge the sequence tables from different chunks
all.seqtab <- mergeSequenceTables(tables = list(all.seqtab_chunk1,
                                                all.seqtab_chunk2,
                                                all.seqtab_chunk3,
                                                all.seqtab_chunk4,
                                                all.seqtab_chunk5),
                                  repeats = "sum")

# Initialise an empty list to store the formatted data frames for each chunk
fasta_tab_list <- list()

# Loop through each chunk to create fasta.tab
for (chunk in list(all.seqtab_chunk1, all.seqtab_chunk2, all.seqtab_chunk3,
                   all.seqtab_chunk4, all.seqtab_chunk5)) {
  fasta_chunk <- as.data.frame(chunk) %>%
  rownames_to_column('sample.names') %>%
  as_tibble() %>%
  pivot_longer(-sample.names, names_to = 'seq', values_to = 'size') %>%
  filter(size > 0) %>%
  ungroup() %>%
  mutate(seq.name = paste0(sample.names, '.fasta', '.',
                           row_number(), ';size=', size)) %>%
  select(seq.name, seq)
  
  fasta_tab_list[[length(fasta_tab_list) + 1]] <- fasta_chunk
}

# Combine the formatted data frames from different chunks
fasta.tab <- do.call(rbind, fasta_tab_list)

# Save the rds file for the merged sequence table
saveRDS(all.seqtab, file.path(path, '04.Denoised_DADA2/all_seqtab.rds'))

# Write and save the fasta file for the merged sequence table
write.fasta(as.list(fasta.tab$seq), fasta.tab$seq.name,
            file.path(path, '06.Chimera_filtered_DADA2/all.fasta'),
            open = 'w', nbchar = 60, as.string = FALSE)