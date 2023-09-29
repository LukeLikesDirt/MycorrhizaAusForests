# Load the required packages
require(dada2, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/envs/R-packages')
require(seqinr, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/envs/R-packages')
require(tidyverse, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/envs/R-packages')

# Set the working directory to the project directory
setwd("/data/group/frankslab/project/LFlorence/MycorrhizaAusForests")

# Define the path to the directory containing the denoised sequence tables
denoised <- "data/AusMicrobiome/16S/04.Denoised"

# Define the number of sequencing runs that will be processed
num_runs <- 43

# Initialise an empty list to store the sequence tables
seqtab_list <- list()

# Loop through the sequencing runs and merge the sequence tables
for (i in 1:num_runs) {
  seqtab_file <- file.path(denoised, paste0('run', i, '_seqtab.rds'))
  seqtab <- readRDS(seqtab_file)
  seqtab_list[[i]] <- seqtab
}

# Merge all sequence tables
all.seqtab <- do.call(mergeSequenceTables, seqtab_list)

## Convert rds to fasta for chimera detection in VSEARCH

## Format data frame
fasta.tab <- as.data.frame(all.seqtab) %>%
    rownames_to_column('sample.names') %>%
    as_tibble() %>%
    pivot_longer(-sample.names, names_to = 'seq', values_to = 'size') %>%
    filter(size > 0) %>%
    group_by(sample.names) %>%
    mutate(seq.name = paste0(sample.names, '.fasta', '.', 1:n(), ';size=', size)) %>%
    ungroup() %>%
    select(seq.name, seq)

## Write the fasta file
write.fasta(as.list(fasta.tab$seq), fasta.tab$seq.name, file.path(denoised, 'all.fasta'), open = 'w', nbchar = 60, as.string = FALSE)
