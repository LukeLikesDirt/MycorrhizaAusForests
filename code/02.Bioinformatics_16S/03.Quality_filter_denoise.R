# This script performs quality filtering, denoising, and sequence table construction for paired-end 16S rRNA gene amplicon data using DADA2.

# Load DADA2
require(dada2, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/envs/R-packages')

# Set the working directory to the project directory
setwd("/data/group/frankslab/project/LFlorence/MycorrhizaAusForests")

# Define the path to the directory containing the runs
path <- "data/AusMicrobiome/16S"

# Define the number of sequencing runs thta will be processed
num_runs <- 43

# Define the parameters for quality filtering
mEE <- c(2,3) # maxEE: maximum number of expected errors per read
tL <- 0       # truncLen: fixed length for truncating reads
tQ <- 0       # truncQ: quality threshold for truncating reads
mNs <- 0      # maxN: maximum number of Ns allowed in a read
# NOTE: tL and tQ are set to 0 to because the reads are already quality truncated then trimmed to the same length with Trimmomatic

# Create subdirectories for quality filtered 'fwd' and 'rev' fastq files
for (i in 1:num_runs) {
  run_dir = file.path(path, '03.Quality_filtered', paste0('run', i))
  dir.create(run_dir, recursive = TRUE) # create parent directories if they don't exist
  dir.create(file.path(run_dir, 'fwd'))
  dir.create(file.path(run_dir, 'rev'))
}
## Subdirectories for merged and denoised sequence tables
for (i in 1:num_runs) {
  run_dir = file.path(path, '04.Denoised', paste0('run', i))
  dir.create(run_dir, recursive = TRUE) # create parent directories if they don't exist
}

# Quality filter each run individually
for (run in 1:num_runs) {
  
  # Define the input and output directories for the current run
  trim.fwd <- file.path(path, "02.Quality_trimmed", paste0("run", run), "fwd")
  trim.rev <- file.path(path, "02.Quality_trimmed", paste0("run", run), "rev")
  qualFilt.fwd <- file.path(path, "03.Quality_filtered", paste0("run", run), "fwd")
  qualFilt.rev <- file.path(path, "03.Quality_filtered", paste0("run", run), "rev")
  dnoise <- file.path(path, "04.Denoised", paste0("run", run))
  
  # List file names for 'fwd' and 'rev' reads for the current run
  fns.fwd <- sort(list.files(trim.fwd, pattern = "R1.fastq.gz", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fns.fwd), "_"), `[`, 1)  ## Assumes filename = samplename_XXX.fastq.gz
  names(fns.fwd) <- sample.names
  fns.rev <- sort(list.files(trim.rev, pattern = "R2.fastq.gz", full.names = TRUE))
  sample.names.rev <- sapply(strsplit(basename(fns.rev), "_"), `[`, 1)
  names(fns.rev) <- sample.names.rev
  
  # Check that forward and reverse files match for the current run
  stopifnot(identical(sample.names, sample.names.rev))
  
  # Quality filtering for the current run
  qualFilt <- filterAndTrim(fwd = fns.fwd, filt = file.path(qualFilt.fwd, basename(fns.fwd)),
                            rev = fns.rev, filt.rev = file.path(qualFilt.rev, basename(fns.rev)),
                            maxEE = mEE, truncLen = tL, truncQ = tQ, maxN = mNs,
                            rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE,
                            matchIDs = TRUE)

  # Track reads through the quality filtering pipeline
  colnames(qualFilt) <- c("input_qualFilter", "output_qualFilter")
  rownames(qualFilt) <- sample.names
  # Save track file to denoised directory
  write.csv(qualFilt, file = file.path(path, "03.Quality_filtered", paste0("run", run), "stats.csv"))

}

# Denoise each run individually
for (run in 1:num_runs) {
  
  # Define the input and output directories for the current run
  qualFilt <- file.path(path, "03.Quality_filtered", paste0("run", run))
  qualFilt.fwd <- file.path(qualFilt, "fwd")
  qualFilt.rev <- file.path(qualFilt, "rev")
  dnoise <- file.path(path, "04.Denoised", paste0("run", run))

  # Rename sample files to samplename.fastq.gz and list names for denoising for the current run
  filts.fwd <- list.files(qualFilt.fwd, pattern = "R1.fastq.gz", full.names = TRUE)
  sample.names <- sapply(strsplit(basename(filts.fwd), "_"), `[`, 1)  ## Assumes filename = samplename_XXX.fastq.gz
  names(filts.fwd) <- sample.names
  filts.rev <- list.files(qualFilt.rev, pattern = "R2.fastq.gz", full.names = TRUE) ## Assumes filename = samplename_XXX.fastq.gz
  sample.names.rev <- sapply(strsplit(basename(filts.rev), "_"), `[`, 1)
  names(filts.rev) <- sample.names.rev
  
  # Check that forward and reverse files match for the current run
  stopifnot(identical(sample.names, sample.names.rev))
  
  # Learn error rates for the current run
  set.seed(1986)
  err.fwd <- learnErrors(filts.fwd, nbases = 1e8, multithread = TRUE, randomize = TRUE)
  err.rev <- learnErrors(filts.rev, nbases = 1e8, multithread = TRUE, randomize = TRUE)
  
  # Denoise and merge paired-end reads for the current run
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  for(sam in sample.names) {
    message("Processing:", sam)
    derep.fwd <- derepFastq(filts.fwd[[sam]])
    dd.fwd <- dada(derep.fwd, err = err.fwd, multithread = TRUE)
    derep.rev <- derepFastq(filts.rev[[sam]])
    dd.rev <- dada(derep.rev, err = err.rev, multithread = TRUE)
    merger <- mergePairs(dd.fwd, derep.fwd, dd.rev, derep.rev)
    mergers[[sam]] <- merger
  }
  rm(derep.fwd, derep.rev)
  
# Build sequence tables

  seqtab <- makeSequenceTable(mergers)
  save(seqtab, file = file.path(dnoise, "seqtab.rds"))

  # Track reads through the DADA2 pipeline
  quality <- read.csv(file.path(qualFilt, "stats.csv"))
  getN <- function(x) sum(getUniques(x))
  track <- cbind(quality, sapply(mergers, getN))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("sample_names", "input_qualFilter", "output_qualFilter", "merged")
  # Save track file to denoised directory
  write.csv(track, file = file.path(dnoise, "stats.csv"), row.names = FALSE)
  
}
