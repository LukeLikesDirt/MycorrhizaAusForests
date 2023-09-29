# This script performs quality filtering, denoising, and sequence table construction for paired-end 16S rRNA gene amplicon data using DADA2.

## Load DADA2
require(dada2, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/envs/R-packages')

# Set the working directory to the project directory
setwd("/data/group/frankslab/project/LFlorence/MycorrhizaAusForests")

# Define the path to the directory containing the runs
path <- "data/AusMicrobiome/ITS"

# Define the number of sequencing runs thta will be processed
num_runs <- 44

# Define the parameters for quality filtering
mEE <- 1      # maxEE: maximum number of expected errors per read
tL <- 0       # truncLen: fixed length for truncating reads
tQ <- 0       # truncQ: quality threshold for truncating reads
mNs <- 0      # maxN: maximum number of Ns allowed in a read
# NOTE: tL and tQ are set to 0 to because the reads are already trimmed and truncated with ITSxpress and Trimmomatic
# NOTE: I am using a maxEE of 1, which is relativley strict when using a "denoising" workflow. However, I am employin strict quality filtering because reads could not be merged.

# Create subdirectories for quality filtered 'fwd' and 'rev' fastq files
for (i in 1:num_runs) {
  run_dir = file.path(path, '04.Quality_filtered', paste0('run', i))
  dir.create(run_dir, recursive = TRUE) # create parent directories if they don't exist
  dir.create(file.path(run_dir, 'fwd'))
  dir.create(file.path(run_dir, 'rev'))
}
## Subdirectories for merged and denoised sequence tables
for (i in 1:num_runs) {
  run_dir = file.path(path, '05.Denoised', paste0('run', i))
  dir.create(run_dir, recursive = TRUE) # create parent directories if they don't exist
}

# Quality filter each run individually
for (run in 1:num_runs) {
    
    ## Input and output directories
    trimmed_dir <- file.path(path, paste0('03.Quality_trimmed/run', run_number, '/'))
    qualFilt_dir <- file.path(path, paste0('04.Quality_filtered/run', run_number))
    denoised_dir <- file.path(path, '05.Denoised')
    
    ## List of files to be filtered
    fns = list.files(trimmed, pattern = 'fastq.gz')
    sample.names <- sapply(strsplit(basename(fns), "_"), `[`, 1)  ## Assumes filename = samplename_XXX.fastq.gz
    names(fns) <- sample.names
    ## Quality filtering
    # Quality filtering for the current run
    qualFilt <- filterAndTrim(fns, file.path(qualFilt_dir, basename(fns)),
                              maxEE = mEE, truncLen = tL, truncQ = tQ, maxN = mNs,
                              rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE,
                              matchIDs = TRUE)
    
    # Track reads through the quality filtering pipeline
    colnames(qualFilt) <- c("input_qualFilter", "output_qualFilter")
    rownames(qualFilt) <- sample.names
    # Save track file to denoised directory
    write.csv(qualFilt, file = file.path(path, "04.Quality_filtered", paste0("run", run), "stats.csv"))

}

# Denoise each run individually
for (run in 1:num_runs) {

    # Define the input and output directories for the current run
    qualFilt <- file.path(path, "03.Quality_filtered", paste0("run", run))
    dnoise <- file.path(path, "04.Denoised", paste0("run", run))

    ## Denoising
    filts = list.files(qualFiltered, pattern = 'fastq.gz', full.names = T)
    sample.names = sapply(strsplit(basename(filts), '_'), `[`, 1)
    names(filts) = sample.names
    
    ## Learn error rates
    set.seed(1986)
    err = learnErrors(filts, nbases = 1e8, multithread = T, randomize = T)
    
    ## Infer sequence variants
    dds = vector('list', length(sample.names))
    names(dds) = sample.names
    for (sam in sample.names) {
        cat('Processing:', sam, '\n')
        derep = derepFastq(filts[[sam]])
        dds[[sam]] = dada(derep, err = err, multithread = TRUE)
    }
    
    ## Build sequence table
    seqtab = makeSequenceTable(dds)
    saveRDS(seqtab, file.path(path, paste0('04.Denoised/run', run_number, '_seqtab.rds')))

    # Track reads through the DADA2 pipeline
    quality <- read.csv(file.path(qualFilt, "stats.csv"))
    getN <- function(x) sum(getUniques(x))
    track <- cbind(quality, sapply(dds, getN))
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("sample_names", "input_qualFilter", "output_qualFilter", "output_denoised")
    # Save track file to denoised directory
    write.csv(track, file = file.path(dnoise, "stats.csv"), row.names = FALSE)
  
}
