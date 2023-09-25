
## Load DADA2
require(dada2, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/envs/R-packages')

## Organise directory
path = '/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/ITS'    ## Main ITS data directory
## Directories for quality filter output
dir.create(file.path(path, '04.Quality_filtered'))   
for (i in 1:44) {
  dir.create(file.path(path, '04.Quality_filtered', paste0('run', i)))
}
## Directories for denoised output
dir.create(file.path(path, '05.Denoised'))    ## Directories for denoised files

# Define the number of runs
num_runs = 44

for (run_number in 1:num_runs) {
    cat('Processing Run', run_number, '\n')
    
    ## Input and output directories
    trimmed = file.path(path, paste0('03.Quality_trimmed/run', run_number, '/'))
    qualFiltered = file.path(path, paste0('04.Quality_filtered/run', run_number))
    denoised = file.path(path, '05.Denoised')
    
    ## List of files to be filtered
    fns = list.files(trimmed, pattern = 'fastq.gz')
    
    ## Quality filtering
    filterAndTrim(file.path(trimmed, fns),
                  file.path(qualFiltered, fns),
                  truncLen = 0,
                  maxEE = 1,
                  minLen = 80,
                  rm.phix = T,
                  compress = T,
                  verbose = T,
                  multithread = T)
    
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
}