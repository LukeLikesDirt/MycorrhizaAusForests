
conda activate R

## Load DADA2
require(dada2, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/envs/R-packages')

## Organise directory
path = '/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/16S/test'    ## Main ITS data directory
dir.create(file.path(path, '02.Quality_filtered')) ## Directories for quality filter output
dir.create(file.path(path, '03.Denoised'))    ## Directories for denoised files
## Define input and output directories
raw_data = file.path(path, '01.Raw_data')
qualFiltered = file.path(path, '02.Quality_filtered')
denoised = file.path(path, '03.Denoised')
## List of files to be filtered
fnsF = sort(list.files(raw_data, pattern="R1.fastq.gz"))
fnsR = sort(list.files(raw_data, pattern="R2.fastq.gz"))
if(length(fnsF) != length(fnsR)) stop("Forward and reverse files do not match.")    

## Quality filtering
r1 = filterAndTrim(file.path(raw_data, fnsF), file.path(qualFiltered, fnsF),
              truncLen=0, maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)
r2 = filterAndTrim(file.path(raw_data, fnsR), file.path(qualFiltered, fnsR),
              truncLen=250, maxEE=3, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)


f1 = filterAndTrim(fwd=file.path(raw_data, fnsF), filt=file.path(qualFiltered, fnsF),
              rev=file.path(raw_data, fnsR), filt.rev=file.path(qualFiltered, fnsR),
              truncLen=c(280,230), maxEE=c(2,3), truncQ=2, maxN=0, rm.phix=TRUE,
              trimLeft = 10,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

    
## Denoising
filtsF = list.files(qualFiltered, pattern = 'R1.fastq.gz', full.names = T)
sample.namesF = sapply(strsplit(basename(filtsF), '_'), `[`, 1)
names(filtsF) = sample.names
    
## Learn error rates
set.seed(1986)
errF = learnErrors(filtsF, nbases = 1e8, multithread = T, randomize = T)
errR = learnErrors(filtsR, nbases = 1e8, multithread = T, randomize = T)
    
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