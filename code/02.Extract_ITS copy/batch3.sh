#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1440:00:00
#SBATCH --partition=long
#SBATCH --mem=100G

echo "Starting at: $(date)"

## Prepare Illumina forward reads for denoising with DADA2 by extracting the ITS1 subregion using ITSxpress and quality truncating reads using Trimmomatic
## Reads cannot be merged due to the high error rates: error rates are routine for Australian Microbiome data, hence only the forward reads are used
## Primers have not been removed, which is also due to the high error rates

## Organise subdirectories
## Organise subdirectories
path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/ITS  # Path to the main 'ITS' data directory
mkdir $path/02.ITS_extracted
mkdir $path/02.ITS_extracted/run15
mkdir $path/02.ITS_extracted/run16
mkdir $path/02.ITS_extracted/run17
mkdir $path/02.ITS_extracted/run18
mkdir $path/02.ITS_extracted/run19
mkdir $path/02.ITS_extracted/run20
mkdir $path/02.ITS_extracted/run21
mkdir $path/03.Quality_trimmed
mkdir $path/03.Quality_trimmed/run15
mkdir $path/03.Quality_trimmed/run16
mkdir $path/03.Quality_trimmed/run17
mkdir $path/03.Quality_trimmed/run18
mkdir $path/03.Quality_trimmed/run19
mkdir $path/03.Quality_trimmed/run20
mkdir $path/03.Quality_trimmed/run21
raw_data=$path/01.Raw_data    # Path to raw data
ITSx=$path/02.ITS_extracted    # Path to the folder for the extracted ITS reads
trimmed=$path/03.Quality_trimmed    # Path to the folder for the quality truncated reads

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda
conda activate shell

echo
echo "Run 15: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run15

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run15/logfile.txt \
    --outfile $ITSx/run15/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run15

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run15/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run15/*.fastq.gz -o $trimmed/run15/
multiqc $trimmed/run15/. -o $trimmed/run15/
rm $trimmed/run15/*fastqc.zip $trimmed/run15/*fastqc.html $trimmed/run15/multiqc_data/*
rmdir $trimmed/run15/multiqc_data

echo
echo "Run 16: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run16

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run16/logfile.txt \
    --outfile $ITSx/run16/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run16

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run16/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run16/*.fastq.gz -o $trimmed/run16/
multiqc $trimmed/run16/. -o $trimmed/run16/
rm $trimmed/run16/*fastqc.zip $trimmed/run16/*fastqc.html $trimmed/run16/multiqc_data/*
rmdir $trimmed/run16/multiqc_data

echo
echo "Run 17: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run17

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run17/logfile.txt \
    --outfile $ITSx/run17/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run17

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run17/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run17/*.fastq.gz -o $trimmed/run17/
multiqc $trimmed/run17/. -o $trimmed/run17/
rm $trimmed/run17/*fastqc.zip $trimmed/run17/*fastqc.html $trimmed/run17/multiqc_data/*
rmdir $trimmed/run17/multiqc_data

echo
echo "Run 18: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run18

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run18/logfile.txt \
    --outfile $ITSx/run18/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run18

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run18/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run18/*.fastq.gz -o $trimmed/run18/
multiqc $trimmed/run18/. -o $trimmed/run18/
rm $trimmed/run18/*fastqc.zip $trimmed/run18/*fastqc.html $trimmed/run18/multiqc_data/*
rmdir $trimmed/run18/multiqc_data

echo
echo "Run 19: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run19

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run19/logfile.txt \
    --outfile $ITSx/run19/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run19

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run19/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run19/*.fastq.gz -o $trimmed/run19/
multiqc $trimmed/run19/. -o $trimmed/run19/
rm $trimmed/run19/*fastqc.zip $trimmed/run19/*fastqc.html $trimmed/run19/multiqc_data/*
rmdir $trimmed/run19/multiqc_data

echo
echo "Run 20: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run20

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run20/logfile.txt \
    --outfile $ITSx/run20/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run20

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run20/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run20/*.fastq.gz -o $trimmed/run20/
multiqc $trimmed/run20/. -o $trimmed/run20/
rm $trimmed/run20/*fastqc.zip $trimmed/run20/*fastqc.html $trimmed/run20/multiqc_data/*
rmdir $trimmed/run20/multiqc_data

echo
echo "Run 21: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run21

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run21/logfile.txt \
    --outfile $ITSx/run21/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run21

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run21/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run21/*.fastq.gz -o $trimmed/run21/
multiqc $trimmed/run21/. -o $trimmed/run21/
rm $trimmed/run21/*fastqc.zip $trimmed/run21/*fastqc.html $trimmed/run21/multiqc_data/*
rmdir $trimmed/run21/multiqc_data

conda deactivate
echo
echo "Finished at: $(date)"

