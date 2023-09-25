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

path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/ITS  # Path to the main 'ITS' data directory
mkdir $path/02.ITS_extracted
mkdir $path/02.ITS_extracted/run29
mkdir $path/02.ITS_extracted/run30
mkdir $path/02.ITS_extracted/run31
mkdir $path/02.ITS_extracted/run32
mkdir $path/02.ITS_extracted/run33
mkdir $path/02.ITS_extracted/run34
mkdir $path/02.ITS_extracted/run35
mkdir $path/03.Quality_trimmed
mkdir $path/03.Quality_trimmed/run29
mkdir $path/03.Quality_trimmed/run30
mkdir $path/03.Quality_trimmed/run31
mkdir $path/03.Quality_trimmed/run32
mkdir $path/03.Quality_trimmed/run33
mkdir $path/03.Quality_trimmed/run34
mkdir $path/03.Quality_trimmed/run35
raw_data=$path/01.Raw_data    # Path to raw data
ITSx=$path/02.ITS_extracted    # Path to the folder for the extracted ITS reads
trimmed=$path/03.Quality_trimmed    # Path to the folder for the quality truncated reads

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda
conda activate shell

echo
echo " Run 29: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run29

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run29/logfile.txt \
    --outfile $ITSx/run29/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run29

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run29/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run29/*.fastq.gz -o $trimmed/run29/
multiqc $trimmed/run29/. -o $trimmed/run29/
rm $trimmed/run29/*fastqc.zip $trimmed/run29/*fastqc.html $trimmed/run29/multiqc_data/*
rmdir $trimmed/run29/multiqc_data

echo
echo " Run 30: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run30

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run30/logfile.txt \
    --outfile $ITSx/run30/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run30

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run30/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run30/*.fastq.gz -o $trimmed/run30/
multiqc $trimmed/run30/. -o $trimmed/run30/
rm $trimmed/run30/*fastqc.zip $trimmed/run30/*fastqc.html $trimmed/run30/multiqc_data/*
rmdir $trimmed/run30/multiqc_data

echo
echo "Run 31: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run31

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run31/logfile.txt \
    --outfile $ITSx/run31/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run31

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run31/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run31/*.fastq.gz -o $trimmed/run31/
multiqc $trimmed/run31/. -o $trimmed/run31/
rm $trimmed/run31/*fastqc.zip $trimmed/run31/*fastqc.html $trimmed/run31/multiqc_data/*
rmdir $trimmed/run31/multiqc_data

echo
echo "Run 32: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run32

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run32/logfile.txt \
    --outfile $ITSx/run32/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run32

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run32/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run32/*.fastq.gz -o $trimmed/run32/
multiqc $trimmed/run32/. -o $trimmed/run32/
rm $trimmed/run32/*fastqc.zip $trimmed/run32/*fastqc.html $trimmed/run32/multiqc_data/*
rmdir $trimmed/run32/multiqc_data

echo
echo "Run 33: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run33

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run33/logfile.txt \
    --outfile $ITSx/run33/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run33

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run33/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run33/*.fastq.gz -o $trimmed/run33/
multiqc $trimmed/run33/. -o $trimmed/run33/
rm $trimmed/run33/*fastqc.zip $trimmed/run33/*fastqc.html $trimmed/run33/multiqc_data/*
rmdir $trimmed/run33/multiqc_data

echo
echo "Run 34: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run34

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run34/logfile.txt \
    --outfile $ITSx/run34/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run34

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run34/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run34/*.fastq.gz -o $trimmed/run34/
multiqc $trimmed/run34/. -o $trimmed/run34/
rm $trimmed/run34/*fastqc.zip $trimmed/run34/*fastqc.html $trimmed/run34/multiqc_data/*
rmdir $trimmed/run34/multiqc_data

echo
echo "Run 35: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run35

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run35/logfile.txt \
    --outfile $ITSx/run35/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run35

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run35/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run35/*.fastq.gz -o $trimmed/run35/
multiqc $trimmed/run35/. -o $trimmed/run35/
rm $trimmed/run35/*fastqc.zip $trimmed/run35/*fastqc.html $trimmed/run35/multiqc_data/*
rmdir $trimmed/run35/multiqc_data

conda deactivate
echo
echo "Finished at: $(date)"

