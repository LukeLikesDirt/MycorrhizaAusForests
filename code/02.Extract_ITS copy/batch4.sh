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
mkdir $path/02.ITS_extracted/run22
mkdir $path/02.ITS_extracted/run23
mkdir $path/02.ITS_extracted/run24
mkdir $path/02.ITS_extracted/run25
mkdir $path/02.ITS_extracted/run26
mkdir $path/02.ITS_extracted/run27
mkdir $path/02.ITS_extracted/run28
mkdir $path/03.Quality_trimmed
mkdir $path/03.Quality_trimmed/run22
mkdir $path/03.Quality_trimmed/run23
mkdir $path/03.Quality_trimmed/run24
mkdir $path/03.Quality_trimmed/run25
mkdir $path/03.Quality_trimmed/run26
mkdir $path/03.Quality_trimmed/run27
mkdir $path/03.Quality_trimmed/run28
raw_data=$path/01.Raw_data    # Path to raw data
ITSx=$path/02.ITS_extracted    # Path to the folder for the extracted ITS reads
trimmed=$path/03.Quality_trimmed    # Path to the folder for the quality truncated reads

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda
conda activate shell

echo
echo " Run 22: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run22

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run22/logfile.txt \
    --outfile $ITSx/run22/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run22

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run22/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run22/*.fastq.gz -o $trimmed/run22/
multiqc $trimmed/run22/. -o $trimmed/run22/
rm $trimmed/run22/*fastqc.zip $trimmed/run22/*fastqc.html $trimmed/run22/multiqc_data/*
rmdir $trimmed/run22/multiqc_data

echo
echo " Run 23: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run23

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run23/logfile.txt \
    --outfile $ITSx/run23/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run23

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run23/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run23/*.fastq.gz -o $trimmed/run23/
multiqc $trimmed/run23/. -o $trimmed/run23/
rm $trimmed/run23/*fastqc.zip $trimmed/run23/*fastqc.html $trimmed/run23/multiqc_data/*
rmdir $trimmed/run23/multiqc_data

echo
echo "Run 24: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run24

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run24/logfile.txt \
    --outfile $ITSx/run24/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run24

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run24/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run24/*.fastq.gz -o $trimmed/run24/
multiqc $trimmed/run24/. -o $trimmed/run24/
rm $trimmed/run24/*fastqc.zip $trimmed/run24/*fastqc.html $trimmed/run24/multiqc_data/*
rmdir $trimmed/run24/multiqc_data

echo
echo "Run 25: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run25

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run25/logfile.txt \
    --outfile $ITSx/run25/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run25

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run25/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run25/*.fastq.gz -o $trimmed/run25/
multiqc $trimmed/run25/. -o $trimmed/run25/
rm $trimmed/run25/*fastqc.zip $trimmed/run25/*fastqc.html $trimmed/run25/multiqc_data/*
rmdir $trimmed/run25/multiqc_data

echo
echo "Run 26: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run26

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run26/logfile.txt \
    --outfile $ITSx/run26/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run26

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run26/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run26/*.fastq.gz -o $trimmed/run26/
multiqc $trimmed/run26/. -o $trimmed/run26/
rm $trimmed/run26/*fastqc.zip $trimmed/run26/*fastqc.html $trimmed/run26/multiqc_data/*
rmdir $trimmed/run26/multiqc_data

echo
echo "Run 27: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run27

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run27/logfile.txt \
    --outfile $ITSx/run27/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run27

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run27/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run27/*.fastq.gz -o $trimmed/run27/
multiqc $trimmed/run27/. -o $trimmed/run27/
rm $trimmed/run27/*fastqc.zip $trimmed/run27/*fastqc.html $trimmed/run27/multiqc_data/*
rmdir $trimmed/run27/multiqc_data

echo
echo "Run 28: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run28

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run28/logfile.txt \
    --outfile $ITSx/run28/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run28

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run28/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run28/*.fastq.gz -o $trimmed/run28/
multiqc $trimmed/run28/. -o $trimmed/run28/
rm $trimmed/run28/*fastqc.zip $trimmed/run28/*fastqc.html $trimmed/run28/multiqc_data/*
rmdir $trimmed/run28/multiqc_data

conda deactivate
echo
echo "Finished at: $(date)"


