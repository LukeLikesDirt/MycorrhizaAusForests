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
path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/ITS  # Path to the main 'ITS' data directory
mkdir $path/02.ITS_extracted
mkdir $path/02.ITS_extracted/run8
mkdir $path/02.ITS_extracted/run9
mkdir $path/02.ITS_extracted/run10
mkdir $path/02.ITS_extracted/run11
mkdir $path/02.ITS_extracted/run12
mkdir $path/02.ITS_extracted/run13
mkdir $path/02.ITS_extracted/run14
mkdir $path/03.Quality_trimmed
mkdir $path/03.Quality_trimmed/run8
mkdir $path/03.Quality_trimmed/run9
mkdir $path/03.Quality_trimmed/run10
mkdir $path/03.Quality_trimmed/run11
mkdir $path/03.Quality_trimmed/run12
mkdir $path/03.Quality_trimmed/run13
mkdir $path/03.Quality_trimmed/run14
raw_data=$path/01.Raw_data    # Path to raw data
ITSx=$path/02.ITS_extracted    # Path to the folder for the extracted ITS reads
trimmed=$path/03.Quality_trimmed    # Path to the folder for the quality truncated reads

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda
conda activate shell

echo
echo "Run 8: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run8

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run8/logfile.txt \
    --outfile $ITSx/run8/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run8

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run8/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run8/*.fastq.gz -o $trimmed/run8/
multiqc $trimmed/run8/. -o $trimmed/run8/
rm $trimmed/run8/*fastqc.zip $trimmed/run8/*fastqc.html $trimmed/run8/multiqc_data/*
rmdir $trimmed/run8/multiqc_data

echo
echo "Run 9: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run9

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run9/logfile.txt \
    --outfile $ITSx/run9/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run9

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run9/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run9/*.fastq.gz -o $trimmed/run9/
multiqc $trimmed/run9/. -o $trimmed/run9/
rm $trimmed/run9/*fastqc.zip $trimmed/run9/*fastqc.html $trimmed/run9/multiqc_data/*
rmdir $trimmed/run9/multiqc_data

echo
echo "Run 10: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run10

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run10/logfile.txt \
    --outfile $ITSx/run10/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run10

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run10/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run10/*.fastq.gz -o $trimmed/run10/
multiqc $trimmed/run10/. -o $trimmed/run10/
rm $trimmed/run10/*fastqc.zip $trimmed/run10/*fastqc.html $trimmed/run10/multiqc_data/*
rmdir $trimmed/run10/multiqc_data

echo
echo "Run 11: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run11

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run11/logfile.txt \
    --outfile $ITSx/run11/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run11

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run11/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run11/*.fastq.gz -o $trimmed/run11/
multiqc $trimmed/run11/. -o $trimmed/run11/
rm $trimmed/run11/*fastqc.zip $trimmed/run11/*fastqc.html $trimmed/run11/multiqc_data/*
rmdir $trimmed/run11/multiqc_data

echo
echo "Run 12: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run12

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run12/logfile.txt \
    --outfile $ITSx/run12/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run12

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run12/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run12/*.fastq.gz -o $trimmed/run12/
multiqc $trimmed/run12/. -o $trimmed/run12/
rm $trimmed/run12/*fastqc.zip $trimmed/run12/*fastqc.html $trimmed/run12/multiqc_data/*
rmdir $trimmed/run12/multiqc_data

echo
echo "Run 13: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run13

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run13/logfile.txt \
    --outfile $ITSx/run13/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run13

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run13/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run13/*.fastq.gz -o $trimmed/run13/
multiqc $trimmed/run13/. -o $trimmed/run13/
rm $trimmed/run13/*fastqc.zip $trimmed/run13/*fastqc.html $trimmed/run13/multiqc_data/*
rmdir $trimmed/run13/multiqc_data

echo
echo "Run 14: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run14

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run14/logfile.txt \
    --outfile $ITSx/run14/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run14

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run14/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run14/*.fastq.gz -o $trimmed/run14/
multiqc $trimmed/run14/. -o $trimmed/run14/
rm $trimmed/run14/*fastqc.zip $trimmed/run14/*fastqc.html $trimmed/run14/multiqc_data/*
rmdir $trimmed/run14/multiqc_data

conda deactivate
echo
echo "Finished at: $(date)"


