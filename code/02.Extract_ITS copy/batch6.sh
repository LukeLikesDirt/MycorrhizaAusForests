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
mkdir $path/02.ITS_extracted/run36
mkdir $path/02.ITS_extracted/run37
mkdir $path/02.ITS_extracted/run38
mkdir $path/02.ITS_extracted/run39
mkdir $path/02.ITS_extracted/run40
mkdir $path/02.ITS_extracted/run41
mkdir $path/02.ITS_extracted/run42
mkdir $path/02.ITS_extracted/run43
mkdir $path/02.ITS_extracted/run44
mkdir $path/03.Quality_trimmed
mkdir $path/03.Quality_trimmed/run36
mkdir $path/03.Quality_trimmed/run37
mkdir $path/03.Quality_trimmed/run38
mkdir $path/03.Quality_trimmed/run39
mkdir $path/03.Quality_trimmed/run40
mkdir $path/03.Quality_trimmed/run41
mkdir $path/03.Quality_trimmed/run42
mkdir $path/03.Quality_trimmed/run43
mkdir $path/03.Quality_trimmed/run44
raw_data=$path/01.Raw_data    # Path to raw data
ITSx=$path/02.ITS_extracted    # Path to the folder for the extracted ITS reads
trimmed=$path/03.Quality_trimmed    # Path to the folder for the quality truncated reads

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda
conda activate shell

echo
echo " Run 36: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run36

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run36/logfile.txt \
    --outfile $ITSx/run36/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run36

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run36/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run36/*.fastq.gz -o $trimmed/run36/
multiqc $trimmed/run36/. -o $trimmed/run36/
rm $trimmed/run36/*fastqc.zip $trimmed/run36/*fastqc.html $trimmed/run36/multiqc_data/*
rmdir $trimmed/run36/multiqc_data

echo
echo " Run 37: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run37

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run37/logfile.txt \
    --outfile $ITSx/run37/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run37

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run37/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run37/*.fastq.gz -o $trimmed/run37/
multiqc $trimmed/run37/. -o $trimmed/run37/
rm $trimmed/run37/*fastqc.zip $trimmed/run37/*fastqc.html $trimmed/run37/multiqc_data/*
rmdir $trimmed/run37/multiqc_data

echo
echo "Run 38: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run38

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run38/logfile.txt \
    --outfile $ITSx/run38/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run38

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run38/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run38/*.fastq.gz -o $trimmed/run38/
multiqc $trimmed/run38/. -o $trimmed/run38/
rm $trimmed/run38/*fastqc.zip $trimmed/run38/*fastqc.html $trimmed/run38/multiqc_data/*
rmdir $trimmed/run38/multiqc_data

echo
echo "Run 39: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run39

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run39/logfile.txt \
    --outfile $ITSx/run39/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run39

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run39/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run39/*.fastq.gz -o $trimmed/run39/
multiqc $trimmed/run39/. -o $trimmed/run39/
rm $trimmed/run39/*fastqc.zip $trimmed/run39/*fastqc.html $trimmed/run39/multiqc_data/*
rmdir $trimmed/run39/multiqc_data

echo
echo "Run 40: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run40

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run40/logfile.txt \
    --outfile $ITSx/run40/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run40

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run40/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run40/*.fastq.gz -o $trimmed/run40/
multiqc $trimmed/run40/. -o $trimmed/run40/
rm $trimmed/run40/*fastqc.zip $trimmed/run40/*fastqc.html $trimmed/run40/multiqc_data/*
rmdir $trimmed/run40/multiqc_data

echo
echo "Run 41: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run41

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run41/logfile.txt \
    --outfile $ITSx/run41/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run41

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run41/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run41/*.fastq.gz -o $trimmed/run41/
multiqc $trimmed/run41/. -o $trimmed/run41/
rm $trimmed/run41/*fastqc.zip $trimmed/run41/*fastqc.html $trimmed/run41/multiqc_data/*
rmdir $trimmed/run41/multiqc_data

echo
echo "Run 42: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run42

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run42/logfile.txt \
    --outfile $ITSx/run42/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run42

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run42/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run42/*.fastq.gz -o $trimmed/run42/
multiqc $trimmed/run42/. -o $trimmed/run42/
rm $trimmed/run42/*fastqc.zip $trimmed/run42/*fastqc.html $trimmed/run42/multiqc_data/*
rmdir $trimmed/run42/multiqc_data

echo
echo "Run 43: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run43

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run43/logfile.txt \
    --outfile $ITSx/run43/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run43

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run43/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run43/*.fastq.gz -o $trimmed/run43/
multiqc $trimmed/run43/. -o $trimmed/run43/
rm $trimmed/run43/*fastqc.zip $trimmed/run43/*fastqc.html $trimmed/run43/multiqc_data/*
rmdir $trimmed/run43/multiqc_data

echo
echo "Run 44: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run44

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run44/logfile.txt \
    --outfile $ITSx/run44/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run44

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run44/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run44/*.fastq.gz -o $trimmed/run44/
multiqc $trimmed/run44/. -o $trimmed/run44/
rm $trimmed/run44/*fastqc.zip $trimmed/run44/*fastqc.html $trimmed/run44/multiqc_data/*
rmdir $trimmed/run44/multiqc_data

conda deactivate
echo
echo "Finished at: $(date)"

