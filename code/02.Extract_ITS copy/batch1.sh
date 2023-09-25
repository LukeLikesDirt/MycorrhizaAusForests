#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1440:00:00
#SBATCH --partition=long
#SBATCH --mem=100G

echo "Starting at: $(date)"

## Prepare Illumina forward reads for denoising with DADA2 by extracting the ITS1 subregion using ITSxpress and quality trun0cating reads using Trimmomatic
## Reads cannot be merged due to the high error rates: error rates are routine for Australian Microbiome data, hence only the forward reads are used
## Primers have not been removed, which is also due to the high error rates

## Organise subdirectories
path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/ITS  # Path to the main 'ITS' data directory
mkdir $path/02.ITS_extracted
mkdir $path/02.ITS_extracted/run1
mkdir $path/02.ITS_extracted/run2
mkdir $path/02.ITS_extracted/run3
mkdir $path/02.ITS_extracted/run4
mkdir $path/02.ITS_extracted/run5
mkdir $path/02.ITS_extracted/run6
mkdir $path/02.ITS_extracted/run7
mkdir $path/03.Quality_trimmed
mkdir $path/03.Quality_trimmed/run1
mkdir $path/03.Quality_trimmed/run2
mkdir $path/03.Quality_trimmed/run3
mkdir $path/03.Quality_trimmed/run4
mkdir $path/03.Quality_trimmed/run5
mkdir $path/03.Quality_trimmed/run6
mkdir $path/03.Quality_trimmed/run7
raw_data=$path/01.Raw_data    # Path to raw data
ITSx=$path/02.ITS_extracted    # Path to the folder for the extracted ITS reads
trimmed=$path/03.Quality_trimmed    # Path to the folder for the quality truncated reads

## Activate conda environment
source /data/group/frankslab/home/21258990/mambaforge/etc/profile.d/conda.sh    # Path to conda
conda activate shell

echo
echo "Run 1: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run1

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run1/logfile.txt \
    --outfile $ITSx/run1/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run1

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run1/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run1/*.fastq.gz -o $trimmed/run1/
multiqc $trimmed/run1/. -o $trimmed/run1/
rm $trimmed/run1/*fastqc.zip $trimmed/run1/*fastqc.html $trimmed/run1/multiqc_data/*
rmdir $trimmed/run1/multiqc_data

echo
echo "Run 2: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run2

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run2/logfile.txt \
    --outfile $ITSx/run2/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run2

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run2/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run2/*.fastq.gz -o $trimmed/run2/
multiqc $trimmed/run2/. -o $trimmed/run2/
rm $trimmed/run2/*fastqc.zip $trimmed/run2/*fastqc.html $trimmed/run2/multiqc_data/*
rmdir $trimmed/run2/multiqc_data

echo
echo "Run 3: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run3

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run3/logfile.txt \
    --outfile $ITSx/run3/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run3

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run3/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run3/*.fastq.gz -o $trimmed/run3/
multiqc $trimmed/run3/. -o $trimmed/run3/
rm $trimmed/run3/*fastqc.zip $trimmed/run3/*fastqc.html $trimmed/run3/multiqc_data/*
rmdir $trimmed/run3/multiqc_data

echo
echo "Run 4: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run4

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run4/logfile.txt \
    --outfile $ITSx/run4/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run4

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run4/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run4/*.fastq.gz -o $trimmed/run4/
multiqc $trimmed/run4/. -o $trimmed/run4/
rm $trimmed/run4/*fastqc.zip $trimmed/run4/*fastqc.html $trimmed/run4/multiqc_data/*
rmdir $trimmed/run4/multiqc_data

echo
echo "Run 5: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run5

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run5/logfile.txt \
    --outfile $ITSx/run5/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run5

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run5/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run5/*.fastq.gz -o $trimmed/run5/
multiqc $trimmed/run5/. -o $trimmed/run5/
rm $trimmed/run5/*fastqc.zip $trimmed/run5/*fastqc.html $trimmed/run5/multiqc_data/*
rmdir $trimmed/run5/multiqc_data

echo
echo "Run 6: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run6

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run6/logfile.txt \
    --outfile $ITSx/run6/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run6

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run6/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run6/*.fastq.gz -o $trimmed/run6/
multiqc $trimmed/run6/. -o $trimmed/run6/
rm $trimmed/run6/*fastqc.zip $trimmed/run6/*fastqc.html $trimmed/run6/multiqc_data/*
rmdir $trimmed/run6/multiqc_data

echo
echo "Run 7: Extract ITS1 and quality truncate reads"

## Extract ITS1
cd $raw_data/run7

for f in *_R1.fastq.gz; do

itsxpress \
    --single_end \
    --fastq ${f} \
    --cluster_id 1.0 \
    --region ITS1 \
    --taxa All \
    --log $ITSx/run7/logfile.txt \
    --outfile $ITSx/run7/${f} \
    --threads 8

done

## Quality trim reads
cd $ITSx/run7

for f in *.fastq.gz; do

    trimmomatic SE $f $trimmed/run7/$f \
    SLIDINGWINDOW:4:15 MINLEN:80 \
    -threads 8

done

## Generate quality report
fastqc $trimmed/run7/*.fastq.gz -o $trimmed/run7/
multiqc $trimmed/run7/. -o $trimmed/run7/
rm $trimmed/run7/*fastqc.zip $trimmed/run7/*fastqc.html $trimmed/run7/multiqc_data/*
rmdir $trimmed/run7/multiqc_data

conda deactivate
echo
echo "Finished at: $(date)"
