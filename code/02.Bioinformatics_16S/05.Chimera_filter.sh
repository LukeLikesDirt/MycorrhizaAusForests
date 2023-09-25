echo
echo Obtaining Gold reference database for chimera detection
path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/16S
mkdir $path/reference_datasets
REF=$path/reference_datasets/gold.fasta
GOLD="https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip"

if [ ! -e gold.fasta ]; then

    if [ ! -e silva.gold.bacteria.zip ]; then
        wget -P$path/reference_datasets/ $GOLD
    fi

    echo Decompressing and reformatting...
    unzip -p silva.gold.bacteria.zip silva.gold.align | \
        sed -e "s/[.-]//g" > gold.fasta

fi
