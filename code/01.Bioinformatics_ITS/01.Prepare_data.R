

## Organise ITS data from the Australian Microbiome Initiative (AMI) dataset
##
## First, I will subet samples by seuqencing run prior to denoising with
## DADA2 because:
##      (1) the AMI dataset spans multiple sequencing runs
##      (2) different runs have unique error profiles
##      (3) and DADA2 denoising is based on estimated error rates.
##
## Second, I will organise subdirectorie by seuqencing run and move fastq files
## to the raw data subdirectories.
##
##
## Third, I will (1) rename control samples to ratain run ID, and (2) append
## 'control type' (i.e. mock community, positive, negative)to the ID so that
## both run ID and control type are retained when renaming all files within the
## DADA2 step of this pipeline.
##
## Forth, I will add 's' to the beginning of sample names so that they are
## ready to be analysed in R, which requires sample names to start with a
## letter.
##
## NOTES:
## Data within the subdiretories 'bpa_2e86107a_20230916T0721' and
## 'bpa_c4ca4509_20230916T0936' and were downloaded from the AMI website on
## 16/09/2023. The 'bpa_2e86107a_20230916T0721' subdirectory contains all soil
## ITS test samples from the AMI website, and the 'bpa_c4ca4509_20230916T0936'
## subdirectory contains the control samples. Data within the subdiretories
## 'bpa_4f999bb9_20230717T0504' and 'bpa_10ad26d2_20230717T0504' are from a
## spcific Govenrment of Victoria project and are not yet avialble on the AMI
## website, yet were provided to me by the state goverment for this project.


###############################################################################
## (1) Subset data by sequencing run ##########################################
###############################################################################

## Activate conda environment
conda activate R
## Activate R terminal
R
## Load required packages
require(tidyverse)
## Define path to main data directory
path = '/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/ITS'

## Filter ID run names
data1 = read.csv(file.path(path, '/bpa_e1f6092f_20230923T1206/package_metadata/package_metadata_bpa_e1f6092f_20230923T1206_amdb-genomics-amplicon.csv'), header = T) %>% 
    # Select sample titles
    select(title) %>%  
    # Rename tile to run ID
    mutate(title = as.character(gsub('Australian Microbiome Amplicons ITS 102.100.100/', '', title))) %>%
    separate(title, c('sample_ID', 'run_ID'), sep = ' ') %>%
    # Reduce to unique run IDs
    select(run_ID) %>%
    unique() %>%
    as_tibble()
data2 = read.csv(file.path(path, '/bpa_e1f6092f_20230923T1206/package_metadata/package_metadata_bpa_e1f6092f_20230923T1206_base-genomics-amplicon.csv'), header = T) %>% 
    select(title) %>%
    mutate(title = as.character(gsub('BASE Amplicons ITS ', '', title))) %>%
    separate(title, c('sample_ID', 'run_ID'), sep = ' ') %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()
data3 = read.csv(file.path(path, '/bpa_4f999bb9_20230717T0504/package_metadata/package_metadata_bpa_4f999bb9_20230717T0504_amdb-genomics-amplicon.csv'), header = T) %>% 
    select(title) %>%
    mutate(title = as.character(gsub('Australian Microbiome Amplicons ITS 102.100.100/', '', title))) %>%
    separate(title, c('sample_ID', 'run_ID'), sep = ' ') %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()
data4 = read.csv(file.path(path, '/bpa_10ad26d2_20230717T0504/package_metadata/package_metadata_bpa_10ad26d2_20230717T0504_amdb-genomics-amplicon.csv'), header = T) %>% 
    select(title) %>%
    mutate(title = as.character(gsub('Australian Microbiome Amplicons ITS 102.100.100/', '', title))) %>%
    separate(title, c('sample_ID', 'run_ID'), sep = ' ') %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()

## Combine datasets
bind_rows(data1, data2, data3, data4) %>%
    as_tibble() %>%
    print(n = Inf)
## There are 44 unique runs

## Most of these data did not come with thier control sample
## Here I match runs with control samples, which were provided as a seperate file

## Load control sample data
data5 = read.csv(file.path(path, '/bpa_07533810_20230925T0954/package_metadata/package_metadata_bpa_07533810_20230925T0954_amdb-genomics-amplicon-control.csv'), header = T) %>% 
    select(title) %>%
    mutate(run_ID = as.character(gsub('Australian Microbiome Amplicons Control ITS ', '', title))) %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()
data6 = read.csv(file.path(path, '/bpa_07533810_20230925T0954/package_metadata/package_metadata_bpa_07533810_20230925T0954_base-genomics-amplicon-control.csv'), header = T) %>% 
    select(title) %>%
    mutate(run_ID = as.character(gsub('BASE Amplicons Control ITS ', '', title))) %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()

## Combine datasets
bind_rows(data5, data6) %>%
    as_tibble() %>%
    print(n = Inf)
## There are 38 unique runs.
## I expected 40 unique runs for controls because 2 of the 44 runs were from
## the specific government of Victoria project. Therefore, therefore some runs
## must not have control samples.

## Quit R
q()
## Do not save workspace image
n
## Deactivate conda environment
conda deactivate

###############################################################################
## (2) Rename control samples #################################################
###############################################################################

## I'll move and rename control samples first so I can manually check them for 
## errors easily, before pairing them with all samples.
## Before I rename the control samples, I organise the subdirectories by run ID.

## Organise subdirectories
path=/data/group/frankslab/project/LFlorence/MycorrhizasAustralianForests/data/AusMicrobiome/ITS    # Path to main data directory
mkdir $path/01.Raw_data    # Path to send raw data
raw_data=$path/01.Raw_data    # Define path raw data subdirectory
mkdir $raw_data/run{1..44}    # Make subdirectories for each run

## Move control samples to the raw data subdirectory
for f in $path/bpa_c4ca4509_20230916T0936/*fastq.gz; do mv $f $raw_data
done
## The control samples from the Victorian Government are in the 'bpa...504'
## subdirectory with test samples, so I'll move them to the raw data
## subdirectory manually...

## Move control samples to run ID subdirectories
for f in $raw_data/*; do
  case "$f" in
    *J9TWH*) mv "$f" "$raw_data/run1" ;;
    *J9GNL*) mv "$f" "$raw_data/run2" ;;
    *J9GGM*) mv "$f" "$raw_data/run3" ;;
    *J9FL8*) mv "$f" "$raw_data/run4" ;;
    *K3W8G*) mv "$f" "$raw_data/run5" ;;
    *KBJ9G*) mv "$f" "$raw_data/run6" ;;
    *JC892*) mv "$f" "$raw_data/run7" ;;
    *JC328*) mv "$f" "$raw_data/run8" ;;
    *JY4PG*) mv "$f" "$raw_data/run9" ;;
    *JRK29*) mv "$f" "$raw_data/run10" ;;
    *JRK2C*) mv "$f" "$raw_data/run11" ;;
    *KCVJJ*) mv "$f" "$raw_data/run12" ;;
    *KG28C*) mv "$f" "$raw_data/run13" ;;
    *AAVYC*) mv "$f" "$raw_data/run14" ;;
    *AC988*) mv "$f" "$raw_data/run15" ;;
    *A64HK*) mv "$f" "$raw_data/run16" ;;
    *A64JJ*) mv "$f" "$raw_data/run17" ;;
    *AC9CM*) mv "$f" "$raw_data/run18" ;;
    *B3C4H*) mv "$f" "$raw_data/run19" ;;
    *ANVM7*) mv "$f" "$raw_data/run20" ;;
    *A5P6W*) mv "$f" "$raw_data/run21" ;;
    *AAUF6*) mv "$f" "$raw_data/run22" ;;
    *AADNU*) mv "$f" "$raw_data/run23" ;;
    *ALBRY*) mv "$f" "$raw_data/run24" ;;
    *AN5TK*) mv "$f" "$raw_data/run25" ;;
    *BKG34*) mv "$f" "$raw_data/run26" ;;
    *AGEDA*) mv "$f" "$raw_data/run27" ;;
    *AGEFG*) mv "$f" "$raw_data/run28" ;;
    *AGECM*) mv "$f" "$raw_data/run29" ;;
    *B39FT*) mv "$f" "$raw_data/run30" ;;
    *BC8BY*) mv "$f" "$raw_data/run31" ;;
    *BC267*) mv "$f" "$raw_data/run32" ;;
    *C43MN*) mv "$f" "$raw_data/run33" ;;
    *B39G7*) mv "$f" "$raw_data/run34" ;;
    *BJTVT*) mv "$f" "$raw_data/run35" ;;
    *BY8DT*) mv "$f" "$raw_data/run36" ;;
    *BY245*) mv "$f" "$raw_data/run37" ;;
    *C5M89*) mv "$f" "$raw_data/run38" ;;
    *C43P8*) mv "$f" "$raw_data/run39" ;;
    *A64JH*) mv "$f" "$raw_data/run40" ;;
    *A79NR*) mv "$f" "$raw_data/run41" ;;
    *A90G8*) mv "$f" "$raw_data/run42" ;;
    *KP69D*) mv "$f" "$raw_data/run43" ;;
    *KP8LM*) mv "$f" "$raw_data/run44" ;;
  esac
done

## Remove control samples that are not in the main dataset that I am working
## with.
for f in $raw_data/*.fastq.gz; do
    rm "$f"
done

## Because we won't be merging these reads, I will first delete the 'R2' reads
## along with the 'I1' and 'I2', which are the indexes for the 'R1' and 'R2'
## reads, respectively.
for run_dir in $raw_data/run{1..44}; do
  rm "$run_dir"/*R2* "$run_dir"/*I1* "$run_dir"/*I2*
done

## File name prefixes change depending on the date of data intake into the
## database. I'll nest loops to rename all control samples

## Rename the prefix of mock community samples
## Start from complex to simple names to avoid renaming the wrong files
for run_dir in $raw_data/run{1..44}; do
  for prefix in 'ATCC1010MOCK_ITS_' 'Fungal_mock_community_ITS_AGRF_ACACTAGATCCG_' 'Fungal_mock_community_ITS_AGRF_GCTCGAAGATCG_' 'Fungal_mock_community_ITS_UNSW_AACCGCGGTCAA_' 'Fungal_mock_community_ITS_UNSW_ACACTAGATCCG_' 'Fungal_mock_community_ITS_UNSW_CCATACATAGCT_' 'Fungal_mock_community_ITS_UNSW_CGAGTTGTAGCG_' 'Fungal-mock-community_ITS_AGRF_CCAAGTCTTACA_' 'Fungal-mock-community_ITS_AGRF_' 'Fungal__mock_Community_ITS_AGRF_' 'Fungal_mock_community_ITS_'; do
    for f in "$run_dir"/*"$prefix"*; do
      if [ -e "$f" ]; then
        new_name=$(basename "$f" | sed "s/$prefix/MockCommunity/")
        mv "$f" "$run_dir/$new_name"
      fi
    done
  done
done

## Fix the messed-up a couple of runs
## Two runs have 'AGRF' before the run ID
for run_dir in $raw_data/run{1..44}; do
  for prefix in 'MockCommunityAGRF'; do
    for f in "$run_dir"/*"$prefix"*; do
      if [ -e "$f" ]; then
        new_name=$(basename "$f" | sed "s/$prefix/MockCommunity/")
        mv "$f" "$run_dir/$new_name"
      fi
    done
  done
done

## There is a couple of samples with an '_' before the run ID
for run_dir in $raw_data/run{1..44}; do
  for prefix in 'MockCommunity_'; do
    for f in "$run_dir"/*"$prefix"*; do
      if [ -e "$f" ]; then
        new_name=$(basename "$f" | sed "s/$prefix/MockCommunity/")
        mv "$f" "$run_dir/$new_name"
      fi
    done
  done
done

## There are two random bacterial mock communities in run 20
## I'll remove them here:
for f in $raw_data/run20/Bac_mock_*; do
    rm $f
done

## Now I will rename positive control samples:
for run_dir in $raw_data/run{1..44}; do
  for prefix in 'Soil_DNA_ITS_AGRF_CCAAGTCTTACA_' 'Soil_DNA_ITS_UNSW_AACCGCGGTCAA_' 'Soil_DNA_ITS_UNSW_CGAGTTGTAGCG_' 'Soil_DNA_ITS_UNSW_TAAGGTAAGGTG_' 'Soil_DNA_ITS_AGRF_AACGAGAACTGA_' 'Soil_DNA_ITS_AGRF_AACCGCGGTCAA_' 'Soil_DNA_ITS_UNSW_CCAAGTCTTACA_' 'Soil_DNA_ITS_UNSW_CCAAGTCTTACA_' 'Soil-DNA_ITS_AGRF_ACACTAGATCCG_' 'Soil_DNA_ITS_AGRF_AACCGCGGTCAA_' 'Soil_DNA_ITS_AGRF_CCATACATAGCT_' 'Soil-DNA_ITS_AGRF_' 'Soil_DNA_ITS_AGRF_' 'Soil_DNA_ITS_'; do
    for f in "$run_dir"/*"$prefix"*; do
      if [ -e "$f" ]; then
        new_name=$(basename "$f" | sed "s/$prefix/PositiveControl/")
        mv "$f" "$run_dir/$new_name"
      fi
    done
  done
done

## Rename negative control samples:
for run_dir in $raw_data/run{1..44}; do
  for prefix in 'No_Template_Control_ITS_' 'Neg1_ITS_AGRF_' 'Neg2_ITS_AGRF_' 'NEG_1_ITS_AGRF_CGAGTTGTAGCG_' 'NEG_1_ITS_AGRF_GAAGAAGCGGTA_' 'NEG_1_ITS_AGRF_' 'Neg_ITS_AGRF_'; do
    for f in "$run_dir"/*"$prefix"*; do
      if [ -e "$f" ]; then
        new_name=$(basename "$f" | sed "s/$prefix/NegativeControl/")
        mv "$f" "$run_dir/$new_name"
      fi
    done
  done
done

## Lastly I will rename the controls that have been run in duplicate or
## triplicate.
## I'm not sure how I could loop them so I will run them one by one.
## Run 1
mv "$raw_data/run1/MockCommunityJ9TWH_GCTTCGGT-ACGACGTG_S1_L001_R1.fastq.gz" "$raw_data/run1/MockCommunity1J9TWH_GCTTCGGT-ACGACGTG_S1_L001_R1.fastq.gz"
mv "$raw_data/run1/MockCommunityJ9TWH_TCTGTTGC-ACGACGTG_S2_L001_R1.fastq.gz" "$raw_data/run1/MockCommunity2J9TWH_TCTGTTGC-ACGACGTG_S2_L001_R1.fastq.gz"
mv "$raw_data/run1/NegativeControlJ9TWH_GCTTCGGT-CGTCGCTA_S5_L001_R1.fastq.gz" "$raw_data/run1/NegativeControl1J9TWH_GCTTCGGT-CGTCGCTA_S5_L001_R1.fastq.gz"
mv "$raw_data/run1/NegativeControlJ9TWH_TCTGTTGC-CGTCGCTA_S6_L001_R1.fastq.gz" "$raw_data/run1/NegativeControl2J9TWH_TCTGTTGC-CGTCGCTA_S6_L001_R1.fastq.gz"
mv "$raw_data/run1/PositiveControlJ9TWH_GCTTCGGT-ATATACAC_S3_L001_R1.fastq.gz" "$raw_data/run1/PositiveControl1J9TWH_GCTTCGGT-ATATACAC_S3_L001_R1.fastq.gz"
mv "$raw_data/run1/PositiveControlJ9TWH_TCTGTTGC-ATATACAC_S4_L001_R1.fastq.gz" "$raw_data/run1/PositiveControl2J9TWH_TCTGTTGC-ATATACAC_S4_L001_R1.fastq.gz"
## Run 5
mv "$raw_data/run5/MockCommunityK3W8G_TCTGTTGC-TGCGTACG_S105_L001_R1.fastq.gz" "$raw_data/run5/MockCommunity1K3W8G_TCTGTTGC-TGCGTACG_S105_L001_R1.fastq.gz"
mv "$raw_data/run5/MockCommunityK3W8G_TGACCTCC-GCTCTAGT_S77_L001_R1.fastq.gz" "$raw_data/run5/MockCommunity2K3W8G_TGACCTCC-GCTCTAGT_S77_L001_R1.fastq.gz"
mv "$raw_data/run5/NegativeControlK3W8G_TACGGTAT-ACGACGTG_S107_L001_R1.fastq.gz" "$raw_data/run5/NegativeControl1K3W8G_TACGGTAT-ACGACGTG_S107_L001_R1.fastq.gz"
mv "$raw_data/run5/NegativeControlK3W8G_TGACCTCC-TGCGTACG_S78_L001_R1.fastq.gz" "$raw_data/run5/NegativeControl2K3W8G_TGACCTCC-TGCGTACG_S78_L001_R1.fastq.gz"
mv "$raw_data/run5/PositiveControlK3W8G_TCTGTTGC-TAGTGTAG_S106_L001_R1.fastq.gz" "$raw_data/run5/PositiveControl1K3W8G_TCTGTTGC-TAGTGTAG_S106_L001_R1.fastq.gz"
mv "$raw_data/run5/PositiveControlK3W8G_TGACCTCC-GACACTGA_S79_L001_R1.fastq.gz" "$raw_data/run5/PositiveControl2K3W8G_TGACCTCC-GACACTGA_S79_L001_R1.fastq.gz"
## Run 6
mv "$raw_data/run6/MockCommunityKBJ9G_TGACCTCC-CTAGAGCT_S79_L001_R1.fastq.gz" "$raw_data/run6/MockCommunity1KBJ9G_TGACCTCC-CTAGAGCT_S79_L001_R1.fastq.gz"
mv "$raw_data/run6/MockCommunityKBJ9G_TGACCTCC-GCTCTAGT_S80_L001_R1.fastq.gz" "$raw_data/run6/MockCommunity2KBJ9G_TGACCTCC-GCTCTAGT_S80_L001_R1.fastq.gz"
## Run 7
mv "$raw_data/run7/MockCommunityJC892_TCTGTTGC-ATATACAC_S2_L001_R1.fastq.gz" "$raw_data/run7/MockCommunity1JC892_TCTGTTGC-ATATACAC_S2_L001_R1.fastq.gz"
mv "$raw_data/run7/MockCommunityJC892_TGACCTCC-ACGACGTG_S1_L001_R1.fastq.gz" "$raw_data/run7/MockCommunity2JC892_TGACCTCC-ACGACGTG_S1_L001_R1.fastq.gz"
mv "$raw_data/run7/NegativeControlJC892_TCTGTTGC-CTAGAGCT_S6_L001_R1.fastq.gz" "$raw_data/run7/NegativeControl1JC892_TCTGTTGC-CTAGAGCT_S6_L001_R1.fastq.gz"
mv "$raw_data/run7/NegativeControlJC892_TGACCTCC-CGTCGCTA_S5_L001_R1.fastq.gz" "$raw_data/run7/NegativeControl2JC892_TGACCTCC-CGTCGCTA_S5_L001_R1.fastq.gz"
mv "$raw_data/run7/PositiveControlJC892_TCTGTTGC-CGTCGCTA_S4_L001_R1.fastq.gz" "$raw_data/run7/PositiveControl1JC892_TCTGTTGC-CGTCGCTA_S4_L001_R1.fastq.gz"
mv "$raw_data/run7/PositiveControlJC892_TGACCTCC-ATATACAC_S3_L001_R1.fastq.gz" "$raw_data/run7/PositiveControl2JC892_TGACCTCC-ATATACAC_S3_L001_R1.fastq.gz"
## Run 9
mv "$raw_data/run8/MockCommunityJC328_ACAAGGAG-CTAGAGCT_S1_L001_R1.fastq.gz" "$raw_data/run8/MockCommunity1JC328_ACAAGGAG-CTAGAGCT_S1_L001_R1.fastq.gz"
mv "$raw_data/run8/MockCommunityJC328_TGTAATTG-TGCGTACG_S2_L001_R1.fastq.gz" "$raw_data/run8/MockCommunity2JC328_TGTAATTG-TGCGTACG_S2_L001_R1.fastq.gz"
mv "$raw_data/run8/NegativeControlJC328_AATGGAGC-ACGACGTG_S6_L001_R1.fastq.gz" "$raw_data/run8/NegativeControl1JC328_AATGGAGC-ACGACGTG_S6_L001_R1.fastq.gz"
mv "$raw_data/run8/NegativeControlJC328_ACAAGGAG-GACACTGA_S5_L001_R1.fastq.gz" "$raw_data/run8/NegativeControl2JC328_ACAAGGAG-GACACTGA_S5_L001_R1.fastq.gz"
mv "$raw_data/run8/PositiveControlJC328_ACAAGGAG-GCTCTAGT_S3_L001_R1.fastq.gz" "$raw_data/run8/PositiveControl1JC328_ACAAGGAG-GCTCTAGT_S3_L001_R1.fastq.gz"
mv "$raw_data/run8/PositiveControlJC328_TGTAATTG-TAGTGTAG_S4_L001_R1.fastq.gz" "$raw_data/run8/PositiveControl2JC328_TGTAATTG-TAGTGTAG_S4_L001_R1.fastq.gz"
## Run 9
mv "$raw_data/run9/MockCommunityJY4PG_AATGTCCG-GACACTGA_S85_L001_R1.fastq.gz" "$raw_data/run9/MockCommunity1JY4PG_AATGTCCG-GACACTGA_S85_L001_R1.fastq.gz"
mv "$raw_data/run9/MockCommunityJY4PG_TATCAGGT-GACACTGA_S86_L001_R1.fastq.gz" "$raw_data/run9/MockCommunity2JY4PG_TATCAGGT-GACACTGA_S86_L001_R1.fastq.gz"
mv "$raw_data/run9/NegativeControlJY4PG_AATGTCCG-TAGTGTAG_S88_L001_R1.fastq.gz" "$raw_data/run9/NegativeControl1JY4PG_AATGTCCG-TAGTGTAG_S88_L001_R1.fastq.gz"
mv "$raw_data/run9/NegativeControlJY4PG_TATCAGGT-TAGTGTAG_S87_L001_R1.fastq.gz" "$raw_data/run9/NegativeControl2JY4PG_TATCAGGT-TAGTGTAG_S87_L001_R1.fastq.gz"
mv "$raw_data/run9/PositiveControlJY4PG_AATGTCCG-TGCGTACG_S89_L001_R1.fastq.gz" "$raw_data/run9/PositiveControl1JY4PG_AATGTCCG-TGCGTACG_S89_L001_R1.fastq.gz"
mv "$raw_data/run9/PositiveControlJY4PG_TATCAGGT-TGCGTACG_S90_L001_R1.fastq.gz" "$raw_data/run9/PositiveControl2JY4PG_TATCAGGT-TGCGTACG_S90_L001_R1.fastq.gz"
## Run 12
mv "$raw_data/run12/MockCommunityKCVJJ_GACTCTTG-GACACTGA_S95_L001_R1.fastq.gz" "$raw_data/run12/MockCommunity1KCVJJ_GACTCTTG-GACACTGA_S95_L001_R1.fastq.gz"
mv "$raw_data/run12/MockCommunityKCVJJ_GACTCTTG-GCTCTAGT_S93_L001_R1.fastq.gz" "$raw_data/run12/MockCommunity2KCVJJ_GACTCTTG-GCTCTAGT_S93_L001_R1.fastq.gz"
mv "$raw_data/run12/MockCommunityKCVJJ_GAGAGAAT-ACGACGTG_S94_L001_R1.fastq.gz" "$raw_data/run12/MockCommunity3KCVJJ_GAGAGAAT-ACGACGTG_S94_L001_R1.fastq.gz"
mv "$raw_data/run12/MockCommunityKCVJJ_GAGAGAAT-ATATACAC_S96_L001_R1.fastq.gz" "$raw_data/run12/MockCommunity4KCVJJ_GAGAGAAT-ATATACAC_S96_L001_R1.fastq.gz"
mv "$raw_data/run12/NegativeControlKCVJJ_GACTCTTG-TAGTGTAG_S97_L001_R1.fastq.gz" "$raw_data/run12/NegativeControl1KCVJJ_GACTCTTG-TAGTGTAG_S97_L001_R1.fastq.gz"
mv "$raw_data/run12/NegativeControlKCVJJ_GAGAGAAT-CTAGAGCT_S98_L001_R1.fastq.gz" "$raw_data/run12/NegativeControl2KCVJJ_GAGAGAAT-CTAGAGCT_S98_L001_R1.fastq.gz"
mv "$raw_data/run12/PositiveControlKCVJJ_GACTCTTG-TGCGTACG_S99_L001_R1.fastq.gz" "$raw_data/run12/PositiveControl1KCVJJ_GACTCTTG-TGCGTACG_S99_L001_R1.fastq.gz"
mv "$raw_data/run12/PositiveControlKCVJJ_GAGAGAAT-CGTCGCTA_S100_L001_R1.fastq.gz" "$raw_data/run12/PositiveControl2KCVJJ_GAGAGAAT-CGTCGCTA_S100_L001_R1.fastq.gz"
## Run 13
mv "$raw_data/run13/MockCommunityKG28C_AATGTCCG-ACGACGTG_S134_L001_R1.fastq.gz" "$raw_data/run13/MockCommunity1KG28C_AATGTCCG-ACGACGTG_S134_L001_R1.fastq.gz"
mv "$raw_data/run13/MockCommunityKG28C_AATGTCCG-ATATACAC_S137_L001_R1.fastq.gz" "$raw_data/run13/MockCommunity2KG28C_AATGTCCG-ATATACAC_S137_L001_R1.fastq.gz"
mv "$raw_data/run13/MockCommunityKG28C_ACAAGGAG-GACACTGA_S139_L001_R1.fastq.gz" "$raw_data/run13/MockCommunity3KG28C_ACAAGGAG-GACACTGA_S139_L001_R1.fastq.gz"
mv "$raw_data/run13/MockCommunityKG28C_ACAAGGAG-GCTCTAGT_S136_L001_R1.fastq.gz" "$raw_data/run13/MockCommunity4KG28C_ACAAGGAG-GCTCTAGT_S136_L001_R1.fastq.gz"
mv "$raw_data/run13/MockCommunityKG28C_GCAGGATA-GACACTGA_S135_L001_R1.fastq.gz" "$raw_data/run13/MockCommunity5KG28C_GCAGGATA-GACACTGA_S135_L001_R1.fastq.gz"
mv "$raw_data/run13/MockCommunityKG28C_GCAGGATA-TGCGTACG_S138_L001_R1.fastq.gz" "$raw_data/run13/MockCommunity6KG28C_GCAGGATA-TGCGTACG_S138_L001_R1.fastq.gz"
mv "$raw_data/run13/NegativeControlKG28C_AATGTCCG-CTAGAGCT_S140_L001_R1.fastq.gz" "$raw_data/run13/NegativeControl1KG28C_AATGTCCG-CTAGAGCT_S140_L001_R1.fastq.gz"
mv "$raw_data/run13/NegativeControlKG28C_ACAAGGAG-TAGTGTAG_S142_L001_R1.fastq.gz" "$raw_data/run13/NegativeControl2KG28C_ACAAGGAG-TAGTGTAG_S142_L001_R1.fastq.gz"
mv "$raw_data/run13/NegativeControlKG28C_GACTCTTG-ACGACGTG_S141_L001_R1.fastq.gz" "$raw_data/run13/NegativeControl3KG28C_GACTCTTG-ACGACGTG_S141_L001_R1.fastq.gz"
mv "$raw_data/run13/PositiveControlKG28C_AATGTCCG-CGTCGCTA_S143_L001_R1.fastq.gz" "$raw_data/run13/PositiveControl1KG28C_AATGTCCG-CGTCGCTA_S143_L001_R1.fastq.gz"
mv "$raw_data/run13/PositiveControlKG28C_ACAAGGAG-TGCGTACG_S145_L001_R1.fastq.gz" "$raw_data/run13/PositiveControl2KG28C_ACAAGGAG-TGCGTACG_S145_L001_R1.fastq.gz"
mv "$raw_data/run13/PositiveControlKG28C_GCAGGATA-TAGTGTAG_S144_L001_R1.fastq.gz" "$raw_data/run13/PositiveControl3KG28C_GCAGGATA-TAGTGTAG_S144_L001_R1.fastq.gz"
## Run 19
mv "$raw_data/run19/NegativeControlB3C4H_AGAGCCTACGTT_L001_R1.fastq.gz" "$raw_data/run19/NegativeControl1B3C4H_AGAGCCTACGTT_L001_R1.fastq.gz"
mv "$raw_data/run19/NegativeControlB3C4H_TATGCACCAGTG_L001_R1.fastq.gz" "$raw_data/run19/NegativeControl2B3C4H_TATGCACCAGTG_L001_R1.fastq.gz"
## Run 34
mv "$raw_data/run34/NegativeControlB39G7_AAGATGGATCAG_L001_R1.fastq.gz" "$raw_data/run34/NegativeControl1B39G7_AAGATGGATCAG_L001_R1.fastq.gz"
mv "$raw_data/run34/NegativeControlB39G7_AGCCGGCACATA_L001_R1.fastq.gz" "$raw_data/run34/NegativeControl2B39G7_AGCCGGCACATA_L001_R1.fastq.gz"
## Run 43
mv "$raw_data/run43/MockCommunityJC328_ACAAGGAG-CTAGAGCT_S1_L001_R1.fastq.gz" "$raw_data/run43/MockCommunity1JC328_ACAAGGAG-CTAGAGCT_S1_L001_R1.fastq.gz"
mv "$raw_data/run43/MockCommunityJC328_TGTAATTG-TGCGTACG_S2_L001_R1.fastq.gz" "$raw_data/run43/MockCommunity2JC328_TGTAATTG-TGCGTACG_S2_L001_R1.fastq.gz"
mv "$raw_data/run43/NegativeControlJC328_AATGGAGC-ACGACGTG_S6_L001_R1.fastq.gz" "$raw_data/run43/NegativeControl1JC328_AATGGAGC-ACGACGTG_S6_L001_R1.fastq.gz"
mv "$raw_data/run43/NegativeControlJC328_ACAAGGAG-GACACTGA_S5_L001_R1.fastq.gz" "$raw_data/run43/NegativeControl2JC328_ACAAGGAG-GACACTGA_S5_L001_R1.fastq.gz"
mv "$raw_data/run43/PositiveControlJC328_ACAAGGAG-GCTCTAGT_S3_L001_R1.fastq.gz" "$raw_data/run43/PositiveControl1JC328_ACAAGGAG-GCTCTAGT_S3_L001_R1.fastq.gz"
mv "$raw_data/run43/PositiveControlJC328_TGTAATTG-TAGTGTAG_S4_L001_R1.fastq.gz" "$raw_data/run43/PositiveControl2JC328_TGTAATTG-TAGTGTAG_S4_L001_R1.fastq.gz"
## Run 44
mv "$raw_data/run44/MockCommunityKP8LM_ACAAGGAG-GACACTGA_S111_L001_R1.fastq.gz" "$raw_data/run44/MockCommunity1KP8LM_ACAAGGAG-GACACTGA_S111_L001_R1.fastq.gz"
mv "$raw_data/run44/MockCommunityKP8LM_CGTCCGAA-CTAGAGCT_S112_L001_R1.fastq.gz" "$raw_data/run44/MockCommunity2KP8LM_CGTCCGAA-CTAGAGCT_S112_L001_R1.fastq.gz"
mv "$raw_data/run44/NegativeControlKP8LM_ACAAGGAG-TAGTGTAG_S113_L001_R1.fastq.gz" "$raw_data/run44/NegativeControl1KP8LM_ACAAGGAG-TAGTGTAG_S113_L001_R1.fastq.gz"
mv "$raw_data/run44/NegativeControlKP8LM_CGTCCGAA-GACACTGA_S114_L001_R1.fastq.gz" "$raw_data/run44/NegativeControl2KP8LM_CGTCCGAA-GACACTGA_S114_L001_R1.fastq.gz"
mv "$raw_data/run44/PositiveControlKP8LM_ACAAGGAG-TGCGTACG_S115_L001_R1.fastq.gz" "$raw_data/run44/PositiveControl1KP8LM_ACAAGGAG-TGCGTACG_S115_L001_R1.fastq.gz"
mv "$raw_data/run44/PositiveControlKP8LM_CGTCCGAA-GCTCTAGT_S116_L001_R1.fastq.gz" "$raw_data/run44/PositiveControl2KP8LM_CGTCCGAA-GCTCTAGT_S116_L001_R1.fastq.gz"

###############################################################################
## (3) Rename samples and move to run ID subdirectories #######################
###############################################################################

## Move fastq files to the raw data subdirectory
for f in $path/bpa*/*fastq.gz; do mv $f $01.Raw_data

## I will add an 's' to the start of each sample name because R doesn't like
## sample names that start with numbers.
for f in $raw_data/*.fastq.gz; do
  new_name=$(basename "$f" | sed 's/^/s/')
  mv "$f" "$raw_data/$new_name"
done

## Move fastq files to run subdirectories
for f in $raw_data/*; do
  case "$f" in
    *J9TWH*) mv "$f" "$raw_data/run1" ;;
    *J9GNL*) mv "$f" "$raw_data/run2" ;;
    *J9GGM*) mv "$f" "$raw_data/run3" ;;
    *J9FL8*) mv "$f" "$raw_data/run4" ;;
    *K3W8G*) mv "$f" "$raw_data/run5" ;;
    *KBJ9G*) mv "$f" "$raw_data/run6" ;;
    *JC892*) mv "$f" "$raw_data/run7" ;;
    *JC328*) mv "$f" "$raw_data/run8" ;;
    *JY4PG*) mv "$f" "$raw_data/run9" ;;
    *JRK29*) mv "$f" "$raw_data/run10" ;;
    *JRK2C*) mv "$f" "$raw_data/run11" ;;
    *KCVJJ*) mv "$f" "$raw_data/run12" ;;
    *KG28C*) mv "$f" "$raw_data/run13" ;;
    *AAVYC*) mv "$f" "$raw_data/run14" ;;
    *AC988*) mv "$f" "$raw_data/run15" ;;
    *A64HK*) mv "$f" "$raw_data/run16" ;;
    *A64JJ*) mv "$f" "$raw_data/run17" ;;
    *AC9CM*) mv "$f" "$raw_data/run18" ;;
    *B3C4H*) mv "$f" "$raw_data/run19" ;;
    *ANVM7*) mv "$f" "$raw_data/run20" ;;
    *A5P6W*) mv "$f" "$raw_data/run21" ;;
    *AAUF6*) mv "$f" "$raw_data/run22" ;;
    *AADNU*) mv "$f" "$raw_data/run23" ;;
    *ALBRY*) mv "$f" "$raw_data/run24" ;;
    *AN5TK*) mv "$f" "$raw_data/run25" ;;
    *BKG34*) mv "$f" "$raw_data/run26" ;;
    *AGEDA*) mv "$f" "$raw_data/run27" ;;
    *AGEFG*) mv "$f" "$raw_data/run28" ;;
    *AGECM*) mv "$f" "$raw_data/run29" ;;
    *B39FT*) mv "$f" "$raw_data/run30" ;;
    *BC8BY*) mv "$f" "$raw_data/run31" ;;
    *BC267*) mv "$f" "$raw_data/run32" ;;
    *C43MN*) mv "$f" "$raw_data/run33" ;;
    *B39G7*) mv "$f" "$raw_data/run34" ;;
    *BJTVT*) mv "$f" "$raw_data/run35" ;;
    *BY8DT*) mv "$f" "$raw_data/run36" ;;
    *BY245*) mv "$f" "$raw_data/run37" ;;
    *C5M89*) mv "$f" "$raw_data/run38" ;;
    *C43P8*) mv "$f" "$raw_data/run39" ;;
    *A64JH*) mv "$f" "$raw_data/run40" ;;
    *A79NR*) mv "$f" "$raw_data/run41" ;;
    *A90G8*) mv "$f" "$raw_data/run42" ;;
    *KP69D*) mv "$f" "$raw_data/run43" ;;
    *KP8LM*) mv "$f" "$raw_data/run44" ;;
  esac
done

## The samples are now ready for ITS extractions with ITSxpress.