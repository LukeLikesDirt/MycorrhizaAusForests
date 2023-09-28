

## Organise 16S data from the Australian Microbiome Initiative (AMI) dataset
##
## First, I will subset samples by seuqencing run prior to denoising to improve
## DADA2 error estimations because:
##      (1) The AMI dataset spans multiple sequencing runs.
##      (2) Different runs have unique error profiles.
##
## Second, I will organise subdirectories by sequencing run and will:
## (1) rename control samples to retain run ID, and (2) append 'control type'
## (i.e. mock community, positive, negative) to the ID so that both run ID
## and control type are retained when renaming all files within the DADA2
## step of this pipeline.
##
## Third, I will move test samples to the run subdirectories and add 's' to the
## beginning of sample names so that they are ready to be analysed in R, which
## requires sample names to start with a letter.
##
## NOTES:## Data within the subdirectories 'bpa_6775ac8e_20230923T1207' and
## 'bpa_de0d5f3a_20230925T0948' were downloaded from the Bioplatforms
## website (https://data.bioplatforms.com/organization/australian-microbiome)
## on 25/09/2023. The subdirectory contains all 'bpa_6775ac8e_20230923T1207'
## soil ITS test samplest from the Bioplatforms website, and the
## 'bpa_de0d5f3a_20230925T0948' subdirectory contains the control samples,
## using the'ITS' and 'Soil' filter tags, as well as the search term
## 'depth:0.0' to restrict samples to the upper (0-10cm) soil layer, and
## excluding sample from the lower (10-30cm) soi layer.
## Data within the subdirectories 'bpa_76ac789f_20230717T0504' and
## 'bpa_761e3497_20230717T0504' are from a spcific Govenrment of Victoria
## and were attained directly from the Govenrment of Victoria as they were
## not yet avialble on the Bioplatforms website. Samples from the lower soil
## layer were mannually removed.
## The 'bpa_8fc89f0d_20230926T0251' is a random riun that was for some reason
## not included in the bulk download whaen it should have been.
##
## NOTE: Samples from flow_id (run ID) K34JG were not found in Bioplatforms
## database when using the above search trems as there was not 'depth' variable
## in thier metadata. Howerver these data were from Antarctic soil samples and
## therefroe are not needed in my downstream analyses so I havnt included them 
## here.


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
path = '/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/data/AusMicrobiome/16S'

## Filter ID run names
data1 = read.csv(file.path(path, '/bpa_6775ac8e_20230923T1207/package_metadata/package_metadata_bpa_6775ac8e_20230923T1207_amdb-genomics-amplicon.csv'), header = T) %>% 
    # Select sample titles
    select(title) %>%  
    # Rename tile to run ID
    mutate(title = as.character(gsub('Australian Microbiome Amplicons 16S 102.100.100/', '', title))) %>%
    separate(title, c('sample_ID', 'run_ID'), sep = ' ') %>%
    # Reduce to unique run IDs
    select(run_ID) %>%
    unique() %>%
    as_tibble()
data2 = read.csv(file.path(path, '/bpa_6775ac8e_20230923T1207/package_metadata/package_metadata_bpa_6775ac8e_20230923T1207_base-genomics-amplicon.csv'), header = T) %>% 
    select(title) %>%
    mutate(title = as.character(gsub('BASE Amplicons 16S ', '', title))) %>%
    separate(title, c('sample_ID', 'run_ID'), sep = ' ') %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()
data3 = read.csv(file.path(path, '/bpa_76ac789f_20230717T0504/package_metadata/package_metadata_bpa_76ac789f_20230717T0504_amdb-genomics-amplicon.csv'), header = T) %>% 
    select(title) %>%
    mutate(title = as.character(gsub('Australian Microbiome Amplicons 16S 102.100.100/', '', title))) %>%
    separate(title, c('sample_ID', 'run_ID'), sep = ' ') %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()
data4 = read.csv(file.path(path, '/bpa_761e3497_20230717T0504/package_metadata/package_metadata_bpa_761e3497_20230717T0504_amdb-genomics-amplicon.csv'), header = T) %>% 
    select(title) %>%
    mutate(title = as.character(gsub('Australian Microbiome Amplicons 16S 102.100.100/', '', title))) %>%
    separate(title, c('sample_ID', 'run_ID'), sep = ' ') %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()
data5 = read.csv(file.path(path, '/bpa_8fc89f0d_20230926T0251/package_metadata/package_metadata_bpa_8fc89f0d_20230926T0251_amdb-genomics-amplicon.csv'), header = T) %>% 
    # Select sample titles
    select(title) %>%  
    # Rename tile to run ID
    mutate(title = as.character(gsub('Australian Microbiome Amplicons 16S 102.100.100/', '', title))) %>%
    separate(title, c('sample_ID', 'run_ID'), sep = ' ') %>%
    # Reduce to unique run IDs
    select(run_ID) %>%
    unique() %>%
    as_tibble()

## Combine datasets
bind_rows(data1, data2, data3, data4) %>%
    as_tibble() %>%
    print(n = Inf)
## There are 44 unique runs

## Here I match runs with control samples, which were provided as a separate file

## Load control sample data
data5 = read.csv(file.path(path, '/bpa_de0d5f3a_20230925T0948/package_metadata/package_metadata_bpa_de0d5f3a_20230925T0948_amdb-genomics-amplicon-control.csv'), header = T) %>% 
    select(title) %>%
    mutate(run_ID = as.character(gsub('Australian Microbiome Amplicons Control 16S ', '', title))) %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()
data6 = read.csv(file.path(path, '/bpa_de0d5f3a_20230925T0948/package_metadata/package_metadata_bpa_de0d5f3a_20230925T0948_base-genomics-amplicon-control.csv'), header = T) %>% 
    select(title) %>%
    mutate(run_ID = as.character(gsub('BASE Amplicons Control 16S ', '', title))) %>%
    select(run_ID) %>%
    unique() %>%
    as_tibble()

## Combine datasets
bind_rows(data5, data6) %>%
    as_tibble() %>%
    print(n = Inf)
## There are 91 unique runs because there are multiple control samples per run
## because control include samlpes from an distinct inicitive to the AMI.

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
path=/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/data/AusMicrobiome/16S   # Path to main data directory
mkdir $path/01.Raw_data         # Path to send raw data
raw_data=$path/01.Raw_data      # Define path raw data subdirectory
mkdir $raw_data/run{1..44}      # Make subdirectories for each run

## Before I start I will delete all files containing 'I1' and 'I2' read
## files, which are the indexes for the 'R1' and 'R2' reads, respectively.
## The 'find' function will search for all '*fastq.gz' files in subdirectories
find "$path" -type f -name '*fastq.gz' | while read -r f; do
    # Check if the file contains any specified patterns and then remove it
    if [[ "$f" == *I1* || "$f" == *I2* ]]; then
        rm "$f"
    fi
done

## Move control samples to the raw data subdirectory
for f in $path/bpa_de0d5f3a_20230925T0948/*fastq.gz; do mv $f $raw_data
done
## NOTE
## The control samples from the Victorian Government are in the
## 'bpa...20230717T0504' subdirectories with test samples, so I'll move
## them to the 'raw data' subdirectory manually.

## Move control samples to run ID subdirectories
for f in $raw_data/*; do
  case "$f" in
    *J6562*) mv "$f" "$raw_data/run1" ;;
    *J655F*) mv "$f" "$raw_data/run2" ;;
    *J654K*) mv "$f" "$raw_data/run3" ;;
    *J6KGL*) mv "$f" "$raw_data/run4" ;;
    *J7C52*) mv "$f" "$raw_data/run5" ;;
    *J7B92*) mv "$f" "$raw_data/run6" ;;
    *J7LRK*) mv "$f" "$raw_data/run7" ;;
    *K39YG*) mv "$f" "$raw_data/run8" ;;
    *KD7T7*) mv "$f" "$raw_data/run9" ;;
    *JC8RL*) mv "$f" "$raw_data/run10" ;;
    *J8GW4*) mv "$f" "$raw_data/run11" ;;
    *JK5PT*) mv "$f" "$raw_data/run12" ;;
    *KG38F*) mv "$f" "$raw_data/run13" ;;
    *KFPWM*) mv "$f" "$raw_data/run14" ;;
    *A815D*) mv "$f" "$raw_data/run15" ;;
    *AC9BY*) mv "$f" "$raw_data/run16" ;;
    *A819B*) mv "$f" "$raw_data/run17" ;;
    *A5YK6*) mv "$f" "$raw_data/run18" ;;
    *AC9E5*) mv "$f" "$raw_data/run19" ;;
    *B3C7L*) mv "$f" "$raw_data/run20" ;;
    *AGBMG*) mv "$f" "$raw_data/run21" ;;
    *A801W*) mv "$f" "$raw_data/run22" ;;
    *AAHKD*) mv "$f" "$raw_data/run23" ;;
    *AAF2V*) mv "$f" "$raw_data/run24" ;;
    *A815N*) mv "$f" "$raw_data/run25" ;;
    *AGEFE*) mv "$f" "$raw_data/run26" ;;
    *AGEK1*) mv "$f" "$raw_data/run27" ;;
    *AGEJK*) mv "$f" "$raw_data/run28" ;;
    *B3BDY*) mv "$f" "$raw_data/run29" ;;
    *B3L2P*) mv "$f" "$raw_data/run30" ;;
    *BDGJV*) mv "$f" "$raw_data/run31" ;;
    *BDC3D*) mv "$f" "$raw_data/run32" ;;
    *BM3J6*) mv "$f" "$raw_data/run33" ;;
    *BWYGW*) mv "$f" "$raw_data/run34" ;;
    *C27WK*) mv "$f" "$raw_data/run35" ;;
    *C444T*) mv "$f" "$raw_data/run36" ;;
    *C5WDL*) mv "$f" "$raw_data/run37" ;;
    *C3WMR*) mv "$f" "$raw_data/run38" ;;
    *A5YJC*) mv "$f" "$raw_data/run39" ;;
    *A8190*) mv "$f" "$raw_data/run40" ;;
    *A5K1H*) mv "$f" "$raw_data/run41" ;;
    *KN482*) mv "$f" "$raw_data/run42" ;;
    *KN48L*) mv "$f" "$raw_data/run43" ;;
    *K39YG*) mv "$f" "$raw_data/run44" ;;
  esac
done

## Remove control samples that are not in the AMI dataset that I am working
## with.
for f in $raw_data/*.fastq.gz; do
    rm "$f"
done

## There are a few 'Undetermined' control samples that I will remove.
## There are also multiple runs where file size for control samples = 15kb,
## which indicates failed sequencing: I will remove these files using a
## threshold of 100kb
find "$raw_data" -type f -name '*fastq.gz' | while read -r f; do
    # Check if the file contains the "Undetermined" pattern or files that are smaller than 100KB (where 1KB = 1024 bytes)
    if [[ "$f" == *Undetermined* || $(stat -c %s "$f") -lt 102400 ]]; then
        rm "$f"
    fi
done

## File name prefixes change depending on the date of data intake into the
## database. I'll nest loops to rename all control samples

## Rename the prefix of mock community samples
## Start from complex to simple names to avoid renaming the wrong files
for run_dir in $raw_data/run{1..44}; do
  for prefix in 'ATCC1002MOCK_16S_' 'Bac_mock_community_16S_AGRF_GAATAGAGCCAA_' 'Bac_mock_community_16S_AGRF_GCGGCAATTACG_' 'Bac_mock_community_16S_AGRF_GTACGTGGGATC_' 'Bac_mock_community_16S_AGRF_' 'Bac_mock_community_16S_UNSW_GCACGACAACAC_' 'Bac_mock_community_16S_UNSW_GCGGCAATTACG_' 'Bac_mock_community_16S_UNSW_CGAGAAGAGAAC_' 'Bac_mock_community_16S_UNSW_GAATAGAGCCAA_' 'Bac_mock_community_16S_UNSW_'; do
    for f in "$run_dir"/*"$prefix"*; do
      if [ -e "$f" ]; then
        new_name=$(basename "$f" | sed "s/$prefix/MockCommunity/")
        mv "$f" "$run_dir/$new_name"
      fi
    done
  done
done

## Now I will rename negative control samples:
for run_dir in $raw_data/run{1..44}; do
  for prefix in 'NEG_1_16S_AGRF_GGAGACAAGGGA_' 'NEG1_16S_AGRF_' 'NEG2_16S_AGRF_' 'Neg1_16S_AGRF_' 'Neg2_16S_AGRF_' 'No_Template_Control_16S_', 'blank_16S_'; do
    for f in "$run_dir"/*"$prefix"*; do
      if [ -e "$f" ]; then
        new_name=$(basename "$f" | sed "s/$prefix/NegativeControl/")
        mv "$f" "$run_dir/$new_name"
      fi
    done
  done
done

## Now I will rename positive control samples:
for run_dir in $raw_data/run{1..44}; do
  for prefix in 'Soil_DNA_16S_AGRF_CCACCTACTCCA_' 'Soil_DNA_16S_AGRF_GGTGTCTATTGT_' 'Soil_DNA_16S_AGRF_CGAGAAGAGAAC_' 'Soil_DNA_16S_AGRF_' 'Soil_DNA_16S_UNSW_CGAGAAGAGAAC_' 'Soil_DNA_16S_UNSW_GCGGCAATTACG_' 'Soil_DNA_16S_UNSW_GGAAACCACCAC_' 'Soil_DNA_16S_UNSW_CCACCTACTCCA_' 'Soil_DNA_16S_'; do
    for f in "$run_dir"/*"$prefix"*; do
      if [ -e "$f" ]; then
        new_name=$(basename "$f" | sed "s/$prefix/PositiveControl/")
        mv "$f" "$run_dir/$new_name"
      fi
    done
  done
done

## Lastly I will rename the controls that have been run in duplicate,
## triplicate or quadruplicate.
## I'm not sure how I could loop this step so I will run these one by one.
## Run 2
find "$raw_data/run2/" -type f -name 'NegativeControlJ655F_AAGAGATG-TCTACACT_S29_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJ655F_AAGAGATG-TCTACACT_S29_L001_/NegativeControl1J655F_AAGAGATG-TCTACACT_S29_L001_}"' {} \;
find "$raw_data/run2/" -type f -name 'NegativeControlJ655F_CGATCCGT-CGATCTAC_S3_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJ655F_CGATCCGT-CGATCTAC_S3_L001_/NegativeControl2J655F_CGATCCGT-CGATCTAC_S3_L001_}"' {} \;
find "$raw_data/run2/" -type f -name 'NegativeControlJ655F_TCGAGGAC-CTAGTATG_S9_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJ655F_TCGAGGAC-CTAGTATG_S9_L001_/NegativeControl3J655F_TCGAGGAC-CTAGTATG_S9_L001_}"' {} \;
## Run 3
find "$raw_data/run3/" -type f -name 'NegativeControlJ654K_CACGCCAT-TCTACACT_S27_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJ654K_CACGCCAT-TCTACACT_S27_/NegativeControl1J654K_CACGCCAT-TCTACACT_S27_}"' {} \;
find "$raw_data/run3/" -type f -name 'NegativeControlJ654K_CAGGCGTA-ACGCGTGA_S29_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJ654K_CAGGCGTA-ACGCGTGA_S29_L001_/NegativeControl2J654K_CAGGCGTA-ACGCGTGA_S29_L001_}"' {} \;
find "$raw_data/run3/" -type f -name 'NegativeControlJ654K_CGATCCGT-CGATCTAC_S3_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJ654K_CGATCCGT-CGATCTAC_S3_L001_/NegativeControl3J654K_CGATCCGT-CGATCTAC_S3_L001_}"' {} \;
## Run 8
find "$raw_data/run8/" -type f -name 'MockCommunityK39YG_AATGCCTC-GTCTAGTG_S104_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityK39YG_AATGCCTC-GTCTAGTG_S104_L001_/MockCommunity1K39YG_AATGCCTC-GTCTAGTG_S104_L001_}"' {} \;
find "$raw_data/run8/" -type f -name 'MockCommunityK39YG_ATGAGACT-TCTACACT_S102_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityK39YG_ATGAGACT-TCTACACT_S102_L001_/MockCommunity2K39YG_ATGAGACT-TCTACACT_S102_L001_}"' {} \;
find "$raw_data/run8/" -type f -name 'MockCommunityK39YG_CAAGCATG-CGATCTAC_S103_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityK39YG_CAAGCATG-CGATCTAC_S103_L001_/MockCommunity3K39YG_CAAGCATG-CGATCTAC_S103_L001_}"' {} \;
find "$raw_data/run8/" -type f -name 'NegativeControlK39YG_AATGCCTC-GATAGCGT_S107_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlK39YG_AATGCCTC-GATAGCGT_S107_L001_/NegativeControl1K39YG_AATGCCTC-GATAGCGT_S107_L001_}"' {} \;
find "$raw_data/run8/" -type f -name 'NegativeControlK39YG_CAAGCATG-GTCTAGTG_S106_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlK39YG_CAAGCATG-GTCTAGTG_S106_L001_/NegativeControl2K39YG_CAAGCATG-GTCTAGTG_S106_L001_}"' {} \;
find "$raw_data/run8/" -type f -name 'NegativeControlK39YG_GAATCTTC-ACGCGTGA_S105_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlK39YG_GAATCTTC-ACGCGTGA_S105_L001_/NegativeControl3K39YG_GAATCTTC-ACGCGTGA_S105_L001_}"' {} \;
find "$raw_data/run8/" -type f -name 'PositiveControlK39YG_AATGCCTC-CTAGTATG_S109_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlK39YG_AATGCCTC-CTAGTATG_S109_L001_/PositiveControl1K39YG_AATGCCTC-CTAGTATG_S109_L001_}"' {} \;
find "$raw_data/run8/" -type f -name 'PositiveControlK39YG_CAAGCATG-ACGACGTG_S110_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlK39YG_CAAGCATG-ACGACGTG_S110_L001_/PositiveControl2K39YG_CAAGCATG-ACGACGTG_S110_L001_}"' {} \;
find "$raw_data/run8/" -type f -name 'PositiveControlK39YG_GAATCTTC-AAGCAGCA_S108_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlK39YG_GAATCTTC-AAGCAGCA_S108_L001_/PositiveControl3K39YG_GAATCTTC-AAGCAGCA_S108_L001_}"' {} \;
## Run 10
find "$raw_data/run10/" -type f -name 'MockCommunityJC8RL_TACACGAT-ACGCGTGA_S4_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityJC8RL_TACACGAT-ACGCGTGA_S4_L001_/MockCommunity1JC8RL_TACACGAT-ACGCGTGA_S4_L001_}"' {} \;
find "$raw_data/run10/" -type f -name 'MockCommunityJC8RL_TCCGAATT-AAGCAGCA_S3_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityJC8RL_TCCGAATT-AAGCAGCA_S3_L001_/MockCommunity2JC8RL_TCCGAATT-AAGCAGCA_S3_L001_}"' {} \;
find "$raw_data/run10/" -type f -name 'NegativeControlJC8RL_TACACGAT-ACGACGTG_S6_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJC8RL_TACACGAT-ACGACGTG_S6_L001_/NegativeControl1JC8RL_TACACGAT-ACGACGTG_S6_L001_}"' {} \;
find "$raw_data/run10/" -type f -name 'NegativeControlJC8RL_TCCGAATT-CGATCTAC_S5_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJC8RL_TCCGAATT-CGATCTAC_S5_L001_/NegativeControl2JC8RL_TCCGAATT-CGATCTAC_S5_L001_}"' {} \;
find "$raw_data/run10/" -type f -name 'PositiveControlJC8RL_TACACGAT-CGATCTAC_S2_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlJC8RL_TACACGAT-CGATCTAC_S2_L001_/PositiveControl1JC8RL_TACACGAT-CGATCTAC_S2_L001_}"' {} \;
find "$raw_data/run10/" -type f -name 'PositiveControlJC8RL_TCCGAATT-ACGCGTGA_S1_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlJC8RL_TCCGAATT-ACGCGTGA_S1_L001_/PositiveControl2JC8RL_TCCGAATT-ACGCGTGA_S1_L001_}"' {} \;
## Run 12
find "$raw_data/run12/" -type f -name 'MockCommunityJK5PT_ACTGACAG-ACGCGTGA_S8_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityJK5PT_ACTGACAG-ACGCGTGA_S8_L001_/MockCommunity1JK5PT_ACTGACAG-ACGCGTGA_S8_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'MockCommunityJK5PT_AGCTGTTG-CTAGTATG_S5_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityJK5PT_AGCTGTTG-CTAGTATG_S5_L001_/MockCommunity2JK5PT_AGCTGTTG-CTAGTATG_S5_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'MockCommunityJK5PT_GTCAATTG-CTAGTATG_S7_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityJK5PT_GTCAATTG-CTAGTATG_S7_L001_/MockCommunity3JK5PT_GTCAATTG-CTAGTATG_S7_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'MockCommunityJK5PT_TAGGAACT-ACGCGTGA_S6_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityJK5PT_TAGGAACT-ACGCGTGA_S6_L001_/MockCommunity4JK5PT_TAGGAACT-ACGCGTGA_S6_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'NegativeControlJK5PT_ACTGACAG-ACGACGTG_S12_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJK5PT_ACTGACAG-ACGACGTG_S12_L001_/NegativeControl1JK5PT_ACTGACAG-ACGACGTG_S12_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'NegativeControlJK5PT_AGCTGTTG-TCTACACT_S9_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJK5PT_AGCTGTTG-TCTACACT_S9_L001_/NegativeControl2JK5PT_AGCTGTTG-TCTACACT_S9_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'NegativeControlJK5PT_GTCAATTG-TCTACACT_S11_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJK5PT_GTCAATTG-TCTACACT_S11_L001_/NegativeControl3JK5PT_GTCAATTG-TCTACACT_S11_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'NegativeControlJK5PT_TAGGAACT-ACGACGTG_S10_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlJK5PT_TAGGAACT-ACGACGTG_S10_L001_/NegativeControl4JK5PT_TAGGAACT-ACGACGTG_S10_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'PositiveControlJK5PT_ACTGACAG-CGATCTAC_S4_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlJK5PT_ACTGACAG-CGATCTAC_S4_L001_/PositiveControl1JK5PT_ACTGACAG-CGATCTAC_S4_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'PositiveControlJK5PT_AGCTGTTG-GATAGCGT_S1_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlJK5PT_AGCTGTTG-GATAGCGT_S1_L001_/PositiveControl2JK5PT_AGCTGTTG-GATAGCGT_S1_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'PositiveControlJK5PT_GTCAATTG-GATAGCGT_S3_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlJK5PT_GTCAATTG-GATAGCGT_S3_L001_/PositiveControl3JK5PT_GTCAATTG-GATAGCGT_S3_L001_}"' {} \;
find "$raw_data/run12/" -type f -name 'PositiveControlJK5PT_TAGGAACT-CGATCTAC_S2_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlJK5PT_TAGGAACT-CGATCTAC_S2_L001_/PositiveControl4JK5PT_TAGGAACT-CGATCTAC_S2_L001_}"' {} \;
## Run 14
find "$raw_data/run14/" -type f -name 'MockCommunityKFPWM_AATGCCTC-CTAGTATG_S134_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityKFPWM_AATGCCTC-CTAGTATG_S134_L001_/MockCommunity1KFPWM_AATGCCTC-CTAGTATG_S134_L001_}"' {} \;
find "$raw_data/run14/" -type f -name 'MockCommunityKFPWM_CGAGAAGA-AAGCAGCA_S135_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityKFPWM_CGAGAAGA-AAGCAGCA_S135_L001_/MockCommunity2KFPWM_CGAGAAGA-AAGCAGCA_S135_L001_}"' {} \;
find "$raw_data/run14/" -type f -name 'NegativeControlKFPWM_AATGCCTC-TCTACACT_S136_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlKFPWM_AATGCCTC-TCTACACT_S136_L001_/NegativeControl1KFPWM_AATGCCTC-TCTACACT_S136_L001_}"' {} \;
find "$raw_data/run14/" -type f -name 'NegativeControlKFPWM_CGAGAAGA-CGATCTAC_S137_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlKFPWM_CGAGAAGA-CGATCTAC_S137_L001_/NegativeControl2KFPWM_CGAGAAGA-CGATCTAC_S137_L001_}"' {} \;
find "$raw_data/run14/" -type f -name 'PositiveControlKFPWM_AATGCCTC-GATAGCGT_S138_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlKFPWM_AATGCCTC-GATAGCGT_S138_L001_/PositiveControl1KFPWM_AATGCCTC-GATAGCGT_S138_L001_}"' {} \;
find "$raw_data/run14/" -type f -name 'PositiveControlKFPWM_CGAGAAGA-ACGCGTGA_S139_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlKFPWM_CGAGAAGA-ACGCGTGA_S139_L001_/PositiveControl2KFPWM_CGAGAAGA-ACGCGTGA_S139_L001_}"' {} \;
## Run 20
find "$raw_data/run20/" -type f -name 'NegativeControlB3C7L_ACGGGACATGCT_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlB3C7L_ACGGGACATGCT_L001_/NegativeControl1B3C7L_ACGGGACATGCT_L001_}"' {} \;
find "$raw_data/run20/" -type f -name 'NegativeControlB3C7L_CTTCGGCAGAAT_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlB3C7L_CTTCGGCAGAAT_L001_/NegativeControl2B3C7L_CTTCGGCAGAAT_L001_}"' {} \;
## Run 29
find "$raw_data/run29/" -type f -name 'NegativeControlB3BDY_GAATAGAGCCAA_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlB3BDY_GAATAGAGCCAA_L001_/NegativeControl1B3BDY_GAATAGAGCCAA_L001_}"' {} \;
find "$raw_data/run29/" -type f -name 'NegativeControlB3BDY_GTACGTGGGATC_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlB3BDY_GTACGTGGGATC_L001_/NegativeControl2B3BDY_GTACGTGGGATC_L001_}"' {} \;
## Run 30
find "$raw_data/run30/" -type f -name 'NegativeControlB3L2P_ATCGGCGTTACA_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlB3L2P_ATCGGCGTTACA_L001_/NegativeControl1B3L2P_ATCGGCGTTACA_L001_}"' {} \;
find "$raw_data/run30/" -type f -name 'NegativeControlB3L2P_CAGCGGTGACAT_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlB3L2P_CAGCGGTGACAT_L001_/NegativeControl2B3L2P_CAGCGGTGACAT_L001_}"' {} \;
## Run 42
find "$raw_data/run42/" -type f -name 'MockCommunityKN482_CTAGCGAA-GTCTAGTG_S109_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityKN482_CTAGCGAA-GTCTAGTG_S109_L001_/MockCommunity1KN482_CTAGCGAA-GTCTAGTG_S109_L001_}"' {} \;
find "$raw_data/run42/" -type f -name 'MockCommunityKN482_GCGGCAAT-CTAGTATG_S110_L001_R1*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityKN482_GCGGCAAT-CTAGTATG_S110_L001_R1/MockCommunity2KN482_GCGGCAAT-CTAGTATG_S110_L001_R1}"' {} \;
find "$raw_data/run42/" -type f -name 'MockCommunityKN482_TACACGAT-AAGCAGCA_S111_L001_*' \
    -exec bash -c 'mv "$0" "${0/MockCommunityKN482_TACACGAT-AAGCAGCA_S111_L001_/MockCommunity3KN482_TACACGAT-AAGCAGCA_S111_L001_}"' {} \;
find "$raw_data/run42/" -type f -name 'NegativeControlKN482_CTAGCGAA-GATAGCGT_S112_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlKN482_CTAGCGAA-GATAGCGT_S112_L001_/NegativeControl1KN482_CTAGCGAA-GATAGCGT_S112_L001_}"' {} \;
find "$raw_data/run42/" -type f -name 'NegativeControlKN482_GCGGCAAT-TCTACACT_S113_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlKN482_GCGGCAAT-TCTACACT_S113_L001_/NegativeControl2KN482_GCGGCAAT-TCTACACT_S113_L001_}"' {} \;
find "$raw_data/run42/" -type f -name 'NegativeControlKN482_TACACGAT-ACGACGTG_S114_L001_*' \
    -exec bash -c 'mv "$0" "${0/NegativeControlKN482_TACACGAT-ACGACGTG_S114_L001_/NegativeControl3KN482_TACACGAT-ACGACGTG_S114_L001_}"' {} \;
find "$raw_data/run42/" -type f -name 'PositiveControlKN482_CTAGCGAA-CTAGTATG_S115_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlKN482_CTAGCGAA-CTAGTATG_S115_L001_/PositiveControl1KN482_CTAGCGAA-CTAGTATG_S115_L001_}"' {} \;
find "$raw_data/run42/" -type f -name 'PositiveControlKN482_GCGGCAAT-GATAGCGT_S116_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlKN482_GCGGCAAT-GATAGCGT_S116_L001_/PositiveControl2KN482_GCGGCAAT-GATAGCGT_S116_L001_}"' {} \;
find "$raw_data/run42/" -type f -name 'PositiveControlKN482_TACACGAT-ACGCGTGA_S117_L001_*' \
    -exec bash -c 'mv "$0" "${0/PositiveControlKN482_TACACGAT-ACGCGTGA_S117_L001_/PositiveControl3KN482_TACACGAT-ACGCGTGA_S117_L001_}"' {} \;

###############################################################################
## (3) Rename samples and move to run ID subdirectories #######################
###############################################################################

## Move fastq files to the raw data subdirectory
# Use 'find' to search for all '*fastq.gz' files in multiple subdirectories
for f in $path/bpa_*/*fastq.gz; do mv $f $raw_data
done

## I will add an 's' to the start of each sample name because R doesn't like
## sample names that start with numbers.
for f in $raw_data/*.fastq.gz; do
  new_name=$(basename "$f" | sed 's/^/s/')
  mv "$f" "$raw_data/$new_name"
done

## Move fastq files to run subdirectories
for f in $raw_data/*; do
  case "$f" in
    *J6562*) mv "$f" "$raw_data/run1" ;;
    *J655F*) mv "$f" "$raw_data/run2" ;;
    *J654K*) mv "$f" "$raw_data/run3" ;;
    *J6KGL*) mv "$f" "$raw_data/run4" ;;
    *J7C52*) mv "$f" "$raw_data/run5" ;;
    *J7B92*) mv "$f" "$raw_data/run6" ;;
    *J7LRK*) mv "$f" "$raw_data/run7" ;;
    *K39YG*) mv "$f" "$raw_data/run8" ;;
    *KD7T7*) mv "$f" "$raw_data/run9" ;;
    *JC8RL*) mv "$f" "$raw_data/run10" ;;
    *J8GW4*) mv "$f" "$raw_data/run11" ;;
    *JK5PT*) mv "$f" "$raw_data/run12" ;;
    *KG38F*) mv "$f" "$raw_data/run13" ;;
    *KFPWM*) mv "$f" "$raw_data/run14" ;;
    *A815D*) mv "$f" "$raw_data/run15" ;;
    *AC9BY*) mv "$f" "$raw_data/run16" ;;
    *A819B*) mv "$f" "$raw_data/run17" ;;
    *A5YK6*) mv "$f" "$raw_data/run18" ;;
    *AC9E5*) mv "$f" "$raw_data/run19" ;;
    *B3C7L*) mv "$f" "$raw_data/run20" ;;
    *AGBMG*) mv "$f" "$raw_data/run21" ;;
    *A801W*) mv "$f" "$raw_data/run22" ;;
    *AAHKD*) mv "$f" "$raw_data/run23" ;;
    *AAF2V*) mv "$f" "$raw_data/run24" ;;
    *A815N*) mv "$f" "$raw_data/run25" ;;
    *AGEFE*) mv "$f" "$raw_data/run26" ;;
    *AGEK1*) mv "$f" "$raw_data/run27" ;;
    *AGEJK*) mv "$f" "$raw_data/run28" ;;
    *B3BDY*) mv "$f" "$raw_data/run29" ;;
    *B3L2P*) mv "$f" "$raw_data/run30" ;;
    *BDGJV*) mv "$f" "$raw_data/run31" ;;
    *BDC3D*) mv "$f" "$raw_data/run32" ;;
    *BM3J6*) mv "$f" "$raw_data/run33" ;;
    *BWYGW*) mv "$f" "$raw_data/run34" ;;
    *C27WK*) mv "$f" "$raw_data/run35" ;;
    *C444T*) mv "$f" "$raw_data/run36" ;;
    *C5WDL*) mv "$f" "$raw_data/run37" ;;
    *C3WMR*) mv "$f" "$raw_data/run38" ;;
    *A5YJC*) mv "$f" "$raw_data/run39" ;;
    *A8190*) mv "$f" "$raw_data/run40" ;;
    *A5K1H*) mv "$f" "$raw_data/run41" ;;
    *KN482*) mv "$f" "$raw_data/run42" ;;
    *KN48L*) mv "$f" "$raw_data/run43" ;;
  esac
done

## The samples are now ready for qulity trimming.