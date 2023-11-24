
# Script: Harmonise plant databases according to the World Flora Online

# Table of contents:
#   (1) Prepare the World Flora Online backbone database
#   (2) Prepare the root trait and tree classification databases
#       (2a) NodDB
#       (2b) FungalRoot
#       (2c) GlobalTreeSearch
#   (3) Prepare the plant distribution databases
#       (3a) BiomassPlotLib 
#       (3b) HAVPlot

# Required packages
require(readxl)
require(WorldFlora)
require(tidyverse)

#### (1) The World Flora Online backbone database ##############################

# Download the WFO backbone
WFO.download(WFO.url =
               paste0("https://files.worldfloraonline.org/files/WFO_Backbone/",
                      "_WFOCompleteBackbone/WFO_Backbone.zip"),
             save.dir = "data/WorldFloraOnline", WFO.remember = TRUE)
# Remove the zip file
file.remove("data/WorldFloraOnline/WFO_Backbone.zip")

# Read in the WFO backbone
WFO.remember()
nrow(WFO.data) # 1,497,586

# Check taxon ranks
unique(WFO.data$taxonRank)

# Filter the WFO to genus and species subset.
# Replace " × " with " ×" for scientific names of hybrid species, and "× " with
# "×" for scientific names of hybrid genera in the WFO – this already needed to
# be done for the World Checklist of Vascular Plants as shown in this Rpub: 
# https://rpubs.com/Roeland-KINDT/812716 - see Section 3.1.

# Species dataset
WFO_species <- WFO.data %>%
  filter(taxonRank == "species") %>%
  mutate(scientificName = gsub(" × ", " ×", scientificName)) %>%
  mutate(scientificName = gsub("× ", "×", scientificName)) %>%
  glimpse()
nrow(WFO_species) # 1,141,953
# Check the mutation
WFO_species %>%
  filter(grepl(" × ", scientificName))
WFO_species %>%
  filter(grepl("× ", scientificName))
WFO_species %>%
  filter(grepl(" ×", scientificName)) %>%
  select(scientificName) %>%
  as_tibble() %>%
  print(n = Inf)
WFO_species %>%
  filter(grepl("(?<!\\s)×", scientificName, perl = TRUE)) %>%
  select(scientificName) %>%
  as_tibble() %>%
  print(n = Inf)

# Genus dataset
WFO_genus <- WFO.data %>%
  filter(taxonRank == "genus") %>%
  mutate(scientificName = gsub("× ", "×", scientificName)) %>%
  glimpse()
nrow(WFO_genus) # 33,090
# Check the mutation
WFO_genus %>%
  filter(grepl("× ", scientificName)) %>%
  glimpse()
WFO_genus %>%
  filter(grepl("(?<!\\s)×", scientificName, perl = TRUE)) %>%
  select(scientificName) %>%
  as_tibble() %>%
  print(n = Inf)

# Save the datasets
write_excel_csv(WFO_species, "data/WorldFloraOnline/WFO_species.csv")
write_excel_csv(WFO_genus, "data/WorldFloraOnline/WFO_genus.csv")

# Remove all objects from the global environment
rm(list = ls())

#### (2) Root trait and tree classification databases ##########################

# Read in the WFO datasets
WFO_species <- read.csv("data/WorldFloraOnline/WFO_species.csv")
WFO_genus <- read.csv("data/WorldFloraOnline/WFO_genus.csv")

##### (2a) NodDB #####

# Download the NodDB database:
URL = "https://files.plutof.ut.ee/doi/FA/CE/FACE809C07AC1CD72932E878E10316DB9B03F8FD5B88C721FFB3F7EAC3F1252F.xlsx" 
file_destination <- "data/NodDB/NodDB.xlsx"
download.file(URL, file_destination)

# Read in the NodDB database:
nod_database <- 
  read_excel("data/NodDB/NodDB.xlsx", sheet = "Table S1 NodDB v1.3a",
             skip = 1, col_names = TRUE) %>%
  select(family_original = family, genus_original = genus,
         nodulation_type = "Consensus estimate") %>%
  mutate(
    nodulation_status = case_when(
      nodulation_type == "None" ~ "non_N_fixer",
      TRUE ~ "N_fixer" # Default value for nodulation_status
    ),
    nodulation_type = case_when(
      nodulation_type == "None" ~ "none",
      nodulation_type == "likely_Rhizobia" ~ "Rhizobia",
      nodulation_type == "unlikely_Rhizobia" ~ "undetermined",
      nodulation_type == "unlikely_Frankia" ~ "undetermined",
      nodulation_type == "likely_present" ~ "undetermined",
      nodulation_type == "Present" ~ "undetermined",
      TRUE ~ nodulation_type
    )
  ) %>%
  na.omit(genus_original) %>%
  as.data.frame() %>%
  glimpse()
# There are 824 genera in the NodDB dataset:
unique(nod_database$genus) %>%
  length()

# Match the tree genera in NodBD to WFO:
nod_database_WFO <- 
  WFO.one(WFO.match.fuzzyjoin(spec.data = nod_database,
                              WFO.data = WFO_genus,
                              spec.name = "genus_original"),
          verbose = FALSE) %>%
  glimpse()

# Direct matches: 823 out of 824
nod_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 1 out of 824
nod_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(genus_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble()

# Save the harmonised NodDB
nod_database_WFO %>%
  select(family, genus = scientificName, nodulation_status, nodulation_type) %>%
  write_excel_csv("data/NodDB/harmonised_NodDB.csv")

##### (2b) FungalRoot #####

# Download the FungalRoot database:
URL = "https://files.plutof.ut.ee/public/orig/77/EE/77EEAF156EA8B35CAE58E18E644ADBDCFC90EB439DE9BA012D4FCAB3ECD3FCAB.zip"
file_destination <- "data/FungalRoot/FungalRoot.zip"
download.file(URL, file_destination)

# Unzip the downdload
unzip(destfile, exdir = "data/FungalRoot/")

# Remove the zip file
file.remove(file_destination)

# Read in the FungalRoot mycorrhizal type database:
fungal_root_database <- read.csv("data/FungalRoot/FungalRoot.csv") %>%
  select(genus_original = Genus, mycorrhizal_type = Mycorrhizal.type) %>%
  # Extract the first word from genus values with multiple words
  mutate(genus_original = str_extract(genus_original, "^\\S+")) %>%
  distinct(genus_original, .keep_all = TRUE) %>%
  glimpse()

# There are 14,527 unique genera in the FungalRoot mycorrhizal type database:
nrow(fungal_root_database)

# Match the tree genera in FungalRoot mycorrhizal type to WFO:
fungal_root_database_WFO <- 
  WFO.one(WFO.match.fuzzyjoin(spec.data = fungal_root_database,
                              WFO.data = WFO_genus,
                              spec.name = "genus_original"),
          verbose = FALSE) %>%
  glimpse()

          # Direct matches: 14,290 out of 14,541
fungal_root_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 229 out of 14,541
fungal_root_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(genus_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# 8 unmatched
fungal_root_database_WFO %>% 
  filter(Matched == FALSE) %>%
  select(genus_original, taxonRank) %>%
  as_tibble() %>%
  print(n = Inf)
fungal_root_database_WFO %>%
  select(genus_original)
# Save the harmonised FungalRoot mycorrhizal type database
fungal_root_database_WFO %>%
  # Remove the unmatched genera
  filter(Matched == TRUE) %>%
  select(family, genus = scientificName, mycorrhizal_type) %>%
  write_excel_csv("data/FungalRoot/harmonised_fungal_root.csv")

# Read in the FungalRoot occurrences database:
fungal_root_occurrences <- read.csv("data/FungalRoot/occurrences.csv") %>%
  select(ID, family_original = family, genus_original = genus,
         species_original = scientificName, taxonRank) %>%
  ###
  ### This code will handle variation in hybrid annotations in the
  ### FungalRoot database but I will stick with genus for now:
  ### #Handle variations in the use of '×' and 'x' to denote hybrid species.
  ### mutate(species_original = gsub(" x ", " × ", species_original)) %>%
  ### mutate(species_original = gsub(" × ", " ×", species_original)) %>%
  ###
  # Drop authorities from species names
  #mutate(species_original = word(species_original, 1, 2)) %>%
  # Drop authorities from genus names
  mutate(genus_original = str_extract(genus_original, "^\\S+")) %>%
  # Extract the genus names from species names if genus names are blank
  mutate(genus_original = ifelse(
    is.na(genus_original) | genus_original == "", 
    str_extract(species_original, "^\\S+"),
    genus_original
  )) %>%
  # Remove occurrences that are only annotated to family or kingdom
  filter(!(taxonRank %in% c("Family", "Kingdom"))) %>%
  select(ID, family_original, genus_original) %>%
  glimpse()

# There are 36,505 occurrences
nrow(fungal_root_occurrences)

# Match the FungalRoot occurrences to WFO
fungal_root_occurrences_WFO <-
  WFO.one(WFO.match.fuzzyjoin(spec.data = fungal_root_occurrences,
                              WFO.data = WFO_genus,
                              spec.name = "genus_original"),
          verbose = FALSE) %>%
  glimpse()

# Direct matches: 36,459 out of 36,505
fungal_root_occurrences_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 46 out of 36,505
fungal_root_occurrences_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(genus_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Save the harmonised FungalRoot occurrences database
fungal_root_occurrences_WFO %>%
  select(ID, family, genus = scientificName) %>%
  write_excel_csv("data/FungalRoot/harmonised_occurrences.csv")

##### (2c) GlobalTreeSearch #####

# Download the GlobalTreeSearch database:
URL <- "https://tools.bgci.org/global_tree_search_trees_1_7.csv"
file_destination <- "data/GlobalTreeSearch/global_tree_search_trees_1_7.csv"
download.file(URL, file_destination)

# Read in the database:
global_tree_database <-
  read.csv("data/GlobalTreeSearch/global_tree_search_trees_1_7.csv") %>%
  select(species_original = TaxonName) %>%
  distinct(species_original, .keep_all = TRUE) %>%
  as.data.frame() %>%
  glimpse()

# There are 57,922 unique tree species in the GlobalTreeSearch database:
nrow(global_tree_database)

# !!! NOTE !!!  I've ran into memory issues so i've had to increase my R memeory
#               limit by accessing my ".Renviron" file and saving "".

# Match the GlobalTreeSearch to WFO
global_tree_database_WFO <-
  WFO.one(WFO.match.fuzzyjoin(spec.data = global_tree_database,
                              WFO.data = WFO_species,
                              spec.name = "species_original"),
          verbose = FALSE) %>%
  glimpse()

# Direct matches: 57,330 out of 57,922
global_tree_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 583 out of 57,922
global_tree_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched
global_tree_database_WFO %>%
  filter(Matched == "FALSE") %>%
  select(species_original)

# Save the harmonised GlobalTreeSearch database
global_tree_database_WFO %>%
  # Remove the unmatched species
  filter(Matched == TRUE) %>%
  select(family, genus, species = scientificName) %>%
  write_excel_csv("data/GlobalTreeSearch/harmonised_global_tree.csv")

# Remove all objects from the global environment
rm(list = ls())

#### (3) Plant distribution databases #########################################

# Read in the WFO datasets
WFO_species <- read.csv("data/WorldFloraOnline/WFO_species.csv")
WFO_genus <- read.csv("data/WorldFloraOnline/WFO_genus.csv")

##### (3a) BiomassPlotLib #####

# Download the BiomassPlotLib dataset:
URL <- "https://field-geoserver.jrsrp.com/geoserver/aus/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=aus%3Abiolib_treelist&outputFormat=csv&srs=EPSG%3A4326&cql_filter=(obs_time%20BETWEEN%201980-01-01%20AND%202023-10-16)%20AND%20INTERSECTS(geom%2CMULTIPOLYGON%20(((130.341796875%20-9.79567758282973%2C%20111.884765625%20-21.861498734372553%2C%20114.43359375%20-36.38591277287651%2C%20131.220703125%20-32.990235559651055%2C%20142.998046875%20-39.97712009843962%2C%20145.107421875%20-43.96119063892025%2C%20149.4140625%20-44.46515101351962%2C%20154.51171875%20-31.05293398570514%2C%20153.80859375%20-23.24134610238612%2C%20146.865234375%20-16.383391123608387%2C%20142.20703125%20-7.536764322084078%2C%20139.5703125%20-14.093957177836224%2C%20137.197265625%20-9.709057068618208%2C%20134.6484375%20-10.487811882056695%2C%20130.341796875%20-9.79567758282973))))"
file_destination <- "data/BiomassPlotLib/biolib_treelist.csv"
download.file(URL, file_destination)

# Read in the database:
tree_list <- read.csv("data/BiomassPlotLib/biolib_treelist.csv") %>%
  rename(species_original = species) %>%
  # Remove 'subspecies', 'variety', and 'authority' annotations by extracting 
  # the first two words from species names
  mutate(species_original = str_split(species_original, "\\s+") %>% 
           sapply(function(x) paste(x[1:2], collapse = " "))) %>%
  # Mutate all non-descriptive species names to genus and remove rows ending with "NA"
  mutate(species_original = ifelse(
    grepl("NA$", species_original), "", gsub(
      "\\bTree\\b|\\bShrub\\b|\\bspecies\\b|\\bsp\\.|\\bspp\\.|\\baff\\.|\\bor\\b",
      "",
      species_original)
  )) %>%
  # Remove observations without any species annotations
  filter(
    !is.na(species_original) & species_original != ""
  ) %>%
  # Create the genus_original column and mutate genus names based on the first
  # species name
  mutate(genus_original = ifelse(
    !is.na(str_extract(species_original, "^\\S+")), str_extract(species_original, "^\\S+"), NA)
  ) %>%
  # Remove individuals with arbitrary genus names
  filter(!genus_original %in% c("Rainforest", "Unknown")) %>%
  # Remove trailing blank spaces from observations annotated to genus
  mutate(species_original = str_trim(species_original)) %>%
  as_tibble() %>%
  glimpse()

# There are 258,466 tree individuals in the BiomassPlotLib dataset:
nrow(tree_list)

# Subset unique genus annotations
tree_genus <- tree_list %>%
  filter(!grepl("\\s", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# There are 50 unique genus annotations
nrow(tree_genus)

# Subset species
tree_species <- tree_list %>%
  filter(grepl("\\s", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# There are 973 unique species annotations
nrow(tree_species)

# Match genera in the BiomassPlotLib to the WFO
tree_genus_WFO <-
  WFO.one(WFO.match.fuzzyjoin(spec.data = tree_genus,
                              WFO.data = WFO_genus,
                              spec.name = "species_original"),
          verbose = FALSE) %>%
  glimpse()

# Direct matches: 50 out of 50
tree_genus_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Match species in the BiomassPlotLib to the WFO
tree_species_WFO <-
  WFO.one(WFO.match.fuzzyjoin(spec.data = tree_species,
                              WFO.data = WFO_species,
                              spec.name = "species_original"),
          verbose = FALSE) %>%
  glimpse()

# Direct matches: 912 out of 979
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 60 out of 979
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched
tree_species_WFO %>%
  filter(Matched == "FALSE") %>%
  select(species_original)

bind_rows(tree_species_WFO, tree_genus_WFO) %>%
  select(species_original, family, genus,
         scientific_name = scientificName, taxon_rank = taxonRank) %>%
  distinct(species_original, .keep_all = TRUE) %>%
  inner_join(tree_list, by = "species_original") %>%
  select(-c(species_original, genus_original)) %>%
  write_excel_csv("data/BiomassPlotLib/harmonised_treelist.csv")

##### (3b) HAVPlot #####

# Download the HAVPlot dataset: https://data.csiro.au/collection/csiro%3A54461v4
URL = "https://s3.data.csiro.au/dapprd/000054461v001/data/speciesAttributes.csv?response-content-disposition=attachment%3B%20filename%3DspeciesAttributes.csv&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20231017T062507Z&X-Amz-SignedHeaders=host&X-Amz-Expires=172800&X-Amz-Credential=2K6P74GE0OUPTMHV3N3A%2F20231017%2FCDC%2Fs3%2Faws4_request&X-Amz-Signature=99cb596aa32cfb8179df6777591bf65ef9ed1eb9d53295f28f84b7a54e4f07fb"
destfile <- "data/HAVPlot/speciesAttributes.csv"
download.file(URL, destfile)

# Read in the database:
HAVPlot <- read.csv("data/HAVPlot/speciesAttributes.csv") %>%
  mutate(unique_id = row_number()) %>%
  glimpse()
# There are 19275 unique taxa in the HAVPlot dataset
nrow(HAVPlot)

HAVPlot_aggregateOrganismObservation <- HAVPlot %>%
  select(unique_id, scientificName) %>%
  inner_join(
    read.csv("data/HAVPlot/aggregateOrganismObservation.csv"),
    by = "scientificName") %>%
  select(-c(scientificName, verbatimScientificName)) %>%
  glimpse()
# There are 5,802,682 individual records in the HAVPlot dataset
nrow(HAVPlot_aggregateOrganismObservation)

# Check taxon ranks
unique(HAVPlot$taxonRank)

# There is 19,204 taxa annotated to spacies and genus
HAVPlot %>%
  filter(is.na(taxonRank) | taxonRank == "family")

# Genus subset
HAVPlot_genus <- HAVPlot %>%
  filter(taxonRank == "genus") %>%
  select(unique_id, scientificName_original = scientificName) %>%
  as.data.frame()
# There are 1,475 unique plant genera
nrow(HAVPlot_genus)

# Species subset
HAVPlot_species <- HAVPlot %>%
  filter(taxonRank == "species") %>%
  select(unique_id, scientificName_original = scientificName) %>%
  as.data.frame()
# There are 17,729 unique plant species
nrow(HAVPlot_species)

# Match genera in the HAVPlot to the WFO
HAVPlot_genus_WFO <-
  WFO.one(WFO.match.fuzzyjoin(spec.data = HAVPlot_genus,
                              WFO.data = WFO_genus,
                              spec.name = "scientificName_original"),
          verbose = FALSE) %>%
  glimpse()

# Direct matches: 1,470 out of 1,475
HAVPlot_genus_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 4 out of 1,475
HAVPlot_genus_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(scientificName_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched
HAVPlot_genus_WFO %>%
  filter(Matched == "FALSE") %>%
  select(scientificName_original)
# That's a species not genus!

# Put "Cardamine tenuifolia" into the species dataset
HAVPlot_species <- HAVPlot_genus_WFO %>%
  filter(Matched == "FALSE") %>%
  select(unique_id, scientificName_original) %>%
  bind_rows(HAVPlot_species) %>%
  glimpse()
# There are 17,730 unique plant species
nrow(HAVPlot_species)

# Match species in the HAVPlot to the WFO
HAVPlot_species_WFO <-
  WFO.one(WFO.match.fuzzyjoin(spec.data = HAVPlot_species,
                              WFO.data = WFO_species,
                              spec.name = "scientificName_original"),
          verbose = FALSE) %>%
  glimpse()

# Direct matches: 16,870 out of 17,729
HAVPlot_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 372 out of 17,729
HAVPlot_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(scientificName_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched: 
# These are mostly genus names with morphological or geographical descriptions.
# I'll grab the first names and match the against the WFO generic dataset
HAVPlot_species_WFO %>%
  filter(Matched == "FALSE") %>%
  select(scientificName_original)

# New subset of unmatched names reduced to genus
HAVPlot_species_to_genus <- HAVPlot_species_WFO %>%
  filter(Matched == FALSE) %>%
  mutate(
    scientificName_unmatched = str_extract(scientificName_original, "^\\S+")
    ) %>%
  select(unique_id, scientificName_original, scientificName_unmatched)
# 488 in the new match subset
nrow(HAVPlot_species_to_genus)

# Try to match unmatched "spices" at to genus
HAVPlot_species_to_genus_WFO <-
  WFO.one(WFO.match.fuzzyjoin(spec.data = HAVPlot_species_to_genus,
                              WFO.data = WFO_genus,
                              spec.name = "scientificName_unmatched"),
          verbose = FALSE) %>%
  glimpse()

# New direct matches: 405 out of 488
HAVPlot_species_to_genus_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 21 out of 488
HAVPlot_species_to_genus_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(scientificName_original, scientificName_unmatched, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched: These look like genus names with part of the authority names
# attached. I'll manually fix the names and then try again against the WFO
HAVPlot_species_to_genus_WFO %>%
  filter(Matched == FALSE) %>%
  select(scientificName_unmatched)
glimpse(HAVPlot_species_to_genus_WFO)

HAVPlot_species_to_genus_2 <- HAVPlot_species_to_genus_WFO %>%
  filter(Matched == FALSE) %>%
  select(unique_id, scientificName_original, scientificName_unmatched) %>%
mutate(
  scientificName_unmatched = case_when(
    scientificName_unmatched == "Trianthemaoorabulka" ~ "Trianthema",
    scientificName_unmatched == "Arthropodiumlbury" ~ "Arthropodium",
    scientificName_unmatched == "Oleariaarnarvon" ~ "Olearia",
    scientificName_unmatched == "Burmanniaathurst" ~ "Burmannia",
    scientificName_unmatched == "Tripterococcusrachylobus" ~ "Tripterococcus",
    scientificName_unmatched == "Dansiealtanmoui" ~ "Dansiea",
    scientificName_unmatched == "Ackamaellenden" ~ "Ackama",
    scientificName_unmatched == "Cyperusoleman" ~ "Cyperus",
    scientificName_unmatched == "Fimbristylisharles" ~ "Fimbristylis",
    scientificName_unmatched == "Lepidospermarcher" ~ "Lepidosperma",
    scientificName_unmatched == "Lepidospermaandalup" ~ "Lepidosperma",
    scientificName_unmatched == "Lepidospermaarracarrup" ~ "Lepidosperma",
    scientificName_unmatched == "Tetrarialackwood" ~ "Tetraria",
    scientificName_unmatched == "Hibbertialackdown" ~ "Hibbertia",
    scientificName_unmatched == "Leucopogonoolbunda" ~ "Leucopogon",
    scientificName_unmatched == "Leucopogonoorabbin" ~ "Leucopogon",
    scientificName_unmatched == "Leucopogonurrum" ~ "Leucopogon",
    scientificName_unmatched == "Leucopogonoolgardie" ~ "Leucopogon",
    scientificName_unmatched == "Leucopogonorrigin" ~ "Leucopogon",
    scientificName_unmatched == "Leucopogonoujinup" ~ "Leucopogon",
    scientificName_unmatched == "Erythroxylumrewer" ~ "Erythroxylum",
    scientificName_unmatched == "Erythroxylumholmondely" ~ "Erythroxylum",
    scientificName_unmatched == "Acaciaharters" ~ "Acacia",
    scientificName_unmatched == "Austrodolichosrnhem" ~ "Austrodolichos",
    scientificName_unmatched == "Derrislaudie" ~ "Derris",
    scientificName_unmatched == "Dillwyniaarren" ~ "Dillwynia",
    scientificName_unmatched == "Dillwyniaoolgardie" ~ "Dillwynia",
    scientificName_unmatched == "Galactiandoom" ~ "Galactia",
    scientificName_unmatched == "Indigoferaungaroo" ~ "Indigofera",
    scientificName_unmatched == "Mirbeliaursarioides" ~ "Mirbelia",
    scientificName_unmatched == "Tephrosiaarkly" ~ "Tephrosia",
    scientificName_unmatched == "Tephrosiaungaroo" ~ "Tephrosia",
    scientificName_unmatched == "Tephrosiaarnarvon" ~ "Tephrosia",
    scientificName_unmatched == "Tephrosiaonduplicate" ~ "Tephrosia",
    scientificName_unmatched == "Tephrosiaopperfield" ~ "Tephrosia",
    scientificName_unmatched == "Geraniumlpine" ~ "Geranium",
    scientificName_unmatched == "Goodeniaarnarvon" ~ "Goodenia",
    scientificName_unmatched == "Prostantheraaking" ~ "Prostanthera",
    scientificName_unmatched == "Amyemalligator" ~ "Amyema",
    scientificName_unmatched == "Brachychitonltanmoui" ~ "Brachychiton",
    scientificName_unmatched == "Dicarpidiumrnhem" ~ "Dicarpidium",
    scientificName_unmatched == "Helicteresachsten" ~ "Helicteres",
    scientificName_unmatched == "Helicteresulimba" ~ "Helicteres",
    scientificName_unmatched == "Hibiscusarnarvon" ~ "Hibiscus",
    scientificName_unmatched == "Lasiopetalumordate-leaved" ~ "Lasiopetalum",
    scientificName_unmatched == "Sidambalindum" ~ "Sida",
    scientificName_unmatched == "Sidarticulation" ~ "Sida",
    scientificName_unmatched == "Baeckeaarbalin" ~ "Baeckea",
    scientificName_unmatched == "Baeckeaurakin" ~ "Baeckea",
    scientificName_unmatched == "Chamelauciumendering" ~ "Chamelaucium",
    scientificName_unmatched == "Chamelauciumoolcalalya" ~ "Chamelaucium",
    scientificName_unmatched == "Leptospermumandalup" ~ "Leptospermum",
    scientificName_unmatched == "Thryptomenearrarang" ~ "Thryptomenea",
    scientificName_unmatched == "Phyllanthuslexandra" ~ "Phyllanthus",
    scientificName_unmatched == "Aristidaradshaw" ~ "Aristida",
    scientificName_unmatched == "Austrostipaarlingup" ~ "Austrostipa",
    scientificName_unmatched == "Calandriniaungalbin" ~ "Calandrinia",
    scientificName_unmatched == "Sedopsisulimba" ~ "Portulaca",
    scientificName_unmatched == "Alangiumlaudie" ~ "Alangium",
    scientificName_unmatched == "Samaderaarong" ~ "Samadera",
    scientificName_unmatched == "Argyrodendronoonjie" ~ "Argyrodendron",
    scientificName_unmatched == "Garcinialaudie" ~ "Garcinia",
    TRUE ~ scientificName_unmatched
  )
)

# Match unmatched genera
HAVPlot_species_to_genus_WFO_2 <-
  WFO.one(WFO.match.fuzzyjoin(spec.data = HAVPlot_species_to_genus_2,
                              WFO.data = WFO_genus,
                              spec.name = "scientificName_unmatched"),
          verbose = FALSE) %>%
  glimpse()

# Direct matches for unmatched genera: 61 out of 62
HAVPlot_species_to_genus_WFO_2 %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 1 out of 62
HAVPlot_species_to_genus_WFO_2 %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(scientificName_original, scientificName_unmatched, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)


# Join the matched subsets
HAVPlot_harmonised <- HAVPlot_species_to_genus_WFO %>%
  bind_rows(HAVPlot_species_to_genus_WFO_2) %>%
  select(-scientificName_unmatched) %>%
  bind_rows(HAVPlot_species_WFO, HAVPlot_genus_WFO) %>%
  filter(Matched == TRUE & scientificName_original != "Gen.") %>%
  select(unique_id, family, genus, 
         scientificName, taxonRank)
# 19,203 out of 19,275 were matched
nrow(HAVPlot_harmonised)
glimpse(HAVPlot_harmonised)

# Harmonised the HAVPlots aggregate organism dataset
HAVPlot_aggregateOrganismObservation_harmonised <-
  HAVPlot_harmonised %>%
  select(unique_id, scientificName) %>%
  inner_join(HAVPlot_aggregateOrganismObservation,
             by = c("unique_id")) %>%
  select(-unique_id)
# 5,787,940 out of 5,799,617 individual records were matched
nrow(HAVPlot_aggregateOrganismObservation_harmonised)
glimpse(HAVPlot_aggregateOrganismObservation_harmonised)

# Save the full harmonised datasets
HAVPlot %>%
  select(unique_id, AusNativeStatus) %>%
  inner_join(HAVPlot_harmonised, by = "unique_id") %>% 
  select(-unique_id) %>%
  write_excel_csv("data/HAVPlot/harmonised_speciesAttributes.csv")
HAVPlot_aggregateOrganismObservation_harmonised %>%
  write_excel_csv("data/HAVPlot/harmonised_aggregateOrganismObservation.csv")
