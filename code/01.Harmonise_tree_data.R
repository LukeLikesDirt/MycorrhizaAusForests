
# Script: Harmonise plant databases according to the World Flora Online

# Table of contents:
#   (1) Prepare the World Flora Online backbone database
#   (2) Download and harmonise GlobalTreeSearch
#   (3) Download and harmonise FungalRoot
#   (4) Download and harmonise HAVPlot
#       (4a) Harmonise tree names
#       (4b) Remove plots dominated by exotics
#       (4c) Retain most recent survey
#       (4d) Generate abundance dataset
#       (4e) Generate presence dataset
#   (5) Harmonise the Australian Plant Census (APC)
#   (6) Download and harmonise GBIF
#       (6a) Clean GBIF occurrences
#       (6b) Harmonise GBIF tree list
#   (7) Prepare occurrence data

# Required packages
require(terra)
require(data.table)
require(WorldFlora)
require(tidyverse)
require(jsonlite)

# (1) World Flora Online backbone ##############################################

# Download the WFO backbone database
WFO.download(
  WFO.url = paste0("https://zenodo.org/records/14538251/files/_DwC_backbone_R.zip?download=1"),
  save.dir = "data/WorldFloraOnline",
  WFO.remember = TRUE
)
# Remove the zip file
file.remove("data/WorldFloraOnline/WFO_Backbone.zip")

# Read in the WFO backbone
WFO.remember()
nrow(WFO.data) # 1,638,552

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
nrow(WFO_species) # 1,183,139
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
nrow(WFO_genus) # 44,450
# Check the mutation
WFO_genus %>%
  filter(grepl("× ", scientificName)) %>%
  glimpse()
WFO_genus %>%
  filter(grepl("(?<!\\s)×", scientificName, perl = TRUE)) %>%
  select(scientificName) %>%
  as_tibble() %>%
  print(n = Inf)

# Family dataset
WFO_family <- WFO.data %>%
  filter(taxonRank == "family") %>%
  distinct() %>%
  glimpse()

# Save the datasets
fwrite(WFO_species, "data/WorldFloraOnline/WFO_species.txt", sep = "\t")
fwrite(WFO_genus, "data/WorldFloraOnline/WFO_genus.txt", sep = "\t")
fwrite(WFO_family, "data/WorldFloraOnline/WFO_family.txt", sep = "\t")

# Remove all objects from the global environment
rm(list = ls())

# (2) Harmonise the FungalRoot databases #######################################

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Download the FungalRoot database:
dir.create("data/FungalRoot")
URL = "https://files.plutof.ut.ee/public/orig/77/EE/77EEAF156EA8B35CAE58E18E644ADBDCFC90EB439DE9BA012D4FCAB3ECD3FCAB.zip"
file_destination <- "data/FungalRoot/FungalRoot.zip"
download.file(URL, file_destination)

# Unzip the downdload
unzip(file_destination, exdir = "data/FungalRoot/")

# Remove the zip file
file.remove(file_destination)

##### (2a) Harmonise genera #####

# Read in the FungalRoot mycorrhizal type database:
fungal_root_database <- fread("data/FungalRoot/FungalRoot.csv") %>%
  select(genus_original = Genus, mycorrhizal_type = `Mycorrhizal type`) %>%
  # Extract the first word from genus values with multiple words
  mutate(genus_original = str_extract(genus_original, "^\\S+")) %>%
  distinct(genus_original, .keep_all = TRUE) %>%
  glimpse()

# There are 14,527 unique genera in the FungalRoot mycorrhizal type database:
nrow(fungal_root_database)

# Match the tree genera in FungalRoot mycorrhizal type to WFO:
fungal_root_database_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = fungal_root_database,
    WFO.data = WFO_genus,
    spec.name = "genus_original"
  ),
  verbose = FALSE) %>%
  glimpse()

# Direct matches: 14,367 out of 14,541
fungal_root_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 158 out of 14,541
fungal_root_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(genus_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# 2 unmatched
fungal_root_database_WFO %>% 
  filter(Matched == FALSE) %>%
  select(genus_original, taxonRank) %>%
  as_tibble() %>%
  print(n = Inf)

# Save the harmonised FungalRoot database
fungal_root_database_WFO %>%
  # Remove the unmatched genera
  filter(Matched == TRUE) %>%
  select(family, genus = scientificName, mycorrhizal_type) %>%
  fwrite("data/FungalRoot/harmonised_fungal_root.txt", sep = "\t")

# Read in the FungalRoot occurrences database:
fungal_root_occurrences <- fread("data/FungalRoot/occurrences.csv") %>%
  select(
    ID, family_original = family, genus_original = genus,
    species_original = scientificName, taxonRank
  ) %>%
  # Drop authorities names
  mutate(genus_original = str_extract(genus_original, "^\\S+")) %>%
  # Fill in blank genus names by extracting the genus names from the species
  mutate(
    genus_original = ifelse(
      is.na(genus_original) | genus_original == "", 
      str_extract(species_original, "^\\S+"),
      genus_original
    )
  ) %>%
  # Remove occurrences that are only annotated to family or kingdom
  filter(!(taxonRank %in% c("Family", "Kingdom"))) %>%
  select(ID, family_original, genus_original) %>%
  glimpse()

# There are 36,505 occurrences
nrow(fungal_root_occurrences)

# Match the FungalRoot occurrences to WFO
fungal_root_occurrences_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = fungal_root_occurrences,
    WFO.data = WFO_genus,
    spec.name = "genus_original"
  ),
  verbose = FALSE) %>%
  glimpse()

# Direct matches: 36,486 out of 36,505
fungal_root_occurrences_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 19 out of 36,505
fungal_root_occurrences_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(genus_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Save the harmonised FungalRoot occurrences database
fungal_root_occurrences_WFO %>%
  select(ID, family, genus = scientificName) %>%
  fwrite("data/FungalRoot/harmonised_occurrences.txt", sep = "\t")

# Fungal root measurements
fungal_root_measurements <- inner_join(
  fread("data/FungalRoot/occurrences.csv") %>%
    filter(taxonRank %in% c("Variety", "Subspecies", "Species")) %>%
    select(ID, family,	genus, scientificName),
  fread("data/FungalRoot/measurements.csv") %>%
    select(ID = `Core ID`, measurementType, measurementValue) %>%
    pivot_wider(names_from = measurementType, values_from = measurementValue) %>%
    rename(
      source = Source,
      name = Name,
      mycorrhizal_type = `Mycorrhiza type`,
      AM_intensity = `AM intensity`,
      AM_method = `AM method`,
      elevation_min = `Elevation min`,
      ECM_method = `EcM method`,
      host_age = `Host age`,
      ECM_intensity = `EcM intensity`,
      remark = `Remark: type`,
      ECM_frequency = `EcM frequency`,
      AM_frequency = `AM frequency`,
      ERM_intensity = `ErM intensity`,
      ERM_method = `ErM method`,
      ORM_method = `OrM method`,
      ORM_intensity = `OrM intensity`,
      elevation_max = `Elevation max`
    ) %>%
    mutate(
      mycorrhizal_type = case_when(
        mycorrhizal_type == "EcM, AM undetermined" ~ "EcM",
        mycorrhizal_type == "non-ectomycorrhizal (AM undetermined)" ~ "uncertain",
        mycorrhizal_type == "EcM, no AM colonization" ~ "EcM",
        mycorrhizal_type == "non-mycorrhizal" ~ "NM",
        mycorrhizal_type == "Other" ~ "uncertain",
        mycorrhizal_type == "EcM,AM" ~ "EcM-AM",
        mycorrhizal_type == "AM-like (non-vascular plants)" ~ "AM-like_non-vascular_plant",
        mycorrhizal_type == "ErM,EcM" ~ "ErM",
        mycorrhizal_type == "ErM,AM" ~ "ErM",
        TRUE ~ mycorrhizal_type
      )),
  by = "ID"
) %>%
  select(ID, species_original = scientificName, mycorrhizal_type) %>%
  # Remove 'subspecies', 'variety', and 'authority' annotations by extracting 
  # the first two words from species names
  mutate(
    species_original = str_split(species_original, "\\s+") %>% 
      sapply(function(x) paste(x[1:2], collapse = " "))
  ) %>%
  glimpse()

# Species list
species_list <- fungal_root_measurements %>%
  select(species_original) %>%
  distinct() %>%
  glimpse()

# Match the GlobalTreeSearch to WFO
fungal_root_measurements_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = species_list,
    WFO.data = WFO_species,
    spec.name = "species_original"
  ),
  verbose = FALSE) %>%
  glimpse()

# Direct matches: 13,566 out of 13,770
fungal_root_measurements_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 192 out of 13,770
fungal_root_measurements_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched
fungal_root_measurements_WFO %>%
  filter(Matched == "FALSE") %>%
  select(species_original)

# Save the harmonised GlobalTreeSearch database
fungal_root_measurements_WFO %>%
  # Remove the unmatched species
  filter(Matched == TRUE) %>%
  left_join(
    fungal_root_measurements %>%
      select(ID, species_original, mycorrhizal_type),
    by = "species_original"
  ) %>%
  select(ID, family, genus, scientific_name = scientificName, mycorrhizal_type) %>%
  fwrite("data/FungalRoot/harmonised_measurements.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

##### (2b) Harmonise species #####

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Read in the database with IDs
fungal_root_ids <- fread(
  "data/FungalRoot/occurrences.csv",
  header = TRUE
) %>%
  select(ID, species_original = scientificName) %>%
  # Grab the first two words from species names to remove 'subspecies', 'variety',
  # and 'authority' annotations
  mutate(
    species_original = str_split(species_original, "\\s+") %>% 
      sapply(function(x) paste(x[1:2], collapse = " "))
  ) %>%
  distinct()

# Read in unique species:
fungal_root_species <- fread(
  "data/FungalRoot/occurrences.csv",
  header = TRUE
) %>%
  select(species_original = scientificName) %>%
  mutate(species_original = gsub(" × ", " ×", species_original)) %>%
  mutate(species_original = gsub("× ", "×", species_original)) %>%
  mutate(species_original = gsub(" x ", " ×", species_original)) %>%
  # Grab the first two words from species names to remove 'subspecies', 'variety',
  # and 'authority' annotations
  mutate(
    species_original = str_split(species_original, "\\s+") %>% 
      sapply(function(x) paste(x[1:2], collapse = " "))
  ) %>%
  distinct(species_original, .keep_all = TRUE) %>%
  as.data.frame() %>%
  glimpse()

# Check the mutation
fungal_root_species %>%
  filter(grepl(" × ", species_original))
fungal_root_species %>%
  filter(grepl(" x ", species_original))
fungal_root_species %>%
  filter(grepl("× ", species_original))
fungal_root_species %>%
  filter(grepl(" ×", species_original)) %>%
  select(species_original) %>%
  as_tibble() %>%
  print(n = Inf)

# How many unique species are there?
nrow(fungal_root_species) # 14,356

# Split species list into 1,000 row chunks
chunk_size <- 1000
n_rows <- nrow(fungal_root_species)
n_chunks <- ceiling(n_rows / chunk_size)

# Create an empty list to store results
fungal_root_WFO_list <- list()

# Loop through chunks and match to WFO
for (i in 1:n_chunks) {
  # Calculate start and end rows for this chunk
  start_row <- (i - 1) * chunk_size + 1
  end_row <- min(i * chunk_size, n_rows)
  
  # Extract chunk
  fungal_root_chunk <- fungal_root_species[start_row:end_row, , drop = FALSE]
  
  # Print progress
  cat(sprintf("Processing chunk %d of %d (rows %d to %d)\n", 
              i, n_chunks, start_row, end_row))
  
  # Match the chunk to WFO
  fungal_root_WFO_chunk <- WFO.one(
    WFO.match.fuzzyjoin(
      spec.data = fungal_root_chunk,
      WFO.data = WFO_species,
      spec.name = "species_original"
    ),
    verbose = FALSE
  )
  
  # Store result in list
  fungal_root_WFO_list[[i]] <- fungal_root_WFO_chunk
}

# Combine all chunks into one dataframe
fungal_root_WFO <- bind_rows(fungal_root_WFO_list) %>%
  glimpse()

# Clean up temporary objects
rm(fungal_root_WFO_list, fungal_root_chunk, fungal_root_WFO_chunk)

# Direct matches: 13,572 out of 14,356
fungal_root_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 555 out of 14,356
fungal_root_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched 222 out of 14,356
fungal_root_WFO %>%
  filter(Matched == "FALSE") %>%
  select(species_original)

# Some fuzzy matches for the species epithet have been based on "NA" or 
# authority names, so I will remove matches based on having "NA" or species 
# epithet starting with a cpital letter or finishing in a period.
fungal_root_WFO_cleaned <- fungal_root_WFO %>%
  mutate(
    valid_match = case_when(
      Matched == FALSE ~ FALSE,
      Fuzzy == TRUE & (grepl("\\bNA\\b", species_original) |
                         str_detect(word(species_original, 2), "^[A-Z]") |
                         str_detect(word(species_original, 2), "\\.$")) ~ FALSE,
      TRUE ~ TRUE
    )
  ) %>%
  filter(
    valid_match == TRUE
  ) %>%
  select(-valid_match)

# How many species pass the cleaning?
nrow(fungal_root_WFO_cleaned) # 13,763 out of 14,356

# Save the harmonised FungalRoot species list
fungal_root_WFO_cleaned %>%
  # Remove the unmatched species
  filter(Matched == TRUE) %>%
  left_join(
    fungal_root_ids,
    by = "species_original"
  ) %>%
  select(ID, family, genus, scientific_name = scientificName) %>%
  fwrite("data/FungalRoot/harmonised_occurrences.txt", sep = "\t")

# Clean up
rm(list = ls())
gc()

# (3) Harmonise GlobalTreeSearch ###############################################

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Download the GlobalTreeSearch database:
dir.create("data/GlobalTreeSearch")
URL <- "https://tools.bgci.org/global_tree_search_trees_1_9.csv"
file_destination <- "data/GlobalTreeSearch/global_tree_search_trees_1_9.csv"
download.file(URL, file_destination)

# Read in the genus-level estimates from FungalRoot:
fungal_root_genus <- fread("data/FungalRoot/harmonised_fungal_root.txt") %>%
  select(genus, mycorrhizal_type) %>%
  # Remove uncertain mycorrhizal types
  filter(mycorrhizal_type != "uncertain") %>%
  # Mutate NM-AM to AM (NM-AM is a facilitative status not type)
  mutate(
    mycorrhizal_type = ifelse(mycorrhizal_type == "NM-AM", "AM", mycorrhizal_type),
    mycorrhizal_type = ifelse(mycorrhizal_type == "nM", "NM", mycorrhizal_type),
    mycorrhizal_type = ifelse(str_detect(mycorrhizal_type, "species-specific") == TRUE, "uncertain", mycorrhizal_type)
  ) %>%
  unique()

# Check for duplicated genera
fungal_root_genus %>%
  group_by(genus) %>%
  filter(n() > 1) %>%
  arrange(genus) %>%
  as_tibble() %>%
  print(n = Inf)

# Resolve duplicated genera
# Create a function to resolve duplicated mycorrhizal types for a genus
resolve_mycorrhizal_types <- function(types) {
  # Priority order (higher number = higher priority)
  priority <- c(
    "NM" = 1,
    "OM" = 2,
    "AM" = 3,
    "EcM" = 4,
    "EcM-AM" = 5,
    "ErM" = 6,
    "Thysanothus" = 7
  )
  
  # Get the highest priority type present
  types_present <- types[types %in% names(priority)]
  
  if (length(types_present) > 0) {
    priorities <- priority[types_present]
    return(types_present[which.max(priorities)])
  } else {
    # If type not in our priority list, just return the first one
    return(types[1])
  }
}

# Apply to deduplicate genera
fungal_root_genus <- fungal_root_genus %>%
  group_by(genus) %>%
  summarise(mycorrhizal_type = resolve_mycorrhizal_types(mycorrhizal_type)) %>%
  ungroup()

# Check if duplicates remain
fungal_root_genus %>%
  group_by(genus) %>%
  filter(n() > 1)

# Unique mycorrhizal types:
unique(fungal_root_genus$mycorrhizal_type)

# Read in the database for Australian trees. Unlike the global export above,
# GlobalTreeSearch has no direct-download URL for a country-filtered list, so
# this is a manual step: go to https://tools.bgci.org/global_tree_search.php,
# filter to Australia, and save the export to the path below. This file is
# the sole source of the native_status field used throughout the rest of
# this pipeline, so getting the current export matters.
australian_tree_database <- fread(
  "data/GlobalTreeSearch/global_tree_search_trees_1_9_australia.csv",
  header = TRUE
) %>%
  select(species_original = taxon) %>%
  distinct(species_original, .keep_all = TRUE) %>%
  as.data.frame() %>%
  glimpse()

# Read in the global database and add missing Australian species
global_tree_database <- fread(
  "data/GlobalTreeSearch/global_tree_search_trees_1_9.csv",
  header = TRUE
) %>%
  select(species_original = TaxonName) %>%
  distinct(species_original, .keep_all = TRUE) %>%
  # Add Australian species that aren't in the global database
  bind_rows(
    anti_join(australian_tree_database, ., by = "species_original")
  ) %>%
  as.data.frame() %>%
  # Add native status for Australian trees
  mutate(
    native_status = ifelse(
      species_original %in% australian_tree_database$species_original,
      "native", "non-native"
    )
  ) %>%
  glimpse()

# There are 57,609 unique tree species in the GlobalTreeSearch database:
nrow(global_tree_database)

# There are 3,249 native tree species in the GlobalTreeSearch database:
nrow(australian_tree_database)
global_tree_database %>%
  filter(native_status == "native") %>%
  nrow()

# Match the GlobalTreeSearch to WFO
global_tree_database_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = global_tree_database,
    WFO.data = WFO_species,
    spec.name = "species_original"
  ),
  verbose = FALSE) %>%
  glimpse()

# Direct matches: 57,050 out of 57,609
global_tree_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 550 out of 57,609
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
  filter(Matched == TRUE) %>%
  select(family, genus, scientific_name = scientificName, species_original) %>%
  # Add native status for Australian trees
  left_join(
    global_tree_database %>%
      select(species_original, native_status),
    by = "species_original"
  ) %>%
  # Extract species epithet from original scientific name
  mutate(species_epithet = word(scientific_name, 2)) %>%
  # Resolve genus names for unaccepted genera and synonym genera
  mutate(
    genus = case_when(
      genus == "Afromorus" ~ "Morus",
      genus == "Archidasyphyllum" ~ "Dasyphyllum",
      genus == "Ceodes" ~ "Pisonia",
      genus == "Lychnophorella" ~ "Lychnophora",
      genus == "Macrolearia" ~ "Olearia",
      genus == "Ardisia" & family == "Ericaceae" ~ "Leptecophylla",
      genus == "Esenbeckia" & family == "Ptychomniaceae" ~ "Garovaglia",
      genus == "Schizocalyx" & family == "Salvadoraceae" ~ "Dobera",
      genus == "Spiranthera" & family == "Pittosporaceae" ~ "Billardiera",
      genus == "Volkameria" & family == "Clethraceae" ~ "Clethra",
      TRUE ~ genus
    ),
    # Reconstruct full scientific name with updated genus
    scientific_name = paste(genus, species_epithet),
    # Add family information for genera with missing family
    family = case_when(
      genus == "Andea" ~ "Lauraceae",
      genus == "Bahiana" ~ "Euphorbiaceae",
      genus == "Morus" ~ "Moraceae",
      genus == "Dasyphyllum" ~ "Asteraceae",
      genus == "Hoffmannanthus" ~ "Asteraceae",
      genus == "Lachanodes" ~ "Asteraceae",
      genus == "Leptogonum" ~ "Polygonaceae",
      genus == "Lundinia" ~ "Asteraceae",
      genus == "Lychnophora" ~ "Asteraceae",
      genus == "Niemeyera" ~ "Sapotaceae",
      genus == "Olearia" ~ "Asteraceae",
      genus == "Maschalostachys" ~ "Asteraceae",
      genus == "Melanodendron" ~ "Asteraceae",
      genus == "Nahuatlea" ~ "Asteraceae",
      genus == "Neoarytera" ~ "Sapindaceae",
      genus == "Nototrichium" ~ "Amaranthaceae",
      genus == "Pladaroxylon" ~ "Asteraceae",
      genus == "Rockia" ~ "Nyctaginaceae",
      genus == "Scyphostegia" ~ "Salicaceae",
      TRUE ~ family
    )
  ) %>%
  select(family, genus, scientific_name, native_status) %>%
  unique() %>%
  # Add mycorrhizal type for genera in FungalRoot
  left_join(
    fungal_root_genus,
    by = "genus"
  ) %>%
  # Mutate NA  mycorrhizal type to "uncertain"
  mutate(mycorrhizal_type = ifelse(is.na(mycorrhizal_type), "uncertain", mycorrhizal_type)) %>%
  unique() %>%
  fwrite("generated_data/global_tree_mycorrhizal_types.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

# (4) Harmonise HAVPlot ########################################################

#### (4a) Harmonise tree names ####

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Download the HAVPlot dataset: https://doi.org/10.25919/5cex-4s70
# The CSIRO Data Access Portal signs S3 download links per request and they
# expire after 48 hours, so a fixed URL cannot be committed to source. Instead,
# query the portal's public files API for the current signed links and
# download the four files this script uses.
dir.create("data/HAVPlot")
havplot_files <- fromJSON(
  "https://data.csiro.au/dap/ws/v2/collections/csiro%3A54461v4/files?folder=/"
)$file
for (havplot_filename in c(
  "speciesAttributes.csv", "aggregateOrganismObservation.csv",
  "plotObservation.csv", "plot.csv"
)) {
  download.file(
    havplot_files$link$href[havplot_files$filename == havplot_filename],
    file.path("data/HAVPlot", havplot_filename)
  )
}
rm(havplot_files, havplot_filename)

# Read in the database:
tree_list <- fread("data/HAVPlot/speciesAttributes.csv") %>%
  filter(
    taxonRank %in% c("species", "genus"),
  ) %>%
  mutate(
    scientificName = gsub(" x ", " ×", scientificName),
    # Remove 'subspecies', 'variety', and 'authority' annotations by extracting
    species_original = str_split(scientificName, "\\s+") %>%
      sapply(function(x) paste(x[1:2], collapse = " ")),
    # Create the genus_original column and mutate genus names based on the first
    # species name
    genus_original = ifelse(
      !is.na(str_extract(species_original, "^\\S+")),
      str_extract(species_original, "^\\S+"), NA)
  ) %>%
  select(
    species_original, genus_original, aus_native_status = AusNativeStatus
  ) %>%
  unique(.) %>%
  glimpse(.)

HAVPlot_aggregateOrganismObservation <- fread(
  "data/HAVPlot/aggregateOrganismObservation.csv"
) %>%
  select(
    plot_observation = plotObservationID, species_original = scientificName,
    abundance_value = abundanceValue, abundance_unit = abundanceUnits
  ) %>%
  mutate(
    species_original = gsub(" x ", " ×", species_original),
    species_original = str_split(species_original, "\\s+") %>%
      sapply(function(x) paste(x[1:2], collapse = " ")) 
  ) %>%
  filter(
    species_original %in% tree_list$species_original
  )

# There are 5,784,885 tree species occurrences in the HAVPlot dataset:
nrow(HAVPlot_aggregateOrganismObservation)
# There are 18,651 unique species in the HAVPlot dataset:
unique(tree_list$species_original) %>%
  length()

# Subset genera
tree_genus <- tree_list %>%
  # Names without a space (i.e., one word) or with 'sp.'
  filter(!grepl("\\s", species_original) | grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# 1,488 with only genus annotations
nrow(tree_genus)

# Subset species
tree_species <- tree_list %>%
  filter(grepl("\\s", species_original) & !grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# There are 17,163 unique species annotations
unique(tree_species$species_original) %>%
  length()

# Match species
tree_species_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = tree_species,
    WFO.data = WFO_species,
    spec.name = "species_original"),
  verbose = FALSE
) %>%
  glimpse()

# Check matching
# Direct matches: 16,904 out of 17,148
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 243 out of 17,148
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, species_original, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched: 16 unmatched species
tree_species_WFO %>%
  filter(Matched == "FALSE") %>%
  select(species_original) %>%
  as_tibble() %>%
  print(n = Inf)

# Match genera
tree_genus_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = tree_genus,
    WFO.data = WFO_genus,
    spec.name = "species_original"),
  verbose = FALSE
) %>%
  glimpse()

# Check matching
# Direct matches: 1486 out of 1,488
tree_genus_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 2 out of 1,488
tree_genus_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, species_original, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Harmonised tree list
harmonised_tree_list <- tree_list %>%
  inner_join(
    .,
    bind_rows(tree_species_WFO, tree_genus_WFO),
    by = "species_original"
  ) %>%
  select(
    family,	genus,	scientific_name = scientificName,	taxon_rank = taxonRank,
    species_original, aus_native_status
  ) %>%
  filter(!is.na(scientific_name), !is.na(genus))

#### (4b) Remove plots dominated by exotics ####

harmonised_tree_list_native <- HAVPlot_aggregateOrganismObservation %>%
  inner_join(
    .,
    harmonised_tree_list,
    by = "species_original",
    relationship = "many-to-many"
  ) %>%
  # Filter to tree species
  filter(
    scientific_name %in% fread("generated_data/global_tree_mycorrhizal_types.txt")$scientific_name
  ) %>%
  # Calculate the proportion of native trees and the proportion of trees
  # identified to species-level at each plot observation
  mutate(
    abundance = 1,
    abundance_native = case_when(
      aus_native_status == "native" ~ 1,
      aus_native_status == "non-native" ~ 0
    )) %>%
  group_by(plot_observation) %>%
  mutate(
    proportion_native = sum(abundance_native) / sum(abundance)
  ) %>%
  ungroup() %>%
  # Remove to plot where native species are < 70% of all species and where
  # observations to species-level is < 90%, and then reduce the dataset to 
  # to observations with species level annotations. I want species level 
  # annotations because endemism and vulnerability will be estimated at the
  # species level.
  filter(
    proportion_native >= 0.7
  ) %>%
  select(-c(species_original, abundance, abundance_native, proportion_native))

# Number of plots before filtering: 212,968
HAVPlot_aggregateOrganismObservation %>%
  distinct(plot_observation) %>%
  nrow()
# Number of plots after filtering: 174867
harmonised_tree_list_native %>%
  distinct(plot_observation) %>%
  nrow()

#### (4c) Retain most recent survey ####

plot_info <- fread("data/HAVPlot/plotObservation.csv") %>%
  mutate(plot_observation = plotObservationID) %>%
  # Grab surveys that should capture complete tree diversity:
  filter(
    taxonomicScope %in% c("vascular plants", "perennial species","woody species")
  ) %>%
  # Grab the most recent survey date
  arrange(desc(obsStartDate)) %>%
  group_by(plotID) %>%
  slice(1) %>%
  ungroup(.) %>%
  inner_join(
    .,
    fread("data/HAVPlot/plot.csv"),
    by = c("plotID", "projectID")
  ) %>%
  select(
    site = plotID, plot_observation = plotObservationID, longitude = decimalLongitude,
    latitude = decimalLatitude, date = obsStartDate, project = projectID, area
  ) %>%
  # Remove very small plots
  filter(
    area >= 400
  )

#### (4d) Generate abundance dataset ####

abundance_data <- harmonised_tree_list_native %>%
  filter(
    abundance_unit == "individuals"
  ) %>%
  inner_join(
    .,
    plot_info,
    by = "plot_observation"
  ) %>%
  mutate(
    date = format(ymd(date), format = "%d/%m/%Y"),
    abundance = as.numeric(abundance_value),
    basal_area = NA_integer_,
    biomass = NA_integer_
  ) %>%
  # Remove plots with less than 5 individuals
  group_by(site) %>%
  filter(sum(abundance) >= 5) %>%
  ungroup() %>%
  select(
    family, genus, scientific_name, taxon_rank, site, longitude, latitude,
    date, project, area, abundance, basal_area, biomass
  ) %>%
  fwrite("data/HAVPlot/harmonised_tree_list_abundance.txt", sep = "\t")

#### (4e) Generate presence dataset ####
harmonised_tree_list_native %>%
  inner_join(
    .,
    plot_info,
    by = "plot_observation"
  ) %>%
  mutate(
    abundance = 1,
    date = format(ymd(date), format = "%d/%m/%Y")
  ) %>%
  select(
    family, genus, scientific_name, taxon_rank, site, longitude, latitude,
    date, project, area, abundance
  ) %>%
  fwrite("data/HAVPlot/harmonised_tree_list_presence.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

# (5) Harmonise the Australian Plant Census (APC) ##############################

# The Australian Plant Census (APC) has no stable direct-download URL, so this
# step is manual. Export the current taxon checklist from the Australian
# Plant Census section of https://biodiversity.org.au/nsl/services/search
# (or the APC-specific search at https://biodiversity.org.au/nsl/services/APC)
# and save the CSV as "data/APC/APC-taxon-<export-date>.csv", updating the
# filename read below to match.

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# The APC has no clean native/naturalised column. Distribution is recorded as
# free text in taxonDistribution, e.g. "WA (naturalised), Qld, NSW" or
# "NSW (native and naturalised)" -- one entry per state/territory, with an
# optional parenthetical qualifier where the species isn't simply native
# there. taxonDistribution is only ever populated on taxonomicStatus ==
# "accepted" rows (verified against the full export: every synonym-type row
# is blank), so no separate status filter is needed here -- a blank string
# already resolves to NA below. Classify a species as native to Australia if
# it's native (or "native and ...") in at least one state/territory, and as
# non-native only if every state/territory it's recorded in is naturalised,
# doubtfully naturalised, or formerly naturalised. "(uncertain origin)" on
# its own is treated as no evidence either way, since that's what it means.
# Validated against known cases: Eucalyptus globulus, Banksia integrifolia
# and Grevillea robusta correctly classify as native (all naturalised in
# parts of their own range); Lantana camara, Pinus radiata and Schinus
# terebinthifolia correctly classify as non-native.
classify_apc_distribution <- function(dist) {
  if (is.na(dist) || dist == "") return(NA_character_)
  entries <- trimws(strsplit(dist, ",")[[1]])
  native_evidence <- vapply(entries, function(entry) {
    qualifier <- str_match(entry, "\\(([^)]*)\\)")[, 2]
    if (is.na(qualifier)) return(TRUE)
    if (grepl("native", qualifier, fixed = TRUE)) return(TRUE)
    if (identical(qualifier, "uncertain origin")) return(NA)
    if (grepl("naturalised", qualifier, fixed = TRUE)) return(FALSE)
    return(TRUE) # e.g. "presumed extinct" alone -- was native there
  }, logical(1))
  if (any(native_evidence == TRUE, na.rm = TRUE)) return("native")
  if (any(native_evidence == FALSE, na.rm = TRUE)) return("non-native")
  NA_character_
}

# Read in the database:
apc_database <- fread(
  "data/APC/APC-taxon-2025-03-18-0622.csv",
  header = TRUE
) %>%
  filter(taxonRank %in% c(
    "[infraspecies]", "[n/a]", "[unranked]", "Forma", "Nothovarietas",
    "Species", "Subforma", "Subspecies", "Subvarietas", "Varietas"
  )) %>%
  mutate(
    apc_native_status = vapply(taxonDistribution, classify_apc_distribution, character(1))
  ) %>%
  select(species_original = canonicalName, apc_native_status) %>%
  mutate(
    # Handle hybrids
    species_original = gsub(" x ", " ×", species_original),
    # Remove 'subspecies' and 'variety' annotations by extracting the first two
    # words from species names. This collapses infraspecific records (which
    # rarely carry their own taxonDistribution) onto the same species_original
    # as their parent species, so a plain distinct() risks arbitrarily keeping
    # a blank-status infraspecific row over the real species-level one --
    # aggregate to prefer whichever row has a non-NA status instead.
    species_original = map_chr(
      str_split(species_original, "\\s+"),
      ~ paste(head(.x, 2), collapse = " ")
    )
  ) %>%
  group_by(species_original) %>%
  summarise(
    apc_native_status = if (all(is.na(apc_native_status))) NA_character_ else apc_native_status[!is.na(apc_native_status)][1],
    .groups = "drop"
  ) %>%
  as.data.frame() %>%
  glimpse()

# Split data into chunks
apc_database_split <- split(
  apc_database,
  ceiling(seq_along(apc_database$species_original) / (nrow(apc_database) / 20))
)

# Apply the matching function to each chunk
processed_chunks <- lapply(apc_database_split, function(chunk) {
  WFO.one(
    WFO.match.fuzzyjoin(
      spec.data = chunk,
      WFO.data = WFO_species,
      spec.name = "species_original"
    ),
    verbose = FALSE
  )
})

# Combine the processed chunks
apc_database_WFO_combined <- do.call(rbind, processed_chunks)

# Glimpse the final collapsed dataset
glimpse(apc_database_WFO_combined)

# Check matching
apc_database_WFO_combined %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()
apc_database_WFO_combined %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)
apc_database_WFO_combined %>%
  filter(Matched == "FALSE") %>%
  as_tibble() %>%
  print()

# Filter to GlobalTreeSearch species
apc_tree_database <- apc_database_WFO_combined %>%
  select(family, genus, scientific_name = scientificName) %>%
  # Make amendments consistent with GlobalTreeSearch
  # Extract species epithet from original scientific name
  mutate(species_epithet = word(scientific_name, 2)) %>%
  # Resolve genus names for unaccepted genera and synonym genera
  mutate(
    genus = case_when(
      genus == "Afromorus" ~ "Morus",
      genus == "Archidasyphyllum" ~ "Dasyphyllum",
      genus == "Ceodes" ~ "Pisonia",
      genus == "Lychnophorella" ~ "Lychnophora",
      genus == "Macrolearia" ~ "Olearia",
      genus == "Ardisia" & family == "Ericaceae" ~ "Leptecophylla",
      genus == "Esenbeckia" & family == "Ptychomniaceae" ~ "Garovaglia",
      genus == "Schizocalyx" & family == "Salvadoraceae" ~ "Dobera",
      genus == "Spiranthera" & family == "Pittosporaceae" ~ "Billardiera",
      genus == "Volkameria" & family == "Clethraceae" ~ "Clethra",
      TRUE ~ genus
    ),
    # Reconstruct full scientific name with updated genus
    scientific_name = paste(genus, species_epithet),
    # Add family information for genera with missing family
    family = case_when(
      genus == "Morus" ~ "Moraceae",
      genus == "Dasyphyllum" ~ "Asteraceae",
      genus == "Hoffmannanthus" ~ "Asteraceae",
      genus == "Lachanodes" ~ "Asteraceae",
      genus == "Leptogonum" ~ "Polygonaceae",
      genus == "Lundinia" ~ "Asteraceae",
      genus == "Lychnophora" ~ "Asteraceae",
      genus == "Niemeyera" ~ "Sapotaceae",
      genus == "Olearia" ~ "Asteraceae",
      genus == "Maschalostachys" ~ "Asteraceae",
      genus == "Melanodendron" ~ "Asteraceae",
      genus == "Nahuatlea" ~ "Asteraceae",
      genus == "Neoarytera" ~ "Sapindaceae",
      genus == "Nototrichium" ~ "Amaranthaceae",
      genus == "Pladaroxylon" ~ "Asteraceae",
      genus == "Rockia" ~ "Nyctaginaceae",
      genus == "Scyphostegia" ~ "Salicaceae",
      TRUE ~ family
    )
  ) %>%
  select(family, genus, scientific_name) %>%
  filter(
    # Keep native Australian tree species only, according to GlobalTreeSearch
    scientific_name %in% (
      fread("generated_data/global_tree_mycorrhizal_types.txt") %>%
        filter(native_status == "native") %>%
        pull(scientific_name)
    )
  ) %>%
  unique(.)

# Save the harmonised APC flora and tree list.
# harmonised_apc_flora_list.txt now carries native_status, derived above from
# taxonDistribution, so Figure_S20.R's full-flora "plants" object can be
# filtered to native_status == "native" the same way its tree-only "trees"
# object already is via global_tree_mycorrhizal_types.txt. Where a species
# has multiple matching rows, take whichever native_status is non-NA (there
# should be at most one real value per resolved name, since taxonDistribution
# is only ever populated on the accepted-status row).
apc_database_WFO_combined %>%
  select(family, genus, scientific_name = scientificName, native_status = apc_native_status) %>%
  group_by(family, genus, scientific_name) %>%
  summarise(
    native_status = if (all(is.na(native_status))) NA_character_ else native_status[!is.na(native_status)][1],
    .groups = "drop"
  ) %>%
  fwrite("data/APC/harmonised_apc_flora_list.txt", sep = "\t")
apc_tree_database %>%
  fwrite("data/APC/harmonised_aus_tree_list.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

# (6) Harmonise GBIF data ######################################################

#### (6a) Clean GBIF occurrences ####

# Two manual inputs are required here; neither is downloaded automatically:
# - data/gbif/0071622-250525065834625.csv: the raw occurrence export for this
#   GBIF download, cited at https://doi.org/10.15468/DL.ATKQFB. Re-download
#   from that DOI (or run a fresh GBIF occurrence search under the same
#   filters) to reproduce or update this file.
# - data/gbif/normalized.csv: the Australian native tree checklist matched
#   against the GBIF backbone taxonomy, providing the "species" column used
#   to filter the raw download below. How the shipped copy was produced is
#   not recorded in this script -- likely GBIF's species-match tool
#   (https://www.gbif.org/tools/species-lookup) or the rgbif package's
#   name_backbone_checklist() run against the native tree species list from
#   generated_data/global_tree_mycorrhizal_types.txt, but this needs
#   confirming rather than assumed.

# Target species: Australian native tree species normalised to GBIF backbone
aus_tree_species <- fread("data/gbif/normalized.csv")$species

# Initialise empty result
gbif_trees_clean <- data.frame()

# Set chunk parameters
chunk_size <- 1000000
current_row <- 1
chunk_number <- 1

# Read header first
header <- names(fread("data/gbif/0071622-250525065834625.csv", quote = "", nrows = 0))

# Process in chunks
repeat {
  cat("Processing chunk", chunk_number, "\n")
  
  # Read chunk with header
  if (chunk_number == 1) {
    # First chunk includes header
    chunk <- fread("data/gbif/0071622-250525065834625.csv", 
                   quote = "", 
                   nrows = chunk_size)
    current_row <- chunk_size + 1
  } else {
    # Subsequent chunks need header added
    chunk <- fread("data/gbif/0071622-250525065834625.csv", 
                   quote = "", 
                   skip = current_row,
                   nrows = chunk_size,
                   col.names = header)
    current_row <- current_row + chunk_size
  }
  
  # Break if no more data
  if (nrow(chunk) == 0) break
  
  # Filter to target species
  chunk_filtered <- chunk %>%
    filter(
      species %in% aus_tree_species,
      !is.na(species) & species != "",
      !is.na(decimalLongitude) & !is.na(decimalLatitude)
    )
  
  # Only process if we have matching species
  if (nrow(chunk_filtered) > 0) {
    # Clean coordinates
    chunk_clean <- chunk_filtered %>%
      CoordinateCleaner::clean_coordinates(
        lon = "decimalLongitude",
        lat = "decimalLatitude", 
        species = "species",
        tests = c("capitals", "centroids", "institutions", "zeros", "seas", "duplicates"),
        verbose = FALSE
      ) %>%
      select(
        family, genus, species, 
        longitude = decimalLongitude, latitude = decimalLatitude
      )
    
    # Add to results
    gbif_trees_clean <- rbind(gbif_trees_clean, chunk_clean)
    
    cat("Added", nrow(chunk_clean), "clean records. Total so far:", nrow(gbif_trees_clean), "\n")
  }
  
  # Clean up memory
  rm(chunk, chunk_filtered)
  if (exists("chunk_clean")) rm(chunk_clean)
  gc()
  
  # Update for next chunk
  chunk_number <- chunk_number + 1
}

# Final result
cat("Final dataset has", nrow(gbif_trees_clean), "records\n")
glimpse(gbif_trees_clean)

#### (6b) Harmonise GBIF tree list ####

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Read in normalised names and APC/WFO verbatim names
species_list <- fread("data/gbif/normalized.csv") %>%
  select(species_original = species) %>%
  unique(.)

# Match species
gbif_species_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = species_list,
    WFO.data = WFO_species,
    spec.name = "species_original"
  ),
  verbose = FALSE
) %>%
  glimpse()

# Direct matches: 3,744 out of 3,772
gbif_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 27 out of 3,772
gbif_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched: 1 out of 3,772
gbif_species_WFO %>%
  filter(Matched == "FALSE") %>%
  select(species_original) %>%
  as_tibble() %>%
  print(n = Inf)

# Retain Australian native tree species according to APC tree list
gbif_tree_species <- gbif_species_WFO %>%
  select(
    family, genus, 
    scientific_name = scientificName, 
    species_original
  ) %>%
  # Make amendments consistent with GlobalTreeSearch
  # Extract species epithet from original scientific name
  mutate(species_epithet = word(scientific_name, 2)) %>%
  # Resolve genus names for unaccepted genera and synonym genera
  mutate(
    genus = case_when(
      genus == "Afromorus" ~ "Morus",
      genus == "Archidasyphyllum" ~ "Dasyphyllum",
      genus == "Ceodes" ~ "Pisonia",
      genus == "Lychnophorella" ~ "Lychnophora",
      genus == "Macrolearia" ~ "Olearia",
      genus == "Ardisia" & family == "Ericaceae" ~ "Leptecophylla",
      genus == "Esenbeckia" & family == "Ptychomniaceae" ~ "Garovaglia",
      genus == "Schizocalyx" & family == "Salvadoraceae" ~ "Dobera",
      genus == "Spiranthera" & family == "Pittosporaceae" ~ "Billardiera",
      genus == "Volkameria" & family == "Clethraceae" ~ "Clethra",
      TRUE ~ genus
    ),
    # Reconstruct full scientific name with updated genus
    scientific_name = paste(genus, species_epithet)
  ) %>%
  select(scientific_name, species_original) %>%
  # Filter to APC and add genus and species according to APC harmonised names
  inner_join(
    fread("data/APC/harmonised_aus_tree_list.txt") %>%
      select(scientific_name, family, genus),
    by = "scientific_name"
  ) %>%
  unique(.) %>%
  select(
    family, genus, scientific_name, species_original
  ) %>%
  glimpse()

# Update taxonomy in the cleaned GBIF occurrence data
gbif_trees_final <- gbif_trees_clean %>%
  select(species_original = species, longitude, latitude) %>%
  inner_join(
    .,
    gbif_tree_species,
    by = "species_original",
    relationship = "many-to-many"
  ) %>%
  select(family, genus, scientific_name, longitude, latitude) %>%
  unique(.) %>%
  glimpse()

# Save the harmonised GBIF tree list
gbif_trees_final %>%
  fwrite("data/GBIF/harmonised_gbif_tree_list.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

# (7) Prepare occurrence data ##################################################

#### (7a) Load data ####

# Download the Forests of Australia (2023) raster (see covariates.txt: aus_for23).
# ABARES does not expose a stable direct-download link for the raw GeoTIFF, so
# this fetches the same file from the Department of Agriculture, Fisheries and
# Forestry's spatial data page:
# https://www.agriculture.gov.au/abares/forestsaustralia/forest-data-maps-and-tools/spatial-data/forest-cover
dir.create("data/aus_forests_23", recursive = TRUE)
URL <- "https://www.agriculture.gov.au/sites/default/files/documents/Geo_Tiff.zip"
destfile <- "data/aus_forests_23/Geo_Tiff.zip"
download.file(URL, destfile)

# Unzip the download (contains aus_for23.tif plus its metadata siblings)
unzip(destfile, exdir = "data/aus_forests_23/")

# Remove the zip file
file.remove(destfile)

# Forest raster
forest_rast <- rast(
  "data/aus_forests_23/aus_for23.tif"
)

# Read in the harmonised HAVPlot data
HAVPlot_data <- fread("data/HAVPlot/harmonised_tree_list_presence.txt") %>%
  select(scientific_name, longitude, latitude)

# Read in the harmonised GBIF data
gbif_data <- fread("data/GBIF/harmonised_gbif_tree_list.txt") %>%
  select(scientific_name, longitude, latitude)

# Read in the tree list with mycorrhizal types
aus_tree_list <- fread("generated_data/global_tree_mycorrhizal_types.txt") %>%
  # Filter to native trees
  filter(native_status == "native") %>%
  # Remove "uncertain" mycorrhizal types
  filter(mycorrhizal_type != "uncertain") %>%
  # Mutate eucalypts to EcM
  mutate(
    mycorrhizal_type = case_when(
      genus == "Eucalyptus" ~ "EcM",
      genus == "Corymbia" ~ "EcM",
      genus == "Angophora" ~ "EcM",
      TRUE ~ mycorrhizal_type
    )
  ) %>%
  select(family, genus, scientific_name, mycorrhizal_type)

#### (7b) Occurrence data #####

# Combine the occurrence data with tree list of Australian native trees
trees <- bind_rows(
  HAVPlot_data,
  gbif_data
) %>%
  # Remove duplicates
  unique(.) %>%
  # Add family and genus information from the harmonised APC tree list and retain
  # only unique combinations of Australian native trees
  inner_join(
    .,
    aus_tree_list,
    by = "scientific_name",
    relationship = "many-to-many"
  ) %>%
  select(family, genus, scientific_name, mycorrhizal_type, longitude, latitude) %>%
  unique(.) %>%
  # Generate site IDs based on unique geographic coordinates
  group_by(longitude, latitude) %>%
  mutate(
    site = paste0("site_", cur_group_id()) 
  ) %>%
  ungroup() %>%
  glimpse()

# Number of tree species (printed at runtime)
unique(trees$scientific_name) %>%
  length() %>%
  message("Number of tree species: ", .)

# Number of sites (printed at runtime)
unique(trees$site) %>%
  length() %>%
  message("Number of sites: ", .)

# Check for unassigned mycorrhizal types
trees %>%
  filter(is.na(mycorrhizal_type) | mycorrhizal_type == "" | mycorrhizal_type == "uncertain")

#### (7c) Site data ####

# Filter sites to forests

# Extract unique sites with coordinates
all_sites <- trees %>%
  select(site, longitude, latitude) %>%
  distinct() %>%
  glimpse(.)

# Project longlat according to the forest_rast
coords <- all_sites %>%
  select(x = longitude, y = latitude) %>%
  as.matrix() %>%
  # Longitude and latitude as a spatial vector
  vect(., crs = '+proj=longlat') %>%
  # Project coordinates according to forest_rast
  project(., forest_rast)

# Forest structure values for each plot
structure_values <- terra::extract(forest_rast, coords) %>%
  select(FOR_CATEGO)

# Forest sites
forest_sites <- all_sites %>%
  bind_cols(structure_values) %>%
  filter(
    FOR_CATEGO == "Native forest"
  ) %>%
  select(-FOR_CATEGO)

# Forest trees
forest_trees <- trees %>%
  filter(
    site %in% forest_sites$site
  ) %>%
  distinct() %>%
  glimpse(.)

# Number of forest tree species (printed at runtime)
forest_trees %>%
  select(scientific_name) %>%
  distinct() %>%
  nrow() %>%
  message("Number of tree species: ", .)

# Number of forest sites (printed at runtime)
forest_sites %>%
  select(site) %>%
  distinct() %>%
  nrow() %>%
  message("Number of sites: ", .)

### (7d) Save forest tree data ####

# Compute the number of observations per site as a proxy to sampling effort
sample_effort <- forest_trees %>%
  group_by(site) %>%
  summarise(
    n_obs = n()
  ) %>%
  ungroup() %>%
  glimpse(.)

# Save the forest tree data
forest_trees %>%
  distinct() %>%
  fwrite(
    "data/presence/trees.txt",
    sep = "\t"
  )

# Save the forest site data
forest_sites %>%
  inner_join(
    sample_effort,
    by = "site"
  ) %>%
  distinct() %>%
  fwrite(
    "data/presence/sites.txt",
    sep = "\t"
  )

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()
