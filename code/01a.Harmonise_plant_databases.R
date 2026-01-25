
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
require(data.table)
require(WorldFlora)
require(tidyverse)

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

# (2) Harmonise GlobalTreeSearch ###############################################

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Download the GlobalTreeSearch database:
dir.create("data/GlobalTreeSearch")
URL <- "https://tools.bgci.org/global_tree_search_trees_1_8.csv"
file_destination <- "data/GlobalTreeSearch/global_tree_search_trees_1_8.csv"
download.file(URL, file_destination)

# Read in the database:
global_tree_database <- fread(
  "data/GlobalTreeSearch/global_tree_search_trees_1_8.csv",
  header = TRUE
) %>%
  select(species_original = TaxonName) %>%
  distinct(species_original, .keep_all = TRUE) %>%
  as.data.frame() %>%
  glimpse()

# There are 57,681 unique tree species in the GlobalTreeSearch database:
nrow(global_tree_database)

# Match the GlobalTreeSearch to WFO
global_tree_database_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = global_tree_database,
    WFO.data = WFO_species,
    spec.name = "species_original"
  ),
  verbose = FALSE) %>%
  glimpse()

# Direct matches: 57,130 out of 57,681
global_tree_database_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 545 out of 57,681
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
  select(family, genus, scientific_name = scientificName) %>%
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
  unique() %>%
  fwrite("data/GlobalTreeSearch/harmonised_global_tree.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

# (3) Harmonise the FungalRoot databases #######################################

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

##### (3a) Harmonise genera #####

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
fungal_root_measurments <- inner_join(
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
species_list <- fungal_root_measurments %>%
  select(species_original) %>%
  distinct() %>%
  glimpse()

# Match the GlobalTreeSearch to WFO
fungal_root_measurments_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = species_list,
    WFO.data = WFO_species,
    spec.name = "species_original"
  ),
  verbose = FALSE) %>%
  glimpse()

# Direct matches: 13,566 out of 13,770
fungal_root_measurments_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 192 out of 13,770
fungal_root_measurments_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched
fungal_root_measurments_WFO %>%
  filter(Matched == "FALSE") %>%
  select(species_original)

# Save the harmonised GlobalTreeSearch database
fungal_root_measurments_WFO %>%
  # Remove the unmatched species
  filter(Matched == TRUE) %>%
  left_join(
    fungal_root_measurments %>%
      select(ID, species_original, mycorrhizal_type),
    by = "species_original"
  ) %>%
  select(ID, family, genus, scientific_name = scientificName, mycorrhizal_type) %>%
  fwrite("data/FungalRoot/harmonised_measurements.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

##### (3a) Harmonise species #####

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

# Save a harmonised FungalRoot tree species list
fread("data/GlobalTreeSearch/harmonised_global_tree.txt") %>%
  inner_join(
    fread("data/FungalRoot/harmonised_species_list.txt"),
    by = c("family", "genus", "scientific_name")
  ) %>%
  fwrite("data/FungalRoot/harmonised_tree_occurrences.txt", sep = "\t")

# Clean up
rm(list = ls())
gc()

# (4) Harmonise BiomassPlotLib #################################################

# Note: The BiomassPlotLib dataset has been supplemented with data from LTERN
# Tropical Forest Plots and the Natural Values Atlas. Data for 
# LTERN and Natural Values Atlas need to be sourced separately.
# LTERN: https://www.ltern.org.au/ltern-plot-networks/tropical-rainforest
# Natural Values Atlas: https://www.environment.gov.au/land/nva

##### (4a) BiomassPlotLib #####

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Download the BiomassPlotLib dataset:
dir.create("data/BiomassPlotLib/")
URL <- "https://field-geoserver.jrsrp.com/geoserver/aus/wfs?service=WFS&version=1.0.0&request=GetFeature&typeName=aus%3Abiolib_treelist&outputFormat=csv&srs=EPSG%3A4326&cql_filter=(obs_time%20BETWEEN%201980-01-01%20AND%202025-01-09)%20AND%20INTERSECTS(geom%2CMULTIPOLYGON%20(((130.341796875%20-9.79567758282973%2C%20111.884765625%20-21.861498734372553%2C%20114.43359375%20-36.38591277287651%2C%20131.220703125%20-32.990235559651055%2C%20142.998046875%20-39.97712009843962%2C%20145.107421875%20-43.96119063892025%2C%20149.4140625%20-44.46515101351962%2C%20154.51171875%20-31.05293398570514%2C%20153.80859375%20-23.24134610238612%2C%20146.865234375%20-16.383391123608387%2C%20142.20703125%20-7.536764322084078%2C%20139.5703125%20-14.093957177836224%2C%20137.197265625%20-9.709057068618208%2C%20134.6484375%20-10.487811882056695%2C%20130.341796875%20-9.79567758282973))))"
file_destination <- "data/BiomassPlotLib/biolib_treelist.csv"
download.file(URL, file_destination)

# Read in the database:
tree_list <- fread("data/BiomassPlotLib/biolib_treelist.csv") %>%
  rename(species_original = species) %>%
  mutate(
    # Remove 'subspecies', 'variety', and 'authority' annotations by extracting 
    # the first two words from species names
    species_original = str_split(species_original, "\\s+") %>% 
      sapply(function(x) paste(x[1:2], collapse = " ")),
    # Mutate all non-descriptive species names to genus and remove rows ending with "NA"
    species_original = ifelse(
      grepl("NA$", species_original), "", gsub(
        "\\bTree\\b|\\bShrub\\b|\\bspecies\\b|\\bsp\\.|\\bspp\\.|\\baff\\.|\\bor\\b",
        "",
        species_original)),
    # Create the genus_original column and mutate genus names based on the first
    # species name
    genus_original = ifelse(
      !is.na(str_extract(species_original, "^\\S+")),
      str_extract(species_original, "^\\S+"), NA),
    # Remove trailing blank spaces from observations annotated to genus
    species_original = str_trim(species_original),
    # Mutate tree ID to include dataset prefix
    tree_id = paste0("bpl_", tree),
    # Mutate area in ha to m2
    area = sitearea_ha * 10000
  ) %>%
  filter(
    # Remove observations without any species annotations
    !is.na(species_original) & species_original != "" & species_original != "NA",
    # Remove individuals with arbitrary genus names
    !genus_original %in% c("Rainforest", "Unknown")
  ) %>%
  # Organise names for merging purposes
  select(
    tree_id, species_original, genus_original, site, longitude, latitude,
    date = estdate, project, area, diameter_cm = diameter, height_m = ht,
    total_biomass_kg = tb_drymass, agb_drymass_kg = agb_drymass, bgb_drymass_kg = bgb_drymass
  ) %>%
  as_tibble() %>%
  glimpse()

# Check the data:
# Biomass estimates were previously misaligned with tree individuals. This was
# obvious because DBH and biomass had no relationship. Further, allometric models
# were not consistent for individual species. For example, biomass estimates for 
# Eucalyptus marginata were made using various agb allometric models such as
# “Eucalypt trees”, “Single stemmed acacia trees”, “Multi-stemmed acacias and
# mallees”, “Other trees - high wood density”, “Other trees - low wood density”,
# and “Shrubs”. Stephen Roxburgh (CSIRO) supplied me with the correct biomass
# estimates for the dataset.
plot(tree_list$bgb_drymass_kg, tree_list$agb_drymass_kg)
plot(tree_list$total_biomass_kg, tree_list$agb_drymass_kg)
plot(tree_list$total_biomass_kg, tree_list$bgb_drymass_kg)
plot(tree_list$diameter_cm, tree_list$total_biomass_kg)
plot(tree_list$diameter_cm, tree_list$height_m)

# There are 258,466 tree individuals in the BiomassPlotLib dataset:
nrow(tree_list)

# Tree genera
tree_genus <- tree_list %>%
  filter(!grepl("\\s", species_original) | grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# There are 50 tree genera without species annotations
nrow(tree_genus)

# Subset species
tree_species <- tree_list %>%
  filter(grepl("\\s", species_original) & !grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# There are 973 unique species
nrow(tree_species)

# Match genera in the BiomassPlotLib to the WFO
tree_genus_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = tree_genus,
    WFO.data = WFO_genus,
    spec.name = "species_original"
  ),
  verbose = FALSE
) %>%
  glimpse()

# Direct matches: 50 out of 50
tree_genus_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Match species in the BiomassPlotLib to the WFO
tree_species_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = tree_species,
    WFO.data = WFO_species,
    spec.name = "species_original"
  ),
  verbose = FALSE) %>%
  glimpse()

# Direct matches: 915 out of 979
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 57 out of 979
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# One unmatched
tree_species_WFO %>%
  filter(Matched == "FALSE") %>%
  select(species_original)

# Number of sites before filtering
tree_list %>%
  select(site) %>%
  distinct() %>%
  nrow()

# Number of sites after filtering
tree_list %>%
  filter(!is.na(diameter_cm), diameter_cm != 0) %>%
  mutate(abundance = 1) %>%
  group_by(site) %>%
  mutate(tree_abundance = sum(abundance)) %>%
  ungroup() %>%
  filter(tree_abundance >= 5) %>%
  select(site) %>%
  distinct() %>%
  nrow()

# Save the harmonised BiomassPlotLib database
bind_rows(tree_species_WFO, tree_genus_WFO) %>%
  select(species_original, family, genus,
         scientific_name = scientificName, taxon_rank = taxonRank) %>%
  inner_join(tree_list, by = "species_original") %>%
  select(-c(species_original, genus_original)) %>%
  # Remove individuals without diameter measurements
  filter(!is.na(diameter_cm), diameter_cm != 0) %>%
  # Limit to plots with at least 5 tree individuals
  mutate(abundance = 1) %>%
  group_by(site) %>%
  mutate(tree_abundance = sum(abundance)) %>%
  ungroup() %>%
  filter(tree_abundance >= 5) %>%
  select(-c(abundance, tree_abundance)) %>%
  fwrite("data/BiomassPlotLib/harmonised_treelist.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

#### (4b) LTERN ####

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Download the LTERN Permanent Rainforest Plots of North Queensland dataset:
URL <- "https://s3.data.csiro.au/dapprd/000006638v003/data/CSIRO_permanent_plots.zip?response-content-disposition=attachment%3B%20filename%3D%22CSIRO_permanent_plots.zip%22&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250109T060445Z&X-Amz-SignedHeaders=host&X-Amz-Expires=172800&X-Amz-Credential=BF3Z3R4Y8ZD2SPFRU2DQ%2F20250109%2FCDC%2Fs3%2Faws4_request&X-Amz-Signature=5f99eaff2291f6ef23da88523eb4eb9cd3d9583fc607799820f6f5a13eb535fa"
file_destination <- "data/CSIRO_permanent_plots.zip"
download.file(URL, file_destination)

# Unzip the downdload
unzip(file_destination, exdir = "data/")

# Remove the zip file
file.remove(file_destination)

# Site metadata is in html format and needs to be manually downloaded and
# copied to a data file:
# https://s3.data.csiro.au/dapprd/000006638v003/data/metadata.htm?response-content-disposition=attachment%3B%20filename%3D%22metadata.htm%22&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250109T222839Z&X-Amz-SignedHeaders=host&X-Amz-Expires=172800&X-Amz-Credential=BS0OVAQZ8VKKF7BYBKLX%2F20250109%2FCDC%2Fs3%2Faws4_request&X-Amz-Signature=4c3fa99d2a6d25707600b64a9760bc1558145e930a41386a8823014b29570e74

# Grab the most recent surveys for each plot
plots_2008 <- inner_join(
  fread("data/CSIRO_PermanentPlots_Data/data/CSIRO_PermanentPlots_TreeMeasurementData.csv"),
  fread("data/CSIRO_PermanentPlots_Data/data/site_metadata.csv"),
  by = "epNumber") %>%
  filter(year == 2008) %>%
  filter(site == "LTERN_McllwraithRange")
plots_2011 <- inner_join(
  fread("data/CSIRO_PermanentPlots_Data/data/CSIRO_PermanentPlots_TreeMeasurementData.csv"),
  fread("data/CSIRO_PermanentPlots_Data/data/site_metadata.csv"),
  by = "epNumber") %>%
  filter(year == 2011)
plots_2012 <- inner_join(
  fread("data/CSIRO_PermanentPlots_Data/data/CSIRO_PermanentPlots_TreeMeasurementData.csv"),
  fread("data/CSIRO_PermanentPlots_Data/data/site_metadata.csv"),
  by = "epNumber") %>%
  filter(year == 2012)
plots_2013 <- inner_join(
  fread("data/CSIRO_PermanentPlots_Data/data/CSIRO_PermanentPlots_TreeMeasurementData.csv"),
  fread("data/CSIRO_PermanentPlots_Data/data/site_metadata.csv"),
  by = "epNumber") %>%
  filter(year == 2013)

tree_list <- bind_rows(plots_2008, plots_2011, plots_2012, plots_2013) %>%
  # Remove individuals without diameter measurements
  filter(!is.na(dbh_centimetres), dbh_centimetres != 0) %>%
  mutate(abundance = 1) %>%
  group_by(site) %>%
  mutate(tree_abundance = sum(abundance)) %>%
  ungroup() %>%
  filter(tree_abundance >= 5) %>%
  select(-c(abundance, tree_abundance)) %>%
  mutate(
    # Remove 'subspecies', 'variety', and 'authority' annotations by extracting 
    # the first two words from species names
    species_original = str_split(taxon, "\\s+") %>% 
      sapply(function(x) paste(x[1:2], collapse = " ")),
    # Create the genus_original column and mutate genus names based on the first
    # species name
    genus_original = ifelse(
      !is.na(str_extract(species_original, "^\\S+")),
      str_extract(species_original, "^\\S+"), NA),
    # Remove trailing blank spaces from observations annotated to genus
    species_original = str_trim(species_original),
    # Mutate genus names based on the first species name
    genus_original = ifelse(
      !is.na(str_extract(species_original, "^\\S+")),
      str_extract(species_original, "^\\S+"), NA),
    # Mutate tree ID to include dataset prefix
    tree_id = paste0("ltern_", stemNumber),
  ) %>%
  # Remove individuals with arbitrary genus names
  filter(species_original != "Gen.(Aq520454) sp.") %>%
  # Remove observations without any species annotations
  filter(
    !is.na(species_original) & species_original != "" & species_original != "NA"
  ) %>%
  # Organise names for merging purposes
  select(
    tree_id, species_original, genus_original, site, longitude, latitude,
    date = date_last_survey, project, area = area,
    diameter_cm = dbh_centimetres, height_m = height_metres
  ) %>%
  as_tibble() %>%
  glimpse()

# There are 7605 tree individuals in the LTERN dataset:
nrow(tree_list)
# There are 437 unique species in the LTERN dataset
unique(tree_list$species_original) %>%
  length()

# Subset genera
tree_genus <- tree_list %>%
  # Names without a space (i.e., one word) or with 'sp.'
  filter(!grepl("\\s", species_original) | grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# 16 genus annotations
nrow(tree_genus)

# Subset species
tree_species <- tree_list %>%
  filter(grepl("\\s", species_original) & !grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# There are 421 unique species annotations
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
# Direct matches: 410 out of 421
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 11 out of 421
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, species_original, Old.name) %>%
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
# Direct matches: 16 out of 16
tree_genus_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

tree_list %>%
  inner_join(
    .,
    bind_rows(tree_species_WFO, tree_genus_WFO),
    by = "species_original"
  ) %>%
  select(
    family,	genus,	scientific_name = scientificName,	taxon_rank = taxonRank,
    tree_id, site,	longitude,	latitude,	date,	project, area,
    diameter_cm,	height_m
  ) %>%
  # Add blank columns for biomass measurements
  mutate(
    total_biomass_kg = NA, agb_drymass_kg = NA, bgb_drymass_kg = NA
  ) %>%
  filter(!is.na(scientific_name), !is.na(genus)) %>%
  fwrite("data/CSIRO_PermanentPlots_Data/harmonised_treelist.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

#### (4c) Natural Values Atlas ####

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Read in the Natural Values Atlas dataset: Supplied by Dr Adelina Latinovic
# from the Department of Natural Resources and Environment Tasmania
tree_list <- fread("data/NaturalValuesAtlas/treelist.csv") %>%
  select(
    tree_id = ID, species_original = SPECIES_NAME, genus_original = GENUS,
    site = OBSERVATION_FOREIGN_ID, longitude = LONGITUDE, latitude = LATITUDE,
    date = OBSERVATION_DATE, project = PROJECT_CODE, PLOT_LENGTH, PLOT_WIDTH,
    DIAMETER_AT_BREAST_HEIGHT_SINGLE_STEM, DIAMETER_AT_BREAST_HEIGHT_MULTISTEM
  ) %>%
  mutate(
    area = PLOT_LENGTH * PLOT_WIDTH,
    diameter_cm = ifelse(
      !is.na(DIAMETER_AT_BREAST_HEIGHT_SINGLE_STEM),
      DIAMETER_AT_BREAST_HEIGHT_SINGLE_STEM,
      DIAMETER_AT_BREAST_HEIGHT_MULTISTEM
    ),
    height_m = NA,
  ) %>%
  select(-c(
    PLOT_LENGTH, PLOT_WIDTH, DIAMETER_AT_BREAST_HEIGHT_SINGLE_STEM,
    DIAMETER_AT_BREAST_HEIGHT_MULTISTEM
  )) %>%
  # Remove individuals without diameter measurements
  filter(!is.na(diameter_cm), diameter_cm != 0) %>%
  mutate(abundance = 1) %>%
  group_by(site) %>%
  mutate(tree_abundance = sum(abundance)) %>%
  ungroup() %>%
  filter(tree_abundance >= 5) %>%
  select(-c(abundance, tree_abundance)) %>%
  mutate(
    # Remove 'subspecies', 'variety', and 'authority' annotations by extracting 
    # the first two words from species names
    species_original = str_split(species_original, "\\s+") %>% 
      sapply(function(x) paste(x[1:2], collapse = " ")),
    # Create the genus_original column and mutate genus names based on the first
    # species name
    genus_original = ifelse(
      !is.na(str_extract(species_original, "^\\S+")),
      str_extract(species_original, "^\\S+"), NA),
    # Remove trailing blank spaces from observations annotated to genus
    species_original = str_trim(species_original),
    # Mutate genus names based on the first species name
    genus_original = ifelse(
      !is.na(str_extract(species_original, "^\\S+")),
      str_extract(species_original, "^\\S+"), NA),
    # Mutate tree ID to include dataset prefix
    tree_id = paste0("nva_", tree_id),
    # Mutate site names
    site = str_remove(site, "^\\d+\\.\\d+\\s*"),
    site = str_replace_all(site, "[^\\w\\d]+", "_"),
    site = str_replace_all(site, "_{2,}", "_"),
    # Mutate date format
    date = format(dmy(date), format = "%d/%m/%Y")
  ) %>%
  # Remove observations without any species annotations
  filter(
    !is.na(species_original) & species_original != "" & species_original != "NA"
  ) %>%
  as_tibble() %>%
  glimpse()

# There are 238 tree individuals in the NVA dataset:
nrow(tree_list)
# There are 14 unique species in the NVA dataset
unique(tree_list$species_original) %>%
  length()

# Subset genera
tree_genus <- tree_list %>%
  # Names without a space (i.e., one word)
  filter(!grepl("\\s", species_original) | grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# No genus only annotations
nrow(tree_genus)

# Subset species
tree_species <- tree_list %>%
  filter(grepl("\\s", species_original) & !grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# There are 14 unique species annotations
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
# Direct matches: 14 out of 14
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

tree_list %>%
  inner_join(
    .,
    bind_rows(tree_species_WFO),
    by = "species_original"
  ) %>%
  select(
    family,	genus,	scientific_name = scientificName,	taxon_rank = taxonRank,
    tree_id, site,	longitude,	latitude,	date,	project, area,
    diameter_cm,	height_m
  ) %>%
  # Add blank columns for biomass measurements
  mutate(
    total_biomass_kg = NA, agb_drymass_kg = NA, bgb_drymass_kg = NA
  ) %>%
  filter(!is.na(scientific_name), !is.na(genus)) %>%
  fwrite("data/NaturalValuesAtlas/harmonised_treelist.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

#### (4d) FORESTCHECK #####

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Latitudinal and longitudinal coordinates for FORESTCHECK sites:
site_info <- fread("data/FORESTCHECK/sites.csv") %>%
  select(
    plot = locationID, longitude = decimalLongitude, latitude = decimalLatitude
  ) %>%
  unique() %>%
  filter(grepl("FC", plot)) %>%
  as_tibble() %>%
  glimpse()

# Read in the FORESTCHECK dataset (supplied by Dr Adrian Pinder, Biodiversity 
# and Conservation Science,Department of Biodiversity, Conservation and 
# Attractions):
tree_list <- fread("data/FORESTCHECK/StandStructureAssessment.csv") %>%
  inner_join(
    .,
    site_info,
    by = "plot") %>%
  mutate(
    tree_id = paste0("fc_", row_number()),
    species = str_to_sentence(species),
    species_original = case_when(
      species == "Jarrah" ~ "Eucalyptus marginata",
      species == "Marri" ~ "Corymbia calophylla",
      TRUE ~ species
    ),
    # Create the genus_original column and mutate genus names based on the first
    # species name
    genus_original = ifelse(
      !is.na(str_extract(species_original, "^\\S+")),
      str_extract(species_original, "^\\S+"), NA),
    site = location,
    project = "FORESTCHECK",
    date = format(ymd(paste(year, month, day, sep = "/")), format = "%d/%m/%Y"),
    height_m = NA,
    total_biomass_kg = NA, agb_drymass_kg = NA, bgb_drymass_kg = NA
  ) %>%
  # Remove individuals without diameter measurements
  filter(!is.na(diameter_cm), diameter_cm != 0) %>%
  mutate(abundance = 1) %>%
  group_by(site) %>%
  mutate(tree_abundance = sum(abundance)) %>%
  ungroup() %>%
  filter(tree_abundance >= 5) %>%
  select(-c(abundance, tree_abundance)) %>%
  # Remove observations without any species annotations
  filter(
    !is.na(species_original) & species_original != "" & species_original != "NA"
  ) %>%
  as_tibble() %>%
  glimpse()

# There are 4635 tree individuals in the FC dataset:
nrow(tree_list)
# There are 16 unique species in the FC dataset
unique(tree_list$species_original) %>%
  length()

# Subset genera
tree_genus <- tree_list %>%
  # Names without a space (i.e., one word)
  filter(!grepl("\\s", species_original) | grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# No genus only annotations
nrow(tree_genus)

# Subset species
tree_species <- tree_list %>%
  filter(grepl("\\s", species_original) & !grepl("sp\\.", species_original)) %>%
  distinct(species_original) %>%
  as.data.frame() %>%
  glimpse()
# There are 16 unique species annotations
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
# Direct matches: 11 out of 16
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 5 out of 16
tree_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

tree_list %>%
  inner_join(
    .,
    bind_rows(tree_species_WFO),
    by = "species_original"
  ) %>%
  select(
    family,	genus,	scientific_name = scientificName,	taxon_rank = taxonRank,
    tree_id, site,	longitude,	latitude,	date,	project, area = plot_area,
    diameter_cm,	height_m, total_biomass_kg, agb_drymass_kg, bgb_drymass_kg
  ) %>%
  filter(!is.na(scientific_name), !is.na(genus)) %>%
  fwrite("data/FORESTCHECK/harmonised_treelist.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

# (5) Harmonise HAVPlot ########################################################

#### (5a) Harmonise tree names ####

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Download the HAVPlot dataset:
dir.create("data/HAVPlot")
URL = "https://s3.data.csiro.au/dapprd/000054461v001/data/speciesAttributes.csv?response-content-disposition=attachment%3B%20filename%3DspeciesAttributes.csv&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20231017T062507Z&X-Amz-SignedHeaders=host&X-Amz-Expires=172800&X-Amz-Credential=2K6P74GE0OUPTMHV3N3A%2F20231017%2FCDC%2Fs3%2Faws4_request&X-Amz-Signature=99cb596aa32cfb8179df6777591bf65ef9ed1eb9d53295f28f84b7a54e4f07fb"
destfile <- "data/HAVPlot/Harmonised_Australian_Vegetation_Plot_dataset_(HAVPlot)-QEzDvqEq-.zip"
download.file(URL, destfile)

# Unzip the download
unzip(destfile, exdir = "data/HAVPlot/")

# Remove the zip file
file.remove(destfile)

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

#### (5b) Remove plots dominated by exotics ####

harmonised_tree_list_native <- HAVPlot_aggregateOrganismObservation %>%
  inner_join(
    .,
    harmonised_tree_list,
    by = "species_original",
    relationship = "many-to-many"
  ) %>%
  # Filter to tree species
  filter(
    scientific_name %in% fread("data/GlobalTreeSearch/harmonised_global_tree.txt")$scientific_name
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

#### (5c) Retain most recent survey ####

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

#### (5d) Generate abundance dataset ####

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

#### (5e) Generate presence dataset ####
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

# (6) Harmonise Global Register of Introduced and Invasive Species (GRIIS) #####

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Download the HAVPlot dataset:
dir.create("data/GRIIS")
URL = "https://cloud.gbif.org/griis/archive.do?r=griis-australia&v=1.10"
destfile <- "data/GRIIS/dwca-griis-australia-v1.10.zip"
download.file(URL, destfile)

# Unzip the download
unzip(destfile, exdir = "data/GRIIS/")

# Remove the zip file
file.remove(destfile)

# Read in the database:
griis_database <- fread(
  "data/GRIIS/taxon.txt"
) %>%
  filter(
    kingdom == "Plantae",
    taxonRank %in% c("SPECIES", "SUBSPECIES", "VARIETY")
  ) %>%
  select(species_original = scientificName) %>%
  mutate(
    # Handle hybrids
    species_original = gsub(" x ", " ×", species_original),
    # Remove 'subspecies' and 'variety' annotations by extracting the first two
    # words from species names
    species_original = map_chr(
      str_split(species_original, "\\s+"), 
      ~ paste(head(.x, 2), collapse = " ")
    )
  ) %>%
  distinct(species_original, .keep_all = TRUE) %>%
  as.data.frame() %>%
  glimpse()

# There are 2,641 unique species in the GRIIS database:
nrow(griis_database)

# Match species
griis_species_WFO <- WFO.one(
  WFO.match.fuzzyjoin(
    spec.data = griis_database,
    WFO.data = WFO_species,
    spec.name = "species_original"
  ),
  verbose = FALSE
) %>%
  glimpse()

# Check matching
# Direct matches: 2,597 out of 2,641
griis_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 29 out of 2,641
griis_species_WFO %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched: 15 out of 2,641
griis_species_WFO %>%
  filter(Matched == "FALSE") %>%
  select(species_original) %>%
  as_tibble() %>%
  print(n = Inf)

# Harmonised tree list
griis_tree_list <- griis_species_WFO %>%
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
  unique(.) %>%
  inner_join(
    fread("data/GlobalTreeSearch/harmonised_global_tree.txt") %>%
      select("scientific_name"),
    by = c("scientific_name")
  )

# Save the harmonised GRIIS tree list
griis_tree_list %>%
  fwrite("data/GRIIS/harmonised_griis_tree_list.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

# (7) Harmonise Australian Plant Census ########################################

# Download the Australian Plant Census dataset: https://biodiversity.org.au/

# Read in the WFO datasets
WFO_species <- fread("data/WorldFloraOnline/WFO_species.txt")
WFO_genus <- fread("data/WorldFloraOnline/WFO_genus.txt")

# Read in the database:
apc_database <- fread(
  "data/APC/APC-taxon-2025-03-18-0622.csv",
  header = TRUE
) %>%
  filter(taxonRank %in% c(
    "[infraspecies]", "[n/a]", "[unranked]", "Forma", "Nothovarietas",
    "Species", "Subforma", "Subspecies", "Subvarietas", "Varietas"
  )) %>%
  select(species_original = canonicalName) %>%
  mutate(
    # Handle hybrids
    species_original = gsub(" x ", " ×", species_original),
    # Remove 'subspecies' and 'variety' annotations by extracting the first two
    # words from species names
    species_original = map_chr(
      str_split(species_original, "\\s+"), 
      ~ paste(head(.x, 2), collapse = " ")
    )
  ) %>%
  distinct(species_original, .keep_all = TRUE) %>%
  as.data.frame() %>%
  glimpse()

# There are 62,556 unique species in the APC database:
nrow(apc_database)

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

# Direct matches: 54,415 out of 62,556
apc_database_WFO_combined %>%
  filter(Matched == TRUE & Fuzzy == FALSE) %>%
  nrow()

# Fuzzy matches: 6,951 out of 62,556
apc_database_WFO_combined %>%
  filter(Matched == TRUE & Fuzzy == TRUE) %>%
  select(species_original, Fuzzy.dist, scientificName, Old.name) %>%
  as_tibble() %>%
  print(n = Inf)

# Unmatched: 1,190 out of 62,556
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
    # Keep tree species only
    scientific_name %in% fread("data/GlobalTreeSearch/harmonised_global_tree.txt")$scientific_name,
    # Remove non-native species
    !scientific_name %in% fread("data/GRIIS/harmonised_griis_tree_list.txt")$scientific_name
  ) %>%
  unique(.)

# How many tree species are there: 3,802
n_distinct(apc_tree_database$scientific_name)

# Save the harmonised APC flora and tree list
apc_database_WFO_combined %>%
  select(family, genus, scientific_name = scientificName) %>%
  unique(.) %>%
  fwrite("data/APC/harmonised_apc_flora_list.txt", sep = "\t")
apc_tree_database %>%
  fwrite("data/APC/harmonised_apc_tree_list.txt", sep = "\t")
apc_tree_database %>%
  fwrite("output/generated_data/apc_tree_list.txt", sep = "\t")

# Remove all objects from the global environment and free unused memory
rm(list = ls())
gc()

# (8) Harmonise GBIF data ######################################################

#### (8a) Clean GBIF data ####

# Download the GBIF dataset used in this analysis: https://doi.org/10.15468/DL.ATKQFB

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

#### (8b) Harmonise GBIF tree list ####

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
    fread("data/APC/harmonised_apc_tree_list.txt") %>%
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
