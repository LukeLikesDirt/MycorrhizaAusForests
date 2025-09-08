
# Load libraries
require(terra)
require(data.table)
require(tidyverse)

# (1) Load data ################################################################

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

# Read in the harmonised abundance data from multiple sources
abundance_data <- bind_rows(
  fread("data/BiomassPlotLib/harmonised_treelist.txt"),
  fread("data/CSIRO_PermanentPlots_Data/harmonised_treelist.txt"),
  fread("data/FORESTCHECK/harmonised_treelist.txt"),
  fread("data/NaturalValuesAtlas/harmonised_treelist.txt")
) %>%
  select(
    scientific_name, longitude, latitude
  ) %>%
  unique()

# Read in the APC tree list with mycorrhizal types
apc_tree_list <- fread("output/generated_data/australian_tree_mycorrhizal_types.txt") %>%
  select(family, genus, scientific_name, mycorrhizal_type)

# (2) Occurence data ################################################################

# Combine the occurrence data with APC tree list of Australian native trees
trees <- bind_rows(
  HAVPlot_data,
  gbif_data,
  abundance_data
  ) %>%
  # Remove duplicates
  unique(.) %>%
  # Add family and genus information from the harmonised APC tree list and retain
  # only unique combinations of Australian native trees
  inner_join(
    .,
    apc_tree_list,
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

# Number of tree species: 3,427
unique(trees$scientific_name) %>%
  length() %>%
  message("Number of tree species: ", .)

# Number of sites: 1,374,342
unique(trees$site) %>%
  length() %>%
  message("Number of sites: ", .)

# Check for unassigned mycorrhizal types
trees %>%
  filter(is.na(mycorrhizal_type) | mycorrhizal_type == "" | mycorrhizal_type == "uncertain")

# (3) Site data ################################################################

#### Filter sites to forests ####

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

# Number of forest tree species: 3,287
forest_trees %>%
  select(scientific_name) %>%
  distinct() %>%
  nrow() %>%
  message("Number of tree species: ", .)

# Number of forest sites: 723,868
forest_sites %>%
  select(site) %>%
  distinct() %>%
  nrow() %>%
  message("Number of sites: ", .)

# Save forest tree data ########################################################

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
