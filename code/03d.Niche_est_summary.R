# Load libraries
library(fst)
library(terra)
library(data.table)
library(tidyverse)
source("code/functions.R")

# Load and prepare data ########################################################

# Load forest raster with relevant layers
forest_rast <- rast("data/aus_forests_23/predictors_10.tif") %>%
  tidyterra::select(climate_zone, biome, ecoregion, habitat_type)

# Load species and niche data
tree_species <- fread("data/presence/trees_10.txt") %>%
  select(family, genus, species = scientific_name, mycorrhizal_type, latitude) %>%
  filter(mycorrhizal_type != "ErM") %>%
  # Compute latitudinal range size based on 95% quantiles
  group_by(family, genus, species, mycorrhizal_type) %>%
  summarise(
    lat_min = quantile(abs(latitude), 0.025, na.rm = TRUE),
    lat_max = quantile(abs(latitude), 0.975, na.rm = TRUE),
    lat_range = lat_max - lat_min
  ) %>%
  select(family, genus, species, mycorrhizal_type, lat_range)

niche_breadth_est <- read_fst("data/niche_estimates_enmeval/niche_model_results_enmeval.fst", as.data.table = TRUE) %>%
  filter(pass_validation == TRUE, auc >= 0.6) %>%
  inner_join(tree_species, by = "species")

cat("Median AUC:", median(niche_breadth_est$auc, na.rm = TRUE), "± SD:", sd(niche_breadth_est$auc, na.rm = TRUE), "\n")
cat("Median omission:", median(niche_breadth_est$omission, na.rm = TRUE), "± SD:", sd(niche_breadth_est$omission, na.rm = TRUE), "\n")

# Compute ecoregions and biomes

# Load all tree occurrences and extract forest data
tree_occurrences <- fread("data/presence/trees_10.txt") %>%
  filter(mycorrhizal_type != "ErM", !is.na(latitude), !is.na(longitude)) %>%
  rename(species = scientific_name)

# Generate points from tree occurrences
tree_points <- vect(tree_occurrences, geom = c("longitude", "latitude"), crs = "EPSG:4326") %>%
  project(forest_rast)

# Compute species biomes and ecoregions
species_biome <- terra::extract(forest_rast, tree_points, xy = TRUE) %>%
  bind_cols(tree_occurrences[, .(family, genus, species, mycorrhizal_type)]) %>%
  # Filter to species in niche_breadth_est
  filter(species %in% niche_breadth_est$species) %>%
  fill_missing_knn(x = "x", y = "y", variable = "biome") %>%
  fill_missing_knn(x = "x", y = "y", variable = "ecoregion") %>%
  fill_missing_knn(x = "x", y = "y", variable = "habitat_type") %>%
  fill_missing_knn(x = "x", y = "y", variable = "climate_zone") %>%
  mutate(
    biome_group = case_when(
      biome %in% c("Tropical & Subtropical Grasslands, Savannas & Shrublands", 
                   "Tropical & Subtropical Moist Broadleaf Forests") ~ "Tropical",
      TRUE ~ "Non-tropical"
    )
  )

# Check the distribution of biomes
dominant_biome_group <- species_biome %>%
  group_by(species, biome_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(desc(prop)) %>%
  slice(1) %>%
  select(species, biome_group)
dominant_biome_group %>%
  group_by(biome_group) %>%
  count(biome_group)

# Count species by biome
species_biome %>%
  group_by(species, biome) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(desc(prop)) %>%
  slice(1) %>%
  select(species, dominant_biome = biome) %>%
  mutate(biome_group = case_when(
    str_detect(dominant_biome, "Tropical") ~ "Tropical",
    str_detect(dominant_biome, "Temperate") ~ "Temperate",
    str_detect(dominant_biome, "Mediterranean") ~ "Mediterranean",
    str_detect(dominant_biome, "Montane") ~ "Montane",
    str_detect(dominant_biome, "Deserts") ~ "Desert",
    TRUE ~ "Other"
  )) %>%
  group_by(biome_group) %>%
  count(biome_group) %>%
  ungroup() %>%
  arrange(desc(n))

# Save the niche estimates
niche_breadth_est %>%
  inner_join(
    dominant_biome_group, 
    by = c("species") 
    ) %>%
  select(
    family, genus, species, mycorrhizal_type,
    AOO, EOO,
    env_B2, geo_B2, ex_dent, 
    env_B2_corrected, geo_B2_corrected,
    RC1_position, RC2_position, RC3_position,
    lat_position,
    lat_range,
    biome = biome_group
  ) %>%
  # Compute "Tropical" and "Non-tropical" climate zones based on lat_position
  mutate(
    climate_zone = case_when(
      lat_position > -23.45 ~ "Tropical",
      lat_position < -23.45 ~ "Non-tropical"
    )
  ) %>%
  fwrite("data/niche_estimates_enmeval/niche_estimates.txt", sep = "\t")

# Also save to generated_data for github repository                 
niche_breadth_est %>%
  inner_join(
    dominant_biome_group, 
    by = c("species") 
  ) %>%
  select(
    family, genus, species, mycorrhizal_type,
    AOO, EOO,
    env_B2, geo_B2, ex_dent, 
    env_B2_corrected, geo_B2_corrected,
    RC1_position, RC2_position, RC3_position,
    lat_position,
    lat_range,
    biome = biome_group
  ) %>%
  # Compute "Tropical" and "Non-tropical" climate zones based on lat_position
  mutate(
    climate_zone = case_when(
      lat_position > -23.45 ~ "Tropical",
      lat_position < -23.45 ~ "Non-tropical"
    )
  ) %>%
  fwrite("output/generated_data/niche_estimates.txt", sep = "\t")                    
                    