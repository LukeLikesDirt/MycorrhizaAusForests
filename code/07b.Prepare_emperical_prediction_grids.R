require(terra)
require(data.table)
require(tidyverse)
source("code/functions.R")
source("code/map_australia.R")

# (1) Tree data ################################################################

# Load the tree data
trees <- data.table::fread(
  "data/presence/trees.txt"
) %>%
  mutate(
    mycorrhizal_type = case_when(
      mycorrhizal_type == 'NM-AM' ~ "AM",
      TRUE ~ mycorrhizal_type
    )
  ) %>%
  # Filter to trees in the emperical dataset
  filter(scientific_name %in% fread(
    "generated_data/emperical_mycorrhizal_type.txt")$species
    )

# Load the site data
sites <- fread("data/presence/sites.txt") %>%
  select(site, longitude, latitude, n_obs) %>%
  # Filter to the sites in the tree dataset
  filter(site %in% trees$site)

# Load the combined raster with predictors
combined_raster <- rast(
  "data/aus_forests_23/predictors_10.tif"
)

# Vectorise the site data
sites_vect <- sites %>% 
  vect(geom = c("longitude", "latitude"), crs = "EPSG:4326") %>%
  project(crs(combined_raster)) %>%
  # Crop to the native forest raster extent
  crop(combined_raster)

# Extract the predictors
sites_predictors <- terra::extract(
  combined_raster,
  sites_vect,
  cells = TRUE,
  xy = TRUE
) %>%
  inner_join(
    sites %>% rownames_to_column("ID") %>% mutate(ID = as.numeric(ID)),
    by = "ID"
  ) %>%
  select(site, n_obs, x, y, climate_zone, biome, ecoregion, habitat_type, everything()) %>%
  as_tibble() %>%
  # Fill missing biomes based on nearest neighbours
  fill_missing_knn(
    x = "x",
    y = "y",
    variable = "biome"
  ) %>%
  # Fill missing ecoregions based on nearest neighbours
  fill_missing_knn(
    x = "x",
    y = "y",
    variable = "ecoregion"
  ) %>%
  # Fill missing climate zones based on nearest neighbours
  fill_missing_knn(
    x = "x",
    y = "y",
    variable = "climate_zone"
  ) %>%
  # Fill missing habitat types based on nearest neighbours
  fill_missing_knn(
    x = "x",
    y = "y",
    variable = "habitat_type"
  ) %>%
  print()

# Get GDA 94 longitude and latitude for each grid
sites_predictors <- sites_predictors %>%
  select(site, x, y) %>%
  vect(geom = c("x", "y"), crs = crs(combined_raster)) %>%
  # Project to WGS84 (EPSG:4326)
  project(crs(aus_map$data)) %>%
  # Extract coordinates
  crds() %>%
  # Rename coordinates
  as_tibble() %>%
  rename(longitude = x, latitude = y) %>%
  bind_cols(
    sites_predictors %>% select(-c(ID, longitude, latitude))
  ) %>%
  select(
    site, n_obs, longitude, latitude, x_albers = x, y_albers = y,
    climate_zone, biome, ecoregion, habitat_type, everything()
  )

#### (1.1) Sites for relative richness ####

# Upscale sites to grid cells
sites_upscaled <- sites_predictors %>%
  # Prior to adding tree information, determine n_obs per grid and select the 
  # most common climate_zone, biome, ecoregion and habitat_type
  # Sum n_obs per grid cell
  group_by(cell) %>%
  mutate(
    n_obs = sum(n_obs)
  ) %>%
  ungroup() %>%
  # Compute most common climate_zone
  group_by(cell, climate_zone) %>%
  mutate(
    most_common_climate_zone = n()
  ) %>%
  ungroup() %>%
  # Compute most common biome
  group_by(cell, biome) %>%
  mutate(
    most_common_biome = n()
  ) %>%
  ungroup() %>%
  # Compute most common ecoregion
  group_by(cell, ecoregion) %>%
  mutate(
    most_common_ecoregion = n()
  ) %>%
  ungroup() %>%
  # Compute most common habitat type
  group_by(cell, habitat_type) %>%
  mutate(
    most_common_habitat_type = n()
  ) %>%
  # Add tree information for richness estimates
  left_join(
    trees %>% select(site, scientific_name, mycorrhizal_type),
    by = "site"
  ) %>%
  # Average across grids and estimate richness
  group_by(cell) %>%
  summarise(
    across(
      where(is.numeric) &
        !contains("richness") &
        !contains("n_obs") &
        !contains("most_common_"),
      \(x) mean(x, na.rm = TRUE)
    ),
    # Total species richness
    richness = n_distinct(scientific_name),
    # Species richness for each mycorrhizal type
    richness_AM = n_distinct(scientific_name[mycorrhizal_type == "AM"], na.rm = TRUE),
    richness_EcM = n_distinct(scientific_name[mycorrhizal_type == "EcM"], na.rm = TRUE),
    richness_EcM_AM = n_distinct(scientific_name[mycorrhizal_type == "EcM-AM"], na.rm = TRUE),
    richness_NM = n_distinct(scientific_name[mycorrhizal_type == "NM"], na.rm = TRUE),
    richness_ErM = n_distinct(scientific_name[mycorrhizal_type == "ErM"], na.rm = TRUE),
    n_sites = n_distinct(site),
    # n_obs is summed above in the mutate() function, so take the first value
    n_obs = first(n_obs),
    forest_cover = first(forest_cover),
    climate_zone = climate_zone[which.max(most_common_climate_zone)],
    biome = biome[which.max(most_common_ecoregion)],
    ecoregion = ecoregion[which.max(most_common_ecoregion)],
    habitat_type = habitat_type[which.max(most_common_habitat_type)],
    .groups = "drop"
  ) %>%
  select(
    cell, longitude, latitude, x_albers, y_albers, forest_cover, n_obs, n_sites,
    starts_with("richness"), climate_zone, biome, ecoregion, habitat_type,
    everything()
  ) %>%
  print()

# Cells without 'forest' represent forest sites that occurred in cells with
# forest cover less than 10%
nrow(sites_upscaled)
nrow(sites_upscaled %>% filter(!is.na(RC1), forest_cover >= 10))
# 13,605 cells have forest cover > 10% and are paired with predictor variables
# Cells with at least 10 observations
nrow(sites_upscaled %>% filter(!is.na(RC1), forest_cover >= 10, n_obs >= 10))
# 9,760 cells have forest cover > 10% and at least 10 observations

# Save the data for richness estimation
sites_upscaled %>% 
  filter(!is.na(RC1), forest_cover >= 10, n_obs >= 10) %>%
  select(x_albers, y_albers) %>%
  # Inner join with the upscaled sites
  inner_join(
    sites_upscaled,
    by = c("x_albers", "y_albers")
  ) %>%
  # Mutate the cell names
  mutate(cell = paste0("cell_", cell)) %>%
  # Order the columns
  select(
    cell, longitude, latitude, x_albers, y_albers, forest_cover, n_obs,
    starts_with("richness"), climate_zone, biome, ecoregion, habitat_type,
    everything()
  ) %>%
  fwrite("data/presence/sites_relative_richness_est_10_emp.txt", sep = "\t")

# Save the richness prediction data
fread("data/aus_forests_23/predictors_10.txt") %>%
  # Fill missing biomes based on nearest neighbours
  fill_missing_knn(
    x = "x_albers",
    y = "y_albers",
    variable = "biome"
  ) %>%
  # Fill missing ecoregions based on nearest neighbours
  fill_missing_knn(
    x = "x_albers",
    y = "y_albers",
    variable = "ecoregion"
  ) %>%
  # Fill missing climate zones based on nearest neighbours
  fill_missing_knn(
    x = "x_albers",
    y = "y_albers",
    variable = "climate_zone"
  ) %>%
  # Fill missing habitat types based on nearest neighbours
  fill_missing_knn(
    x = "x_albers",
    y = "y_albers",
    variable = "habitat_type"
  ) %>%
  anti_join(
    .,
    fread("data/presence/sites_relative_richness_est_10_emp.txt"),
    by = c("x_albers", "y_albers")
  ) %>%
  # Order the columns
  select(
    longitude, latitude, x_albers, y_albers, forest_cover,
    starts_with("richness"), climate_zone, biome, ecoregion,
    everything()
  ) %>%
  fwrite("data/presence/sites_relative_richness_pred_10_emp.txt", sep = "\t")

# Save the small sites for niche and biographical modelling
sites_upscaled %>%
  # Filter out sites without forest cover and missing RCs values
  filter(!is.na(RC1), forest_cover >= 10) %>%
  # Remove sites in the estimation data
  anti_join(
    .,
    fread("data/presence/sites_relative_richness_est_10_emp.txt") %>%
      select(x_albers, y_albers),
    by = c("x_albers", "y_albers")
  ) %>%
  # Mutate the cell names
  mutate(cell = paste0("cell_", cell)) %>%
  # Select the columns
  select(
    cell, longitude, latitude, x_albers, y_albers, forest_cover, n_obs, n_sites,
    starts_with("richness"), climate_zone, biome, ecoregion,
    everything()
  ) %>%
  fwrite("data/presence/sites_10_under_sampled_emp.txt", sep = "\t")

#### (1.2) Trees data for niche models ####

# Upscale tree data to 10km x 10km grid cells
trees_upscaled <- sites_predictors %>%
  # Remove non-forest sites and sites without RCs values
  filter(forest_cover >= 10, !is.na(RC1)) %>%
  mutate(cell = paste0("cell_", cell)) %>%
  group_by(cell) %>%
  mutate(
    n_obs = sum(n_obs)
  ) %>%
  # Add the most common climate zone, biome, ecoregion and habitat type
  group_by(cell, climate_zone) %>%
  mutate(
    most_common_climate_zone = n()
  ) %>%
  ungroup() %>%
  group_by(cell, biome) %>%
  mutate(
    most_common_biome = n()
  ) %>%
  ungroup() %>%
  group_by(cell, ecoregion) %>%
  mutate(
    most_common_ecoregion = n()
  ) %>%
  ungroup() %>%
  group_by(cell, habitat_type) %>%
  mutate(
    most_common_habitat_type = n()
  ) %>%
  ungroup() %>%
  # Add tree information for richness estimates
  left_join(
    trees %>% select(-longitude, -latitude),
    by = "site"
  ) %>%
  # Average across grids
  group_by(cell, scientific_name) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    family = first(family),
    genus = first(genus),
    mycorrhizal_type = first(mycorrhizal_type),
    climate_zone = climate_zone[which.max(most_common_climate_zone)],
    biome = biome[which.max(most_common_biome)],
    ecoregion = ecoregion[which.max(most_common_ecoregion)],
    habitat_type = habitat_type[which.max(most_common_habitat_type)],
    # n_obs is summed above in the mutate() function, so take the first value
    n_obs = first(n_obs),
    n_sites = n_distinct(site),
    .groups = "drop"
  ) %>%
  select(
    cell, family, genus, scientific_name, mycorrhizal_type,
    longitude, latitude, x_albers, y_albers, forest_cover, n_obs,
    climate_zone, biome, ecoregion, habitat_type, everything()
  ) %>%
  print()

# Save the upscale tree data
fwrite(
  trees_upscaled,
  "data/presence/trees_10_emp.txt",
  sep = "\t"
)

# Count the number of unique species
n_distinct(trees_upscaled$scientific_name)

# Count the number of unique grid cells
n_distinct(trees_upscaled$cell)

# Count the number of unique species in at least 10 grid cells
trees_upscaled %>%
  group_by(scientific_name) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 10) %>%
  nrow()
