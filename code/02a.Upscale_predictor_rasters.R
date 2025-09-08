
# Load libraries
require(terra)
require(data.table)
require(tidyverse)
source("code/map_australia.R")

# (1) Prepare forest layers ####################################################

# Mask the forest raster:
# I will eventually work with a 10km x 10km raster for downstream analysies.
# To upsacle predictor variable to 10km x 10kmn I need focus on values at 
# 100m x 100m resolution to focus on 'forest-related' values during up-scaling.
# Becuase the 100m x 100m raster is very big I need to mask grid cells to those
# covered by the final 10km x 10km raster for computation efficiency.

# Load the forest raster
forest_rast_unmasked <- rast(
  "data/aus_forests_23/aus_for23.tif"
)

# Extract the "Native forest" values from the raster
rat <- levels(forest_rast_unmasked)[[1]]
native_values <- rat$VALUE[rat$FOR_CATEGO == "Native forest"]

# Create a binary mask for "Native forest"
native_forest_rast <- classify(
  forest_rast_unmasked,
  rcl = cbind(native_values, 1),
  others = 0
) %>%
  rename(forest = FOR_CATEGO)

# Aggregate the raster to 10km x 10km grids using the sum function
native_forest_rast_10 <- aggregate(
  native_forest_rast,
  fact = 100,
  fun = sum,
  na.rm = TRUE
) %>%
  # Mutate forest into a percentage of cover
  mutate(forest = forest / 100)

# Mask 10km cells with less than 10% forest cover
native_forest_rast_10[native_forest_rast_10 < 10] <- NA

# Resample to match 100m resolution (without altering the original raster)
native_forest_rast_resampled <- resample(
  native_forest_rast_10, native_forest_rast,
  method = "near"
)

# Apply the mask to the original 100m raster
native_forest_rast_masked <- mask(
  native_forest_rast, native_forest_rast_resampled
)

# Finalise the masking: Convert 0s (non-native forest) to NA 
native_forest_rast_masked[native_forest_rast_masked == 0] <- NA

#### * 50km grids * ####

# Aggregate to 50 km (i.e. 500 × 100 m)
native_forest_rast_50 <- aggregate(
  native_forest_rast_masked,
  fact = 500,
  fun = sum,
  na.rm = TRUE
)

# Convert to percentage of native forest cover
# 500 x 500 = 250,000 cells per 50km grid cell
native_forest_rast_50_perc <- native_forest_rast_50 / 250000 * 100

# Mask out cells with less than 10% native forest cover
native_forest_rast_50_perc[native_forest_rast_50_perc < 10] <- NA

#### * 100km grids * ####

# Aggregate to 100 km (i.e. 1000 × 100 m)
native_forest_rast_100 <- aggregate(
  native_forest_rast_masked,
  fact = 1000,
  fun = sum,
  na.rm = TRUE
)

# Convert to percentage of native forest cover
# 1,000 x 1,000 = 1,000,000 cells per 100km grid cell
native_forest_rast_100_perc <- native_forest_rast_100 / 1000000 * 100

# Mask out cells with less than 10% native forest cover
native_forest_rast_100_perc[native_forest_rast_100_perc < 5] <- NA

# Save the masked rasters
writeRaster(
  native_forest_rast_masked,
  "data/aus_forests_23/aus_for23_masked.tif"
)
writeRaster(
  native_forest_rast_10,
  "data/aus_forests_23/aus_for23_masked_10.tif"
)
writeRaster(
  native_forest_rast_50_perc,
  "data/aus_forests_23/aus_for23_masked_50.tif"
)
writeRaster(
  native_forest_rast_100_perc,
  "data/aus_forests_23/aus_for23_masked_100.tif"
)

# Clean up the environment and free unused memory
rm(list = ls())
gc()

# (2) Extract predictors #######################################################

# Load libraries
require(terra)
require(data.table)
require(tidyverse)
source("code/map_australia.R")

#### (2.1) Forest rasters ####

# Read in the masked rasters: 
# 100m resolution is for up-scaling to 10km resolution and the 10km resolution
# is for efficiency in projecting the predictors to Albers.
native_forest_rast <- rast(
  "data/aus_forests_23/aus_for23_masked.tif"
)
native_forest_rast_10 <- rast(
  "data/aus_forests_23/aus_for23_masked_10.tif"
)

# Generate a vector for data extraction for continuous predictors.
# Use only values extracted from forested areas because I don't want values 
# from non-forested areas to skew the predictors.
native_forest_vect <- native_forest_rast %>%
  as.points()

# Generate a vector for data extraction for categorical predictors and for 
# habitat condition.
# Use values from all areas within a cell (i.e. not just forested areas) so to 
# gauge the surrounding habitat condition. This is because poor habitat 
# condition in areas with high human activity influences the vunerability of
# trees.
native_forest_vect_10 <- native_forest_rast_10 %>%
  as.points()

#### (2.2) BIOCLIM rasters ####

# List of BIOCLIM raster files
bio_files <- c(
  "data/georeferenced_predictors/bioclim/bio1.tif",
  "data/georeferenced_predictors/bioclim/bio4.tif",
  "data/georeferenced_predictors/bioclim/bio8.tif",
  "data/georeferenced_predictors/bioclim/bio9.tif",
  "data/georeferenced_predictors/bioclim/bio12.tif",
  "data/georeferenced_predictors/bioclim/bio15.tif",
  "data/georeferenced_predictors/bioclim/bio18.tif",
  "data/georeferenced_predictors/bioclim/bio19.tif",
  "data/georeferenced_predictors/bioclim/bio28.tif",
  "data/georeferenced_predictors/bioclim/bio20.tif",
  "data/georeferenced_predictors/bioclim/bio23.tif",
  "data/georeferenced_predictors/bioclim/bio31.tif",
  "data/georeferenced_predictors/bioclim/bio34.tif",
  "data/georeferenced_predictors/bioclim/bio35.tif"
)

# Function to process each raster file
process_raster <- function(file_path, native_forest_vect, native_forest_rast_10) {
  
  # Extract bioclim values
  predictor_values <-  rast(file_path) %>%
    project("EPSG:3577", threads = TRUE) %>%
    terra::extract(
      native_forest_vect,
      xy = TRUE
    )
  
  # Vectorise extracted soil values
  predictor_vect <- vect(
    predictor_values,
    geom = c("x", "y"),
    crs = crs(native_forest_rast)
  )
  
  # Upscale the data to 10km resolution
  predictor_rast <- rasterize(
    predictor_vect,
    native_forest_rast_10,
    field = tools::file_path_sans_ext(basename(file_path)),
    fun = mean,
    na.rm = TRUE
  )
  
  # Save the upscaled raster
  writeRaster(
    predictor_rast,
    paste0("data/aus_forests_23/", tools::file_path_sans_ext(basename(file_path)), "_10.tif"),
    overwrite = TRUE
  )
  
}

# Apply the function to each file in bio_files
lapply(bio_files, process_raster, native_forest_vect, native_forest_rast_10)

#### (2.3) Soil Grid rasters ####

# List of soil raster files
soil_files <- c(
  "data/georeferenced_predictors/soil_grid/bulk_density_15cm.tif",
  "data/georeferenced_predictors/soil_grid/bulk_density_5cm.tif",
  "data/georeferenced_predictors/soil_grid/CEC_15cm.tif",
  "data/georeferenced_predictors/soil_grid/CEC_5cm.tif",
  "data/georeferenced_predictors/soil_grid/ECE_15cm.tif",
  "data/georeferenced_predictors/soil_grid/ECE_5cm.tif",
  "data/georeferenced_predictors/soil_grid/N_total_15cm.tif",
  "data/georeferenced_predictors/soil_grid/N_total_5cm.tif",
  "data/georeferenced_predictors/soil_grid/P_total_15cm.tif",
  "data/georeferenced_predictors/soil_grid/P_total_5cm.tif",
  "data/georeferenced_predictors/soil_grid/pH_CaCl2_15cm.tif",
  "data/georeferenced_predictors/soil_grid/pH_CaCl2_5cm.tif",
  "data/georeferenced_predictors/soil_grid/SOC_15cm.tif",
  "data/georeferenced_predictors/soil_grid/SOC_5cm.tif",
  "data/georeferenced_predictors/soil_grid/MAOC_proportion_5cm.tif",
  "data/georeferenced_predictors/soil_grid/MAOC_proportion_15cm.tif",
  "data/georeferenced_predictors/soil_grid/POC_proportion_5cm.tif",
  "data/georeferenced_predictors/soil_grid/POC_proportion_15cm.tif",
  "data/georeferenced_predictors/soil_grid/PyOC_proportion_5cm.tif",
  "data/georeferenced_predictors/soil_grid/PyOC_proportion_15cm.tif",
  "data/georeferenced_predictors/soil_grid/P_available_15cm.tif",     
  "data/georeferenced_predictors/soil_grid/P_available_5cm.tif" 
)

# Lookup table for renaming the raster files
rename_lookup <- list(
  "bulk_density_15cm.tif" = "bulk_density_15cm",
  "bulk_density_5cm.tif" = "bulk_density_5cm",
  "CEC_15cm.tif" = "CEC_15cm",
  "CEC_5cm.tif" = "CEC_5cm",
  "ECE_15cm.tif" = "ECE_15cm",
  "ECE_5cm.tif" = "ECE_5cm",
  "N_total_15cm.tif" = "N_total_15cm",
  "N_total_5cm.tif" = "N_total_5cm",
  "P_total_15cm.tif" = "P_total_15cm",
  "P_total_5cm.tif" = "P_total_5cm",
  "pH_CaCl2_15cm.tif" = "pH_CaCl2_15cm",
  "pH_CaCl2_5cm.tif" = "pH_CaCl2_5cm",
  "SOC_15cm.tif" = "SOC_15cm",
  "SOC_5cm.tif" = "SOC_5cm",
  "MAOC_proportion_5cm.tif" = "MAOC_proportion_5cm",
  "MAOC_proportion_15cm.tif" = "MAOC_proportion_15cm",
  "POC_proportion_5cm.tif" = "POC_proportion_5cm",
  "POC_proportion_15cm.tif" = "POC_proportion_15cm",
  "PyOC_proportion_5cm.tif" = "PyOC_proportion_5cm",
  "PyOC_proportion_15cm.tif" = "PyOC_proportion_15cm",
  "P_available_15cm.tif" = "P_available_15cm",
  "P_available_5cm.tif" = "P_available_5cm"
)

# Function to process each raster file
process_raster <- function(file_path, native_forest_vect, native_forest_rast_10) {
  # Extract the base name of the file
  file_name <- basename(file_path)
  # Get the corresponding field name from the lookup table
  field_name <- rename_lookup[[file_name]]
  
  # Extract soil grid values
  predictor_values <- rast(file_path) %>%
    # Project to Albers equal area
    project("EPSG:3577", threads = TRUE) %>%
    terra::extract(
      native_forest_vect,
      xy = TRUE
    )
  
  # Rename the second column to the value in 'field_name' (column 1 is ID)
  names(predictor_values)[2] <- field_name
  
  # Vectorise extracted soil grid values
  predictor_vect <- vect(
    predictor_values,
    geom = c("x", "y"),
    crs = crs(native_forest_rast)
  )
  
  # Upscale the data to 10km resolution
  predictor_rast <- rasterize(
    predictor_vect,
    native_forest_rast_10,
    field = field_name,
    fun = mean,
    na.rm = TRUE
  )
  
  # Save the upscaled raster
  writeRaster(
    predictor_rast,
    paste0("data/aus_forests_23/", tools::file_path_sans_ext(file_name), "_10.tif"),
    overwrite = TRUE
  )
}

# Apply the function to each file in soil_files
lapply(soil_files, process_raster, native_forest_vect, native_forest_rast_10)

#### (2.4) Elevation ####

# Extract elevation values
elevation_values <- rast("data/georeferenced_predictors/elevation.tif") %>%
  project("EPSG:3577", threads = TRUE) %>%
  terra::extract(
    native_forest_vect,
    xy = TRUE
  )

# Vectorise extracted elevation values
elevation_vect <- vect(
  elevation_values %>% filter(!is.na(elevation)),
  geom = c("x", "y"),
  crs = crs(native_forest_rast)
)

# Upscale the data to 10km resolution
elevation_rast <- rasterize(
  elevation_vect,
  native_forest_rast_10,
  field = "elevation",
  fun = mean,
  na.rm = TRUE
) %>%
  rename(elevation = mean)

# Save the upscaled raster
writeRaster(
  elevation_rast,
  "data/aus_forests_23/elevation_10.tif",
  overwrite = TRUE
)

#### (2.5) Woody veg cover ####

# Woody veg cover values
wood_veg_values <- rast(
  "data/georeferenced_predictors/vegetation_structure/vegetation_cover.tif"
) %>%
  project("EPSG:3577", threads = TRUE) %>%
  terra::extract(
    native_forest_vect,
    xy = TRUE
  ) %>% filter(
    !is.na(fcov_totalc1000)
  ) %>%
  mutate(
    wood_veg_cover = (fcov_totalc1000 * 0.01)
  ) %>%
  select(x, y, wood_veg_cover)

# Vectorise extracted elevation values
wood_veg_vect <- vect(
  wood_veg_values,
  geom = c("x", "y"),
  crs = crs(native_forest_rast)
)

# Upscale the data to 10km resolution
wood_veg_rast <- rasterize(
  wood_veg_vect,
  native_forest_rast_10,
  field = "wood_veg_cover",
  fun = mean,
  na.rm = TRUE
) %>%
  rename(wood_veg_cover = mean)

# Save the upscaled raster
writeRaster(
  wood_veg_rast,
  "data/aus_forests_23/wood_veg_10.tif",
  overwrite = TRUE
)

#### (2.6) Climate zones ####

# Read and project the climate zone raster
climate_zone_rast <- rast(
  "data/georeferenced_predictors/kpn/kpnall.txt"
)
crs(climate_zone_rast) <- "EPSG:4283"

# Check the unique values
unique(climate_zone_rast[["kpnall"]])

# Create a template raster at 10km resolution
climate_zone_rast_10 <- native_forest_rast_10

# Get cell centre points for the 10km raster
climate_zone_vect_10 <- as_tibble(
  xyFromCell(climate_zone_rast_10, 1:ncell(climate_zone_rast_10))
) %>%
  vect(
    geom = c("x", "y"),
    crs = crs(climate_zone_rast_10)
  )

# Project the climate zone raster
climate_zone_rast_projected <- project(
  climate_zone_rast,
  crs(climate_zone_rast_10),
  method = "near"
)

# Resample to match rasters
climate_zone_rast_resampled <- resample(
  climate_zone_rast_projected,
  native_forest_rast_10,
  method = "near"
)

# Mask the climate zone raster to match the forest raster
climate_zone_rast_masked <- mask(
  climate_zone_rast_resampled,
  native_forest_rast_10
)

# Check the plot
plot(climate_zone_rast_masked)

# Check the unique values
unique(climate_zone_rast_masked[["kpnall"]])

# Mutate the climate zone values
climate_zone_10 <- climate_zone_rast_masked %>%
  tidyterra::mutate(
    kpnall = case_when(
      # Equatorial
      kpnall %in% c(41, 42) ~ "Equatorial",
      # Tropical
      between(as.numeric(kpnall), 35, 37) ~ "Tropical",
      # Subtropical
      between(as.numeric(kpnall), 31, 34) ~ "Subtropical",
      # Desert
      between(as.numeric(kpnall), 21, 24) ~ "Desert",
      # Grassland
      between(as.numeric(kpnall), 11, 15) ~ "Grassland",
      # Temperate
      between(as.numeric(kpnall), 1, 9) ~ "Temperate",
      TRUE ~ NA_character_
    ),
    # Level the factor
    kpnall = factor(kpnall, levels = c(
      "Equatorial", "Tropical", "Subtropical", "Desert", "Grassland", "Temperate"
    ))
  ) %>%
  rename(climate_zone = kpnall)
  
# Check the plot
plot(climate_zone_10)

# Save the final climate zone raster
writeRaster(
  climate_zone_10,
  "data/aus_forests_23/climate_zones_10.tif",
  overwrite = TRUE
)

#### (2.7) Habitat condition raster ####

# Read and project the habitat condition raster
habitat_condition_raster <- rast(
  "data/georeferenced_predictors/habitat_condition/HCAS30_HCB_1988_2022.tif"
) %>%
  project("EPSG:3577", threads = TRUE)

# Direct aggregation (most efficient)
# This aggregates the 250m pixels to 10km resolution by taking the mean
# First, we need to determine the aggregation factor
# 10km / 250m = 40, so we aggregate by factor of 40

habitat_condition_10km <- habitat_condition_raster %>%
  aggregate(fact = 40, fun = mean, na.rm = TRUE)
plot(habitat_condition_10km)

# Now resample to match exactly the native_forest_rast_10 grid
habitat_condition_aligned <- habitat_condition_10km %>%
  resample(native_forest_rast_10, method = "bilinear")

# Mask to only forest areas
habitat_condition_forest <- mask(
  habitat_condition_aligned, 
  native_forest_rast_10
  ) %>%
  rename(habitat_condition = HCAS30_HCB_1988_2022)

# Plot the final habitat condition raster
plot(habitat_condition_forest)

# Save the final habitat condition raster
writeRaster(
  habitat_condition_forest,
  "data/aus_forests_23/habitat_condition_10.tif",
  overwrite = TRUE
)

#### (2.8) Habitat type ####

# Load habitat type rasters with keyword filtering
habitat_type_rast <- list.files(
  "data/habitat/lvl2_frac_1km_ver004",
  pattern = "iucn_habitatclassification.*(?i)(Forest|Grassland|Shrubland|Savanna).*ver004\\.tif$",
  full.names = TRUE
) %>%
  rast()

# Project the forest raster to match the habitat raster CRS (more efficient)
native_forest_rast_projected <- project(
  native_forest_rast_10,
  crs(habitat_type_rast),
  threads = TRUE
)

# Crop habitat raster to the projected forest extent to reduce processing area
habitat_type_rast_cropped <- crop(habitat_type_rast, native_forest_rast_projected)

# Sum habitat values within each 10km native forest cell
habitat_type_rast_10 <- resample(
  habitat_type_rast_cropped, 
  native_forest_rast_projected, 
  method = "sum"
  ) %>%
  # Project back to the original forest raster CRS and grid
  project(
    .,
    native_forest_rast_10,
    threads = TRUE
    )

# Crop and mask to match native forest extent exactly
habitat_type_rast_10 <- crop(habitat_type_rast_10, native_forest_rast_10)
habitat_type_rast_10 <- mask(habitat_type_rast_10, native_forest_rast_10)

# Get habitat names from the raster
habitat_names <- names(habitat_type_rast_10)

# Extract clean habitat names from layer names (between __ and __ver004)
clean_names <- gsub(".*__\\d+_(.+)__ver004.*", "\\1", habitat_names)

# Keep only minor classes (those with " – " or " - ")
minor_mask <- grepl(" – | - ", clean_names)
minor_classes <- habitat_names[minor_mask]
minor_clean_names <- clean_names[minor_mask]

# Create raster with only minor classes
minor_habitat_rast <- habitat_type_rast_10[[minor_classes]]

# Find dominant minor class in each cell
dominant_minor <- which.max(minor_habitat_rast)
dominant_minor_names <- minor_clean_names[values(dominant_minor)]

# Get maximum values to identify where minor classes have data
minor_max_values <- app(minor_habitat_rast, max, na.rm = TRUE)

# Create final raster with only cells that have minor habitat data
dominant_habitat <- habitat_type_rast_10[[1]] * 0  # Create template raster
values(dominant_habitat) <- ifelse(
  values(minor_max_values) > 0 & !is.na(dominant_minor_names),
  dominant_minor_names,
  NA
)
names(dominant_habitat) <- "habitat_type"

# Check the plot
plot(dominant_habitat)

# Save the final climate zone raster
writeRaster(
  dominant_habitat,
  "data/aus_forests_23/habitat_type_10.tif",
  overwrite = TRUE
)

# Clean up the environment and free unused memory
rm(list = ls())
gc()
