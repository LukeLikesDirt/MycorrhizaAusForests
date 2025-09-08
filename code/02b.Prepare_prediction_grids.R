
# Load libraries
require(psych)
require(terra)
require(viridis)
require(data.table)
require(ggtext)
require(tidyverse)
source("code/functions.R")
source("code/map_australia.R")

# (1) Rotated component analysis ###############################################

#### (1.1) Forest cover ####

# Read in the masked rasters: 
native_forest_rast_10 <- rast(
  "data/aus_forests_23/aus_for23_masked_10.tif"
)

# Convert to a tibble with cell IDs
native_forest_values <- native_forest_rast_10 %>%
  as_tibble(xy = TRUE) %>%
  na.omit() %>%
  rename(forest_cover = forest)

# Check predictor distributions:
native_forest_values %>% 
  ggplot(aes(x = forest_cover)) +
  geom_density()
native_forest_values %>%
  ggplot(aes(x = log10(forest_cover))) +
  geom_density()

#### (1.2) Biomes, Ecoregions and Bioregions ####

# Read and project the ecoregions vector
ecoregions_vect <- vect(
  "data/georeferenced_predictors/Ecoregions2017/Ecoregions2017.shp"
) %>%
  select(
    biome = BIOME_NAME,
    ecoregion = ECO_NAME
  ) %>%
  project(crs(native_forest_rast_10))

# Rasterise the ecoregion vector
ecoregions_raster <- terra::rasterize(
  x = ecoregions_vect,
  y = native_forest_rast_10,
  field = "ecoregion",
)

# Rasterise the biome vector
biomes_raster <- terra::rasterize(
  x = ecoregions_vect,
  y = native_forest_rast_10,
  field = "biome",
)

# Read and project the ecoregions vector
bioregions_vect <- vect(
  "data/georeferenced_predictors/IBRA7_regions/ibra7_regions.shp"
) %>%
  select(
    bioregion = REG_NAME_7
  ) %>%
  project(crs(native_forest_rast_10))

# Rasterise the ecoregion vector
bioregions_raster <- terra::rasterize(
  x = bioregions_vect,
  y = native_forest_rast_10,
  field = "bioregion",
)

#### (1.3) Climate zones ####

# Read climate zones vector
climate_zone_rast <- rast(
  "data/aus_forests_23/climate_zones_10.tif"
)

#### (1.4) Habitat type ####

# Read in the habitat type raster
habitat_type_rast <- rast(
  "data/aus_forests_23/habitat_type_10.tif"
)

#### (1.5) BIOCLIM ####

# Load the BIOCLIM rasters
bio_rasters <- list.files(
  "data/aus_forests_23/",
  pattern = "bio\\d+\\_10.tif$",
  full.names = TRUE
) %>%
  rast()

# Assign variable names
names(bio_rasters) <- sub(
  "(bio\\d+)_.*", "\\1",
  basename(list.files(
    "data/aus_forests_23/",
    pattern = "bio\\d+\\_10.tif$",
    full.names = TRUE
  ))
)

# Convert to a tibble with cell IDs
bio_values <- bio_rasters %>%
  as_tibble(xy = TRUE) %>%
  na.omit()

# Check predictor distributions:
bio_values %>% 
  select(starts_with("bio")) %>%
  mutate(decomp = exp(0.095 * bio1 + -0.00014 * bio1^2) * (1 - exp(-1.21 * bio12))) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free")
bio_values %>%
  select(starts_with("bio")) %>%
  mutate(decomp = exp(0.095 * bio1 + -0.00014 * bio1^2) * (1 - exp(-1.21 * bio12))) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = log(value))) +
  geom_density() +
  facet_wrap(~name, scales = "free")

#### (1.4) Soil Grid #####

# Load the Soil Grids rasters
soil_rasters <- list.files(
  "data/aus_forests_23/",
  pattern = "\\_5cm_10.tif$",
  full.names = TRUE
) %>%
  rast()

# Assign variable names
names(soil_rasters) <- sub(
  "_5cm_10\\.tif$", "",
  basename(list.files(
    "data/aus_forests_23/",
    pattern = "_5cm_10\\.tif$",
    full.names = TRUE
  ))
)

# Convert to a tibble
soil_values <- soil_rasters %>%
  as_tibble(xy = TRUE) %>%
  na.omit() %>%
  mutate(
    CN_ratio = SOC / N_total,
    CP_ratio = SOC / P_total
  )

# Check normality of predictors %>%
soil_values %>% select(-x, -y) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~name, scales = "free")
soil_values %>% select(-x, -y) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = log10(value))) +
  geom_density() +
  facet_wrap(~name, scales = "free")

#### (1.5) Elevation Grid #####

# Load the elevation raster
elevation_raster <- rast(
  "data/aus_forests_23/elevation_10.tif"
)

# Convert to a tibble with cell IDs
elevation_values <- elevation_raster %>%
  as_tibble(xy = TRUE) %>%
  na.omit()

# Check predictor distributions:
elevation_values %>% 
  ggplot(aes(x = elevation)) +
  geom_density()
elevation_values %>%
  mutate(
    elevation = case_when(
      elevation < 1 ~ 1,
      TRUE ~ elevation
    )
  ) %>%
  ggplot(aes(x = log10(elevation))) +
  geom_density()

#### (1.6) Woody vegetation cover Grid #####

# Load the elevation raster
woody_veg_raster <- rast(
  "data/aus_forests_23/wood_veg_10.tif"
)

# Convert to a tibble with cell IDs
woody_veg_values <- woody_veg_raster %>%
  as_tibble(xy = TRUE) %>%
  na.omit()

# Check predictor distributions:
woody_veg_values %>% 
  ggplot(aes(x = wood_veg_cover)) +
  geom_density()

#### (1.7) Habitat condition ####

# Read the habitat condition raster
habitat_rast <- rast(
  "data/aus_forests_23/habitat_condition_10.tif"
)

# Convert to a tibble with cell IDs
habitat_values <- habitat_rast %>%
  as_tibble(xy = TRUE) %>%
  na.omit()

# Check distribution:
habitat_values %>% 
  ggplot(aes(x = habitat_condition)) +
  geom_density()

#### (1.7) Perform RCA ####

variables <- inner_join(native_forest_values, bio_values, by = c("x", "y")) %>%
  inner_join(soil_values, by = c("x", "y")) %>%
  inner_join(elevation_values, by = c("x", "y")) %>%
  select(
    bio1, bio4, bio12, bio15, bio28, bio31, SOC, N_total, P_available,
    pH_CaCl2
  ) %>%
  mutate(
    decomp = exp(0.095 * bio1 + -0.00014 * bio1^2) * (1 - exp(-1.21 * bio12)),
    bio12 = log10(bio12),
    N_total = log10(N_total),
    P_available = log10(P_available),
    SOC = log10(SOC)
  ) %>%
  # Rename variables
  rename(
    `MAT (°C)` = bio1, `MAP (mm, logarithm)` = bio12, MAMI = bio28, 
    `Temp. seasonality` = bio4, `Precip. seasonality` = bio15, 
    `MI seasonality` = bio31, MAD = decomp,
    `Soil OC (%, logarithm)` = SOC, `Soil pH` = pH_CaCl2, 
    `Soil total N (%, logarithm)` = N_total, `Soil available P (mg/kg, logarithm)` = P_available
  ) %>%
  na.omit()

# RCA
rca <- variables %>%
  mutate_all(scale) %>%
  psych::principal(
    .,
    nfactors = 3,
    rotate = "varimax"
  )

# Extract the scores
scores <- bind_cols(
  # Add x and y coordinates
  inner_join(native_forest_values, bio_values, by = c("x", "y")) %>%
    inner_join(soil_values, by = c("x", "y")) %>%
    inner_join(elevation_values, by = c("x", "y")) %>%
    inner_join(woody_veg_values, by = c("x", "y")) %>%
    select(x, y),
  # Add RCA scores
  rca$scores[,1:3]
)

# Create a vector
predictors_rca_vect <- scores %>%
  vect(
    geom = c("x", "y"),
    crs = crs(native_forest_rast_10)
  )

# Create a raster stack for all RCA predictors
rca_rasters <- list()

# Loop through predictor names and rasterize each one
for (predictor in names(predictors_rca_vect)) {
  raster <- rasterize(
    predictors_rca_vect,
    native_forest_rast_10,
    field = predictor,
    fun = "mean",
    background = NA
  )
  rca_rasters[[predictor]] <- raster
}

# Stack the rasters
rca_raster_stack <- rast(rca_rasters)

# Combine forest grid with PCA/RCA raster stack
combined_raster <- c(
  native_forest_rast_10 %>% rename(forest_cover = forest),
  climate_zone_rast,
  biomes_raster,
  ecoregions_raster,
  bioregions_raster,
  habitat_type_rast,
  rca_raster_stack,
  bio_rasters %>% mutate(
    decomp = exp(0.095 * bio1 + -0.00014 * bio1^2) * (1 - exp(-1.21 * bio12)),
  ),
  soil_rasters %>%
    mutate(
      CN_ratio = SOC / N_total,
      CP_ratio = SOC / P_total
    ),
  elevation_raster,
  woody_veg_raster,
  habitat_rast
)

# Save the raster for predictions
writeRaster(
  combined_raster,
  "data/aus_forests_23/predictors_10.tif",
  overwrite = TRUE
)

# Get longitudes and latitudes
projected_points <- combined_raster %>%
  as_tibble(xy = TRUE) %>%
  na.omit() %>%
  mutate(
    x_albers = x,
    y_albers = y
  ) %>%
  vect(geom = c("x", "y"), crs = crs(native_forest_rast_10)) %>%
  project(crs(aus_map$data))

# Extract transformed coordinates
coords <- geom(projected_points)

# Convert back to tibble and include coordinates
predictors_10 <- as_tibble(projected_points) %>%
  mutate(
    longitude = coords[,"x"],
    latitude = coords[,"y"]
  )
glimpse(predictors_10)

# Save a tibble for predictions
fwrite(
  predictors_10,
  "data/aus_forests_23/predictors_10.txt",
  sep = "\t"
)

#### (1.8) Plot RC correlations ####

# Define the factor levels
factor_levels <- c(
  "MAT (°C)", "MAP (mm, logarithm)", "MAMI",
  "Temp. seasonality", "Precip. seasonality", "MI seasonality",
  "MAD", "Soil OC (%, logarithm)", "Soil pH",
  "Soil total N (%, logarithm)", "Soil available P (mg/kg, logarithm)"
)

#### * RC1 ####

# Correlations between raw variables and RCA scores
RC1_correlations <- bind_cols(
  scores,
  variables
) %>%
  select(-x, -y, -RC2, -RC3) %>%
  pivot_longer(
    cols = -c(RC1),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = factor_levels
    )
  ) %>%
  ggplot(aes(x = RC1, y = value)) +
  geom_point(
    shape = 21,
    alpha = 0.2,
    size = 0.5
  ) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(
    aes(label = after_stat(rr.label)),
    color = "red",
    r.accuracy = 0.01
  ) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(aspect.ratio = 1)

# Display the plot
print(RC1_correlations)

# Save the plot
ggsave(
  "output/supplementary_predictors/RC1_correlations.png",
  RC1_correlations,
  width = 15.75,
  height = 22,
  bg = "white",
  units = "cm"
)

#### * RC2 ####

# Correlations between raw variables and RCA scores
RC2_correlations <- bind_cols(
  scores,
  variables
) %>%
  select(
    -x, -y, -RC1, -RC3
  ) %>%
  pivot_longer(
    cols = -c(RC2),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = factor_levels
    )
  ) %>%
  ggplot(aes(x = RC2, y = value)) +
  geom_point(
    shape = 21,
    alpha = 0.2,
    size = 0.5
  ) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(
    aes(label = after_stat(rr.label)),
    color = "red",
    r.accuracy = 0.01
  ) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(aspect.ratio = 1)

# Display the plot
print(RC2_correlations)

# Save the plot
ggsave(
  "output/supplementary_predictors/RC2_correlations.png",
  RC2_correlations,
  width = 15.75,
  height = 22,
  bg = "white",
  units = "cm"
)

#### * RC3 ####

# Correlations between raw variables and RCA scores
RC3_correlations <- bind_cols(
  scores,
  variables
) %>%
  select(
    -x, -y, -RC1, -RC2
  ) %>%
  pivot_longer(
    cols = -c(RC3),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = factor_levels
    )
  ) %>%
  ggplot(aes(x = RC3, y = value)) +
  geom_point(
    shape = 21,
    alpha = 0.2,
    size = 0.5
  ) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(
    aes(label = after_stat(rr.label)),
    color = "red",
    r.accuracy = 0.01
  ) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(aspect.ratio = 1)

# Display the plot
print(RC3_correlations)

# Save the plot
ggsave(
  "output/supplementary_predictors/RC3_correlations.png",
  RC3_correlations,
  width = 15.75,
  height = 22,
  bg = "white",
  units = "cm"
)

#### (1.9) Plot RC maps ####

#### * RC1 ####

RC1_plots <- patchwork::wrap_plots(
  # RC1
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = RC1),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "RC1 (44%)"),
  # decomp
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = decomp),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "MAD"),
  # bio1
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = bio1),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "MAT (°C)"),
  # bio15
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = bio15),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "Precip. seasonality"),
  # bio04
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = bio4),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "Temp. seasonality"),
  ncol = 2
)

# Display the plot
print(RC1_plots)

# Save the plots
ggsave(
  "output/supplementary_predictors/RC1_map.png",
  RC1_plots,
  width = 15.75,
  height = 22,
  units = "cm"
)

#### * RC2 #### 

RC2_plots <- patchwork::wrap_plots(
  # RC2
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = RC2),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "RC2 (37%)"),
  # bio28
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = bio28),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "MAMI"),
  # bio12
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = bio12),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent",
      trans = "log10"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "MAP (mm, log-scaled)"),
  # SOC
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = SOC),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent",
      trans = "log10"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "Soil total OC (%, log-scaled)"),
  # pH_CaCl2
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = pH_CaCl2),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "Soil pH"),

  ncol = 2
)

# Display the plot
print(RC2_plots)

# Save the plots
ggsave(
  "output/supplementary_predictors/RC2_maps.png",
  RC2_plots,
  width = 15.75,
  height = 22,
  units = "cm"
)

#### * RC3 ####

RC3_noRC_plots <- patchwork::wrap_plots(
  # RC3
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = RC3),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank(),
      plot.tag = element_markdown()
    ) +
    labs(x = NULL,y = NULL, tag = "(**a**)", fill = "RC3 (19%)"),
  # P_total
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = P_available),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent",
      trans = "log10"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "Soil avail. P (mg/kg), log-scaled)"),
  # Moisture index seasonality
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = bio31),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank(),
      plot.tag = element_markdown()
    ) +
    labs(x = NULL,y = NULL, tag = "(**b**)", fill = "MI seasonality"),
  # Soil total N
  aus_map +
    tidyterra::geom_spatraster(
      data = combined_raster,
      aes(fill = N_total),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent",
      trans = "log10"
    ) +
    theme(
      legend.text.align = 1,
      legend.position = c(0, 1),
      legend.justification = c(-0.08, 0.97),
      legend.background = element_blank(),
      legend.key.size = unit(8, "pt"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.text = element_blank()
    ) +
    labs(x = NULL,y = NULL, fill = "Soil total N (%), log-scaled)"),
  ncol = 2
)

# Display the plot
print(RC3_noRC_plots)

# Save the plots
ggsave(
  "output/supplementary_predictors/RC3_map.png",
  RC3_noRC_plots,
  width = 15.75,
  height = 16,
  units = "cm"
)

# (2) Tree data ############################################################

# Load the tree data
trees <- data.table::fread(
  "data/presence/trees.txt"
) %>%
  mutate(
    mycorrhizal_type = case_when(
      mycorrhizal_type == 'NM-AM' ~ "AM",
      TRUE ~ mycorrhizal_type
    )
  )

# Load the site data
sites <- fread("data/presence/sites.txt") %>%
  select(site, longitude, latitude, n_obs)

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

#### (2.1) Sites for relative richness ####

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
# 17,787 cells have forest cover > 10% and are paired with predictor variables
# Cells with at least 10 observations
nrow(sites_upscaled %>% filter(!is.na(RC1), forest_cover >= 10, n_obs >= 10))
# 11,811 cells have forest cover > 10% and at least 10 observations

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
  fwrite("data/presence/sites_relative_richness_est_10.txt", sep = "\t")

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
    fread("data/presence/sites_relative_richness_est_10.txt"),
    by = c("x_albers", "y_albers")
  ) %>%
  # Order the columns
  select(
    longitude, latitude, x_albers, y_albers, forest_cover,
    starts_with("richness"), climate_zone, biome, ecoregion,
    everything()
  ) %>%
  fwrite("data/presence/sites_relative_richness_pred_10.txt", sep = "\t")

# Save the small sites for niche and biographical modelling
sites_upscaled %>%
  # Filter out sites without forest cover and missing RCs values
  filter(!is.na(RC1), forest_cover >= 10) %>%
  # Remove sites in the estimation data
  anti_join(
      .,
      fread("data/presence/sites_relative_richness_est_10.txt") %>%
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
  fwrite("data/presence/sites_10_under_sampled.txt", sep = "\t")

#### (2.3) Trees data for niche models ####

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
  "data/presence/trees_10.txt",
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
