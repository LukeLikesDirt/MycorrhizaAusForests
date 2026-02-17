
# Load libraries
require(terra)
require(data.table)
require(tidyverse)
source("code/functions.R")
source("code/map_australia.R")

# (1) Map geographic distribution of forests ###################################

# Read in the prediction raster
native_forest_rast_10 <- rast("data/aus_forests_23/predictors_10.tif")
global(!is.na(native_forest_rast_10), "sum") 

# Create a raster for relative richness estimation cells
relative_richness_df <- fread("data/presence/sites_relative_richness_est_10.txt")

# Vectorise the relative richness data
relative_richness_vect <- relative_richness_df %>%
  vect(geom = c("x_albers", "y_albers"), crs = crs(native_forest_rast_10)) %>%
  mutate(cell = 1) %>%
  select(cell)

# Rasterise the relative richness vector
relative_richness_rast <- rasterize(
  relative_richness_vect,
  native_forest_rast_10,
  field = "cell",
  fun = "mean"
) %>%
  rename(cell = mean)

# Create a raster for additional sites (used in niche models)
additional_sites_df <- fread("data/presence/sites_10_under_sampled.txt")

# Vectorise the additional sites data
additional_sites_vect <- additional_sites_df %>%
  vect(geom = c("x_albers", "y_albers"), crs = crs(native_forest_rast_10)) %>%
  mutate(cell = 1) %>%
  select(cell)

# Rasterise the additional sites vector
additional_sites_rast <- rasterize(
  additional_sites_vect,
  native_forest_rast_10,
  field = "cell",
  fun = "mean"
) %>%
  rename(cell = mean)

# How many sites in total?
nrow(additional_sites_df)
nrow(relative_richness_df)
sum(nrow(additional_sites_df) + nrow(relative_richness_df))

# How many sites in the raster?
sum(relative_richness_df$n_sites) + sum(additional_sites_df$n_sites)

# Join the rasters for mapping
combined_rast <- relative_richness_rast
combined_rast[!is.na(relative_richness_rast)] <- 1  # Value 1 for green sites
combined_rast[!is.na(additional_sites_rast)] <- 2   # Value 2 for red sites

# Take a glimpse at the rasters
native_forest_rast_10
combined_rast

#### * Geographic distribution * ####
forest_plots <- patchwork::wrap_plots(
  aus_map_2 +
    tidyterra::geom_spatraster(
      data = native_forest_rast_10,
      aes(fill = forest_cover),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      na.value = "transparent",
      breaks = c(20, 40, 60, 80),
    ) +
    theme(
      legend.text.align = 1,
      legend.text = element_text(size = rel(0.7)),
      legend.title = element_text(size = rel(0.8)),
      legend.position = c(0, 1),
      legend.justification = c(0.025, 0.975),
      legend.background = element_blank(),
      legend.frame = element_rect(colour = "grey30", linewidth = 0.1),
      legend.key.size = unit(7, "pt"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      plot.tag = ggtext::element_markdown()
    ) +
    labs(x = NULL,y = NULL, fill = "Forest extent (%)", tag = "**a**"),
  # Our samples
  aus_map_2 +
    tidyterra::geom_spatraster(
      data = combined_rast,
      aes(fill = cell),
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    scale_fill_gradient(
      low = "blue", high = "red",
      na.value = "transparent"
    ) +
    theme(
      axis.text.y = element_blank(),
      plot.tag = ggtext::element_markdown()
    ) +
    labs(x = NULL, y = NULL, tag = "**b**"),
  ncol = 2
)

# Save plots
ggsave(
  "output/figure_S1.png", forest_plots, width = 15, height = 8, units = "cm"
  )
