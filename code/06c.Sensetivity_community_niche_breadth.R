
# Modelling community-level environemtal breadth of Australian trees

# Predict community-level environmental breadth of Australian trees under
# the assumption that niche breadth is the product of competition (proxied 
# by species richness) and climate stability (proxied by bio4 and bio31).
# bio4 is the annual temperature variability and bio31 is the annual soil 
# moisture variability.

# Packages
require(INLA) # inla.list.models()
require(fmesher)
require(ggtext)
require(terra)
require(patchwork)
require(tidyverse)
source("code/map_australia.R")
set.seed(1986)

# Helper functions and constants ###############################################

# Common theme
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 9
common_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    axis.ticks = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
    plot.tag = element_markdown(size = tag_size),
    strip.text = element_text(face = "bold", size = strip_size),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    aspect.ratio = 1
  )

# Helper function create covariates for estimation and prediction
create_covariates <- function(data_est, data_pred, n_marg_pred = 100) {
  
  n_est <- nrow(data_est)
  n_pred <- nrow(data_pred)
  
  bind_rows(
    # Data for estimation and prediciton
    tibble(
      richness = c(data_est[["richness"]], data_pred[["richness"]]),
      bio4 = c(data_est[["bio4"]], data_pred[["bio4"]]),
      bio31 = c(data_est[["bio31"]], data_pred[["bio31"]]),
      tag = c(rep("est", n_est), rep("pred", n_pred))
    ),
    # Marginal effects of richness
    tibble(
      richness = seq(min(data_est[["richness"]]), max(data_est[["richness"]]), length.out = n_marg_pred),
      bio4 = median(data_est[["bio4"]]),
      bio31 = median(data_est[["bio31"]]),
      tag = "pred_richness"
    ),
    # Marginal effects of bio4
    tibble(
      richness = median(data_est[["richness"]]),
      bio4 = seq(min(data_est[["bio4"]]), max(data_est[["bio4"]]), length.out = n_marg_pred),
      bio31 = median(data_est[["bio31"]]),
      tag = "pred_bio4"
    ),
    # Marginal effects of bio31
    tibble(
      richness = median(data_est[["richness"]]),
      bio4 = median(data_est[["bio4"]]),
      bio31 = seq(min(data_est[["bio31"]]), max(data_est[["bio31"]]), length.out = n_marg_pred),
      tag = "pred_bio31"
    )
  )
}

# Helper functions to create estimateion stack
create_estimation_stack <- function(tag, data, covariates, A_matrix, spatial_index) {
  inla.stack(
    tag = tag,
    data = list(env_breadth = data[["env_breadth"]]),
    A = list(1, 1, A_matrix),
    effects = list(
      intercept = rep(1, nrow(data)),
      X = covariates %>% filter(tag == "est") %>% select(-tag),
      spatial = spatial_index
    )
  )
}

# Helper function for spatial prediction stack
create_prediction_stack <- function(tag, covariates, A_matrix, spatial_index) {
  
  covars <- covariates %>% filter(tag == !!tag) %>% select(-tag)
  n_pred <- nrow(covars)
  
  inla.stack(
    tag = tag,
    data = list(env_breadth = rep(NA_real_, n_pred)),
    A = list(1, 1, A_matrix),
    effects = list(
      intercept = rep(1, n_pred),
      X = covars,
      spatial = spatial_index
    )
  )
}

# Helper functions to marginal effects stack
create_marginal_stack <- function(tag, covariates) {
  
  covars <- covariates %>% filter(tag == !!tag) %>% select(-tag)
  n_pred <- nrow(covars)
  
  inla.stack(
    tag = tag,
    data = list(env_breadth = rep(NA_real_, n_pred)),
    A = list(1, 1),
    effects = list(
      intercept = rep(1, n_pred),
      X = covars
    )
  )
}

# Helper function to extract and back-transformation marginal effects
extract_marginal_effects <- function(
    model, covariates, 
    index_richness, index_bio4, index_bio31,
    original_data, model_name
) {
  
  # Get scaling parameters from original data
  env_breadth_center <- attr(original_data$env_breadth, "scaled:center")
  env_breadth_scale <- attr(original_data$env_breadth, "scaled:scale")
  richness_center <- attr(original_data$richness, "scaled:center")
  richness_scale <- attr(original_data$richness, "scaled:scale")
  bio4_center <- attr(original_data$bio4, "scaled:center")
  bio4_scale <- attr(original_data$bio4, "scaled:scale")
  bio31_center <- attr(original_data$bio31, "scaled:center")
  bio31_scale <- attr(original_data$bio31, "scaled:scale")
  
  # Extract predictions and back-transform
  richness_effects <- tibble(
    response_scaled = model$summary.linear.predictor$mean[index_richness],
    predictor_scaled = covariates %>% filter(tag == "pred_richness") %>% pull(richness),
    upper_scaled = model$summary.linear.predictor$`0.975quant`[index_richness],
    lower_scaled = model$summary.linear.predictor$`0.025quant`[index_richness],
    variable = "richness",
    model = model_name
  ) %>%
    mutate(
      response = response_scaled * env_breadth_scale + env_breadth_center,
      predictor = predictor_scaled * richness_scale + richness_center,
      upper = upper_scaled * env_breadth_scale + env_breadth_center,
      lower = lower_scaled * env_breadth_scale + env_breadth_center,
      predictor_std = predictor_scaled,
      upper_std = upper_scaled,
      lower_std = lower_scaled
    )
  
  bio4_effects <- tibble(
    response_scaled = model$summary.linear.predictor$mean[index_bio4],
    predictor_scaled = covariates %>% filter(tag == "pred_bio4") %>% pull(bio4),
    upper_scaled = model$summary.linear.predictor$`0.975quant`[index_bio4],
    lower_scaled = model$summary.linear.predictor$`0.025quant`[index_bio4],
    variable = "bio4",
    model = model_name
  ) %>%
    mutate(
      response = response_scaled * env_breadth_scale + env_breadth_center,
      predictor = predictor_scaled * bio4_scale + bio4_center,
      upper = upper_scaled * env_breadth_scale + env_breadth_center,
      lower = lower_scaled * env_breadth_scale + env_breadth_center,
      predictor_std = predictor_scaled,
      upper_std = upper_scaled,
      lower_std = lower_scaled
    )
  
  bio31_effects <- tibble(
    response_scaled = model$summary.linear.predictor$mean[index_bio31],
    predictor_scaled = covariates %>% filter(tag == "pred_bio31") %>% pull(bio31),
    upper_scaled = model$summary.linear.predictor$`0.975quant`[index_bio31],
    lower_scaled = model$summary.linear.predictor$`0.025quant`[index_bio31],
    variable = "bio31",
    model = model_name
  ) %>%
    mutate(
      response = response_scaled * env_breadth_scale + env_breadth_center,
      predictor = predictor_scaled * bio31_scale + bio31_center,
      upper = upper_scaled * env_breadth_scale + env_breadth_center,
      lower = lower_scaled * env_breadth_scale + env_breadth_center,
      predictor_std = predictor_scaled,
      upper_std = upper_scaled,
      lower_std = lower_scaled
    )
  
  bind_rows(richness_effects, bio4_effects, bio31_effects)
}

# Helper function to extract model coefficients for saving
extract_model_coefficients <- function(model, myc_type) {
  bind_rows(
    model[["summary.fixed"]] %>%
      rownames_to_column("parameter") %>%
      as_tibble() %>%
      select(
        parameter, median = `0.5quant`, lower = `0.025quant`, upper = `0.975quant`
      ),
    model[["summary.hyperpar"]] %>%
      rownames_to_column("parameter") %>%
      as_tibble() %>%
      select(
        parameter, median = `0.5quant`, lower = `0.025quant`, upper = `0.975quant`
      ) 
  ) %>%
    mutate(
      parameter = recode(parameter, 
                         "intercept" = "Intercept",
                         "richness" = "Tree richness",
                         "bio4" = "Temperature variability (Bio4)",
                         "bio31" = "Soil moisture variability (Bio31)",
                         "Precision for the Gaussian observations" = "Observation Precision",
                         "Range for spatial" = "Spatial Range",
                         "Stdev for spatial" = "Spatial SD"
      )
    )
}

# Helper function to extract fitted values
extract_fitted_values <- function(mod, idx, data_est_type, data_pred_type, type, mesh, forest_raster) {
  
  # Extract the fitted values
  fitted_values <- tibble(
    mycorrhizal_type = type,
    predicted = c(
      mod$summary.fitted[["mean"]][idx[["est"]]],
      mod$summary.fitted[["mean"]][idx[["pred"]]]
    ),
    upper = c(
      mod$summary.fitted[["0.975quant"]][idx[["est"]]],
      mod$summary.fitted[["0.975quant"]][idx[["pred"]]]
    ),
    lower = c(
      mod$summary.fitted[["0.025quant"]][idx[["est"]]],
      mod$summary.fitted[["0.025quant"]][idx[["pred"]]]
    ),
    sd = c(
      mod$summary.fitted[["sd"]][idx[["est"]]],
      mod$summary.fitted[["sd"]][idx[["pred"]]]
    ),
    observed = c(
      data_est_type[["env_breadth_orig"]],
      rep(NA_real_, length(idx[["pred"]]))
    ),
    latitude  = c(
      data_est_type[["latitude"]],
      data_pred_type[["latitude"]]
    ),
    longitude = c(
      data_est_type[["longitude"]],
      data_pred_type[["longitude"]]
    ),
    x_albers = c(
      data_est_type[["x_albers"]],
      data_pred_type[["x_albers"]]
    ),
    y_albers = c(
      data_est_type[["y_albers"]],
      data_pred_type[["y_albers"]]
    ),
    richness = c(
      data_est_type[["richness_orig"]],
      data_pred_type[["richness_orig"]]
    ),
    bio4 = c(
      data_est_type[["bio4_orig"]],
      data_pred_type[["bio4_orig"]]
    ),
    bio31 = c(
      data_est_type[["bio31_orig"]],
      data_pred_type[["bio31_orig"]]
    ),
    source = c(
      rep("est", length(idx[["est"]])),
      rep("pred", length(idx[["pred"]]))
    )) %>%
    mutate(
      x_albers          = as.numeric(x_albers),
      y_albers          = as.numeric(y_albers)
    )
  
  # Add spatial field
  A_posterior <- fmesher::fm_basis(
    mesh,
    loc = matrix(
      c(fitted_values$x_albers / 1000,
        fitted_values$y_albers / 1000),
      ncol = 2
    )
  )
  
  fitted_values <- fitted_values %>%
    mutate(
      spatial_field = as.vector(
        A_posterior %*% mod$summary.random$spatial[, "mean"]
      )
    )
  
  # Back-standardise the predicted, upper and lower bounds
  env_breadth_center <- attr(data_est_type$env_breadth, "scaled:center")
  env_breadth_scale <- attr(data_est_type$env_breadth, "scaled:scale")
  
  fitted_values <- fitted_values %>%
    mutate(
      predicted = predicted * env_breadth_scale + env_breadth_center,
      upper = upper * env_breadth_scale + env_breadth_center,
      lower = lower * env_breadth_scale + env_breadth_center,
      sd = sd * env_breadth_scale,
      uncertainty = upper - lower,
      residuals = observed - predicted,
      pearsons_residuals = residuals / sd
    )
  
  # Rasterise fitted values
  fitted_rast <- rasterize(
    vect(
      fitted_values,
      geom = c("x_albers", "y_albers"),
      crs = crs(forest_raster)
    ),
    forest_raster,
    field = names(fitted_values %>% select(-latitude, -longitude, -x_albers, -y_albers)),
    background = NA
  ) %>%
    project(crs(aus_map$data))
  
  # Standardised to percentiles
  standardised_rast <- standardise_raster_percentiles(fitted_rast[["predicted"]])
  
  return(list(
    fitted_values = fitted_values,
    fitted_raster = fitted_rast,
    fitted_raster_standardised = standardised_rast
  ))
}

# Function to standardise richness values scale for each map using 1% brackets
standardise_raster_percentiles <- function(raster_data) {
  # Create a copy of the raster
  result <- raster_data
  
  # Get the original values and remove NAs for quantile calculation
  values <- terra::values(raster_data)
  non_na_values <- values[!is.na(values)]
  
  # Calculate quantile breaks for 1% intervals (100 brackets)
  breaks <- quantile(non_na_values, probs = seq(0, 1, by = 0.01))
  
  # Create a new vector to hold the quantized values
  new_values <- values
  
  # Assign percentile categories (1-100)
  for (i in 1:100) {
    if (i == 1) {
      # First percentile includes the minimum value
      mask <- values <= breaks[2] & !is.na(values)
    } else if (i == 100) {
      # Last percentile includes the maximum value
      mask <- values > breaks[100] & !is.na(values)
    } else {
      # Middle percentiles
      mask <- values > breaks[i] & values <= breaks[i+1] & !is.na(values)
    }
    
    # Assign the percentile number to matching cells
    new_values[mask] <- i
  }
  
  # Update the raster with the new values
  terra::values(result) <- new_values
  return(result)
}

# Supplementary spatial plots to evaluate model performance
create_spatial_plot <- function(raster, data_column, color_scale = "Spectral", 
                                y_label = NULL, title = NULL, viridis_option = NULL, 
                                transformation = "identity") {
  
  p <- aus_map +
    tidyterra::geom_spatraster(
      data = raster,
      aes(fill = .data[[data_column]]),
      na.rm = TRUE
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = if (is.null(title) || title == "") element_blank() else
        element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = y_label, fill = NULL, title = title)
  
  if (!is.null(viridis_option)) {
    # Viridis scale
    p <- p + scale_fill_viridis_c(
      option = viridis_option,
      na.value = "transparent",
      trans = transformation
    )
  } else if (!is.null(color_scale)) {
    # Only apply brewer scale if a palette name is given
    p <- p + scale_fill_distiller(
      palette = color_scale,
      direction = -1,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3),
      trans = transformation
    )
  }
  
  return(p)
}

# Get the spatial parameters of an INLA model
spatial_params <- function(model, spde, name) {
  spatial_field <- inla.spde2.result(
    inla = model, name = name, spde = spde, do.transfer = TRUE
  )
  
  compute_marginal <- function(marginal, transform = identity) {
    inla.emarginal(transform, marginal[[1]])
  }
  
  result <- c(
    kappa = compute_marginal(spatial_field$marginals.kappa),
    sigma_spatial = compute_marginal(spatial_field$marginals.variance.nominal, sqrt),
    range = compute_marginal(spatial_field$marginals.range.nominal)
  )
  
  return(result)
}

# Helper function to plot the Matern covariance function
plot_matern_correlation <- function(
    spatial_stats, 
    x, y, 
    range = NULL, 
    x_limits = c(0, 1000), 
    x_breaks = seq(0, 1000, by = 100),
    y_limits = c(0, 1),
    y_breaks = seq(0, 1, by = 0.1)
) {
  dist_vect <- seq(0, max(as.matrix(dist(tibble(x, y)))), length.out = 100)
  mat_cor <- (spatial_stats["kappa"] * dist_vect) * besselK(spatial_stats["kappa"] * dist_vect, 1)
  
  plot_data <- tibble(dist_vect = dist_vect, mat_cor = mat_cor)
  
  if (!is.null(range)) {
    p <- ggplot(plot_data, aes(dist_vect, mat_cor)) +
      geom_hline(yintercept = 0.15, linetype = "dotted", color = "red") +
      geom_vline(xintercept = range, linetype = "dotted", color = "red") +
      geom_line() +
      scale_x_continuous(limits = x_limits, breaks = x_breaks) +
      scale_y_continuous(limits = y_limits, breaks = y_breaks) +
      theme_minimal() +
      theme(
        panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
        aspect.ratio = 1
      ) +
      labs(x = "Distance (km)", y = "Matérn correlation")
  } else {
    p <- ggplot(plot_data, aes(dist_vect, mat_cor)) +
      geom_hline(yintercept = 0.15, linetype = "dotted", color = "red") +
      geom_line() +
      scale_x_continuous(limits = x_limits, breaks = x_breaks) +
      scale_y_continuous(limits = y_limits, breaks = y_breaks) +
      theme_minimal() +
      theme(
        panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
        aspect.ratio = 1
      ) +
      labs(x = "Distance (km)", y = "Matérn correlation")
  }
  
  return(p)
}

# Function to compute Moran's I for a given model type
compute_moran_i <- function(model_type) {
  # Get fitted values for estimation data only
  fitted_data <- model_info[[model_type]]$fitted_values %>%
    filter(source == "est")
  
  # Create coordinates matrix
  coords <- as.matrix(fitted_data[, c("x_albers", "y_albers")])
  
  # Create spatial weights using delaunay triangulation approach
  # Use a distance threshold based on 10th percentile of distances
  distances <- as.matrix(dist(coords))
  threshold <- quantile(distances[distances > 0], 0.1)
  
  # Create neighborhood structure
  nb <- spdep::dnearneigh(coords, d1 = 0, d2 = threshold)
  
  # Convert to listw object
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Extract residuals
  residuals <- fitted_data$residuals
  
  # Compute Moran's I (handle zero policy for isolated observations)
  moran_test <- spdep::moran.test(residuals, lw, zero.policy = TRUE)
  
  return(list(
    Moran_I = moran_test$estimate[1],
    p_value = moran_test$p.value,
    model_type = model_type
  ))
}

# Function to standardise niche breadth values 100 percentiles
quantise_tibble_percentiles <- function(df, columns_to_quantise) {
  result <- df
  
  for (col in columns_to_quantise) {
    # Get the original values for quantile calculation
    values <- df[[col]]
    non_na_values <- values[!is.na(values)]
    
    # Calculate quantile breaks for 100 percentiles
    breaks <- quantile(non_na_values, probs = seq(0, 1, by = 0.01))
    
    # Create a new vector to hold the quantized values
    new_values <- values
    
    # Assign percentile categories (1-100)
    for (i in 1:100) {
      if (i == 1) {
        # First percentile includes the minimum value
        mask <- values <= breaks[2] & !is.na(values)
      } else if (i == 100) {
        # Last percentile includes the maximum value
        mask <- values > breaks[100] & !is.na(values)
      } else {
        # Middle percentiles
        mask <- values > breaks[i] & values <= breaks[i+1] & !is.na(values)
      }
      
      # Assign the percentile number to matching cells
      new_values[mask] <- i
    }
    
    # Update the dataframe with the new values
    result[[col]] <- new_values
  }
  
  return(result)
}

# (1) Organise data ############################################################

# Note: Georeferenced predictors were organised and in an earlier step. Richness
# for all trees was estimated in a subsequent step by correcting for sample
# effort. I therefore compile the data here to prepare for the modelling. 
# I limit modelling to cells that have at least 10 observations for each 
# mycorrhizal type to ensure community-level estimates are not skewed by a 
# random single observation of an individual. Therefore, a single species may
# be considered only when supported multiple unique observations.

#### * Niche and richness data * ####

# Read in the niche estimate data
species_niche_data <- data.table::fread(
  "data/niche_estimates_enmeval/niche_estimates.txt"
) %>%
  filter(mycorrhizal_type != "ErM") %>%
  mutate(
    # Sqrt-root transformation to stabilise env breadth
    env_breadth = sqrt(env_B2_corrected)
  )

# Sample effort-corrected richness for prediction
richness_data <- data.table::fread(
  "output/supplimentary_absolute_richness/effort_adjusted_richness.txt"
)

#### * Number of observations * ####

# Count the number of observations for each mycorrhizal type in each cell:
# Read in the forest raster
forest_raster <- terra::rast("data/aus_forests_23/aus_for23_masked_10.tif")

# Trees presence data
presence_vect <- data.table::fread(
  "data/presence/trees.txt",
  stringsAsFactors = TRUE
) %>%
  # Remove ErM mycorrhizal type
  filter(mycorrhizal_type != "ErM") %>%
  # Mutate NM-AM to AM and EcM-AM to Dual
  mutate(
    mycorrhizal_type = recode(
      mycorrhizal_type,
      "NM-AM" = "AM"
    ),
    # Re-level the factors
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "EcM", "EcM-AM", "NM")
    )
  ) %>%
  # Convert to vector
  terra::vect(
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  ) %>%
  # Project to match forest raster CRS
  terra::project(terra::crs(forest_raster))

# Count observations per mycorrhizal type using your forest grid
count_rasters <- list()

for(type in unique(presence_vect$mycorrhizal_type)) {
  # Filter for this mycorrhizal type
  type_data <- presence_vect[presence_vect$mycorrhizal_type == type]
  
  # Rasterize using your forest raster as template
  count_rasters[[type]] <- terra::rasterize(
    type_data, 
    forest_raster,  # Use your existing grid
    fun = "length",  # Count number of points per cell
    background = 0
  )
  names(count_rasters[[type]]) <- type
}

# Combine into single raster stack
n_obs_stack <- terra::rast(count_rasters)

# Mask to only include cells that have forest data (non-NA in forest_raster)
n_obs_stack <- terra::mask(n_obs_stack, forest_raster)

# Vectorise, keeping cell names and coordinates
n_obs_per_cell <- as.data.frame(n_obs_stack, xy = TRUE, cell = TRUE) %>%
  # Add "cell_" prefix
  mutate(cell = paste0("cell_", cell)) %>%
  pivot_longer(
    cols = -c(x, y, cell),
    names_to = "mycorrhizal_type",
    values_to = "n_obs"
  ) %>%
  # Remove cells with < 10 observations
  filter(n_obs >= 10)

#### * Estimation data * ####

# Summarise the niche data for estimation
data_est <- data.table::fread(
  "data/presence/trees_10.txt",
  stringsAsFactors = TRUE
) %>%
  # Remove cells for mycorrhizal types with < 10 observations
  inner_join(
    n_obs_per_cell %>%
      select(mycorrhizal_type, cell),
    by = c("mycorrhizal_type", "cell")
  ) %>%
  select(
    cell, longitude, latitude, 
    species = scientific_name, mycorrhizal_type
  ) %>%
  inner_join(
    species_niche_data %>%
      select(species, env_breadth),
    by = c("species")
  ) %>%
  group_by(mycorrhizal_type, cell, longitude, latitude) %>%
  summarise(
    env_breadth = mean(env_breadth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Recode dual-mycorrhizal type
  mutate(
    mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
    # Level the mycorrhizal types
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "EcM", "Dual", "NM")
    ),
  ) %>%
  # Remove environmental breadth outliers
  filter(
    env_breadth < 0.7
  )

#### * Prediction data * ####

# Combine estimation and prediction data from across Austrlia and standardise
# the predictors
data_pred <- bind_rows(
  data.table::fread("data/presence/sites_relative_richness_est_10.txt") %>%
    select(longitude, latitude, x_albers, y_albers, bio4, bio31),
  data.table::fread("data/presence/sites_relative_richness_pred_10.txt") %>%
    select(longitude, latitude, x_albers, y_albers, bio4, bio31)
) %>%
  inner_join(
    richness_data %>%
      select(richness, longitude, latitude),
    by = c("longitude", "latitude")
  ) %>%
  # Scale the predictors and convetrt albers to km
  mutate(
    bio4_orig = bio4,
    bio31_orig = bio31,
    richness_orig = sqrt(richness),
    bio4 = scale(bio4),
    bio31 = scale(bio31),
    richness = scale(sqrt(richness)),
    x_km = x_albers / 1000,
    y_km = y_albers / 1000
  )

#### * Split data * ####

# Split the estimation data by mycorrhizal type
data_est_list <- data_est %>%
  # Add the prodictors
  inner_join(
    data_pred,
    by = c("longitude", "latitude")
  ) %>%
  split(
    data_est %>% 
      inner_join(data_pred, by = c("longitude", "latitude")) %>% 
      pull(mycorrhizal_type)
  ) %>%
  # Scale the response across all mycorrhizal types
  lapply(function(df) {
    df %>%
      mutate(
        env_breadth_orig = env_breadth,
        env_breadth = scale(env_breadth)
      )
  })

# Split the prediction data by mycorrhizal type
data_pred_list <- list(
  "AM" = data_pred %>% anti_join(
    data_est_list[["AM"]] %>% select(longitude, latitude),
    by = c("longitude", "latitude")
  ),
  "EcM" = data_pred %>% anti_join(
    data_est_list[["EcM"]] %>% select(longitude, latitude),
    by = c("longitude", "latitude")
  ),
  "Dual" = data_pred %>% anti_join(
    data_est_list[["Dual"]] %>% select(longitude, latitude),
    by = c("longitude", "latitude")
  ),
  "NM" = data_pred %>% anti_join(
    data_est_list[["NM"]] %>% select(longitude, latitude),
    by = c("longitude", "latitude")
  )
)

# Niche Breadth Sensitivity Analysis
# Add this after the data organization section and before model fitting

# (2) Create sensitivity datasets for niche breadth ###########################

# Function to create sensitivity datasets by resampling mycorrhizal types
create_niche_sensitivity_datasets <- function(species_niche_data, tree_data, n_sensitivity = 5, resample_prop = 0.1) {
  
  # Calculate mycorrhizal type probabilities from unique species in tree data
  species_types <- tree_data %>%
    distinct(scientific_name, mycorrhizal_type) %>%
    count(mycorrhizal_type)
  myco_counts <- setNames(species_types$n, species_types$mycorrhizal_type)
  myco_probs <- myco_counts / sum(myco_counts)
  
  cat("Mycorrhizal type probabilities calculated from tree data (per species):\n")
  print(round(myco_probs, 3))
  cat("\n")
  
  # Get species per type from tree data
  species_per_type <- tree_data %>%
    group_by(mycorrhizal_type) %>%
    distinct(scientific_name) %>%
    group_by(mycorrhizal_type) %>%
    nest() %>%
    mutate(species = map(data, ~pull(., scientific_name))) %>%
    select(-data) %>%
    deframe()
  
  sensitivity_datasets <- list()
  
  for (i in 1:n_sensitivity) {
    cat("Creating niche breadth sensitivity dataset", i, "...\n")
    
    # Create a copy of the original species niche data
    sens_niche_data <- species_niche_data
    
    # Get all species that appear in both datasets
    common_species <- intersect(unique(sens_niche_data$species), unique(tree_data$scientific_name))
    
    # Determine number of species to resample
    n_resample <- ceiling(length(common_species) * resample_prop)
    
    # Get current mycorrhizal types for common species from tree data
    current_types <- tree_data %>%
      filter(scientific_name %in% common_species) %>%
      distinct(scientific_name, mycorrhizal_type) %>%
      deframe()
    
    # Stratified selection by mycorrhizal type
    expected <- n_resample * myco_probs
    n_resample_per_group <- floor(expected)
    remaining <- n_resample - sum(n_resample_per_group)
    if (remaining > 0) {
      frac <- expected - n_resample_per_group
      add_idx <- order(-frac)[1:remaining]
      n_resample_per_group[add_idx] <- n_resample_per_group[add_idx] + 1
    }
    
    resample_species <- c()
    group_names <- names(myco_probs)
    for (j in seq_along(species_per_type)) {
      group_name <- group_names[j]
      n_res_g <- n_resample_per_group[group_name]
      if (n_res_g > 0 && group_name %in% names(species_per_type)) {
        # Only sample from species that have niche data
        available_species <- intersect(species_per_type[[group_name]], common_species)
        if (length(available_species) > 0) {
          sampled <- sample(available_species, min(n_res_g, length(available_species)))
          resample_species <- c(resample_species, sampled)
        }
      }
    }
    
    # Resample mycorrhizal types for selected species
    new_myco_types <- sample(
      names(myco_probs),
      size = length(resample_species),
      replace = TRUE,
      prob = myco_probs
    )
    
    # Create a mapping for new mycorrhizal types
    species_remap <- setNames(new_myco_types, resample_species)
    
    # Store the sensitivity dataset with species remapping
    sensitivity_datasets[[paste0("sensitivity_", i)]] <- list(
      species_niche_data = sens_niche_data,
      species_remap = species_remap
    )
  }
  
  return(sensitivity_datasets)
}

# Read in the tree data for sensitivity analysis
tree_data_for_sensitivity <- data.table::fread("data/presence/trees_10.txt") %>%
  filter(
    mycorrhizal_type %in% c("AM", "EcM", "EcM-AM", "NM"),
  ) %>%
  # Use all cells that have prediction data
  inner_join(
    data_pred %>% select(longitude, latitude),
    by = c("longitude", "latitude")
  )

# Create sensitivity datasets with 10% resampling
sensitivity_niche_datasets <- create_niche_sensitivity_datasets(
  species_niche_data, tree_data_for_sensitivity, n_sensitivity = 5, resample_prop = 0.1
)

# Process sensitivity datasets to create estimation data
sensitivity_est_niche_list <- list()

for (i in 1:5) {
  cat("Processing niche breadth sensitivity dataset", i, "...\n")
  
  # Get the sensitivity dataset
  sens_data <- sensitivity_niche_datasets[[paste0("sensitivity_", i)]]
  sens_niche_data <- sens_data$species_niche_data
  species_remap <- sens_data$species_remap
  
  # Apply species remapping to tree data
  sens_tree_data <- tree_data_for_sensitivity %>%
    mutate(
      mycorrhizal_type = ifelse(
        scientific_name %in% names(species_remap),
        species_remap[scientific_name],
        as.character(mycorrhizal_type)
      ),
      mycorrhizal_type = factor(
        mycorrhizal_type,
        levels = c("AM", "EcM", "EcM-AM", "NM")
      )
    )
  
  # Create community-level niche breadth estimates using remapped tree data
  sens_processed <- sens_tree_data %>%
    # Only keep cells with minimum observations per type
    inner_join(
      n_obs_per_cell %>% select(mycorrhizal_type, cell),
      by = c("mycorrhizal_type", "cell")
    ) %>%
    select(
      cell, longitude, latitude,
      species = scientific_name, mycorrhizal_type
    ) %>%
    inner_join(
      sens_niche_data %>% select(species, env_breadth),
      by = "species"
    ) %>%
    group_by(mycorrhizal_type, cell, longitude, latitude) %>%
    summarise(
      env_breadth = mean(env_breadth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Recode dual-mycorrhizal type and filter outliers
    mutate(
      mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
      mycorrhizal_type = factor(
        mycorrhizal_type,
        levels = c("AM", "EcM", "Dual", "NM")
      )
    ) %>%
    filter(env_breadth < 0.7) %>%
    # Add predictors and split by mycorrhizal type
    inner_join(data_pred, by = c("longitude", "latitude")) %>%
    split(.$mycorrhizal_type) %>%
    # Scale the response using the same scaling as original data
    lapply(function(df) {
      # Find corresponding original data to get scaling parameters
      original_type <- unique(df$mycorrhizal_type)
      if (original_type %in% names(data_est_list)) {
        original_center <- attr(data_est_list[[original_type]]$env_breadth, "scaled:center")
        original_scale <- attr(data_est_list[[original_type]]$env_breadth, "scaled:scale")
        
        df %>%
          mutate(
            env_breadth_orig = env_breadth,
            env_breadth = scale(env_breadth, center = original_center, scale = original_scale)
          )
      } else {
        df %>%
          mutate(
            env_breadth_orig = env_breadth,
            env_breadth = scale(env_breadth)
          )
      }
    })
  
  sensitivity_est_niche_list[[i]] <- sens_processed
}

# Verify data alignment for first sensitivity dataset
if (length(sensitivity_est_niche_list) > 0 && "AM" %in% names(sensitivity_est_niche_list[[1]])) {
  alignment_check <- plot(
    data_est_list[["AM"]] %>% arrange(cell) %>% pull(bio4),
    sensitivity_est_niche_list[[1]][["AM"]] %>% arrange(cell) %>% pull(bio4),
    main = "Bio4 alignment check",
    xlab = "Original data Bio4",
    ylab = "Sensitivity 1 Bio4"
  )
}

# (3) Sensitivity analysis model fitting #######################################

# Fit models for each sensitivity dataset
sensitivity_model_info <- list()

for (sens_i in 1:5) {
  cat("\n=== Processing Sensitivity Dataset", sens_i, "===\n")
  
  current_est_list <- sensitivity_est_niche_list[[sens_i]]
  
  sensitivity_model_info[[sens_i]] <- list()
  
  for (type in c("AM", "EcM", "Dual", "NM")) {
    if (!type %in% names(current_est_list)) {
      cat("Skipping", type, "- not available in sensitivity dataset", sens_i, "\n")
      next
    }
    
    cat("Processing sensitivity", sens_i, "- type:", type, "\n")
    
    tryCatch({
      # Use the same data structure as original analysis
      data_est_type <- current_est_list[[type]]
      data_pred_type <- data_pred_list[[type]]  # Use original prediction data
      
      cat("Sample size for", type, ":", nrow(data_est_type), "\n")
      
      if (nrow(data_est_type) < 10) {
        cat("Insufficient data for", type, "in sensitivity", sens_i, "\n")
        next
      }
      
      # Build mesh
      loc_est <- as.matrix(data_est_type %>% select(x_km, y_km))
      loc_pred <- as.matrix(data_pred_type %>% select(x_km, y_km))
      
      mesh <- fm_mesh_2d_inla(
        loc = loc_est,
        max.edge = c(1, 5) * max_edge,
        cutoff = max_edge / 5
      )
      fm_crs(mesh) <- crs(aus_map_albers$data)
      
      spde <- inla.spde2.pcmatern(
        mesh = mesh,
        prior.range = c(range_est, 0.01),
        prior.sigma = c(sigma_est, 0.01)
      )
      index_spatial <- inla.spde.make.index("spatial", n.spde = spde$n.spde)
      
      # Projector matrices
      A_est <- inla.spde.make.A(mesh, loc = loc_est)
      A_pred <- inla.spde.make.A(mesh, loc = loc_pred)
      
      # Covariates
      covariates <- create_covariates(data_est_type, data_pred_type)
      
      # Stacks
      stack_est <- create_estimation_stack(
        "est", data_est_type, covariates, A_est, index_spatial
      )
      stack_pred <- create_prediction_stack(
        "pred", covariates, A_pred, index_spatial
      )
      stack_pred_richness <- create_marginal_stack("pred_richness", covariates)
      stack_pred_bio4 <- create_marginal_stack("pred_bio4", covariates)
      stack_pred_bio31 <- create_marginal_stack("pred_bio31", covariates)
      
      stack_type <- inla.stack.join(
        stack_est, stack_pred, stack_pred_richness, stack_pred_bio4, stack_pred_bio31
      )
      
      pred_indexes <- list(
        richness = inla.stack.index(stack_type, "pred_richness")$data,
        bio4 = inla.stack.index(stack_type, "pred_bio4")$data,
        bio31 = inla.stack.index(stack_type, "pred_bio31")$data
      )
      
      # Fit model
      form <- env_breadth ~ -1 + intercept + richness + bio4 + bio31 + f(spatial, model = spde)
      
      if (type == "AM") {
        mod <- inla(
          form,
          family = "gaussian",
          data = inla.stack.data(stack_type),
          control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stack_type)),
          control.fixed = list(
            mean = list(bio31 = 0),
            prec = list(bio31 = 250)
          )
        )
      } else {
        mod <- inla(
          form,
          family = "gaussian",
          data = inla.stack.data(stack_type),
          control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stack_type))
        )
      }
      
      if (is.null(mod) || is.null(mod$summary.fixed)) {
        cat("Model fitting failed for", type, "in sensitivity", sens_i, "\n")
        next
      }
      
      # Extract marginal effects and coefficients
      marginal_effects <- extract_marginal_effects(
        mod, covariates,
        pred_indexes[["richness"]], pred_indexes[["bio4"]], pred_indexes[["bio31"]],
        data_est_type, type
      ) %>%
        mutate(predictor = ifelse(variable == "richness", predictor^2, predictor))
      
      plot_coefficients <- tibble(
        parameter = c("richness", "bio4", "bio31"),
        richness_median = mod$summary.fixed["richness", "0.5quant"],
        richness_lower = mod$summary.fixed["richness", "0.025quant"],
        richness_upper = mod$summary.fixed["richness", "0.975quant"],
        bio4_median = mod$summary.fixed["bio4", "0.5quant"],
        bio4_lower = mod$summary.fixed["bio4", "0.025quant"],
        bio4_upper = mod$summary.fixed["bio4", "0.975quant"],
        bio31_median = mod$summary.fixed["bio31", "0.5quant"],
        bio31_lower = mod$summary.fixed["bio31", "0.025quant"],
        bio31_upper = mod$summary.fixed["bio31", "0.975quant"]
      ) %>%
        mutate(
          annotation = case_when(
            parameter == "richness" ~ sprintf("β = %.3f [%.3f, %.3f]", richness_median, richness_lower, richness_upper),
            parameter == "bio4" ~ sprintf("β = %.3f [%.3f, %.3f]", bio4_median, bio4_lower, bio4_upper),
            parameter == "bio31" ~ sprintf("β = %.3f [%.3f, %.3f]", bio31_median, bio31_lower, bio31_upper),
          ),
          x_pos = Inf, y_pos = Inf,
          mycorrhizal_type = type
        )
      
      # Prepare observed data for plotting
      raw_data <- data_est_type %>%
        select(
          mycorrhizal_type,
          response = env_breadth_orig,
          richness = richness_orig,
          bio4 = bio4_orig,
          bio31 = bio31_orig,
          richness_std = richness,
          bio4_std = bio4,
          bio31_std = bio31
        ) %>%
        mutate(richness = richness^2) %>%
        pivot_longer(
          cols = c("richness", "bio4", "bio31", "richness_std", "bio4_std", "bio31_std"),
          names_to = "variable",
          values_to = "predictor"
        )
      
      # Store results
      sensitivity_model_info[[sens_i]][[type]] <- list(
        plot_marginal_effects = marginal_effects,
        plot_coefficients = plot_coefficients,
        plot_data = raw_data
      )
      
      cat("Successfully completed", type, "for sensitivity", sens_i, "\n")
      
    }, error = function(e) {
      cat("Error processing", type, "in sensitivity", sens_i, ":", e$message, "\n")
    })
  }
}

# (4) Create sensitivity plots #################################################

# Prepare data for all analyses (main + sensitivity)
analysis_list_niche <- list()

# Main analysis
analysis_list_niche[["main"]] <- list(
  marginal = bind_rows(
    model_info[["AM"]]$plot_marginal_effects,
    model_info[["EcM"]]$plot_marginal_effects,
    model_info[["Dual"]]$plot_marginal_effects,
    model_info[["NM"]]$plot_marginal_effects
  ),
  coefficients = bind_rows(
    model_info[["AM"]]$plot_coefficients,
    model_info[["EcM"]]$plot_coefficients,
    model_info[["Dual"]]$plot_coefficients,
    model_info[["NM"]]$plot_coefficients
  ),
  obs = bind_rows(
    model_info[["AM"]]$plot_data,
    model_info[["EcM"]]$plot_data,
    model_info[["Dual"]]$plot_data,
    model_info[["NM"]]$plot_data
  )
)

# Sensitivity analyses
for (i in 1:5) {
  analysis_name <- paste0("sens", i)
  
  if (length(sensitivity_model_info) >= i && !is.null(sensitivity_model_info[[i]])) {
    sens_data <- sensitivity_model_info[[i]]
    
    # Check if we have data for all types
    available_types <- names(sens_data)
    
    if (length(available_types) > 0) {
      analysis_list_niche[[analysis_name]] <- list(
        marginal = bind_rows(lapply(sens_data, function(x) x$plot_marginal_effects)),
        coefficients = bind_rows(lapply(sens_data, function(x) x$plot_coefficients)),
        obs = bind_rows(lapply(sens_data, function(x) x$plot_data))
      )
    }
  }
}

# Create sensitivity analysis plots for each variable
niche_sensitivity_plots <- list()

for (var in c("richness", "bio4", "bio31")) {
  p_list <- list()
  
  for (anal_name in names(analysis_list_niche)) {
    anal <- analysis_list_niche[[anal_name]]
    
    if (is.null(anal)) next
    
    # Filter data for standardized predictors but original scale response
    df_current <- anal$obs %>%
      filter(str_detect(variable, paste0(var, "_std"))) %>%  # Regular scale predictors
      mutate(
        predictor_std = predictor,
        mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
      )
    
    # Filter marginal effects
    names(anal$marginal)
    filtered_effects <- anal$marginal %>%
      filter(variable == var) %>%
      mutate(
        mycorrhizal_type = if (anal_name == "main") mycorrhizal_type else model,
        mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
      )
    
    # Filter coefficients
    filtered_coefficients <- anal$coefficients %>%
      filter(parameter == var) %>%
      mutate(
        mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
      )
    
    # Create plot
    p <- ggplot(data = df_current, aes(x = predictor_std, y = response)) +
      geom_hex(bins = 30, aes(
        alpha = after_stat(ndensity),
        fill = after_stat(ndensity)
      )) +
      scale_alpha(range = c(0.8, 1)) +
      scale_fill_gradient(
        low = "#EEEEEE",
        high = "#212121"
      ) +
      geom_ribbon(
        data = filtered_effects,
        aes(ymin = lower, ymax = upper), fill = "#3366FF", alpha = 0.2
      ) +
      geom_line(
        data = filtered_effects,
        aes(y = response), colour = "#3366FF", linewidth = 0.5
      ) +
      geom_text(
        data = filtered_coefficients,
        aes(x = x_pos, y = y_pos, label = annotation),
        size = 2.25,
        hjust = 1.04,
        vjust = 1.5
      ) +
      facet_wrap(~mycorrhizal_type, nrow = 1) +
      scale_y_continuous(
        limits = c(0, 0.7),
        sec.axis = sec_axis(~ .,
                            name = case_when(
                              anal_name == "main" ~ "Main Analysis",
                              anal_name == "sens1" ~ "Sensitivity 1",
                              anal_name == "sens2" ~ "Sensitivity 2",
                              anal_name == "sens3" ~ "Sensitivity 3",
                              anal_name == "sens4" ~ "Sensitivity 4",
                              anal_name == "sens5" ~ "Sensitivity 5",
                              TRUE ~ ""
                            ),
                            breaks = NULL, labels = NULL)
      ) +
      labs(
        x = case_when(
          var == "richness" ~ "Tree richness",
          var == "bio4" ~ "Temperature variability (Bio4)",
          var == "bio31" ~ "Moisture variability (Bio31)",
          TRUE ~ var
        ),
        y = case_when(
          anal_name == "sens3" ~ "Environmental breadth",
          TRUE ~ ""
        )
      ) +
      common_theme +
      theme(
        strip.text = if(anal_name == "main") element_text(face = "bold", size = strip_size) else element_blank(),
        axis.text.x = if(anal_name == "sens5") element_text(size = text_size) else element_blank(),
        axis.title.x = if(anal_name == "sens5") element_text(size = title_size) else element_blank(),
        axis.title.y.left = if(anal_name == "sens3") element_text(size = title_size, hjust = -0.4) else element_blank(),
        axis.title.y.right = element_text(size = title_size, face = "bold")
      )
    
    # Add conditional x-axis scaling for richness
    # if (var == "richness") {
    #   p <- p + scale_x_sqrt()
    # }
    
    p_list[[anal_name]] <- p
  }
  
  # Wrap into one plot for this variable
  niche_sensitivity_plots[[var]] <- wrap_plots(p_list, ncol = 1)
}

# Display the plots
print(niche_sensitivity_plots[["richness"]])
print(niche_sensitivity_plots[["bio4"]])  
print(niche_sensitivity_plots[["bio31"]])

# Save plots
ggsave("output/supplimentary_community_niche_breadth/sensitivity_richness.png", 
       niche_sensitivity_plots[["richness"]],
       width = 14, height = 20, units = "cm", dpi = 300)
ggsave("output/supplimentary_community_niche_breadth/sensitivity_bio4.png", 
       niche_sensitivity_plots[["bio4"]],
       width = 14, height = 20, units = "cm", dpi = 300)
ggsave("output/supplimentary_community_niche_breadth/sensitivity_bio31.png", 
       niche_sensitivity_plots[["bio31"]],
       width = 14, height = 20, units = "cm", dpi = 300)
