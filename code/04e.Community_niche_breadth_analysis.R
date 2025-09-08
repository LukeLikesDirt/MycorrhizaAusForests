
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
      quantiles = sprintf("[%.3f, %.3f]", lower, upper),
      parameter = recode(parameter, 
                         "intercept" = "Intercept",
                         "richness" = "Tree richness",
                         "bio4" = "Temperature variability (Bio4)",
                         "bio31" = "Soil moisture variability (Bio31)",
                         "Precision for the Gaussian observations" = "Observation Precision",
                         "Range for spatial" = "Spatial Range",
                         "Stdev for spatial" = "Spatial SD"
      ),
      mycorrhizal_type = myc_type
    ) %>%
    select(
      parameter, median, quantiles, mycorrhizal_type
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

#### * Evaluate VIF * ####

# Evaluate the variation inflation factor in the predictors
for (myc_type in names(data_est_list)) {
  cat("VIF for", myc_type, ":\n")
  vif_values <- car::vif(
    lm(env_breadth ~ bio4 + bio31 + richness, data = data_est_list[[myc_type]])
    )
  print(vif_values)
}

# (2) Fit models ###############################################################

# Evaluate distances between estimation sites
lapply(data_est_list, function(df) {
  hist(dist(as.matrix(df %>% select(x_km, y_km))), 
       main = "", 
       xlab = "Distance between sites (km)")
})
# Compute percentiles of distances
lapply(data_est_list, function(df) {
  quantile(dist(as.matrix(df %>% select(x_km, y_km))), probs = seq(0, 0.5, 0.05))
})

# Define max edge for mesh and range and sigma for the SPDE model:
# Generating a mesh 1/5 the "small scale" distance is genrally the most 
# ideal for estimating fine-scale autocorrelation in the data, with a minimum
# edge length 1/5 of the "small scale" distance. 
small_scale <- 300 # <-- small scale based on sample distribution
max_edge <- small_scale / 5
range_est <- 300 # <-- range estimte for SPDE priors - check hyperparameters for ran post mod to evaluate misspecification
sigma_est <- 1 # <-- sigma estimate for SPDE priors - check hyperparameters for ran post mod to evaluate misspecification

# Fit the models for each mycorrhizal type
model_info <- list()

for (type in c("AM", "EcM", "Dual", "NM")) {
  message("----- Processing: ", type, " -----")
  
  tryCatch({
    # --- 1. Define estimation and prediction data ---
    data_est_type <- data_est_list[[type]]
    data_pred_type <- data_pred_list[[type]]
    
    message("Sample sizes (", type, ") :\n", nrow(data_est_type))
    
    # --- 2. Build mesh for this model ---
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
    
    # Verify mesh size matches SPDE
    message("Mesh size: ", mesh$n, "\nSPDE size: ", spde$n.spde)
    
    # --- 3. Projector matrices ---
    A_est  <- inla.spde.make.A(mesh, loc = loc_est)
    A_pred <- inla.spde.make.A(mesh, loc = loc_pred)
    
    # --- 4. Covariates ---
    covariates <- create_covariates(data_est_type, data_pred_type)
    
    # --- 5. Stacks ---
    stack_est <- create_estimation_stack(
      "est", 
      data_est_type, 
      covariates, 
      A_est, 
      index_spatial
    )
    
    stack_pred <- create_prediction_stack(
      "pred",
      covariates,
      A_pred,
      index_spatial
    )
    
    stack_pred_richness <- create_marginal_stack("pred_richness", covariates)
    stack_pred_bio4 <- create_marginal_stack("pred_bio4", covariates)
    stack_pred_bio31 <- create_marginal_stack("pred_bio31", covariates)
    
    stack_type <- inla.stack.join(
      stack_est, 
      stack_pred,
      stack_pred_richness, 
      stack_pred_bio4, 
      stack_pred_bio31
    )
    
    idx <- list(
      est = inla.stack.index(stack_type, "est")$data,
      pred = inla.stack.index(stack_type, "pred")$data
    )
    
    pred_indexes <- list(
      richness = inla.stack.index(stack_type, "pred_richness")$data,
      bio4 = inla.stack.index(stack_type, "pred_bio4")$data,
      bio31 = inla.stack.index(stack_type, "pred_bio31")$data
    )
    
    # --- 6. Fit the model ---
    form <- env_breadth ~ -1 + intercept + richness + bio4 + bio31 + f(spatial, model = spde)
    
    if (type == "AM") {
      mod <- inla(
        form,
        family = "gaussian",
        data = inla.stack.data(stack_type),
        control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stack_type)),
        control.fixed = list(
          mean = list(bio31 = 0),
          prec = list(bio31 = 250) # <-- squeeze bio3 towards zero for stability
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
    
    # Check if model fitted successfully
    if (is.null(mod) || is.null(mod$summary.fixed)) {
      warning(paste("Model fitting failed for type", type))
      next
    }
    
    # --- 7. Extract model coefficients ---
    model_coefficients <- extract_model_coefficients(mod, type) %>%
      mutate(mycorrhizal_type = type)
 
    # Format coefficients for plot annotations
    plot_coefficients <- tibble(
      parameter = c("richness", "bio4", "bio31"),
      intercept_median = mod$summary.fixed["intercept", "0.5quant"],
      intercept_lower = mod$summary.fixed["intercept", "0.025quant"],
      intercept_upper = mod$summary.fixed["intercept", "0.975quant"],
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
      ) %>%
      select(parameter, annotation, x_pos, y_pos, mycorrhizal_type)

    # --- 8. Data for marginal effect plots ---
    
    # Prepare raw data on the original scale
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
      # Back-transform richness from the sqrt scale
      mutate(richness = richness^2) %>%
      pivot_longer(
        cols = c("richness", "bio4", "bio31", "richness_std", "bio4_std", "bio31_std"),
        names_to = "variable",
        values_to = "predictor"
      )
    
    # Marginal effects data
    marginal_effects <- extract_marginal_effects(
      mod, covariates, 
      pred_indexes[["richness"]], pred_indexes[["bio4"]], pred_indexes[["bio31"]],
      data_est_type, type
    ) %>%
      # Back-transform richness from the sqrt scale
      mutate(
        predictor = ifelse(variable == "richness", predictor^2, predictor)
      ) %>%
      select(
        mycorrhizal_type = model,
        variable, response, predictor, upper, lower,
        predictor_std, upper_std, lower_std
      )
     
    # --- 9. Extract the fitted values using helper function ---
    fitted_results <- extract_fitted_values(mod, idx, data_est_type, data_pred_type, type, mesh, forest_raster)
    
    # --- 10. Store model info ---
    model_info[[type]] <- list(
      data_est = data_est_type,
      data_pred = data_pred_type,
      covariates = covariates,
      mesh = mesh,
      spde = spde,
      index_spatial = index_spatial,
      index_model = idx,
      indexes_marg = pred_indexes,
      model_coefficients = model_coefficients,
      plot_coefficients = plot_coefficients,
      plot_data = raw_data,
      plot_marginal_effects = marginal_effects,
      fitted_values = fitted_results$fitted_values,
      fitted_raster = fitted_results$fitted_raster,
      fitted_raster_standardised = fitted_results$fitted_raster_standardised
    )
    
    assign(paste0("model_", type), mod, envir = .GlobalEnv)
    
    message("Successfully completed processing for type: ", type)
    
  }, error = function(e) {
    warning(paste("Error processing type", type, ":", e$message))
    # Continue to next iteration instead of stopping
  })
}

# (3) Model diagnostics ########################################################

#### * Marginal effects * ####
marginal_effect_plots <- list()
for (type in names(model_info)) {
  
  rich_min_max <- range(data_est[["richness"]], na.rm = TRUE)
  
  marginal_effects <- model_info[[type]]$plot_marginal_effects
  plot_data <- model_info[[type]]$plot_data
  
  # Split by variable and create individual plots
  variables <- unique(plot_data$variable)
  
  plots_list <- map(variables, function(var) {
    
    # Filter data for this variable
    plot_data_var <- filter(plot_data, variable == var)
    marginal_effects_var <- filter(marginal_effects, variable == var)
    
    p <- ggplot(data = plot_data_var, aes(x = predictor, y = response)) +
      geom_hex(bins = 30, aes(alpha = after_stat(ndensity)), fill = "grey30") +
      scale_alpha(range = c(0.075, 0.95)) +
      geom_smooth(
        method = "loess",
        colour = '#de2d26',
        linewidth = 0.25
      ) +
      geom_ribbon(
        data = marginal_effects_var,
        aes(ymin = lower, ymax = upper), fill = "#3366FF", alpha = 0.2,
      ) +
      geom_line(
        data = marginal_effects_var,
        aes(y = response), colour = "#3366FF", linewidth = 0.5) +
      scale_y_continuous(limits = c(0, 0.7)) +
      labs(title = paste(type, "-", var)) +
      common_theme +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = if (type == "NM") element_text() else element_blank(),
        axis.text.y = if (var == "richness") element_text() else element_blank()
      )
    
    # Add conditional x-axis scaling
    if (var == "richness") {
      p <- p + scale_x_sqrt()
    } else {
      p <- p + scale_x_continuous()
    }
    
    return(p)
  })
  
  names(plots_list) <- variables
  marginal_effect_plots[[type]] <- plots_list
}

# Extract all individual plots from the nested structure
all_plots <- list()

for (type in names(marginal_effect_plots)) {
  for (variable in names(marginal_effect_plots[[type]][1:3])) {
    plot_name <- paste(type, variable, sep = "_")
    all_plots[[plot_name]] <- marginal_effect_plots[[type]][[variable]]
  }
}

# Wrap all plots
patchwork::wrap_plots(
  all_plots,
  ncol = 3, nrow = 4
)

#### * Model performance * ####

# --- Main figure for comparison ---
sup_plot_1 <- patchwork::wrap_plots(
  create_spatial_plot(model_info[["AM"]]$fitted_raster_standardised, "predicted",
                     title = "AM", y_label = "Predicted env. breadth\n(standardised to percentiles)"),
  create_spatial_plot(model_info[["EcM"]]$fitted_raster_standardised, "predicted",
                     title = "EcM"),
  create_spatial_plot(model_info[["Dual"]]$fitted_raster_standardised, "predicted",
                     title = "Dual"),
  create_spatial_plot(model_info[["NM"]]$fitted_raster_standardised, "predicted",
                     title = "NM"),
  nrow = 1
)

# --- Predicted values (observed scale) ---
sup_plot_2 <- patchwork::wrap_plots(
  create_spatial_plot(model_info[["AM"]]$fitted_raster, "predicted", y_label = "Predicted env. breadth"),
  create_spatial_plot(model_info[["EcM"]]$fitted_raster, "predicted"),
  create_spatial_plot(model_info[["Dual"]]$fitted_raster, "predicted"),
  create_spatial_plot(model_info[["NM"]]$fitted_raster, "predicted"),
  nrow = 1
)

# --- Observed values ---
sup_plot_3 <- patchwork::wrap_plots(
  create_spatial_plot(model_info[["AM"]]$fitted_raster, "observed", y_label = "Observed env. breadth"),
  create_spatial_plot(model_info[["EcM"]]$fitted_raster, "observed"),
  create_spatial_plot(model_info[["Dual"]]$fitted_raster, "observed"),
  create_spatial_plot(model_info[["NM"]]$fitted_raster, "observed"),
  nrow = 1
)

# --- Prediction uncertainty ---
sup_plot_4 <- patchwork::wrap_plots(
  create_spatial_plot(model_info[["AM"]]$fitted_raster, "uncertainty", viridis_option = "C",
                      y_label = "Uncertainty\n(95% credible interval)"),
  create_spatial_plot(model_info[["EcM"]]$fitted_raster, "uncertainty", viridis_option = "C"),
  create_spatial_plot(model_info[["Dual"]]$fitted_raster, "uncertainty", viridis_option = "C"),
  create_spatial_plot(model_info[["NM"]]$fitted_raster, "uncertainty", viridis_option = "C"),
  nrow = 1
)

# --- Posterior predictive checks ---
all_fitted <- map(model_info, "fitted_values")
sup_plot_5 <- bind_rows(all_fitted, .id = "mycorrhizal_type") %>%
  mutate(
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "EcM", "Dual", "NM")
    )
  ) %>%
  filter(source == "est") %>%
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(size = 0.5, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  ggpubr::stat_cor(
    aes(label = after_stat(rr.label)), colour = "red"
  ) +
  facet_wrap(~mycorrhizal_type, nrow = 1) +
  xlab("Posterior mean predictions") +
  ylab("Observed values") +
  common_theme +
  theme(strip.text = element_blank())

# --- Save model diagnostics ---

# Combine the plots
model_performance <- patchwork::wrap_plots(
  sup_plot_1,
  sup_plot_2,
  sup_plot_3,
  sup_plot_4,
  sup_plot_5,
  ncol = 1
)

# Display the plots
print(model_performance)

# Save the supplementary plots
ggsave(
  filename = "output/supplimentary_community_niche_breadth/model_performance_plot.png",
  plot = model_performance,
  width = 22,
  height = 26,
  units = "cm",
  dpi = 300
)

#### * Spatial diagnostics * ####

# --- Matérn correlation function ---

# Compute statistics
spatial_stats_AM <- spatial_params(model_AM, spde, name = "spatial")
spatial_stats_EcM <- spatial_params(model_EcM, spde, name = "spatial")
spatial_stats_Dual <- spatial_params(model_Dual, spde, name = "spatial")
spatial_stats_NM <- spatial_params(model_NM, spde, name = "spatial")

# Compare prior and posterior sigma:
sigma_est
spatial_stats_AM["sigma_spatial"]
spatial_stats_EcM["sigma_spatial"]
spatial_stats_Dual["sigma_spatial"]
spatial_stats_NM["sigma_spatial"]

# Compare prior and posterior range:
# Re-fit models if posterior range is smaller than prior range
range_est
spatial_stats_AM["range"]
spatial_stats_EcM["range"]
spatial_stats_Dual["range"]
spatial_stats_NM["range"]

# Matérn correlation function plot
matern_plot <- names(model_info) %>%
  map(~{
    # Get spatial stats - adjust this path based on where they're stored
    stats <- switch(.x,
                    "AM" = spatial_stats_AM,
                    "EcM" = spatial_stats_EcM, 
                    "Dual" = spatial_stats_Dual,
                    "NM" = spatial_stats_NM
    )
    
    plot <- plot_matern_correlation(
      stats,
      x = model_info[[.x]]$data_est$x_km,  # Using data from model_info
      y = model_info[[.x]]$data_est$y_km,  # Using data from model_info
      range = stats["range"],
      y_limits = c(0, 0.8),
      y_breaks = seq(0, 0.8, by = 0.2),
      x_limits = c(0, 600), 
      x_breaks = seq(0, 600, by = 100)
    ) +
      theme(
        plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
        axis.text.x = element_text(angle = 45, hjust = 1)
        )
    
    if (.x != "AM") {
      plot <- plot + theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank()
      )
    }
    
    return(plot)
  }) %>%
  patchwork::wrap_plots(nrow = 1)

# Display the Matérn correlation plot
print(matern_plot)

# --- Moran's I ---

# Compute Moran's I for all model types
moran_results <- list()
for (model_name in names(model_info)) {
  tryCatch({
    moran_results[[model_name]] <- compute_moran_i(model_name)
    message("Computed Moran's I for ", model_name)
  }, error = function(e) {
    warning("Error computing Moran's I for ", model_name, ": ", e$message)
    moran_results[[model_name]] <- list(
      Moran_I = NA,
      p_value = NA,
      model_type = model_name
    )
  })
}

# Spatial residual plots
spatial_residual_plots <- list()

for (type in names(model_info)) {
  
  moran_str <- paste0("Moran's I = ",
                      round(unname(moran_results[[type]][["Moran_I"]]), 3))
  if (round(moran_results[[type]][["p_value"]], 3) == 1) {
    p_val_str <- "p-value = 1.000"
  } else if (round(moran_results[[type]][["p_value"]], 3) == 0) {
    p_val_str <- paste0("*p*-value = ",
                        round(moran_results[[type]][["p_value"]], 3))
  }
  
  spatial_residual_plots[[type]] <- create_spatial_plot(
    raster = model_info[[type]]$fitted_raster,
    data_column = "residuals",
    color_scale = NULL, # skip brewer scale
    title = type,
    y_label = if (type == "AM") "Residuals" else NULL
  ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white", 
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7)
    ) +
    labs(fill = NULL) +
    geom_richtext(
      x = 114,
      y = -Inf,
      label = p_val_str,
      size = 2.5,
      hjust = 0.1,
      vjust = 0,
      label.color = NA,
      fill = NA
    ) +
    geom_richtext(
      x = 114,
      y = -Inf,
      label = moran_str,
      size = 2.5,
      hjust = 0.1,
      vjust = -0.75,
      label.color = NA,
      fill = NA
    )
}

# --- Posterior spatial effect plots ---
spatial_effect_plots <- list()

for (type in names(model_info)) {
  spatial_effect_plots[[type]] <- create_spatial_plot(
    raster = model_info[[type]]$fitted_raster,
    data_column = "spatial_field",
    color_scale = NULL, # skip brewer scale
    title = type,
    y_label = if (type == "AM") "Posterior spatial effect" else NULL
  ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red", 
      midpoint = 0,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold")
    ) +
    labs(fill = NULL)
}

# --- Show spatial diagnostics plots ---

# Wrap the plots
spatial_diagnostics <- patchwork::wrap_plots(
  patchwork::wrap_plots(spatial_effect_plots, nrow = 1),
  patchwork::wrap_plots(spatial_residual_plots, nrow = 1),
  matern_plot,
  ncol = 1
)

# Display the plots
print(spatial_diagnostics)

# Save the plots
ggsave(
  filename = "output/supplimentary_community_niche_breadth/spatial_diagnostics.png",
  plot = spatial_diagnostics,
  width = 16,
  height = 13.75,
  units = "cm",
  dpi = 300
)

# (4) Environmental space ######################################################

# Load the tree data
data_env_space <- bind_rows(
  tibble(
    predicted = model_info[["AM"]][["fitted_values"]] %>% pull(predicted),
    bio4 = c(model_info[["AM"]][["data_est"]][["bio4_orig"]], model_info[["AM"]][["data_pred"]][["bio4_orig"]]),
    bio31 = c(model_info[["AM"]][["data_est"]][["bio31_orig"]], model_info[["AM"]][["data_pred"]][["bio31_orig"]]),
    bio4_std = c(model_info[["AM"]][["data_est"]][["bio4"]], model_info[["AM"]][["data_pred"]][["bio4"]]),
    bio31_std = c(model_info[["AM"]][["data_est"]][["bio31"]], model_info[["AM"]][["data_pred"]][["bio31"]]),
    mycorrhizal_type = "predicted_AM"
  ),
  tibble(
    predicted = model_info[["EcM"]][["fitted_values"]] %>% pull(predicted),
    bio4 = c(model_info[["EcM"]][["data_est"]][["bio4_orig"]], model_info[["EcM"]][["data_pred"]][["bio4_orig"]]),
    bio31 = c(model_info[["EcM"]][["data_est"]][["bio31_orig"]], model_info[["EcM"]][["data_pred"]][["bio31_orig"]]),
    bio4_std = c(model_info[["EcM"]][["data_est"]][["bio4"]], model_info[["EcM"]][["data_pred"]][["bio4"]]),
    bio31_std = c(model_info[["EcM"]][["data_est"]][["bio31"]], model_info[["EcM"]][["data_pred"]][["bio31"]]),
    mycorrhizal_type = "predicted_EcM"
  ),
  tibble(
    predicted = model_info[["Dual"]][["fitted_values"]] %>% pull(predicted),
    bio4 = c(model_info[["Dual"]][["data_est"]][["bio4_orig"]], model_info[["Dual"]][["data_pred"]][["bio4_orig"]]),
    bio31 = c(model_info[["Dual"]][["data_est"]][["bio31_orig"]], model_info[["Dual"]][["data_pred"]][["bio31_orig"]]),
    bio4_std = c(model_info[["Dual"]][["data_est"]][["bio4"]], model_info[["Dual"]][["data_pred"]][["bio4"]]),
    bio31_std = c(model_info[["Dual"]][["data_est"]][["bio31"]], model_info[["Dual"]][["data_pred"]][["bio31"]]),
    mycorrhizal_type = "predicted_Dual"
  ),
  tibble(
    predicted = model_info[["NM"]][["fitted_values"]] %>% pull(predicted),
    bio4 = c(model_info[["NM"]][["data_est"]][["bio4_orig"]], model_info[["NM"]][["data_pred"]][["bio4_orig"]]),
    bio31 = c(model_info[["NM"]][["data_est"]][["bio31_orig"]], model_info[["NM"]][["data_pred"]][["bio31_orig"]]),
    bio4_std = c(model_info[["NM"]][["data_est"]][["bio4"]], model_info[["NM"]][["data_pred"]][["bio4"]]),
    bio31_std = c(model_info[["NM"]][["data_est"]][["bio31"]], model_info[["NM"]][["data_pred"]][["bio31"]]),
    mycorrhizal_type = "predicted_NM"
  )
) %>%
  pivot_wider(
    names_from = mycorrhizal_type,
    values_from = predicted
  )

# --- Create env bins ---

# Find min and max for axis pairs
min_max_limits  <- data_env_space %>%
  summarise(
    bio4_min = min(bio4, na.rm = TRUE),
    bio4_max = max(bio4, na.rm = TRUE),
    bio31_min = min(bio31, na.rm = TRUE),
    bio31_max = max(bio31, na.rm = TRUE),
    bio4_std_min = min(bio4_std, na.rm = TRUE),
    bio4_std_max = max(bio4_std, na.rm = TRUE),
    bio31_std_min = min(bio31_std, na.rm = TRUE),
    bio31_std_max = max(bio31_std, na.rm = TRUE)
  )

# Create the break sequences outside the mutate function
bio4_breaks <- seq(min_max_limits $bio4_min, min_max_limits $bio4_max, length.out = 21)
bio4_std_breaks <- seq(min_max_limits $bio4_std_min, min_max_limits $bio4_std_max, length.out = 21)
bio31_breaks <- seq(min_max_limits $bio31_min, min_max_limits $bio31_max, length.out = 21)
bio31_std_breaks <- seq(min_max_limits $bio31_std_min, min_max_limits $bio31_std_max, length.out = 21)

# Create 2D bins for bio4 and bio31
bio4_bio31_bin_lookup <- expand.grid(
  bio4_bin = 1:20,
  bio31_bin = 1:20
)

# Assign unique bin numbers
bio4_bio31_bin_lookup[["bin_2d"]] <- 1:nrow(bio4_bio31_bin_lookup)

# Add bin numbers to the data
bio4_bio31_2d_bins <- data_env_space %>%
  mutate(
    # Get the bin numbers again (1-40)
    bio4_bin = cut(bio4, breaks = bio4_breaks, labels = FALSE, include.Low = TRUE),
    bio31_bin = cut(bio31, breaks = bio31_breaks, labels = FALSE, include.Low = TRUE)
  ) %>%
  left_join(bio4_bio31_bin_lookup, by = c("bio4_bin", "bio31_bin"))

#### * Summarise the data * ####

# bio4-bio31
data_bio4_bio31 <- bio4_bio31_2d_bins %>%
  group_by(bio4_bin, bio31_bin) %>%
  summarise(
    predicted_AM = mean(predicted_AM, na.rm = TRUE),
    predicted_EcM = mean(predicted_EcM, na.rm = TRUE),
    predicted_Dual = mean(predicted_Dual, na.rm = TRUE),
    predicted_NM = mean(predicted_NM, na.rm = TRUE)
  )

# Function to standardise richness values to percentiles for colour scale
data_bio4_bio31_quantised <- quantise_tibble_percentiles(
  data_bio4_bio31, 
  c("predicted_AM", "predicted_EcM", "predicted_Dual", "predicted_NM")
)

#### * Visualise * ####

# bio4-bio31 environmental space plot
bio42_plot <- data_bio4_bio31_quantised %>%
  mutate(
    bio4_original = min_max_limits$bio4_std_min + (bio4_bin - 0.5) * (min_max_limits$bio4_std_max - min_max_limits$bio4_std_min) / 20,
    bio31_original = min_max_limits$bio31_std_min + (bio31_bin - 0.5) * (min_max_limits$bio31_std_max - min_max_limits$bio31_std_min) / 20
  ) %>%
  pivot_longer(
    cols = -c(bio4_bin, bio31_bin, bio4_original, bio31_original),
    names_to = "mycorrhizal_type",
    values_to = "predicted"
  ) %>%
  ggplot(aes(x = bio4_original, y = bio31_original, fill = predicted)) +
  geom_tile() +
  scale_fill_stepsn(
    colors = rev(paletteer::paletteer_c("grDevices::Spectral", 100)),
    breaks = 1:100,
    labels = c("Low", rep("", 98), "High"),
    na.value = "transparent"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "grey90", fill = NA),
    legend.position = "none",
    strip.text = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.tag = element_markdown(size = 16),
    aspect.ratio = 1,
    plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
  ) +
  labs(
    x = "Temperature seasonality",
    y = "Moisture index seasonality", 
    tag = "(**c**)"
  ) +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# Display the plots
print(bio42_plot)

# (5) Save the results #########################################################

# Save raster data for figure 5a
writeRaster(
  rast(c(
    AM = model_info[["AM"]][["fitted_raster_standardised"]],
    EcM = model_info[["EcM"]][["fitted_raster_standardised"]],
    Dual = model_info[["Dual"]][["fitted_raster_standardised"]],
    NM = model_info[["NM"]][["fitted_raster_standardised"]]
  )),
  overwrite = TRUE,
  filename = "output/generated_data/figure_5a.tif"
)

# Save data for figure 2c
data_figure_5c <- data_bio4_bio31_quantised
data_figure_5c_limits <- min_max_limits
save(
  data_figure_5c,
  data_figure_5c_limits,
  file = "output/generated_data/figure_5c.RData"
)

# Organise and save the data for figure 6

# Observed data
raw_data <- bind_rows(
  model_info[["AM"]][["plot_data"]],
  model_info[["EcM"]][["plot_data"]],
  model_info[["Dual"]][["plot_data"]],
  model_info[["NM"]][["plot_data"]]
)

# Marginal effects data
marginal_effects_data <- bind_rows(
  model_info[["AM"]][["plot_marginal_effects"]],
  model_info[["EcM"]][["plot_marginal_effects"]],
  model_info[["Dual"]][["plot_marginal_effects"]],
  model_info[["NM"]][["plot_marginal_effects"]]
)

# Model coefficients
coeficients_data <- bind_rows(
  model_info[["AM"]][["plot_coefficients"]],
  model_info[["EcM"]][["plot_coefficients"]],
  model_info[["Dual"]][["plot_coefficients"]],
  model_info[["NM"]][["plot_coefficients"]]
)

# Save the data
save(
  raw_data,
  marginal_effects_data,
  coeficients_data,
  file = "output/generated_data/figure_6.RData"
)

# Save the model information
save(
  model_info,
  file = "output/supplimentary_community_niche_breadth/model_info.RData"
)

# Save the model coefficients
model_info[["AM"]][["model_coefficients"]] %>%
  bind_rows(
    model_info[["EcM"]][["model_coefficients"]],
    model_info[["Dual"]][["model_coefficients"]],
    model_info[["NM"]][["model_coefficients"]]
  ) %>%
  data.table::fwrite(
    "output/supplimentary_community_niche_breadth/model_coefficients.txt",
    sep = "\t"
    )
