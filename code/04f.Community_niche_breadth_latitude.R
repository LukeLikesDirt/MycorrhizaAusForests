
# Packages
require(INLA) # inla.list.models()
require(fmesher)
require(ggtext)
require(tidyverse)
source("code/map_australia.R")
set.seed(1986)

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
create_covariates <- function(data, n_pred = 100) {
  bind_rows(
    # Data for estimation
    tibble(
      latitude = data[["latitude"]],
      tag = rep("est", nrow(data))
    ),
    # Marginal effects of latitude
    tibble(
      latitude = seq(min(data[["latitude"]]), max(data[["latitude"]]), length.out = n_pred),
      tag = "pred_latitude"
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

# Helper functions to marginal effects stack
create_prediction_stack <- function(tag, covariates, n_pred = 100) {
  inla.stack(
    tag = tag,
    data = list(env_breadth = rep(NA_real_, n_pred)),
    A = list(1, 1),
    effects = list(
      intercept = rep(1, n_pred),
      X = covariates %>% filter(tag == !!tag) %>% select(-tag)
    )
  )
}

# Helper function to extract and back-transformation marginal effects
extract_marginal_effects <- function(model, covariates, index_latitude, 
                                     original_data, model_name) {
  
  # Get scaling parameters from original data
  env_breadth_center <- attr(original_data$env_breadth, "scaled:center")
  env_breadth_scale <- attr(original_data$env_breadth, "scaled:scale")
  latitude_center <- attr(original_data$latitude, "scaled:center")
  latitude_scale <- attr(original_data$latitude, "scaled:scale")
  
  # Extract predictions and back-transform
  latitude_effects <- tibble(
    response_scaled = model$summary.linear.predictor$mean[index_latitude],
    predictor_scaled = covariates %>% filter(tag == "pred_latitude") %>% pull(latitude),
    upper_scaled = model$summary.linear.predictor$`0.975quant`[index_latitude],
    lower_scaled = model$summary.linear.predictor$`0.025quant`[index_latitude],
    variable = "latitude",
    model = model_name
  ) %>%
    mutate(
      response = response_scaled * env_breadth_scale + env_breadth_center,
      predictor = predictor_scaled * latitude_scale + latitude_center,
      upper = upper_scaled * env_breadth_scale + env_breadth_center,
      lower = lower_scaled * env_breadth_scale + env_breadth_center
    )
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
    cell, longitude, latitude, x_albers, y_albers,
    species = scientific_name, mycorrhizal_type
  ) %>%
  inner_join(
    species_niche_data %>%
      select(species, env_breadth),
    by = c("species")
  ) %>%
  group_by(mycorrhizal_type, cell, longitude, latitude, x_albers, y_albers) %>%
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
  ) %>%

  mutate(
    # Make latitude absolute
    latitude = abs(latitude),
    # Define tropical/non-tropical
    climate_zone = ifelse(latitude < 23.45, "Tropical", "Non-tropical"),
    # Albers to km
    x_km = x_albers / 1000,
    y_km = y_albers / 1000
  ) 

#### * Standardise and split data * ####

# Standardise full dataset
data_std <- data_est %>%
  mutate(
    env_breadth_orig = env_breadth,
    latitude_orig = latitude,
    env_breadth = scale(env_breadth),
    latitude = scale(abs(latitude))
  )

# Standardise tropical subset
data_tropical_std <- data_est %>%
  filter(climate_zone == "Tropical") %>%
  mutate(
    env_breadth_orig = env_breadth,
    latitude_orig = latitude,
    env_breadth = scale(env_breadth),
    latitude = scale(abs(latitude))
  )

# Standardise non-tropical subset
data_nontropical_std <- data_est %>%
  filter(climate_zone == "Non-tropical") %>%
  mutate(
    env_breadth_orig = env_breadth,
    latitude_orig = latitude,
    env_breadth = scale(env_breadth),
    latitude = scale(abs(latitude))
  )

# Split the data by mycorrhizal type
data_list <- split(data_std, data_std$mycorrhizal_type)
data_list_tropical <- split(data_tropical_std, data_tropical_std$mycorrhizal_type)
data_list_nontropical <- split(data_nontropical_std, data_nontropical_std$mycorrhizal_type)


# (2) Fit models ###############################################################

# Evaluate distances between estimation sites
lapply(data_list, function(df) {
  hist(dist(as.matrix(df %>% select(x_km, y_km))), 
       main = "", 
       xlab = "Distance between sites (km)")
})
lapply(data_list_tropical, function(df) {
  hist(dist(as.matrix(df %>% select(x_km, y_km))), 
       main = "", 
       xlab = "Distance between sites (km)")
})
lapply(data_list_nontropical, function(df) {
  hist(dist(as.matrix(df %>% select(x_km, y_km))), 
       main = "", 
       xlab = "Distance between sites (km)")
})
# Compute percentiles of distances
lapply(data_list, function(df) {
  quantile(dist(as.matrix(df %>% select(x_km, y_km))), probs = seq(0, 0.5, 0.05))
})
lapply(data_list_tropical, function(df) {
  quantile(dist(as.matrix(df %>% select(x_km, y_km))), probs = seq(0, 0.5, 0.05))
})
lapply(data_list_nontropical, function(df) {
  quantile(dist(as.matrix(df %>% select(x_km, y_km))), probs = seq(0, 0.5, 0.05))
})

# 300 km seems to be a good estimate for the "small scale" distance without
# needing individually tune for each model

# Define max edge for mesh and range and sigma for the SPDE model:
# Generating a mesh 1/5 the "small scale" distance is generally the most 
# ideal for estimating fine-scale autocorrelation in the data, with a minimum
# edge length 1/5 of the "small scale" distance. 
small_scale <- 300 # <-- small scale based on sample distribution
max_edge <- small_scale / 5
range_est <- 300 # <-- range estimte for SPDE priors - check hyperparameters for ran post mod to evaluate misspecification
sigma_est <- 1 # <-- sigma estimate for SPDE priors - check hyperparameters for ran post mod to evaluate misspecification

# Fit the models for each mycorrhizal type in each climate zone
model_info <- list()
marginal_effects <- list()
coefficients <- list()
raw_data <- list()

# Define types and zones
zones <- c("all", "tropical", "nontropical")
types <- c("AM", "EcM", "Dual", "NM")

# Fit models
for (type in types) {
  message("\n===== Processing type: ", type, " =====")
  model_info[[type]] <- list()
  
  # Prepare raw plotting data for each type
  if (!is.null(data_list[[type]]) && nrow(data_list[[type]]) > 0) {
    raw_data[[type]] <- data_list[[type]] %>% select(mycorrhizal_type, response = env_breadth_orig, predictor = latitude_orig)
  } else raw_data[[type]] <- tibble()
  
  for (zone in zones) {
    message("-- zone: ", zone)
    
    df_est <- switch(zone,
                     all = data_list[[type]],
                     tropical = data_list_tropical[[type]],
                     nontropical = data_list_nontropical[[type]])
    
    # Locations and mesh
    loc_est <- df_est %>% select(x_km, y_km) %>% as.matrix()
    mesh <- fm_mesh_2d_inla(loc = loc_est, max.edge = c(1,5)*max_edge, cutoff = max_edge/5)
    if (exists("aus_map_albers")) fm_crs(mesh) <- terra::crs(aus_map_albers$data)
    
    # Projector matrix
    A_est <- inla.spde.make.A(mesh = mesh, loc = loc_est)
    
    # SPDE and spatial index
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = c(range_est, 0.01), prior.sigma = c(sigma_est, 0.01))
    index_spatial <- inla.spde.make.index(name = "spatial", n.spde = spde$n.spde)
    
    # Covariates
    covariates <- create_covariates(df_est)
    
    # Stack the data
    stack_est <- create_estimation_stack("est", df_est, covariates, A_est, index_spatial)
    stack_pred_lat <- create_prediction_stack("pred_latitude", covariates)
    stack_full <- inla.stack.join(stack_est, stack_pred_lat)
    
    # Indexing for extraction
    index_est <- inla.stack.index(stack_full, "est")$data
    index_pred_latitude <- inla.stack.index(stack_full, "pred_latitude")$data
    
    # Define the formula
    if (zone == "all") {
      form <- as.formula("env_breadth ~ -1 + intercept + latitude + I(latitude^2) + f(spatial, model = spde)")
    } else {
      form <- as.formula("env_breadth ~ -1 + intercept + latitude + f(spatial, model = spde)")
    }
    
    # Fit the model
    mod <- tryCatch({
      inla(form, family = 'gaussian', data = inla.stack.data(stack_full), control.predictor = list(compute = TRUE, A = inla.stack.A(stack_full)))
    }, error = function(e) {
      warning("INLA failed: ", e$message); NULL
    })
    
    if (is.null(mod) || is.null(mod$summary.fixed)) {
      warning("Model failed for ", type, " (", zone, ")")
      next
    }
    
    # Extract marginal effects
    marginal_effects_type <- extract_marginal_effects(mod, covariates, index_pred_latitude, df_est, paste0(type, if(zone == "all") "" else paste0("_", zone)))
    
    # Extract fitted and observed values
    fitted_values <- tibble(
      fitted = mod$summary.fitted[["mean"]][index_est],
      observed = df_est[["env_breadth"]]
    )
    
    # Store
    model_info[[type]][[zone]] <- list(
      model = mod,
      mesh = mesh,
      spde = spde,
      index_spatial = index_spatial,
      covariates = covariates,
      index_est = index_est,
      index_pred_latitude = index_pred_latitude,
      marginal_effects = marginal_effects_type,
      fitted_values = fitted_values
    )
    
    assign(paste0("model_", type, "_", zone), mod, envir = .GlobalEnv)
    
    message(sprintf("-> Completed %s (%s): n=%d", type, zone, nrow(df_est)))
    
  } # zones
  
  # combine marginal effects for this type
  me_objs <- purrr::map(model_info[[type]], "marginal_effects")
  me_objs <- purrr::keep(me_objs, ~ !is.null(.x) && nrow(.x) > 0)
  marginal_effects[[type]] <- if (length(me_objs) > 0) bind_rows(me_objs) else tibble()
  
  # coefficient annotations (match your AM example format)
  model_all <- model_info[[type]][["all"]]$model %||% NULL
  model_trop <- model_info[[type]][["tropical"]]$model %||% NULL
  model_non  <- model_info[[type]][["nontropical"]]$model %||% NULL
  
  coef_tbl <- tibble(
    model = c(paste0(type, "_1"), paste0(type, "_2"), paste0(type, "_tropical"), paste0(type, "_nontropical")),
    intercept_median = c(
      if (!is.null(model_all)) model_all$summary.fixed["intercept","0.5quant"] else NA,
      if (!is.null(model_all)) model_all$summary.fixed["intercept","0.5quant"] else NA,
      if (!is.null(model_trop)) model_trop$summary.fixed["intercept","0.5quant"] else NA,
      if (!is.null(model_non)) model_non$summary.fixed["intercept","0.5quant"] else NA
    ),
    latitude_median = c(
      if (!is.null(model_all)) model_all$summary.fixed["latitude","0.5quant"] else NA,
      if (!is.null(model_all)) model_all$summary.fixed["I(latitude^2)","0.5quant"] else NA,
      if (!is.null(model_trop)) model_trop$summary.fixed["latitude","0.5quant"] else NA,
      if (!is.null(model_non)) model_non$summary.fixed["latitude","0.5quant"] else NA
    ),
    latitude_lower = c(
      if (!is.null(model_all)) model_all$summary.fixed["latitude","0.025quant"] else NA,
      if (!is.null(model_all)) model_all$summary.fixed["I(latitude^2)","0.025quant"] else NA,
      if (!is.null(model_trop)) model_trop$summary.fixed["latitude","0.025quant"] else NA,
      if (!is.null(model_non)) model_non$summary.fixed["latitude","0.025quant"] else NA
    ),
    latitude_upper = c(
      if (!is.null(model_all)) model_all$summary.fixed["latitude","0.975quant"] else NA,
      if (!is.null(model_all)) model_all$summary.fixed["I(latitude^2)","0.975quant"] else NA,
      if (!is.null(model_trop)) model_trop$summary.fixed["latitude","0.975quant"] else NA,
      if (!is.null(model_non)) model_non$summary.fixed["latitude","0.975quant"] else NA
    )
  )
  
  coefficients[[type]] <- coef_tbl %>%
    mutate(annotation = case_when(
      model == paste0(type, "_1") ~ sprintf("β = %.3f [%.3f, %.3f]", latitude_median, latitude_lower, latitude_upper),
      model == paste0(type, "_2") ~ sprintf("β = %.3f [%.3f, %.3f]", latitude_median, latitude_lower, latitude_upper),
      model == paste0(type, "_tropical") ~ sprintf("β = %.3f [%.3f, %.3f]", latitude_median, latitude_lower, latitude_upper),
      model == paste0(type, "_nontropical") ~ sprintf("β = %.3f [%.3f, %.3f]", latitude_median, latitude_lower, latitude_upper),
      TRUE ~ NA_character_
    ), x_pos = Inf, y_pos = Inf, parameter = "latitude")
  
}

# (3) Validate ################################################################

# Collapse all fitted data into one frame with model labels
fitted_combined <- map2(names(model_info), model_info, function(main_group, sublist) {
  map2(names(sublist), sublist, function(subgroup, model_entry) {
    model_entry$fitted_values %>%
      mutate(model = paste(main_group, subgroup, sep = "_"))
  })
}) %>%
  flatten() %>%
  bind_rows()

# Posterior predictive checks
fitted_combined %>%
  ggplot(aes(x = fitted, y = observed)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label)), colour = "red") +
  facet_wrap(~ model, nrow = 4) +
  labs(
    x = "Posterior mean predicitons",
    y = "Observed values",
  ) +
  common_theme

# (4) Save #####################################################################

# Save marginal effect data for plotting
raw_data <- bind_rows(raw_data) %>%
  # Level mycorrhizal types
  mutate(
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "EcM", "Dual", "NM")
    )
  )

# Marginal effects data
marginal_effects <- bind_rows(
  marginal_effects[["AM"]] %>% mutate(mycorrhizal_type = "AM"),
  marginal_effects[["EcM"]] %>% mutate(mycorrhizal_type = "EcM"),
  marginal_effects[["Dual"]] %>% mutate(mycorrhizal_type = "Dual"),
  marginal_effects[["NM"]] %>% mutate(mycorrhizal_type = "NM")
  ) %>%
  # Level mycorrhizal types
  mutate(
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "EcM", "Dual", "NM")
    )
  )

# Coefficients data
coefficients <- bind_rows(
  coefficients[["AM"]] %>%
    mutate(mycorrhizal_type = "AM"),
  coefficients[["EcM"]] %>%
    mutate(mycorrhizal_type = "EcM"),
  coefficients[["Dual"]] %>%
    mutate(mycorrhizal_type = "Dual"),
  coefficients[["NM"]] %>%
    mutate(mycorrhizal_type = "NM")
) %>%
  # Level mycorrhizal types
  mutate(
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "EcM", "Dual", "NM")
    )
  )

# Save the data for plotting
save(
  raw_data,
  marginal_effects,
  coefficients,
  file = "output/generated_data/figure_5b.RData"
)

# Save the model information
save(
  model_info,
  file = "output/supplimentary_community_niche_breadth/model_info_latitude.RData"
)

