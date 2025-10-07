
# Required packages
require(betareg)
require(parameters)
require(performance)
require(parallel)
require(foreach)
require(doParallel)
require(ggeffects)
require(tidyverse)

# Plot constants
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 8
common_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.grid.minor = element_line(linewidth = 0.15),
    axis.ticks = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
    strip.text = element_text(face = "bold", size = strip_size),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    aspect.ratio = 1
  )

# (1) Organise data ------------------------------------------------------------

# Define domain boundaries: North-south latitude limits for Australia
lat_min <- 10
lat_max <- 45
domain_size <- lat_max - lat_min

# Read in the data
data <- data.table::fread(
  "output/generated_data/niche_estimates.txt"
) %>%
  filter(mycorrhizal_type != "ErM") %>%
  mutate(
    env_breadth = env_B2_corrected,
    # Rename EcM-AM to Dual
    mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
    mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")),
    species = str_replace_all(species, " ", "_"),
    # Work with absolute latitude
    latitude = abs(lat_position)
  ) %>%
  select(
    family, genus, species, mycorrhizal_type, latitude, env_breadth, lat_range
  )

# (2) Predict observed values --------------------------------------------------

# Initialise lists
pred_list <- list()
r2_list <- list()
model_list <- list()
model_prams_list <- list()
preds_species_list <- list()

# Fit models for each mycorrhizal type
for (myc_type in unique(data$mycorrhizal_type)) {
  subset_data <- data[data$mycorrhizal_type == myc_type, ]
  
  # Fit GLM with second-order polynomial
  model <- betareg(env_breadth ~ poly(latitude, 2),
                   data = subset_data,
                   link = "logit"
  )
  
  # Store model
  model_list[[myc_type]] <- model
  
  # Store model parameters
  model_prams_list[[myc_type]] <- parameters(model)
  
  # Get model R2
  r2_list[[myc_type]] <- data.frame(
    mycorrhizal_type = myc_type,
    R2 = performance(model)$R2
  )
  # Get predictions with CIs across the entire domain
  preds <- ggpredict(model, terms = list(latitude = seq(lat_min + 1, lat_max - 1, by = 0.01)))  
  preds <- as.data.frame(preds)
  preds$mycorrhizal_type <- myc_type
  
  # Get predicted values for each species
  preds_species_list[[myc_type]] <- data.frame(
    species = subset_data$species,
    latitude = subset_data$latitude,
    observed = subset_data$env_breadth,
    predicted = predict(model, newdata = subset_data),
    mycorrhizal_type = myc_type
  )
  
  pred_list[[myc_type]] <- preds
}

# Combine predictions
pred_combined <- do.call(rbind, pred_list) %>%
  # Level the mycorrhizal_type factor
  mutate(mycorrhizal_type = factor(
    mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")
  ))

# Combine R2 values
r2_combined <- do.call(rbind, r2_list) %>%
  mutate(mycorrhizal_type = factor(
    mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")
  ))

# Combine species-level predictions
preds_species_combined <- do.call(rbind, preds_species_list) %>%
  mutate(mycorrhizal_type = factor(
    mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")
  ))

# (3) Null models for latitude effects -----------------------------------------

set.seed(1986)
n_sims <- 1000

null_pred_list <- list()

for (myc_type in unique(data$mycorrhizal_type)) {
  # Message when starting each mycorrhizal type
  message(paste0("Processing mycorrhizal type: ", myc_type))
  
  subset_data <- data[data$mycorrhizal_type == myc_type, ]
  
  # Storage for null predictions
  lat_range <- seq(lat_min + 1, lat_max - 1, by = 1)
  null_predictions <- matrix(NA, nrow = length(lat_range), ncol = n_sims)
  
  for (sim in 1:n_sims) {
    # Progress message every 100 iterations
    if (sim %% 100 == 0) {
      message(paste0("  - Completed ", sim, " of ", n_sims, " simulations"))
    }
    
    # Shuffle env_breadth values across latitude
    # This breaks the latitude-breadth relationship while keeping
    # the same distribution of breadth values
    sim_data <- subset_data
    sim_data$env_breadth_shuffled <- sample(sim_data$env_breadth)
    
    # Fit GLM with second-order polynomial
    model <- betareg(env_breadth_shuffled ~ poly(latitude, 2), 
                     data = sim_data)
    
    # Get predictions: ggpredict approach
    # preds <- ggpredict(model, terms = list(latitude = seq(lat_min + 1, lat_max - 1, by = 1)))
    # preds <- ggpredict(model, terms="latitude [all]")
    # null_predictions[, sim] <- preds$predicted
    
    # Get predictions
    pred_data <- data.frame(latitude = lat_range)
    preds <- predict(model, newdata = pred_data)
    null_predictions[, sim] <- preds
  }
  
  message(paste0("  - Completed all ", n_sims, " simulations for ", myc_type))
  
  # Calculate empirical mean and 95% CIs across simulations
  null_summary <- data.frame(
    x = lat_range,
    predicted = rowMeans(null_predictions),
    conf.low = apply(null_predictions, 1, quantile, probs = 0.025),
    conf.high = apply(null_predictions, 1, quantile, probs = 0.975),
    mycorrhizal_type = myc_type
  )
  
  null_pred_list[[myc_type]] <- null_summary
  
  message(paste0("Finished processing ", myc_type, "\n"))
}

# Combine null predictions
null_combined <- do.call(rbind, null_pred_list) %>%
  mutate(mycorrhizal_type = factor(
    mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")
  ))

# Create the plot
ggplot() +
  geom_point(data = data, aes(x = latitude, y = env_breadth), 
             alpha = 0.25, size = 0.3) +
  geom_ribbon(data = null_combined, 
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "#de2d26") +
  geom_line(data = null_combined, aes(x = x, y = predicted), 
            color = "#de2d26", linewidth = 0.5) +
  geom_ribbon(data = pred_combined, 
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "#3366FF") +
  geom_line(data = pred_combined, aes(x = x, y = predicted), 
            color = "#3366FF", linewidth = 0.5) +
  labs(
    x = "Latitude",
    y = "Environmental breadth"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1) +
  facet_wrap(~ mycorrhizal_type, nrow = 1)

# (4) Null model mid-domain effect ---------------------------------------------
# This tests if species' latitudinal ranges were randomly placed within the
# domain, given their observed range sizes, what pattern of overlap (and mean
# breadth) would result from geometry alone?

# Set random seed and number of simulations
set.seed(1986)
n_sims <- 1000

# MDE null model constants 
lat_min <- 10 # <-- Minimum latitude in the dataset
lat_max <- 45 # <-- Maximum latitude in the dataset
prop_sample <- 0.2 # <-- Proportion of species to sample (uncomment in the loop if using)
n_species <- 100 # <-- Number of species to sample in each simulation

# Create latitude grid for predictions
lat_range_grid <- seq(lat_min, lat_max, by = 1)

# Set up parallel processing
n_cores <- detectCores() - 1
message(paste0("Using ", n_cores, " cores for parallel processing"))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary objects to cluster
clusterExport(cl, c("data", "lat_min", "lat_max", "prop_sample", 
                    "n_sims", "lat_range_grid"))

# Get unique mycorrhizal types
myc_types <- unique(data$mycorrhizal_type)

# Run MDE null model in parallel across mycorrhizal types
mde_pred_list <- foreach(
  myc_type = myc_types, 
  .packages = c('dplyr')
) %dopar% {
  
  message(paste0("Processing MDE null model for mycorrhizal type: ", myc_type))
  
  subset_data <- data[data$mycorrhizal_type == myc_type, ]
  
  # Alternative: Compute n_species to sample based on proportion
  # # Compute n_species to sample (20% of the original data)
  # n_species <- round(prop_sample * nrow(subset_data))
  
  # Pre-calculate feasible midpoint ranges (outside simulation loop)
  min_midpoints <- lat_min + (subset_data$lat_range / 2)
  max_midpoints <- lat_max - (subset_data$lat_range / 2)
  half_ranges <- subset_data$lat_range / 2
  
  # Storage for null predictions
  null_breadth_sims <- matrix(NA, nrow = length(lat_range_grid), ncol = n_sims)
  
  for (sim in 1:n_sims) {
    if (sim %% 100 == 0) {
      message(paste0("  - Completed ", sim, " of ", n_sims, " simulations"))
    }
    
    # Randomly place midpoint within feasible bounds for each species
    random_midpoints <- runif(nrow(subset_data), min_midpoints, max_midpoints)
    
    # Calculate range endpoints
    range_starts <- random_midpoints - half_ranges
    range_ends <- random_midpoints + half_ranges
    
    # Randomly sample species for this simulation
    sample_idx <- sample(nrow(subset_data), n_species, replace = FALSE)
    sampled_starts <- range_starts[sample_idx]
    sampled_ends <- range_ends[sample_idx]
    sampled_breadths <- subset_data$env_breadth[sample_idx]
    
    # For each latitude point, calculate mean breadth of overlapping species
    breadth_at_lat <- sapply(lat_range_grid, function(lat) {
      # Which species' ranges overlap this latitude?
      overlapping_idx <- which(sampled_starts <= lat & sampled_ends >= lat)
      
      # Return mean breadth of overlapping species
      if (length(overlapping_idx) > 0) {
        mean(sampled_breadths[overlapping_idx])
      } else {
        NA
      }
    })
    
    null_breadth_sims[, sim] <- breadth_at_lat
  }
  
  message(paste0("  - Completed all ", n_sims, " simulations for ", myc_type))
  
  # Calculate empirical mean and 95% CIs across simulations
  mde_summary <- data.frame(
    x = lat_range_grid,
    predicted = rowMeans(null_breadth_sims, na.rm = TRUE),
    conf.low = apply(null_breadth_sims, 1, quantile, probs = 0.025, na.rm = TRUE),
    conf.high = apply(null_breadth_sims, 1, quantile, probs = 0.975, na.rm = TRUE),
    mycorrhizal_type = myc_type
  )
  
  # Apply loess smoothing to reduce choppiness at 95% CI edges
  # Adjust span parameter (0.1-0.3) for more/less smoothing
  smooth_span <- 0.3
  
  # Only smooth if we have enough non-NA values
  if (sum(!is.na(mde_summary$predicted)) > 10) {
    mde_summary <- mde_summary %>%
      mutate(
        predicted = predict(loess(predicted ~ x, data = ., span = smooth_span, na.action = na.exclude)),
        conf.low = predict(loess(conf.low ~ x, data = ., span = smooth_span, na.action = na.exclude)),
        conf.high = predict(loess(conf.high ~ x, data = ., span = smooth_span, na.action = na.exclude))
      )
  }
  
  message(paste0("Finished processing MDE null for ", myc_type, "\n"))
  
  return(mde_summary)
}

# Stop cluster
stopCluster(cl)

# Name the list elements
names(mde_pred_list) <- myc_types

# Combine MDE null predictions
mde_combined <- do.call(rbind, mde_pred_list) %>%
  mutate(mycorrhizal_type = factor(
    mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")
  ))

# (5) Plot MDE -----------------------------------------------------------------

# Compute the proportion of species witihin the 95% CIs of the null models
prop_within_mde <- data.frame(mycorrhizal_type = character(), prop_within = numeric())
for (myc_type in unique(data$mycorrhizal_type)) {
  subset_data <- data[data$mycorrhizal_type == myc_type, ]
  subset_mde <- mde_combined[mde_combined$mycorrhizal_type == myc_type, ]
  
  within_count <- 0
  
  for (i in 1:nrow(subset_data)) {
    lat <- subset_data$latitude[i]
    breadth <- subset_data$env_breadth[i]
    
    # Find the closest latitude in the MDE predictions
    closest_lat <- subset_mde[which.min(abs(subset_mde$x - lat)), ]
    
    # Check if the breadth is within the 95% CI
    if (!is.na(closest_lat$conf.low) && !is.na(closest_lat$conf.high)) {
      if (breadth >= closest_lat$conf.low && breadth <= closest_lat$conf.high) {
        within_count <- within_count + 1
      }
    }
  }
  
  prop_within <- within_count / nrow(subset_data)
  
  prop_within_mde <- rbind(prop_within_mde, data.frame(
    mycorrhizal_type = myc_type,
    prop_within = prop_within
  ))
}

# Annotate R2 values for observed models
r2_summary <- r2_combined %>%
  rowwise() %>%
  mutate(
    # Build a plot math-safe label where the number is a quoted string
    label = paste0('italic(R)^2 == \"', sprintf("%.2f", R2), '\"')
  ) %>%
  ungroup()

# Plot with both null models
ggplot() +
  geom_point(data = data, aes(x = latitude, y = env_breadth), 
             alpha = 0.25, size = 0.3) +
  # MDE null model
  geom_ribbon(data = mde_combined, 
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "#de2d26") +
  geom_line(data = mde_combined, aes(x = x, y = predicted), 
            color = "#de2d26", linewidth = 0.5) +
  # Observed pattern in blue
  geom_ribbon(data = pred_combined, 
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2, fill = "#3366FF") +
  geom_line(data = pred_combined, aes(x = x, y = predicted), 
            color = "#3366FF", linewidth = 0.5) +
  # Add degrees south symbol to x-axis
  scale_x_continuous(
    labels = function(x) paste0(abs(x), "°S")
  ) +
  # Add R2 annotations
  geom_text(
    data = r2_summary,
    aes(
      x = lat_min + 0.99 * domain_size,
      y = max(data$env_breadth, na.rm = TRUE) * 0.99,
      label = label
    ),
    parse = TRUE,
    hjust = 1,
    vjust = 1,
    size = 2.5,
    colour = "black",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Latitude",
    y = "Environmental breadth"
  ) +
  common_theme +
  facet_wrap(~ mycorrhizal_type, nrow = 1)

# Save the plot
ggsave(
  "output/figure6.png",
  width = 16,
  height = 5,
  units = "cm",
  bg = "white",
  dpi = 300
)
ggsave(
  "output/figure6.tif",
  width = 16,
  height = 5,
  units = "cm",
  bg = "white"
)
