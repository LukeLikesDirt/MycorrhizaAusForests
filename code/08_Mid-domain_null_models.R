
# Required packages
require(data.table)
require(betareg)
require(parameters)
require(performance)
require(parallel)
require(foreach)
require(doParallel)
require(ggeffects)
require(tidyverse)

# (1) Organise data --------------------------------------------------------

# Define domain boundaries: North-south latitude limits for Australia
lat_min <- 10.5
lat_max <- 44.5

# Read in the data
data <- fread(
  "generated_data/niche_estimates.txt"
) %>%
  filter(mycorrhizal_type != "ErM") %>%
  mutate(
    # Sqrt-root transformation to normalise env_breadth
    env_breadth = sqrt(env_B2_corrected),
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

# (2) Predict observed values ----------------------------------------------

# Initialise lists
pred_list <- list()
r2_list <- list()
preds_species_list <- list()

# Fit models for each mycorrhizal type
for (myc_type in unique(data$mycorrhizal_type)) {
  subset_data <- data[data$mycorrhizal_type == myc_type, ]

  # Fit GLM with second-order polynomial
  model <- betareg(env_breadth ~ poly(latitude, 2),
                   data = subset_data,
                   link = "logit"
  )

  # Get model R2
  r2_list[[myc_type]] <- data.frame(
    mycorrhizal_type = myc_type,
    R2 = performance(model)$R2
  )
  # Get predictions with CIs across the entire domain
  preds <- ggpredict(model, terms = list(latitude = seq(lat_min, lat_max, by = 0.01)))
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

# (3) Mid-domain effect (MDE) null model ------------------------------------
# This tests if species' latitudinal ranges were randomly placed within the
# domain, given their observed range sizes, what pattern would result from
# fitting models to the randomly placed data?

n_sims <- 1000

# MDE null model constants
n_species <- 100 # <-- Number of species to sample in each simulation

# Create latitude grid for predictions
lat_range_grid <- seq(lat_min, lat_max, by = 0.01)

# Set up parallel processing
n_cores <- detectCores() - 1
message(paste0("Using ", n_cores, " cores for parallel processing"))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Seed the cluster's RNG stream so results are reproducible. set.seed() in the
# main process does not propagate to %dopar% workers, which each start from
# an independent, non-reproducible stream unless seeded explicitly.
clusterSetRNGStream(cl, 1986)

# Export necessary objects to cluster
clusterExport(cl, c("data", "lat_min", "lat_max", "n_species",
                    "n_sims", "lat_range_grid"))

# Get unique mycorrhizal types
myc_types <- unique(data$mycorrhizal_type)

# Run MDE null model in parallel across mycorrhizal types
mde_pred_list <- foreach(
  myc_type = myc_types,
  .packages = c('dplyr', 'betareg')) %dopar% {

    message(paste0("Processing MDE null model for mycorrhizal type: ", myc_type))

    subset_data <- data[data$mycorrhizal_type == myc_type, ]

    # Get species-level data with range sizes
    species_data <- subset_data %>%
      group_by(species) %>%
      summarise(
        lat_range = first(lat_range),
        env_breadth = first(env_breadth)
      ) %>%
      filter(!is.na(lat_range))

    # Pre-calculate feasible midpoint ranges
    species_data <- species_data %>%
      mutate(
        min_midpoint = lat_min + (lat_range / 2),
        max_midpoint = lat_max - (lat_range / 2)
      )

    # Storage for model predictions from each simulation
    null_prediction_sims <- matrix(NA, nrow = length(lat_range_grid), ncol = n_sims)

    for (sim in 1:n_sims) {
      if (sim %% 100 == 0) {
        message(paste0("  - Completed ", sim, " of ", n_sims, " simulations"))
      }

      # Sample n_species randomly
      sample_size <- min(n_species, nrow(species_data))
      sampled_species <- species_data %>%
        slice_sample(n = sample_size, replace = FALSE)

      # Randomly place each species' range midpoint
      sampled_species <- sampled_species %>%
        mutate(
          random_midpoint = runif(n(), min_midpoint, max_midpoint),
          range_start = random_midpoint - (lat_range / 2),
          range_end = random_midpoint + (lat_range / 2)
        )

      # For each latitude bin, identify which species' ranges overlap
      # and create observation rows
      sim_observations <- lapply(lat_range_grid, function(lat) {
        overlapping <- sampled_species %>%
          filter(range_start <= lat & range_end >= lat)

        if (nrow(overlapping) > 0) {
          data.frame(
            latitude = lat,
            env_breadth = overlapping$env_breadth
          )
        } else {
          NULL
        }
      })

      # Combine all observations
      sim_data <- do.call(rbind, sim_observations)

      # Only fit model if we have enough data
      if (!is.null(sim_data) && nrow(sim_data) > 10) {

        # Fit betareg model with second-order polynomial
        tryCatch({
          model <- betareg(env_breadth ~ poly(latitude, 2),
                           data = sim_data,
                           link = "logit")

          # Get predictions across the latitude grid
          pred_data <- data.frame(latitude = lat_range_grid)
          preds <- predict(model, newdata = pred_data, type = "response")

          null_prediction_sims[, sim] <- preds
        }, error = function(e) {
          # If model fails, store NAs
          null_prediction_sims[, sim] <- NA
        })
      } else {
        null_prediction_sims[, sim] <- NA
      }
    }

    message(paste0("  - Completed all ", n_sims, " simulations for ", myc_type))

    # Calculate empirical mean and 95% CIs from model predictions
    mde_summary <- data.frame(
      x = lat_range_grid,
      predicted = apply(null_prediction_sims, 1, mean, na.rm = TRUE),
      conf.low = apply(null_prediction_sims, 1, quantile, probs = 0.025, na.rm = TRUE),
      conf.high = apply(null_prediction_sims, 1, quantile, probs = 0.975, na.rm = TRUE),
      mycorrhizal_type = myc_type
    )

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

# (4) Proportion of species within the MDE null model's 95% CIs -------------

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

# Take a glimpse at the summaries
prop_within_mde %>% print()
r2_combined %>% print()

# (5) Save the data ----------------------------------------------------------

save(
  data, pred_combined, r2_combined, preds_species_combined,
  mde_combined, prop_within_mde,
  file = "generated_data/figure_6.RData"
)
