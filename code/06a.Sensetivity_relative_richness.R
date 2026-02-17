
# Set the seed
set.seed(1986)

# Required packages
require(patchwork)
require(INLA) # inla.list.models()
require(fmesher)
require(terra)
require(ggtext)
require(tidyverse)
source("code/functions.R")
source("code/map_australia.R")

# Plot constants
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 8
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

# (1) Read original data ########################################################

# Read the original estimation data
data_est <- data.table::fread(
  "data/presence/sites_relative_richness_est_10.txt",
  stringsAsFactors = TRUE
) %>%
  mutate(
    x_km = x_albers / 1000, y_km = y_albers / 1000,
    richness_AM_log_ratio = log((richness_AM + 1) / (richness - richness_AM + 1)),
    richness_EcM_log_ratio = log((richness_EcM + 1) / (richness - richness_EcM + 1)),
    richness_EcM_AM_log_ratio = log((richness_EcM_AM + 1) / (richness - richness_EcM_AM + 1)),
    richness_NM_log_ratio = log((richness_NM + 1) / (richness - richness_NM + 1)),
    RC1 = std(RC1),
    RC2 = std(RC2),
    RC3 = std(RC3)
  )

# Count the number of cells for estimation
n_est <- nrow(data_est) %>% print()

# Read prediction data
data_pred <- data.table::fread(
  "data/presence/sites_relative_richness_pred_10.txt",
  stringsAsFactors = TRUE
) %>%
  mutate(x_km = x_albers / 1000, y_km = y_albers / 1000)

# Count the number of cells for prediction
n_pred <- nrow(data_pred) %>% print()

# (2) Setup spatial components (same across all analyses) ######################

# Site locations for the estimation and prediction data
loc_est <- data_est %>%
  select(x_km, y_km) %>%
  as.matrix()
loc_pred <- data_pred %>%
  select(x_km, y_km) %>%
  as.matrix()

# Distances between estimation sites
hist(dist(loc_est), main = "", xlab = "Distance between sites (km)")
small_scale <- 300

# Maximum edge for the triangulation (about 1/3 of the spatial correlation range)
max_edge <- small_scale / 5

# Prepare the mesh
mesh <- fm_mesh_2d_inla(
  loc = loc_est,
  max.edge = c(1, 5) * max_edge,
  cutoff = max_edge / 5
)

# Define the mesh CRS
fm_crs(mesh) <- crs(aus_map_albers$data)

# Visualise the mesh using ggplot
ggplot() +
  theme_minimal() +
  fmesher::geom_fm(data = mesh) +
  geom_point(
    data = data_est,
    aes(x_km, y_km), size = 0.2, alpha = 0.5
  ) +
  geom_sf(
    data = aus_map_albers$data,
    fill = "transparent",
    col = "red"
  ) +
  labs(x = NULL, y = NULL)

# Define projector matrix for the estimation and prediction data
A_est <- inla.spde.make.A(
  mesh = mesh,
  loc = loc_est
)
dim(A_est)
A_pred <- inla.spde.make.A(
  mesh = mesh,
  loc = loc_pred
)
dim(A_pred)

# Define the SPDE
range_est <- 80
sigma_est <- 1

# Define the SPDE for the estimation
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(range_est, 0.01),
  prior.sigma = c(sigma_est, 0.01)
)

# Define the SPDE index
index_spatial <- inla.spde.make.index(
  name = "spatial",
  n.spde = spde$n.spde
)

# The size of the mesh and n.spde should match
mesh$n
spde$n.spde

# (3) Create sensitivity datasets ###############################################

# Function to create sensitivity datasets by resampling mycorrhizal types
create_sensitivity_datasets <- function(original_data, n_sensitivity = 5, resample_prop = 0.2) {
  
  # Calculate mycorrhizal type probabilities from unique species
  species_types <- original_data %>%
    distinct(scientific_name, mycorrhizal_type) %>%
    count(mycorrhizal_type)
  myco_counts <- setNames(species_types$n, species_types$mycorrhizal_type)
  myco_probs <- myco_counts / sum(myco_counts)
  
  cat("Mycorrhizal type probabilities calculated from data (per species):\n")
  print(round(myco_probs, 3))
  cat("\n")
  
  # Get species per type
  species_per_type <- original_data %>%
    group_by(mycorrhizal_type) %>%
    distinct(scientific_name) %>%
    group_by(mycorrhizal_type) %>%
    nest() %>%
    mutate(species = map(data, ~pull(., scientific_name))) %>%
    select(-data) %>%
    deframe()
  
  sensitivity_datasets <- list()
  
  for (i in 1:n_sensitivity) {
    cat("Creating sensitivity dataset", i, "...\n")
    
    # Create a copy of the original data
    sens_data <- original_data
    
    # Determine number of species to resample
    n_resample <- ceiling(length(unique(sens_data$scientific_name)) * resample_prop)
    
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
      if (n_res_g > 0) {
        # Sample without replacement, cap at group size
        sampled <- sample(species_per_type[[group_name]], min(n_res_g, length(species_per_type[[group_name]])))
        resample_species <- c(resample_species, sampled)
      }
    }
    
    # Resample mycorrhizal types for selected species
    new_myco_types <- sample(
      names(myco_probs),
      size = length(resample_species),
      replace = TRUE,
      prob = myco_probs
    )
    
    # Update all rows for the selected species with new mycorrhizal types
    for (k in seq_along(resample_species)) {
      sp <- resample_species[k]
      new_type <- new_myco_types[k]
      sens_data$mycorrhizal_type[sens_data$scientific_name == sp] <- new_type
    }
    
    # Reapply factor levels and contrasts
    sens_data$mycorrhizal_type <- factor(sens_data$mycorrhizal_type, 
                                         levels = c("AM", "EcM", "EcM-AM", "NM"))
    contrasts(sens_data$mycorrhizal_type) <- contr.sum(levels(sens_data$mycorrhizal_type))
    colnames(contrasts(sens_data$mycorrhizal_type)) <- levels(sens_data$mycorrhizal_type)[1:(length(levels(sens_data$mycorrhizal_type))-1)]
    
    sensitivity_datasets[[paste0("sensitivity_", i)]] <- sens_data
  }
  
  return(sensitivity_datasets)
}

# Read in the tree data
tree_data <- data.table::fread("data/presence/trees_10.txt") %>%
  filter(
    mycorrhizal_type %in% c("AM", "EcM", "EcM-AM", "NM"),
    cell %in% data_est$cell
  )

# Create datasets with 10% resampling
sensitivity_tree_data <- create_sensitivity_datasets(
  tree_data, n_sensitivity = 5, resample_prop = 0.2
)

# Process sensitivity datasets to create estimation data
sensitivity_est_list <- list()
for (i in 1:5) {
  cat("Processing sensitivity dataset", i, "...\n")
  sens_raw <- sensitivity_tree_data[[paste0("sensitivity_", i)]]
  
  sens_processed <- sens_raw %>%
    group_by(cell) %>%
    summarise(
      richness = n_distinct(scientific_name, na.rm = TRUE),
      richness_AM = n_distinct(scientific_name[mycorrhizal_type == "AM"], na.rm = TRUE),
      richness_EcM = n_distinct(scientific_name[mycorrhizal_type == "EcM"], na.rm = TRUE),
      richness_EcM_AM = n_distinct(scientific_name[mycorrhizal_type == "EcM-AM"], na.rm = TRUE),
      richness_NM = n_distinct(scientific_name[mycorrhizal_type == "NM"], na.rm = TRUE),
      richness_ErM = n_distinct(scientific_name[mycorrhizal_type == "ErM"], na.rm = TRUE),
      x_albers = mean(x_albers, na.rm = TRUE),
      y_albers = mean(y_albers, na.rm = TRUE),
      RC1 = mean(RC1, na.rm = TRUE),
      RC2 = mean(RC2, na.rm = TRUE),
      RC3 = mean(RC3, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      x_km = x_albers / 1000, y_km = y_albers / 1000,
      richness_AM_log_ratio = log((richness_AM + 1) / (richness - richness_AM + 1)),
      richness_EcM_log_ratio = log((richness_EcM + 1) / (richness - richness_EcM + 1)),
      richness_EcM_AM_log_ratio = log((richness_EcM_AM + 1) / (richness - richness_EcM_AM + 1)),
      richness_NM_log_ratio = log((richness_NM + 1) / (richness - richness_NM + 1)),
      RC1 = std(RC1),
      RC2 = std(RC2),
      RC3 = std(RC3)
    )
  
  sensitivity_est_list[[i]] <- sens_processed
}

# Verify predictor alignment for first sensitivity dataset
alignment_check <- plot(
  data_est %>% arrange(cell) %>% pull(RC1),
  sensitivity_est_list[[1]] %>% arrange(cell) %>% pull(RC1),
  main = "RC1 alignment check",
  xlab = "Original data RC1",
  ylab = "Sensitivity 1 RC1"
)

# (4) Model formulas ###########################################################

# Define the formulas
formula_AM <- as.formula(paste0(
  "y_AM ~ -1 + intercept + ",
  "RC1 + RC2 + RC3 +",
  "f(spatial, model = spde)"
))
formula_EcM <- as.formula(paste0(
  "y_EcM ~ -1 + intercept + ",
  "RC1 + RC2 + RC3 +",
  "f(spatial, model = spde)"
))
formula_EcM_AM <- as.formula(paste0(
  "y_EcM_AM ~ -1 + intercept + ",
  "RC1 + RC2 + RC3 +",
  "f(spatial, model = spde)"
))
formula_NM <- as.formula(paste0(
  "y_NM ~ -1 + intercept + ",
  "RC1 + RC2 + RC3 +",
  "f(spatial, model = spde)"
))

# (5) Stack creation function ##################################################

create_data_stack <- function(data_est_input) {
  n_est_current <- nrow(data_est_input)
  
  # Covariates for all predictions
  covariates <- bind_rows(
    # Estimation and spatial prediction data
    tibble(
      RC1 = c(data_est_input$RC1, data_pred$RC1),
      RC2 = c(data_est_input$RC2, data_pred$RC2),
      RC3 = c(data_est_input$RC3, data_pred$RC3),
      tag = c(rep("est", n_est_current), rep("pred", n_pred))
    ),
    # Marginal effects prediction data
    tibble(
      RC1 = seq(min(data_est_input$RC1), max(data_est_input$RC1), length.out = 100),
      RC2 = median(data_est_input$RC2),
      RC3 = median(data_est_input$RC3),
      tag = "pred_RC1"
    ),
    tibble(
      RC1 = median(data_est_input$RC1),
      RC2 = seq(min(data_est_input$RC2), max(data_est_input$RC2), length.out = 100),
      RC3 = median(data_est_input$RC3),
      tag = "pred_RC2"
    ),
    tibble(
      RC1 = median(data_est_input$RC1),
      RC2 = median(data_est_input$RC2),
      RC3 = seq(min(data_est_input$RC3), max(data_est_input$RC3), length.out = 100),
      tag = "pred_RC3"
    )
  )
  
  # Stack the estimation data
  stack_est <- inla.stack(
    tag = "est",
    data = list(
      y_AM = data_est_input$richness_AM_log_ratio,
      y_EcM = data_est_input$richness_EcM_log_ratio,
      y_EcM_AM = data_est_input$richness_EcM_AM_log_ratio,
      y_NM = data_est_input$richness_NM_log_ratio
    ),
    A = list(1, 1, A_est),
    effects = list(
      intercept = rep(1, n_est_current),
      X = covariates %>% filter(tag == "est") %>% select(-tag),
      spatial = index_spatial
    )
  )
  
  # Stack the prediction data
  stack_pred <- inla.stack(
    tag = "pred",
    data = list(
      y_AM = rep(NA_real_, n_pred),
      y_EcM = rep(NA_real_, n_pred),
      y_EcM_AM = rep(NA_real_, n_pred),
      y_NM = rep(NA_real_, n_pred)
    ),
    A = list(1, 1, A_pred),
    effects = list(
      intercept = rep(1, n_pred),
      X = covariates %>% filter(tag == "pred") %>% select(-tag),
      spatial = index_spatial
    )
  )
  
  # Stack the data for the marginal effects
  stack_RC1 <- inla.stack(
    tag = "pred_RC1",
    data = list(
      y_AM = rep(NA_real_, 100),
      y_EcM = rep(NA_real_, 100),
      y_EcM_AM = rep(NA_real_, 100),
      y_NM = rep(NA_real_, 100)
    ),
    A = list(1, 1),
    effects = list(
      intercept = rep(1, 100),
      X = covariates %>% filter(tag == "pred_RC1") %>% select(-tag)
    )
  )
  stack_RC2 <- inla.stack(
    tag = "pred_RC2",
    data = list(
      y_AM = rep(NA_real_, 100),
      y_EcM = rep(NA_real_, 100),
      y_EcM_AM = rep(NA_real_, 100),
      y_NM = rep(NA_real_, 100)
    ),
    A = list(1, 1),
    effects = list(
      intercept = rep(1, 100),
      X = covariates %>% filter(tag == "pred_RC2") %>% select(-tag)
    )
  )
  stack_RC3 <- inla.stack(
    tag = "pred_RC3",
    data = list(
      y_AM = rep(NA_real_, 100),
      y_EcM = rep(NA_real_, 100),
      y_EcM_AM = rep(NA_real_, 100),
      y_NM = rep(NA_real_, 100)
    ),
    A = list(1, 1),
    effects = list(
      intercept = rep(1, 100),
      X = covariates %>% filter(tag == "pred_RC3") %>% select(-tag)
    )
  )
  
  # Join the stacks
  stack <- inla.stack.join(
    stack_est,
    stack_pred,
    stack_RC1,
    stack_RC2,
    stack_RC3
  )
  
  # Retrieve the stack indexes
  indices <- list(
    est = inla.stack.index(stack, "est")$data,
    pred = inla.stack.index(stack, "pred")$data,
    RC1 = inla.stack.index(stack, "pred_RC1")$data,
    RC2 = inla.stack.index(stack, "pred_RC2")$data,
    RC3 = inla.stack.index(stack, "pred_RC3")$data
  )
  
  return(list(stack = stack, indices = indices, covariates = covariates))
}

# (6) Model fitting function ###################################################

fit_models <- function(stack_data) {
  # Fit models
  models <- list()
  
  # AM
  cat("  Fitting model: AM\n")
  models$AM <- inla(
    formula_AM,
    family = 'gaussian',
    data = inla.stack.data(stack_data$stack),
    control.predictor = list(
      compute = TRUE, link = 1,
      A = inla.stack.A(stack_data$stack)
    ),
    control.compute = control_compute
  )
  
  # EcM
  cat("  Fitting model: EcM\n")
  models$EcM <- inla(
    formula_EcM,
    family = 'gaussian',
    data = inla.stack.data(stack_data$stack),
    control.predictor = list(
      compute = TRUE, link = 1,
      A = inla.stack.A(stack_data$stack)
    ),
    control.compute = control_compute
  )
  
  # EcM_AM
  cat("  Fitting model: EcM_AM\n")
  models$EcM_AM <- inla(
    formula_EcM_AM,
    family = 'gaussian',
    data = inla.stack.data(stack_data$stack),
    control.predictor = list(
      compute = TRUE, link = 1,
      A = inla.stack.A(stack_data$stack)
    ),
    control.compute = control_compute
  )
  
  # NM
  cat("  Fitting model: NM\n")
  models$NM <- inla(
    formula_NM,
    family = 'gaussian',
    data = inla.stack.data(stack_data$stack),
    control.predictor = list(
      compute = TRUE, link = 1,
      A = inla.stack.A(stack_data$stack)
    ),
    control.compute = control_compute
  )
  
  return(models)
}

# (7) Extract marginal effects function ########################################

extract_marginal_effects <- function(models, indices, covariates) {
  marginal_effects_list <- list()
  
  for (model_type in names(models)) {
    current_model <- models[[model_type]]
    
    marginal_effects_list[[model_type]] <- bind_rows(
      tibble(
        response = current_model[["summary.linear.predictor"]]$mean[indices$RC1],
        predictor = covariates %>% filter(tag == "pred_RC1") %>% pull(RC1),
        upper = current_model[["summary.linear.predictor"]]$`0.975quant`[indices$RC1],
        lower = current_model[["summary.linear.predictor"]]$`0.025quant`[indices$RC1],
        variable = "RC1",
        model = model_type
      ),
      tibble(
        response = current_model[["summary.linear.predictor"]]$mean[indices$RC2],
        predictor = covariates %>% filter(tag == "pred_RC2") %>% pull(RC2),
        upper = current_model[["summary.linear.predictor"]]$`0.975quant`[indices$RC2],
        lower = current_model[["summary.linear.predictor"]]$`0.025quant`[indices$RC2],
        variable = "RC2",
        model = model_type
      ),
      tibble(
        response = current_model[["summary.linear.predictor"]]$mean[indices$RC3],
        predictor = covariates %>% filter(tag == "pred_RC3") %>% pull(RC3),
        upper = current_model[["summary.linear.predictor"]]$`0.975quant`[indices$RC3],
        lower = current_model[["summary.linear.predictor"]]$`0.025quant`[indices$RC3],
        variable = "RC3",
        model = model_type
      )
    )
  }
  
  return(bind_rows(marginal_effects_list))
}

# (8) Main sensitivity analysis loop ###########################################

# Initialise storage for results
all_marginal_effects <- list()
all_fixed_effects <- list()

# Define analysis scenarios (main + 5 sensitivity)
analysis_scenarios <- c("Main analysis", paste0("Sensitivity ", 1:5))

# Run analysis loop
for (i in seq_along(analysis_scenarios)) {
  analysis_name <- analysis_scenarios[i]
  cat("\n=== Running", analysis_name, "===\n")
  
  # Step 1: Select appropriate data
  cat("1. Selecting data...\n")
  if (i == 1) {
    current_data <- data_est  # Main analysis uses original data
  } else {
    current_data <- sensitivity_est_list[[i-1]]  # Sensitivity analyses use resampled data
  }
  
  # Step 2: Create data stack
  cat("2. Creating data stack...\n")
  stack_data <- create_data_stack(current_data)
  
  # Step 3: Fit models
  cat("3. Fitting models...\n")
  models <- fit_models(stack_data)
  
  # Step 4: Extract marginal effects
  cat("4. Extracting marginal effects...\n")
  marginal_effects <- extract_marginal_effects(models, stack_data$indices, stack_data$covariates)
  marginal_effects$analysis <- analysis_name
  all_marginal_effects[[analysis_name]] <- marginal_effects
  
  # Step 5: Extract fixed effects
  cat("5. Extracting fixed effects...\n")
  fixed_effects_list <- list()
  for (model_type in names(models)) {
    fixed <- models[[model_type]]$summary.fixed %>%
      rownames_to_column("parameter") %>%
      select(parameter, median = mean, lower = `0.025quant`, upper = `0.975quant`) %>%
      mutate(model = model_type, analysis = analysis_name)
    fixed_effects_list[[model_type]] <- fixed
  }
  all_fixed_effects[[analysis_name]] <- bind_rows(fixed_effects_list)
  
  cat("Completed", analysis_name, "!\n")
}

# Combine all results
combined_marginal_effects <- bind_rows(all_marginal_effects)
combined_fixed_effects <- bind_rows(all_fixed_effects)

# (9) Data preparation for visualization ########################################

# Prepare observed data for all analyses
prepare_observed_data <- function(data_est_input) {
  data_est_input %>%
    mutate(
      AM = log10((richness_AM + 1) / (richness - richness_AM + 1)),
      EcM = log10((richness_EcM + 1) / (richness - richness_EcM + 1)),
      `EcM-AM` = log10((richness_EcM_AM + 1) / (richness - richness_EcM_AM + 1)),
      NM = log10((richness_NM + 1) / (richness - richness_NM + 1))
    ) %>%
    select(RC1, RC2, RC3, AM, EcM, `EcM-AM`, NM) %>%
    pivot_longer(
      cols = -c(RC1, RC2, RC3),
      names_to = "mycorrhizal_type",
      values_to = "response"
    )
}

# Prepare all analysis data
analysis_list <- list()

# Main analysis
analysis_list[["main"]] <- list(
  obs = prepare_observed_data(data_est),
  marginal = all_marginal_effects[["Main analysis"]],
  fixed = all_fixed_effects[["Main analysis"]]
)

# Sensitivity analyses
for (i in 1:5) {
  analysis_name <- paste0("sens", i)
  sensitivity_name <- paste0("Sensitivity ", i)
  
  analysis_list[[analysis_name]] <- list(
    obs = prepare_observed_data(sensitivity_est_list[[i]]),
    marginal = all_marginal_effects[[sensitivity_name]],
    fixed = all_fixed_effects[[sensitivity_name]]
  )
}

# (10) Visualization ###########################################################

# Create sensitivity analysis plots
rc_plots <- list()
for (i in 1:3) {
  rc_param <- paste0("RC", i)
  p_list <- list()
  for (anal_name in names(analysis_list)) {
    anal <- analysis_list[[anal_name]]
    
    # df_current
    df_current <- anal$obs %>%
      select(mycorrhizal_type, response, all_of(rc_param)) %>%
      rename(RC = rc_param) %>%
      mutate(
        mycorrhizal_type = if_else(mycorrhizal_type == "EcM-AM", "Dual", mycorrhizal_type),
        mycorrhizal_type = if_else(mycorrhizal_type == "EcM_AM", "Dual", mycorrhizal_type),
        mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
      )
    
    # filtered_effects
    filtered_effects <- anal$marginal %>%
      filter(variable == rc_param) %>%
      mutate(
        response = log10(exp(response)),
        lower = log10(exp(lower)),
        upper = log10(exp(upper)),
        mycorrhizal_type = model
      ) %>%
      select(variable, response, predictor, lower, upper, mycorrhizal_type) %>%
      mutate(
        mycorrhizal_type = if_else(mycorrhizal_type == "EcM-AM", "Dual", mycorrhizal_type),
        mycorrhizal_type = if_else(mycorrhizal_type == "EcM_AM", "Dual", mycorrhizal_type),
        mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")),
        anal_name = anal_name
      ) %>%
      rename(RC = predictor)
    
    # filtered_coefficients
    filtered_coefficients <- anal$fixed %>%
      filter(parameter == rc_param) %>%
      mutate(
        annotation = sprintf("β = %.3f [%.3f, %.3f]", median, lower, upper),
        x_pos = -Inf,
        y_pos = Inf,
        mycorrhizal_type = model
      ) %>%
      select(parameter, x_pos, y_pos, annotation, mycorrhizal_type) %>%
      mutate(
        mycorrhizal_type = if_else(mycorrhizal_type == "EcM-AM", "Dual", mycorrhizal_type),
        mycorrhizal_type = if_else(mycorrhizal_type == "EcM_AM", "Dual", mycorrhizal_type),
        mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
      )
    
    # Create plot
    p <- ggplot(data = df_current, aes(x = RC, y = response)) +
      geom_hex(bins = 40, aes(
        alpha = after_stat(ndensity),
        fill = after_stat(ndensity)
      )) +
      scale_alpha(range = c(0.8, 1)) +
      scale_fill_gradient(
        low = "#EEEEEE",
        high = "#212121"
      ) +
      geom_hline(yintercept = 0, linetype = 'dotted', colour = '#de2d26', linewidth = 0.5) +
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
        hjust = -0.04,
        vjust = 2.5
      ) +
      facet_wrap(~mycorrhizal_type, nrow = 1) +
      scale_y_continuous(
        limits = c(-2.25, 2),
        breaks = c(-2, -1, 0, 1, 2),
        labels = c("1:100", "1:10", "1:1", "10:1", "100:1"),
        # Add the right-hand side title here
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
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      labs(
        x = case_when(
          rc_param == "RC1" ~ "RC1 (Temperature & decomposition)",
          rc_param == "RC2" ~ "RC2 (Soil moisture)",
          rc_param == "RC3" ~ "RC3 (Soil phosphorus)",
          TRUE ~ rc_param
        ),
        y = case_when(
          anal_name == "sens3" ~ "Relative Richness (log-scale)",
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
    
    p_list[[anal_name]] <- p
  }
  
  # Wrap into one plot for this RC
  rc_plots[[rc_param]] <- wrap_plots(p_list, ncol = 1)
}


# Display the plots
print(rc_plots[["RC1"]])
print(rc_plots[["RC2"]])  
print(rc_plots[["RC3"]])

# Save plots
ggsave("output/figure_S14.png", rc_plots[["RC1"]],
       width = 14, height = 20, units = "cm", dpi = 300)
ggsave("output/figure_S15.png", rc_plots[["RC2"]],
       width = 14, height = 20, units = "cm", dpi = 300)
ggsave("output/figure_S16.png", rc_plots[["RC3"]],
       width = 14, height = 20, units = "cm", dpi = 300)

