
# Required packages
library(boot)
library(emmeans)
library(ape)
library(adephylo)
library(MPSEM)
library(V.PhyloMaker2)
library(tidyverse)
source("code/fast.phylo.eigenvector.select.R")

# (1) Data preparation ########################################################

# Define the response variables for the models
response_vars <- c(
  "env_breadth"
)

# Read in the data
data <- data.table::fread(
  "data/niche_estimates_enmeval/niche_estimates.txt"
) %>%
  filter(mycorrhizal_type != "ErM") %>%
  mutate(
    # Sqrt-root transformation to normalise env_breadth
    env_breadth = sqrt(env_B2_corrected),
    # Rename EcM-AM to Dual
    mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
    mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")),
    species = str_replace_all(species, " ", "_"),
    biome = ifelse(biome == "Tropical", "tropical", "nontropical"),
    biome = factor(
      biome,
      levels = c("tropical", "nontropical")
    )
  ) %>%
  select(
    family, genus, species, mycorrhizal_type, biome, env_breadth
  )

# Set effect coding for mycorrhizal_type
contrasts(data$mycorrhizal_type) <- contr.sum(
  levels(data$mycorrhizal_type)
)

# Rename the headers
colnames(contrasts(data$mycorrhizal_type)) <- levels(data$mycorrhizal_type)[1:(length(levels(data$mycorrhizal_type))-1)]

# Split the data into tropical and non-tropical
data_tropical <- data %>%
  filter(biome == "tropical")
data_nontropical <- data %>%
  filter(biome == "nontropical")

# Check the contrasts
contrasts(data$mycorrhizal_type)
contrasts(data_tropical$mycorrhizal_type)
contrasts(data_nontropical$mycorrhizal_type)

# Number of species: 2,334 species
unique(data$species) %>% length(.)
# Number of tropical species: 1,224 species
unique(data_tropical$species) %>% length(.)
# Number of nontropical species: 1,110 species
unique(data_nontropical$species) %>% length(.)

#### * Phylogenetic tree * ####

# Add species as row names
rownames(data) <- data[["species"]]
rownames(data_tropical) <- data_tropical[["species"]]
rownames(data_nontropical) <- data_nontropical[["species"]]

# Create a species list
species_list <- data %>%
  dplyr::select(species, genus, family) %>%
  # Change family name to suit the phylo.maker backbone dataset
  mutate(
    family = case_when(
      family == "Viburnaceae" ~ "Adoxaceae",
      TRUE ~ family
    )
  ) %>%
  as.data.frame(.)

# Read in the phylogenetic tree
phylo_tree <- read.tree(
  "generated_data/phylo_tree_mycorrhizal_types.tre"
) 

phylo_tree <- phylo_tree %>%
  # Prune the tree to species list
  drop.tip(
    ., 
    setdiff(phylo_tree$tip.label, species_list$species)
  )

# Check the tips are ordered to match the data
identical(phylo_tree[["tip.label"]], rownames(data))

# Match the phylogenetic tree tips with the data
data <- data[match(phylo_tree[["tip.label"]], data[["species"]]), ]
rownames(data) <- data[["species"]]

# Check the tips are ordered to match the data
identical(phylo_tree[["tip.label"]], rownames(data))

# Prune the tree to nontropical, tropical and mediterranean trees
phylo_tree_tropical <- ape::drop.tip(
  phylo_tree, 
  setdiff(phylo_tree$tip.label, data_tropical$species)
)
phylo_tree_nontropical <- ape::drop.tip(
  phylo_tree, 
  setdiff(phylo_tree$tip.label, data_nontropical$species)
)

# Align the tips and data
data_tropical <- data_tropical[match(phylo_tree_tropical$tip.label, data_tropical$species), ]
rownames(data_tropical) <- data_tropical[["species"]]
data_nontropical <- data_nontropical[match(phylo_tree_nontropical$tip.label, data_nontropical$species), ]
rownames(data_nontropical) <- data_nontropical[["species"]]

#### * Phylogenetic data * ####

# Generate proximity matrices:
phylo_prox <- proxTips(phylo_tree, method = "oriAbouheif", normalize = "row")
phylo_prox_tropical <- proxTips(phylo_tree_tropical, method = "oriAbouheif", normalize = "row")
phylo_prox_nontropical <- proxTips(phylo_tree_nontropical, method = "oriAbouheif", normalize = "row")

# Fix the proximity matrix by setting tiny negative values to zero:
# Count total negative values
sum(phylo_prox < 0)
sum(phylo_prox_tropical < 0)
sum(phylo_prox_nontropical < 0)
# Get the minimum value to see how negative it gets
min(phylo_prox) # <- round to 0: numeric precision errors, not true negative values
min(phylo_prox_tropical) # <- round to 0: numeric precision errors, not true negative values
min(phylo_prox_nontropical) # <- round to 0: numeric precision errors, not true negative values
# Set negative values to zero
phylo_prox[phylo_prox < 0] <- 0
phylo_prox_tropical[phylo_prox_tropical < 0] <- 0
phylo_prox_nontropical[phylo_prox_nontropical < 0] <- 0

# Verify the fix
any(phylo_prox < 0)  # Should be FALSE now
any(phylo_prox_tropical < 0)  # Should be FALSE now
any(phylo_prox_nontropical < 0)  # Should be FALSE now
min(phylo_prox)      # Should now be 0
min(phylo_prox_tropical) # Should now be 0
min(phylo_prox_nontropical) # Should now be 0

# Generate eigenvector matrices:
pem_obj <- PEM.build(Phylo2DirectedGraph(phylo_tree))
pem <- as.data.frame(pem_obj)
colnames(pem) <- paste0("PEM", seq_len(ncol(pem)))
attr(pem, "values") <- pem_obj$d
names(attr(pem, "values")) <- paste("PEM", seq_len(ncol(pem)))

pem_tropical_obj <- PEM.build(Phylo2DirectedGraph(phylo_tree_tropical), a = 1)
pem_tropical <- as.data.frame(pem_tropical_obj)
colnames(pem_tropical) <- paste0("PEM", seq_len(ncol(pem_tropical)))
attr(pem_tropical, "values") <- pem_tropical_obj$d
names(attr(pem_tropical, "values")) <- paste("PEM", seq_len(ncol(pem_tropical)))

pem_nontropical_obj <- PEM.build(Phylo2DirectedGraph(phylo_tree_nontropical))
pem_nontropical <- as.data.frame(pem_nontropical_obj)
colnames(pem_nontropical) <- paste0("PEM", seq_len(ncol(pem_nontropical)))
attr(pem_nontropical, "values") <- pem_nontropical_obj$d
names(attr(pem_nontropical, "values")) <- paste("PEM", seq_len(ncol(pem_nontropical)))

# Verify the attributes are set correctly
head(attr(pem, "values"))
head(attr(pem_tropical, "values"))
head(attr(pem_nontropical, "values"))

# Ensure that data and phylogenetic structures are aligned by species
stopifnot(
  identical(rownames(data), rownames(pem)),
  identical(rownames(data), rownames(phylo_prox)),
  identical(rownames(data_tropical), rownames(pem_tropical)),
  identical(rownames(data_tropical), rownames(phylo_prox_tropical)),
  identical(rownames(data_nontropical), rownames(pem_nontropical)),
  identical(rownames(data_nontropical), rownames(phylo_prox_nontropical))
)

#### * Setup data lists * ####

# Data subsets
data_subsets <- list(
  all = data,
  tropical = data_tropical,
  nontropical = data_nontropical
)

# Tree proximity subsets
phylo_prox_list <- list(
  all = phylo_prox,
  tropical = phylo_prox_tropical,
  nontropical = phylo_prox_nontropical
)

# PEM subsets
pem_list <- list(
  all = pem,
  tropical = pem_tropical,
  nontropical = pem_nontropical
)

# (2) Base models ############################################################

#### * Fit base models * ####
base_models <- list()

for (response in response_vars) {
  base_models[[response]] <- list()
  
  for (subset_name in names(data_subsets)) {
    model_name <- paste0("base_model_", response, 
                         if (subset_name != "all") paste0("_", subset_name) else "")
    
    base_models[[response]][[subset_name]] <- lm(
      as.formula(paste(response, "~ mycorrhizal_type")),
      data = data_subsets[[subset_name]]
    )
  }
}

#### * Extract residuals * ####

residuals_list <- list()

for (response in response_vars) {
  residuals_list[[response]] <- list()
  
  for (subset_name in names(data_subsets)) {
    # Extract residuals
    residuals_list[[response]][[subset_name]] <- 
      residuals(base_models[[response]][[subset_name]])
    
  }
}

# (3) Phylogenetic models ######################################################

#### * MEM selection * ####

# Set the seed
set.seed(1986)

pem_select <- list()

for (response in response_vars) {
  pem_select[[response]] <- list()
  
  for (subset_name in names(data_subsets)) {
    cat("\n=== Running MIR selection for", response, subset_name, "===\n")
    
    pem_select[[response]][[subset_name]] <- fast.phylo.eigenvector.select(
      x = residuals_list[[response]][[subset_name]],
      phylo_eigenvectors = pem_list[[subset_name]],
      phylo_matrix = phylo_prox_list[[subset_name]],
      method = "abouheif",
      abouheif_method = "oriAbouheif",
      alpha = 0.05,
      verbose = TRUE
    )
  }
}

#### * Phylogenetic models * ####

pem_models <- list()
formulas <- list()
for (response in response_vars) {
  pem_models[[response]] <- list()
  formulas[[response]] <- list()
  
  for (subset_name in names(data_subsets)) {
    # Create formula
    pem_vars <- colnames(pem_select[[response]][[subset_name]]$PE.select)
    if (length(pem_vars) > 0) {
      formulas[[response]][[subset_name]] <- as.formula(paste(
        response, "~ mycorrhizal_type +", 
        paste(pem_vars, collapse = " + ")
      ))
    } else {
      formulas[[response]][[subset_name]] <- as.formula(paste(
        response, "~ mycorrhizal_type"
      ))
    }
    
    # Fit model
    pem_models[[response]][[subset_name]] <- lm(
      formulas[[response]][[subset_name]],
      data = bind_cols(data_subsets[[subset_name]], 
                       pem_select[[response]][[subset_name]]$PE.select)
    )
  }
}

# Get phylogenetic autocorrelation stats
all_n <- abouheif.moran(residuals_list[["env_breadth"]][["all"]], phylo_prox_list[["all"]], method = "oriAbouheif")
all_p <- abouheif.moran(residuals(pem_models[["env_breadth"]][["all"]]), phylo_prox_list[["all"]], method = "oriAbouheif")
tropical_n <- abouheif.moran(residuals_list[["env_breadth"]][["tropical"]], phylo_prox_list[["tropical"]], method = "oriAbouheif")
tropical_p <- abouheif.moran(residuals(pem_models[["env_breadth"]][["tropical"]]), phylo_prox_list[["tropical"]], method = "oriAbouheif")
nontropical_n <- abouheif.moran(residuals_list[["env_breadth"]][["nontropical"]], phylo_prox_list[["nontropical"]], method = "oriAbouheif")
nontropical_p <- abouheif.moran(residuals(pem_models[["env_breadth"]][["nontropical"]]), phylo_prox_list[["nontropical"]], method = "oriAbouheif")

phylo_auto_stats <- tibble(
  model = c(
    rep(c("naïve", "phylogenetic"), 3) 
  ),
  group = c(
    "all", "all", 
    "tropical", "tropical", 
    "nontropical",  "nontropical"
  ),
  cmean = c(
    all_n$obs, all_p$obs, 
    tropical_n$obs, tropical_p$obs, 
    nontropical_n$obs, nontropical_p$obs
  ),
  pvalue = c(
    all_n$pvalue, all_p$pvalue, 
    tropical_n$pvalue, tropical_p$pvalue, 
    nontropical_n$pvalue, nontropical_p$pvalue
  )
) %>%
  print()

#### * Model parameters * ####

# Create lists to store model parameters
base_model_params <- list()
pem_model_params <- list()

for (response in response_vars) {
  base_model_params[[response]] <- list()
  pem_model_params[[response]] <- list()
  
  for (subset_name in names(data_subsets)) {
    cat("\n=== Results for", response, subset_name, "===\n")
    
    # Extract and store base model parameters
    base_params <- parameters::parameters(base_models[[response]][[subset_name]])
    base_model_params[[response]][[subset_name]] <- base_params
    
    # Extract and store ME model parameters
    pem_params <- parameters::parameters(pem_models[[response]][[subset_name]])
    pem_model_params[[response]][[subset_name]] <- pem_params
    
    # Model comparisons
    cat("\nBase model:\n")
    print(base_params)
    
    cat("\nPhylogenetic model:\n")
    print(pem_params)
  }
}

#### * Bootstrap effects * ####

# Bootstrapping function
bootstrap_effects <- function(data, indices, pems, formula, response_var) {
  d <- data[indices, ]
  pems_subset <- pems[indices, , drop = FALSE]
  d_full <- bind_cols(d, pems_subset)
  
  model <- tryCatch({
    lm(formula, data = d_full)
  }, error = function(e) return(NULL))
  
  if (is.null(model)) return(rep(NA, 3 * length(unique(data$mycorrhizal_type)) + 1))
  
  # Compute estimated marginal means
  emm <- tryCatch({
    suppressMessages(emmeans(model, "mycorrhizal_type", type = "response", data = d_full) %>% as_tibble())
  }, error = function(e) return(NULL))
  
  if (is.null(emm)) return(rep(NA, 3 * length(unique(data$mycorrhizal_type)) + 1))
  
  # Get emmean values and SE values
  emmean_values <- emm$emmean
  se_values <- emm$SE
  
  # Calculate sample sizes for each mycorrhizal type in this bootstrap sample
  n_values <- d %>%
    count(mycorrhizal_type) %>%
    arrange(mycorrhizal_type) %>%
    pull(n)
  
  # Convert SE to SD: SD = SE * sqrt(n)
  sd_values <- se_values * sqrt(n_values)
  
  # Extract coefficients for mycorrhizal types
  coef_names <- names(coef(model))
  myc_coefs <- coef(model)[grep("^mycorrhizal_type", coef_names)]
  
  # Get the unique mycorrhizal types from the data
  myc_levels <- levels(data$mycorrhizal_type)
  n_levels <- length(myc_levels)
  
  # Initialise effect vectors
  effects <- rep(NA, n_levels)
  names(effects) <- myc_levels
  
  # For effects coding, compute effects for each level
  if (length(myc_coefs) == (n_levels - 1)) {
    # Extract effects for explicitly coded levels
    for (i in 1:(n_levels - 1)) {
      level_name <- myc_levels[i]
      coef_name <- paste0("mycorrhizal_type", level_name)
      if (coef_name %in% names(myc_coefs)) {
        effects[level_name] <- myc_coefs[coef_name]
      }
    }
    
    # Compute effect for the reference level (NM) - sum of all other effects with opposite sign
    ref_level <- myc_levels[n_levels]  # Assuming NM is the last level
    effects[ref_level] <- -sum(effects[1:(n_levels-1)], na.rm = TRUE)
  } else {
    return(rep(NA, 3 * n_levels + 1))
  }
  
  # Extract intercept
  intercept <- coef(model)[1]
  
  # Return emmeans, SDs, coefficient-based effects, and intercept
  return(c(emmean_values, sd_values, effects, intercept))
}

# Run bootstrap for each model
boot_results <- list()

for (response in response_vars) {
  boot_results[[response]] <- list()
  
  for (subset_name in names(data_subsets)) {
    cat("\n=== Running bootstrap for", response, subset_name, "===\n")
    
    # Prepare fomulas, data and pems for this model
    current_data <- data_subsets[[subset_name]]
    current_pems <- pem_select[[response]][[subset_name]]$PE.select
    current_formula <- formulas[[response]][[subset_name]]
    
    boot_results[[response]][[subset_name]] <- boot::boot(
      data = current_data,
      strata = current_data$mycorrhizal_type,
      parallel = "multicore",
      ncpus = 10,
      statistic = function(data, indices) {
        bootstrap_effects(
          data = data,
          indices = indices,
          pems = current_pems,
          formula = current_formula,
          response_var = response
        )
      },
      R = 10000
    )
  }
}

# (4) Summarise results ########################################################

#### * Process results * ####
# Function to process bootstrap results
process_boot_results <- function(boot_result, data_subset) {
  n_levels <- length(levels(data_subset$mycorrhizal_type))
  level_names <- levels(data_subset$mycorrhizal_type)
  
  # Split the bootstrap results into mean, SD, effect, and intercept parts
  mean_boot <- boot_result$t[, 1:n_levels, drop = FALSE]
  sd_boot <- boot_result$t[, (n_levels + 1):(2 * n_levels), drop = FALSE]
  effect_boot <- boot_result$t[, (2 * n_levels + 1):(3 * n_levels), drop = FALSE]
  intercept_boot <- boot_result$t[, 3 * n_levels + 1]
  
  # Process intercept results
  intercept_summary <- tibble(
    mean = median(intercept_boot, na.rm = TRUE),
    lower = quantile(intercept_boot, 0.025, na.rm = TRUE),
    upper = quantile(intercept_boot, 0.975, na.rm = TRUE)
  ) %>%
    mutate(
      intercept_formatted = sprintf("%.3f [%.3f, %.3f]", mean, lower, upper)
    )
  
  # Process mean results
  long_boot_mean <- as_tibble(mean_boot, .name_repair = "minimal") %>%
    setNames(level_names) %>%
    mutate(replicate = row_number()) %>%
    pivot_longer(-replicate, names_to = "mycorrhizal_type", values_to = "emmean") %>%
    drop_na(emmean)
  
  # Process SD results
  long_boot_sd <- as_tibble(sd_boot, .name_repair = "minimal") %>%
    setNames(level_names) %>%
    mutate(replicate = row_number()) %>%
    pivot_longer(-replicate, names_to = "mycorrhizal_type", values_to = "sd") %>%
    drop_na(sd)
  
  # Density data for mean
  density_data_mean <- long_boot_mean %>%
    group_by(mycorrhizal_type) %>%
    mutate(
      lower_q = quantile(emmean, 0.005, na.rm = TRUE),
      upper_q = quantile(emmean, 0.995, na.rm = TRUE)
    ) %>%
    filter(emmean >= lower_q, emmean <= upper_q) %>%
    summarise(
      x = list(density(emmean, adjust = 3)$x),
      y = list(density(emmean, adjust = 3)$y),
      .groups = "drop"
    ) %>%
    unnest(c(x, y)) %>%
    mutate(
      mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
      mycorrhizal_type = factor(mycorrhizal_type, levels = rev(c("AM", "EcM", "Dual", "NM")))
    )
  
  # Density data for SD
  density_data_sd <- long_boot_sd %>%
    group_by(mycorrhizal_type) %>%
    mutate(
      lower_q = quantile(sd, 0.005, na.rm = TRUE),
      upper_q = quantile(sd, 0.995, na.rm = TRUE)
    ) %>%
    filter(sd >= lower_q, sd <= upper_q) %>%
    summarise(
      x = list(density(sd, adjust = 3)$x),
      y = list(density(sd, adjust = 3)$y),
      .groups = "drop"
    ) %>%
    unnest(c(x, y)) %>%
    mutate(
      mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
      mycorrhizal_type = factor(mycorrhizal_type, levels = rev(c("AM", "EcM", "Dual", "NM")))
    )
  
  # Summary data for mean (including mean SD)
  summary_data_mean <- long_boot_mean %>%
    left_join(long_boot_sd, by = c("replicate", "mycorrhizal_type")) %>%
    group_by(mycorrhizal_type) %>%
    summarise(
      mean = median(emmean, na.rm = TRUE),
      sd = median(sd, na.rm = TRUE),
      lower = quantile(emmean, 0.025, na.rm = TRUE),
      upper = quantile(emmean, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
      mycorrhizal_type = factor(mycorrhizal_type, levels = rev(c("AM", "EcM", "Dual", "NM")))
    )
  
  # Summary data for SD
  summary_data_sd <- long_boot_sd %>%
    group_by(mycorrhizal_type) %>%
    summarise(
      mean = median(sd, na.rm = TRUE),
      lower = quantile(sd, 0.025, na.rm = TRUE),
      upper = quantile(sd, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
      mycorrhizal_type = factor(mycorrhizal_type, levels = rev(c("AM", "EcM", "Dual", "NM")))
    )
  
  # Process effect results
  long_boot_effect <- as_tibble(effect_boot, .name_repair = "minimal") %>%
    setNames(level_names) %>%
    mutate(replicate = row_number()) %>%
    pivot_longer(-replicate, names_to = "mycorrhizal_type", values_to = "effect") %>%
    drop_na(effect)
  
  # Density data for effect
  density_data_effect <- long_boot_effect %>%
    group_by(mycorrhizal_type) %>%
    mutate(
      lower_q = quantile(effect, 0.005, na.rm = TRUE),
      upper_q = quantile(effect, 0.995, na.rm = TRUE)
    ) %>%
    filter(effect >= lower_q, effect <= upper_q) %>%
    summarise(
      x = list(density(effect, adjust = 3)$x),
      y = list(density(effect, adjust = 3)$y),
      .groups = "drop"
    ) %>%
    unnest(c(x, y)) %>%
    mutate(
      mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
      mycorrhizal_type = factor(mycorrhizal_type, levels = rev(c("AM", "EcM", "Dual", "NM")))
    )
  
  # Summary data for effect
  summary_data_effect <- long_boot_effect %>%
    group_by(mycorrhizal_type) %>%
    summarise(
      mean = median(effect, na.rm = TRUE),
      lower = quantile(effect, 0.025, na.rm = TRUE),
      upper = quantile(effect, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
      mycorrhizal_type = factor(mycorrhizal_type, levels = rev(c("AM", "EcM", "Dual", "NM")))
    )
  
  return(list(
    mean = list(density = density_data_mean, summary = summary_data_mean),
    effect = list(density = density_data_effect, summary = summary_data_effect),
    intercept = intercept_summary
  ))
}

# Process results for each model
processed_results <- list()

for (response in response_vars) {
  processed_results[[response]] <- list()
  
  for (subset_name in names(data_subsets)) {
    processed_results[[response]][[subset_name]] <- 
      process_boot_results(boot_results[[response]][[subset_name]], 
                           data_subsets[[subset_name]])
  }
}

#### * Compute pairwise differences * ####

# Function to compute pairwise differences from marginal means
compute_pairwise_differences <- function(summary_data) {
  # Define the comparisons we want
  comparisons <- list(
    "Dual - AM" = c("Dual", "AM"),
    "Dual - EcM" = c("Dual", "EcM"), 
    "Dual - NM" = c("Dual", "NM"),
    "NM - AM" = c("NM", "AM"),
    "NM - EcM" = c("NM", "EcM"),
    "AM - EcM" = c("AM", "EcM")
  )
  
  # Initialise results list
  pairwise_results <- list()
  
  for (comp_name in names(comparisons)) {
    comp_pair <- comparisons[[comp_name]]
    reference <- comp_pair[2]  # Second element is reference (subtracted from)
    comparison <- comp_pair[1]  # First element is comparison (subtractor)
    
    # Get data for reference and comparison groups
    ref_data <- summary_data %>% filter(mycorrhizal_type == reference)
    comp_data <- summary_data %>% filter(mycorrhizal_type == comparison)
    
    # Check if both groups exist
    if (nrow(ref_data) == 0 || nrow(comp_data) == 0) {
      warning(paste("Missing data for comparison:", comp_name))
      next
    }
    
    # Calculate differences (comparison - reference)
    mean_diff <- comp_data$mean - ref_data$mean
    lower_diff <- comp_data$lower - ref_data$upper  # Conservative: comp lower - ref upper
    upper_diff <- comp_data$upper - ref_data$lower  # Conservative: comp upper - ref lower
    
    # Calculate percent change: ((comparison - reference) / reference) * 100
    percent_change <- (mean_diff / ref_data$mean) * 100
    percent_change_lower <- (lower_diff / ref_data$mean) * 100
    percent_change_upper <- (upper_diff / ref_data$mean) * 100
    
    # Store results
    pairwise_results[[comp_name]] <- tibble(
      comparison = comp_name,
      reference_group = reference,
      comparison_group = comparison,
      reference_mean = 0,  # Centered on 0 as requested
      reference_lower = ref_data$lower - ref_data$mean,  # Relative to mean
      reference_upper = ref_data$upper - ref_data$mean,  # Relative to mean
      difference_mean = mean_diff,
      difference_lower = lower_diff,
      difference_upper = upper_diff,
      percent_change_mean = percent_change,
      percent_change_lower = percent_change_lower,
      percent_change_upper = percent_change_upper
    )
  }
  
  # Combine all results
  result_table <- bind_rows(pairwise_results) %>%
    mutate(
      # Format the results nicely
      reference_formatted = sprintf("0.000 [%.3f, %.3f]", reference_lower, reference_upper),
      difference_formatted = sprintf("%.3f [%.3f, %.3f]", difference_mean, difference_lower, difference_upper),
      percent_change_formatted = sprintf("%.1f%% [%.1f%%, %.1f%%]", 
                                         percent_change_mean, percent_change_lower, percent_change_upper)
    )
  
  return(result_table)
}

# Compute the pairwise differences for each environment breadth result
env_breadth_all_pairwise <- compute_pairwise_differences(processed_results[["env_breadth"]][["all"]]$mean$summary)
env_breadth_tropical_pairwise <- compute_pairwise_differences(processed_results[["env_breadth"]][["tropical"]]$mean$summary)
env_breadth_nontropical_pairwise <- compute_pairwise_differences(processed_results[["env_breadth"]][["nontropical"]]$mean$summary)

# Create a tibble to store all pairwise results
pairwise_results <- bind_rows(
  env_breadth_all_pairwise %>%
    mutate(model = "All") %>%
    select(
      comparison, reference_group = reference_group, 
      reference_level = reference_formatted,
      difference = difference_formatted,
      percent_change = percent_change_formatted
    ),
  env_breadth_tropical_pairwise %>%
    mutate(model = "Tropical") %>%
    select(
      comparison, reference_group = reference_group, 
      reference_level = reference_formatted,
      difference = difference_formatted,
      percent_change = percent_change_formatted
    )
  ) %>%
  bind_rows(
    env_breadth_nontropical_pairwise %>%
      mutate(model = "Temperate") %>%
      select(
        comparison, reference_group = reference_group, 
        reference_level = reference_formatted,
        difference = difference_formatted,
        percent_change = percent_change_formatted
      )
  )

#### * Extract results * ####

# Organise breadth results
env_breadth_results <- processed_results[["env_breadth"]]
data_breadth <- data 
data_breadth_tropical <- data_tropical
data_breadth_nontropical <- data_nontropical
pem_breadth <- pem
pem_breadth_tropical <- pem_tropical
pem_breadth_nontropical <- pem_nontropical
phylo_prox_breadth <- phylo_prox
phylo_prox_breadth_tropical <- phylo_prox_tropical
phylo_prox_breadth_nontropical <- phylo_prox_nontropical
pem_select_breadth <- pem_select[["env_breadth"]]
base_model_params_breadth <- base_model_params[["env_breadth"]]
pem_model_params_breadth <- pem_model_params[["env_breadth"]]
pairwise_results_breadth <- pairwise_results
phylo_auto_stats_breadth <- phylo_auto_stats

#### * Save generated data * ####
save(
  # Raw data
  data_breadth, data_breadth_tropical, data_breadth_nontropical,
  
  # Phylogenetic data
  pem_breadth, pem_breadth_tropical, pem_breadth_nontropical,
  phylo_prox_breadth, phylo_prox_breadth_tropical, phylo_prox_breadth_nontropical,
  pem_select_breadth,
  
  # Results
  env_breadth_results,
  base_model_params_breadth,
  pem_model_params_breadth,
  pairwise_results_breadth,
  
  file = "generated_data/figure_5.RData"
)

# Save Table S3: Phylogenetic autocorrelation statistics for all models
phylo_auto_stats_breadth %>%
  as_tibble() %>%
  fwrite("output/table_S3.txt", sep = "\t")

# Save Table S4: Coefficients for mycorrhizal type effects on environmental niche breadth for all trees
pem_model_params_breadth[["all"]] %>%
  as_tibble() %>%
  fwrite("output/table_S4.txt", sep = "\t")

# Save Table S5: Coefficients for mycorrhizal type effects on environmental breadth for tropical trees
pem_model_params_breadth[["tropical"]] %>%
  as_tibble() %>%
  fwrite("output/table_S5.txt", sep = "\t")

# Save Table S6: Coefficients for mycorrhizal type effects on environmental breadth for temperate trees
pem_model_params_breadth[["nontropical"]] %>%
  as_tibble() %>%
  fwrite("output/table_S6.txt", sep = "\t")

# Save Table S7: Pairwise comparisons of marginal mean environmental niche breadths. 
pairwise_results_breadth %>%
  fwrite("output/table_S7.txt", sep = "\t")
