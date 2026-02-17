
# Required packages
library(phyloMEM)
library(boot)
library(emmeans)
library(ape)
library(adephylo)
library(MPSEM)
library(V.PhyloMaker2)
library(tidyverse)
source("code/fast.phylo.eigenvector.select.R")

# (1) Data preparation ########################################################

# Read in the data
data <- data.table::fread("data/niche_estimates_enmeval/niche_estimates.txt") %>%
  filter(mycorrhizal_type != "ErM") %>%
  mutate(
    # Rename EcM-AM to Dual
    mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
    mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")),
    species = str_replace_all(species, " ", "_"),
    # Scale the RC variables for the full dataset
    RC1 = as.numeric(scale(RC1_position)),
    RC2 = as.numeric(scale(RC2_position)),
    RC3 = as.numeric(scale(RC3_position)),
    biome = ifelse(biome == "Tropical", "tropical", "nontropical"),
    biome = factor(
      biome,
      levels = c("tropical", "nontropical")
    )
  ) %>%
  select(
    family, genus, species, mycorrhizal_type, biome,
    RC1, RC2, RC3
  )

# Set effect coding for mycorrhizal_type
contrasts(data$mycorrhizal_type) <- contr.sum(
  levels(data$mycorrhizal_type)
)

# Rename the headers
colnames(contrasts(data$mycorrhizal_type)) <- levels(data$mycorrhizal_type)[1:(length(levels(data$mycorrhizal_type))-1)]

# Split the data into tropical and non-tropical and re-standardize RC variables
data_tropical <- data %>%
  filter(biome == "tropical")

data_nontropical <- data %>%
  filter(biome == "nontropical")

# Check the contrasts
contrasts(data$mycorrhizal_type)
contrasts(data_tropical$mycorrhizal_type)
contrasts(data_nontropical$mycorrhizal_type)

# Number of species
cat("Total species:", unique(data$species) %>% length(.), "\n")
cat("Tropical species:", unique(data_tropical$species) %>% length(.), "\n")
cat("Non-tropical species:", unique(data_nontropical$species) %>% length(.), "\n")

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

# Prune the tree to nontropical, tropical trees
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
min(phylo_prox) # <- these are numerical precision errors rather than true negative values, so round to 0
min(phylo_prox_nontropical) # <- these are numerical precision errors rather than true negative values, so round to 0
# Set negative values to zero
phylo_prox[phylo_prox < 0] <- 0
phylo_prox_nontropical[phylo_prox_nontropical < 0] <- 0

# Generate eigenvector matrices: Set branch lengths "a" to 1 to focus on topology
pem_obj <- PEM.build(Phylo2DirectedGraph(phylo_tree))
pem <- as.data.frame(pem_obj)
colnames(pem) <- paste0("PEM", seq_len(ncol(pem)))
attr(pem, "values") <- pem_obj$d
names(attr(pem, "values")) <- paste("PEM", seq_len(ncol(pem)))

pem_tropical_obj <- PEM.build(Phylo2DirectedGraph(phylo_tree_tropical))
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

data_subsets <- list(
  all = data,
  tropical = data_tropical,
  nontropical = data_nontropical
)

pem_list <- list(
  all = pem,
  tropical = pem_tropical,
  nontropical = pem_nontropical
)

phylo_prox_list <- list(
  all = phylo_prox,
  tropical = phylo_prox_tropical,
  nontropical = phylo_prox_nontropical
)

# (2) Base models ############################################################

# Updated response variables to include RC1, RC2, RC3
response_vars <- c("RC1", "RC2", "RC3")

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

#### * PEM selection * ####

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

# Get phylogenetic autocorrelation stats for all response variables
phylo_auto_stats <- list()

for (response in response_vars) {
  all_n <- abouheif.moran(residuals_list[[response]][["all"]], phylo_prox_list[["all"]], method = "oriAbouheif")
  all_p <- abouheif.moran(residuals(pem_models[[response]][["all"]]), phylo_prox_list[["all"]], method = "oriAbouheif")
  tropical_n <- abouheif.moran(residuals_list[[response]][["tropical"]], phylo_prox_list[["tropical"]], method = "oriAbouheif")
  tropical_p <- abouheif.moran(residuals(pem_models[[response]][["tropical"]]), phylo_prox_list[["tropical"]], method = "oriAbouheif")
  nontropical_n <- abouheif.moran(residuals_list[[response]][["nontropical"]], phylo_prox_list[["nontropical"]], method = "oriAbouheif")
  nontropical_p <- abouheif.moran(residuals(pem_models[[response]][["nontropical"]]), phylo_prox_list[["nontropical"]], method = "oriAbouheif")
  
  phylo_auto_stats[[response]] <- tibble(
    response = response,
    model = c(rep(c("naïve", "phylogenetic"), 3)),
    group = c("all", "all", "tropical", "tropical", "nontropical", "nontropical"),
    cmean = c(all_n$obs, all_p$obs, tropical_n$obs, tropical_p$obs, nontropical_n$obs, nontropical_p$obs),
    pvalue = c(all_n$pvalue, all_p$pvalue, tropical_n$pvalue, tropical_p$pvalue, nontropical_n$pvalue, nontropical_p$pvalue)
  )
}

# Combine all phylo auto stats
phylo_auto_stats_combined <- bind_rows(phylo_auto_stats) %>%
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

# Bootstrapping function (same as original)
bootstrap_effects <- function(data, indices, pems, formula, response_var) {
  d <- data[indices, ]
  pems_subset <- pems[indices, , drop = FALSE]
  d_full <- bind_cols(d, pems_subset)
  
  model <- tryCatch({
    lm(formula, data = d_full)
  }, error = function(e) return(NULL))
  
  if (is.null(model)) return(rep(NA, 2 * length(unique(data$mycorrhizal_type)) + 1))
  
  # Compute estimated marginal means
  emm <- tryCatch({
    suppressMessages(emmeans(model, "mycorrhizal_type", type = "response", data = d_full) %>% as_tibble())
  }, error = function(e) return(NULL))
  
  if (is.null(emm)) return(rep(NA, 2 * length(unique(data$mycorrhizal_type)) + 1))
  
  # Get emmean values
  emmean_values <- emm$emmean
  
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
    return(rep(NA, 2 * n_levels + 1))
  }
  
  # Extract intercept
  intercept <- coef(model)[1]
  
  # Return emmeans, coefficient-based effects, and intercept
  return(c(emmean_values, effects, intercept))
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

# (4) Process results ##########################################################

#### * Process results * ####
# Function to process bootstrap results (same as original)
process_boot_results <- function(boot_result, data_subset) {
  n_levels <- length(levels(data_subset$mycorrhizal_type))
  level_names <- levels(data_subset$mycorrhizal_type)
  
  # Split the bootstrap results into mean, effect, and intercept parts
  mean_boot <- boot_result$t[, 1:n_levels, drop = FALSE]
  effect_boot <- boot_result$t[, (n_levels + 1):(2 * n_levels), drop = FALSE]
  intercept_boot <- boot_result$t[, 2 * n_levels + 1]
  
  # Process intercept results
  intercept_summary <- tibble(
    mean = median(intercept_boot, na.rm = TRUE),
    lower = quantile(intercept_boot, 0.025, na.rm = TRUE),
    upper = quantile(intercept_boot, 0.975, na.rm = TRUE)
  ) %>%
    mutate(
      intercept_formatted = sprintf("%.3f [%.3f, %.3f]", mean, lower, upper)
    )
  
  long_boot_mean <- as_tibble(mean_boot, .name_repair = "minimal") %>%
    setNames(level_names) %>%
    mutate(replicate = row_number()) %>%
    pivot_longer(-replicate, names_to = "mycorrhizal_type", values_to = "emmean") %>%
    drop_na(emmean)
  
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
  
  # Summary data for mean
  summary_data_mean <- long_boot_mean %>%
    group_by(mycorrhizal_type) %>%
    summarise(
      mean = median(emmean, na.rm = TRUE),
      lower = quantile(emmean, 0.025, na.rm = TRUE),
      upper = quantile(emmean, 0.975, na.rm = TRUE),
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

#### * Extract results for each RC variable * ####

# Function to extract results for a given response variable
extract_results <- function(response_var, processed_results) {
  list(
    # Mean results
    density_mean_all = processed_results[[response_var]][["all"]]$mean$density,
    mean_all = processed_results[[response_var]][["all"]]$mean$summary,
    density_mean_tropical = processed_results[[response_var]][["tropical"]]$mean$density,
    mean_tropical = processed_results[[response_var]][["tropical"]]$mean$summary,
    density_mean_nontropical = processed_results[[response_var]][["nontropical"]]$mean$density,
    mean_nontropical = processed_results[[response_var]][["nontropical"]]$mean$summary,
    
    # Effect results
    density_effect_all = processed_results[[response_var]][["all"]]$effect$density,
    effect_all = processed_results[[response_var]][["all"]]$effect$summary,
    density_effect_tropical = processed_results[[response_var]][["tropical"]]$effect$density,
    effect_tropical = processed_results[[response_var]][["tropical"]]$effect$summary,
    density_effect_nontropical = processed_results[[response_var]][["nontropical"]]$effect$density,
    effect_nontropical = processed_results[[response_var]][["nontropical"]]$effect$summary,
    
    # Intercept results
    intercept_all = processed_results[[response_var]][["all"]]$intercept,
    intercept_tropical = processed_results[[response_var]][["tropical"]]$intercept,
    intercept_nontropical = processed_results[[response_var]][["nontropical"]]$intercept
  )
}

# Organise the results: Unique name to conduct sensitivity along side side breadth
RC1_results <- extract_results("RC1", processed_results)
RC2_results <- extract_results("RC2", processed_results)
RC3_results <- extract_results("RC3", processed_results)
data_position <- data
data_position_tropical <- data_tropical
data_position_nontropical <- data_nontropical
pem_position <- pem
pem_position_tropical <- pem_tropical
pem_position_nontropical <- pem_nontropical
phylo_prox_position <- phylo_prox
phylo_prox_position_tropical <- phylo_prox_tropical
phylo_prox_position_nontropical <- phylo_prox_nontropical
pem_select_position <- pem_select
base_model_params_position <- base_model_params
pem_model_params_position <- pem_model_params
phylo_auto_stats_position <- phylo_auto_stats_combined

#### * Save generated data * ####
save(
  # Raw data
  data_position, data_position_tropical, data_position_nontropical,
  
  # Phylogenetic data
  pem_position, pem_position_tropical, pem_position_nontropical,
  phylo_prox_position, phylo_prox_position_tropical, phylo_prox_position_nontropical,
  pem_select_position,
  
  # Results
  RC1_results, RC2_results, RC3_results,
  base_model_params_position, pem_model_params_position,
  phylo_auto_stats_position,
  
  file = "generated_data/figure_4.RData"
)
