# Set the seed
set.seed(1986)

# Required packages
library(boot)
library(emmeans)
library(tidyverse)
library(ggtext)
library(ggridges)
library(patchwork)
library(ape)
library(adephylo)
library(MPSEM)
library(V.PhyloMaker2)
library(data.table)
source("code/fast.phylo.eigenvector.select.R")

# (1) Data preparation ########################################################

# Load in niche data: Data prepared in "03e.Niche_breadth_diffs.R" and 
# "03f.Niche_position_diffs.R"
load("output/generated_data/figure_4.RData")
load("output/generated_data/figure_5.RData")

# Phylogenetic lists
pem_list <- list(
  breadth = list(
    all = pem_breadth,
    tropical = pem_breadth_tropical,
    nontropical = pem_breadth_nontropical
  ),
  position = list(
    all = pem_position,
    tropical = pem_position_tropical,
    nontropical = pem_position_nontropical
  )
)

phylo_prox_list <- list(
  breadth = list(
    all = phylo_prox_breadth,
    tropical = phylo_prox_breadth_tropical,
    nontropical = phylo_prox_breadth_nontropical
  ),
  position = list(
    all = phylo_prox_position,
    tropical = phylo_prox_position_tropical,
    nontropical = phylo_prox_position_nontropical
  )
)

# Ensure alignment
stopifnot(
  identical(rownames(data_breadth), rownames(pem_breadth)),
  identical(rownames(data_breadth), rownames(phylo_prox_breadth)),
  identical(rownames(data_breadth_tropical), rownames(pem_breadth_tropical)),
  identical(rownames(data_breadth_tropical), rownames(phylo_prox_breadth_tropical)),
  identical(rownames(data_breadth_nontropical), rownames(pem_breadth_nontropical)),
  identical(rownames(data_breadth_nontropical), rownames(phylo_prox_breadth_nontropical)),
  identical(rownames(data_position), rownames(pem_position)),
  identical(rownames(data_position), rownames(phylo_prox_position)),
  identical(rownames(data_position_tropical), rownames(pem_position_tropical)),
  identical(rownames(data_position_tropical), rownames(phylo_prox_position_tropical)),
  identical(rownames(data_position_nontropical), rownames(pem_position_nontropical)),
  identical(rownames(data_position_nontropical), rownames(phylo_prox_position_nontropical))
)

# (2) Sensitivity Analysis Functions ##########################################

create_sensitivity_datasets <- function(original_data, n_sensitivity = 5, resample_prop = 0.1) {
  myco_counts <- table(original_data$mycorrhizal_type)
  myco_probs <- myco_counts / sum(myco_counts)
  cat("Mycorrhizal type probabilities calculated from data:\n")
  print(round(myco_probs, 3))
  cat("\n")
  
  sensitivity_datasets <- list()
  
  for (i in 1:n_sensitivity) {
    sens_data <- original_data
    n_resample <- ceiling(nrow(sens_data) * resample_prop)
    groups <- split(seq_len(nrow(sens_data)), sens_data$mycorrhizal_type)
    expected <- n_resample * myco_probs
    n_resample_per_group <- floor(expected)
    remaining <- n_resample - sum(n_resample_per_group)
    if (remaining > 0) {
      frac <- expected - n_resample_per_group
      add_idx <- order(-frac)[1:remaining]
      n_resample_per_group[add_idx] <- n_resample_per_group[add_idx] + 1
    }
    resample_indices <- c()
    group_names <- names(myco_probs)
    for (j in seq_along(groups)) {
      n_res_g <- n_resample_per_group[group_names[j]]
      if (n_res_g > 0) {
        sampled <- sample(groups[[group_names[j]]], min(n_res_g, length(groups[[group_names[j]]])))
        resample_indices <- c(resample_indices, sampled)
      }
    }
    
    new_myco_types <- sample(names(myco_probs), size = length(resample_indices), replace = TRUE, prob = myco_probs)
    sens_data$mycorrhizal_type[resample_indices] <- new_myco_types
    sens_data$mycorrhizal_type <- factor(sens_data$mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
    contrasts(sens_data$mycorrhizal_type) <- contr.sum(levels(sens_data$mycorrhizal_type))
    colnames(contrasts(sens_data$mycorrhizal_type)) <- levels(sens_data$mycorrhizal_type)[1:(length(levels(sens_data$mycorrhizal_type))-1)]
    rownames(sens_data) <- sens_data[["species"]]
    sensitivity_datasets[[paste0("sensitivity_", i)]] <- sens_data
  }
  
  return(sensitivity_datasets)
}

# Create sensitivity datasets for 10% resampling
cat("=== Creating sensitivity datasets at 10% resampling ===\n")
cat("\nBreadth - All data:\n")
sensitivity_data_breadth_all <- create_sensitivity_datasets(data_breadth, n_sensitivity = 5, resample_prop = 0.1)
cat("\nBreadth - Tropical data:\n")
sensitivity_data_breadth_tropical <- create_sensitivity_datasets(data_breadth_tropical, n_sensitivity = 5, resample_prop = 0.1)
cat("\nBreadth - Non-tropical data:\n")
sensitivity_data_breadth_nontropical <- create_sensitivity_datasets(data_breadth_nontropical, n_sensitivity = 5, resample_prop = 0.1)
cat("\nPosition - All data:\n")
sensitivity_data_position_all <- create_sensitivity_datasets(data_position, n_sensitivity = 5, resample_prop = 0.1)
cat("\nPosition - Tropical data:\n")
sensitivity_data_position_tropical <- create_sensitivity_datasets(data_position_tropical, n_sensitivity = 5, resample_prop = 0.1)
cat("\nPosition - Non-tropical data:\n")
sensitivity_data_position_nontropical <- create_sensitivity_datasets(data_position_nontropical, n_sensitivity = 5, resample_prop = 0.1)

# Combine all datasets
all_datasets <- list(
  breadth = list(
    original = list(
      all = data_breadth,
      tropical = data_breadth_tropical,
      nontropical = data_breadth_nontropical
    )
  ),
  position = list(
    original = list(
      all = data_position,
      tropical = data_position_tropical,
      nontropical = data_position_nontropical
    )
  )
)

# Add sensitivity datasets
for (i in 1:5) {
  sens_name <- paste0("sensitivity_", i)
  all_datasets$breadth[[sens_name]] <- list(
    all = sensitivity_data_breadth_all[[paste0("sensitivity_", i)]],
    tropical = sensitivity_data_breadth_tropical[[paste0("sensitivity_", i)]] %>% filter(biome == "tropical"),
    nontropical = sensitivity_data_breadth_nontropical[[paste0("sensitivity_", i)]] %>% filter(biome == "nontropical")
  )
  all_datasets$position[[sens_name]] <- list(
    all = sensitivity_data_position_all[[paste0("sensitivity_", i)]],
    tropical = sensitivity_data_position_tropical[[paste0("sensitivity_", i)]] %>% filter(biome == "tropical"),
    nontropical = sensitivity_data_position_nontropical[[paste0("sensitivity_", i)]] %>% filter(biome == "nontropical")
  )
}

# (3) Model Fitting and Bootstrap Analysis ####################################

response_vars <- list(
  breadth = c("env_breadth"),
  position = c("RC1", "RC2", "RC3")
)

# Fit base models
base_models_all <- list(breadth = list(), position = list())
for (data_type in c("breadth", "position")) {
  for (dataset_name in names(all_datasets[[data_type]])) {
    base_models_all[[data_type]][[dataset_name]] <- list()
    for (response in response_vars[[data_type]]) {
      base_models_all[[data_type]][[dataset_name]][[response]] <- list()
      for (subset_name in names(all_datasets[[data_type]][[dataset_name]])) {
        base_models_all[[data_type]][[dataset_name]][[response]][[subset_name]] <- lm(
          as.formula(paste(response, "~ mycorrhizal_type")),
          data = all_datasets[[data_type]][[dataset_name]][[subset_name]]
        )
      }
    }
  }
}

# Extract residuals
residuals_list_all <- list(breadth = list(), position = list())
for (data_type in c("breadth", "position")) {
  for (dataset_name in names(all_datasets[[data_type]])) {
    residuals_list_all[[data_type]][[dataset_name]] <- list()
    for (response in response_vars[[data_type]]) {
      residuals_list_all[[data_type]][[dataset_name]][[response]] <- list()
      for (subset_name in names(all_datasets[[data_type]][[dataset_name]])) {
        residuals_list_all[[data_type]][[dataset_name]][[response]][[subset_name]] <- 
          residuals(base_models_all[[data_type]][[dataset_name]][[response]][[subset_name]])
      }
    }
  }
}

# MEM selection
pem_select_all <- list(breadth = list(), position = list())
for (data_type in c("breadth", "position")) {
  for (dataset_name in names(all_datasets[[data_type]])) {
    pem_select_all[[data_type]][[dataset_name]] <- list()
    for (response in response_vars[[data_type]]) {
      pem_select_all[[data_type]][[dataset_name]][[response]] <- list()
      for (subset_name in names(all_datasets[[data_type]][[dataset_name]])) {
        cat("\n=== Running MEM selection for", data_type, dataset_name, response, subset_name, "===\n")
        pem_select_all[[data_type]][[dataset_name]][[response]][[subset_name]] <- fast.phylo.eigenvector.select(
          x = residuals_list_all[[data_type]][[dataset_name]][[response]][[subset_name]],
          phylo_eigenvectors = pem_list[[data_type]][[subset_name]],
          phylo_matrix = phylo_prox_list[[data_type]][[subset_name]],
          method = "abouheif",
          abouheif_method = "oriAbouheif",
          alpha = 0.05,
          verbose = TRUE
        )
      }
    }
  }
}

# Fit phylogenetic models
formulas_all <- list(breadth = list(), position = list())
pem_models_all <- list(breadth = list(), position = list())
for (data_type in c("breadth", "position")) {
  for (dataset_name in names(all_datasets[[data_type]])) {
    formulas_all[[data_type]][[dataset_name]] <- list()
    pem_models_all[[data_type]][[dataset_name]] <- list()
    for (response in response_vars[[data_type]]) {
      formulas_all[[data_type]][[dataset_name]][[response]] <- list()
      pem_models_all[[data_type]][[dataset_name]][[response]] <- list()
      for (subset_name in names(all_datasets[[data_type]][[dataset_name]])) {
        pem_vars <- colnames(pem_select_all[[data_type]][["original"]][[response]][[subset_name]]$PE.select)
        if (length(pem_vars) > 0) {
          formulas_all[[data_type]][[dataset_name]][[response]][[subset_name]] <- as.formula(paste(
            response, "~ mycorrhizal_type +", paste(pem_vars, collapse = " + ")
          ))
          cat("Formula for", data_type, dataset_name, response, subset_name, ":", 
              as.character(formulas_all[[data_type]][[dataset_name]][[response]][[subset_name]])[3], "\n")
        } else {
          formulas_all[[data_type]][[dataset_name]][[response]][[subset_name]] <- as.formula(paste(
            response, "~ mycorrhizal_type"
          ))
          cat("No pems selected for", data_type, dataset_name, response, subset_name, "- using base formula\n")
        }
        
        # Bind the selected PEMs to the data 
        pem_data <- bind_cols(
          all_datasets[[data_type]][[dataset_name]][[subset_name]], 
          pem_select_all[[data_type]][[dataset_name]][[response]][[subset_name]]$PE.select
        )
        pem_models_all[[data_type]][[dataset_name]][[response]][[subset_name]] <- lm(
          formulas_all[[data_type]][[dataset_name]][[response]][[subset_name]],
          data = pem_data
        )
      }
    }
  }
}

# Bootstrapping function (fixed to properly handle pems)
bootstrap_effects <- function(data, indices, pems, formula, response_var) {
  # Debug: Check inputs
  if (length(indices) == 0) {
    cat("Warning: Empty indices\n")
    return(rep(NA, 2 * length(unique(data$mycorrhizal_type)) + 1))
  }
  
  d <- data[indices, ]
  
  # Handle case where pems might be empty or NULL
  if (is.null(pems) || ncol(pems) == 0) {
    d_full <- d
  } else {
    pems_subset <- pems[indices, , drop = FALSE]
    # Check for column name conflicts
    common_names <- intersect(names(d), names(pems_subset))
    if (length(common_names) > 0) {
      cat("Warning: Column name conflicts:", paste(common_names, collapse = ", "), "\n")
    }
    d_full <- tryCatch({
      bind_cols(d, pems_subset)
    }, error = function(e) {
      cat("Error in bind_cols:", e$message, "\n")
      return(NULL)
    })
  }
  
  if (is.null(d_full)) {
    return(rep(NA, 2 * length(unique(data$mycorrhizal_type)) + 1))
  }
  
  # Check if all required variables are present
  formula_vars <- all.vars(formula)
  missing_vars <- setdiff(formula_vars, names(d_full))
  if (length(missing_vars) > 0) {
    cat("Missing variables in data:", paste(missing_vars, collapse = ", "), "\n")
    return(rep(NA, 2 * length(unique(data$mycorrhizal_type)) + 1))
  }
  
  model <- tryCatch({
    lm(formula, data = d_full)
  }, error = function(e) {
    cat("Error in lm:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(model)) {
    return(rep(NA, 2 * length(unique(data$mycorrhizal_type)) + 1))
  }
  
  # Check for model convergence issues
  if (any(is.na(coef(model)))) {
    cat("Warning: Model has NA coefficients\n")
    return(rep(NA, 2 * length(unique(data$mycorrhizal_type)) + 1))
  }
  
  emm <- tryCatch({
    suppressMessages(emmeans(model, "mycorrhizal_type", type = "response", data = d_full) %>% as_tibble())
  }, error = function(e) {
    cat("Error in emmeans:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(emm) || nrow(emm) == 0) {
    return(rep(NA, 2 * length(unique(data$mycorrhizal_type)) + 1))
  }
  
  emmean_values <- emm$emmean
  coef_names <- names(coef(model))
  myc_coefs <- coef(model)[grep("^mycorrhizal_type", coef_names)]
  myc_levels <- levels(data$mycorrhizal_type)
  n_levels <- length(myc_levels)
  effects <- rep(NA, n_levels)
  names(effects) <- myc_levels
  
  if (length(myc_coefs) == (n_levels - 1)) {
    for (i in 1:(n_levels - 1)) {
      level_name <- myc_levels[i]
      coef_name <- paste0("mycorrhizal_type", level_name)
      if (coef_name %in% names(myc_coefs)) {
        effects[level_name] <- myc_coefs[coef_name]
      }
    }
    ref_level <- myc_levels[n_levels]
    effects[ref_level] <- -sum(effects[1:(n_levels-1)], na.rm = TRUE)
  } else {
    cat("Warning: Unexpected number of mycorrhizal coefficients:", length(myc_coefs), "expected:", n_levels-1, "\n")
    return(rep(NA, 2 * n_levels + 1))
  }
  
  intercept <- coef(model)[1]
  return(c(emmean_values, effects, intercept))
}

# Run bootstrap
boot_results_all <- list(breadth = list(), position = list())
for (data_type in c("breadth", "position")) {
  for (dataset_name in names(all_datasets[[data_type]])) {
    cat("\n=== Processing", data_type, dataset_name, "===\n")
    boot_results_all[[data_type]][[dataset_name]] <- list()
    for (response in response_vars[[data_type]]) {
      boot_results_all[[data_type]][[dataset_name]][[response]] <- list()
      for (subset_name in names(all_datasets[[data_type]][[dataset_name]])) {
        cat("Running bootstrap for", response, subset_name, "\n")
        current_data <- all_datasets[[data_type]][[dataset_name]][[subset_name]]
        # PEMs from the selection with proper structure
        current_pems <- pem_select_all[[data_type]][[dataset_name]][[response]][[subset_name]]$PE.select
        current_formula <- formulas_all[[data_type]][[dataset_name]][[response]][[subset_name]]
        
        # Check for potential issues before bootstrap
        if (nrow(current_data) == 0) {
          cat("  ERROR: Empty dataset\n")
          next
        }
        
        if (nrow(current_pems) != nrow(current_data)) {
          cat("  ERROR: Row mismatch between data and pems\n")
          next
        }
        
        # Test the bootstrap function once before running full bootstrap
        cat("  Testing bootstrap function...\n")
        test_result <- tryCatch({
          bootstrap_effects(data = current_data, indices = 1:nrow(current_data), 
                            pems = current_pems, formula = current_formula, response_var = response)
        }, error = function(e) {
          cat("  ERROR in test bootstrap:", e$message, "\n")
          return(NULL)
        })
        
        if (is.null(test_result) || all(is.na(test_result))) {
          cat("  ERROR: Bootstrap function test failed\n")
          next
        } else {
          cat("  Bootstrap function test successful\n")
        }
        
        # Run the actual bootstrap
        boot_results_all[[data_type]][[dataset_name]][[response]][[subset_name]] <- boot::boot(
          data = current_data,
          strata = current_data$mycorrhizal_type,
          parallel = "multicore",
          ncpus = 10,
          statistic = function(data, indices) {
            bootstrap_effects(data = data, indices = indices, pems = current_pems, formula = current_formula, response_var = response)
          },
          R = 1000
        )
        
        # Check bootstrap results
        n_na <- sum(is.na(boot_results_all[[data_type]][[dataset_name]][[response]][[subset_name]]$t))
        total_vals <- length(boot_results_all[[data_type]][[dataset_name]][[response]][[subset_name]]$t)
        cat("  Bootstrap completed:", n_na, "NAs out of", total_vals, "total values\n")
      }
    }
  }
}

# (4) Results Processing #######################################################

process_boot_results <- function(boot_result, data_subset) {
  n_levels <- length(levels(data_subset$mycorrhizal_type))
  level_names <- levels(data_subset$mycorrhizal_type)
  mean_boot <- boot_result$t[, 1:n_levels, drop = FALSE]
  effect_boot <- boot_result$t[, (n_levels + 1):(2 * n_levels), drop = FALSE]
  intercept_boot <- boot_result$t[, 2 * n_levels + 1]
  intercept_summary <- tibble(
    mean = median(intercept_boot, na.rm = TRUE),
    lower = quantile(intercept_boot, 0.025, na.rm = TRUE),
    upper = quantile(intercept_boot, 0.975, na.rm = TRUE)
  ) %>%
    mutate(intercept_formatted = sprintf("%.3f [%.3f, %.3f]", mean, lower, upper))
  long_boot_effect <- as_tibble(effect_boot, .name_repair = "minimal") %>%
    setNames(level_names) %>%
    mutate(replicate = row_number()) %>%
    pivot_longer(-replicate, names_to = "mycorrhizal_type", values_to = "effect") %>%
    drop_na(effect)
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
    effect = list(density = density_data_effect, summary = summary_data_effect),
    intercept = intercept_summary
  ))
}

# Process results
processed_results_all <- list(breadth = list(), position = list())
for (data_type in c("breadth", "position")) {
  for (dataset_name in names(all_datasets[[data_type]])) {
    processed_results_all[[data_type]][[dataset_name]] <- list()
    for (response in response_vars[[data_type]]) {
      processed_results_all[[data_type]][[dataset_name]][[response]] <- list()
      for (subset_name in names(all_datasets[[data_type]][[dataset_name]])) {
        processed_results_all[[data_type]][[dataset_name]][[response]][[subset_name]] <- 
          process_boot_results(boot_results_all[[data_type]][[dataset_name]][[response]][[subset_name]], 
                               all_datasets[[data_type]][[dataset_name]][[subset_name]])
      }
    }
  }
}

# (5) Visualisation ############################################################

myco_types <- rev(c("AM", "EcM", "Dual", "NM"))
my_myco_colours <- c('NM' = '#d470a2', "AM" = "#E69F00", "EcM" = "#56B4E9", "Dual" = "#009E73")
tag_size <- 16
strip_size <- 10
title_size <- 8
text_size <- 7

common_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    axis.ticks = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
    strip.text = element_text(face = "bold", size = strip_size),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    aspect.ratio = 1.2
  )

# Calculate axis limits
get_all_x_limits <- function(data_type, response, subset_name) {
  all_x_values <- c()
  for (dataset_name in names(processed_results_all[[data_type]])) {
    density_data <- processed_results_all[[data_type]][[dataset_name]][[response]][[subset_name]]$effect$density
    all_x_values <- c(all_x_values, density_data$x)
  }
  return(c(min(all_x_values, na.rm = TRUE), max(all_x_values, na.rm = TRUE)))
}

x_limits <- list()
for (data_type in c("breadth", "position")) {
  x_limits[[data_type]] <- list()
  for (response in response_vars[[data_type]]) {
    x_limits[[data_type]][[response]] <- list()
    for (subset_name in c("all", "tropical", "nontropical")) {
      x_limits[[data_type]][[response]][[subset_name]] <- get_all_x_limits(data_type, response, subset_name)
    }
  }
}

# Effects plot function
effects_plot <- function(density_data, effect_data, intercept_data, y_title, show_plot_title = FALSE, plot_title = "", show_x_axis = TRUE, y_axis_title = FALSE, x_limits) {
  intercept_annotation <- tibble(
    x = x_limits[1],
    y = "AM",
    label = paste0("β<sub>0</sub> = ", intercept_data$intercept_formatted)
  )
  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_density_ridges(
      data = density_data,
      aes(x = x, y = mycorrhizal_type, height = y, fill = mycorrhizal_type),
      stat = "identity", scale = 0.9, colour = alpha("grey30", 0.5), alpha = 0.9
    ) +
    geom_errorbarh(
      data = effect_data,
      aes(y = mycorrhizal_type, xmin = lower, xmax = upper),
      height = 0, size = 0.8
    ) +
    geom_point(
      data = effect_data,
      aes(x = mean, y = mycorrhizal_type),
      shape = 21, fill = "white", stroke = 0.5, size = 1.2
    ) +
    geom_richtext(
      data = intercept_annotation,
      aes(x = x, y = y, label = label),
      hjust = 0.05, vjust = 1, size = 1.85,
      fill = NA, label.color = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = my_myco_colours, limits = myco_types) +
    scale_y_discrete(limits = rev(myco_types), labels = NULL) +
    scale_x_continuous(limits = x_limits, breaks = scales::pretty_breaks()) +
    common_theme +
    theme(
      axis.text.x = if(show_x_axis) element_text(angle = 45, hjust = 1) else element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = if(y_axis_title) element_text() else element_blank(),
      plot.title = if(show_plot_title) element_text() else element_blank(),
      plot.margin = if(show_plot_title) margin(3, 1, 1, 1, "pt") else margin(1, 1, 1, 1, "pt")
    ) +
    labs(title = if(show_plot_title) plot_title else NULL, y = if(y_axis_title) y_title else "")
  return(p)
}

# Create sensitivity panel for a specific subset
create_sensitivity_panel <- function(subset_name) {
  dataset_names <- c("original", paste0("sensitivity_", 1:5))
  y_titles <- c("Main model", paste0("Sensitivity ", 1:5))
  response_vars_ordered <- c("env_breadth", "RC1", "RC2", "RC3")
  plot_titles <- c("Env Breadth", "RC1", "RC2", "RC3")
  data_types <- c("breadth", rep("position", 3))
  plot_list <- list()
  
  for (i in 1:length(dataset_names)) {
    dataset_name <- dataset_names[i]
    y_title <- y_titles[i]
    for (j in 1:length(response_vars_ordered)) {
      response <- response_vars_ordered[j]
      data_type <- data_types[j]
      plot_title_val <- plot_titles[j]
      show_plot_title <- (i == 1)
      show_x_axis <- (i == length(dataset_names))
      y_axis_title <- (j == 1)
      
      results <- processed_results_all[[data_type]][[dataset_name]][[response]][[subset_name]]
      plot_list[[(i-1)*length(response_vars_ordered) + j]] <- effects_plot(
        density_data = results$effect$density,
        effect_data = results$effect$summary,
        intercept_data = results$intercept,
        y_title = y_title,
        show_plot_title = show_plot_title,
        plot_title = plot_title_val,
        show_x_axis = show_x_axis,
        y_axis_title = y_axis_title,
        x_limits = x_limits[[data_type]][[response]][[subset_name]]
      )
    }
  }
  
  return(wrap_plots(plot_list, nrow = 6, ncol = 4, byrow = TRUE))
}

# Create and save plots for each subset
subset_names <- c("all","tropical", "nontropical")
for (i in 1:length(subset_names)) {
  subset_name <- subset_names[i]
  combined_plot <- create_sensitivity_panel(subset_name)
  print(combined_plot)
  ggsave(
    filename = paste0("output/supplimentary_niche/sensitivity_analysis_", subset_name, ".png"),
    plot = combined_plot,
    width = 14,
    height = 21,
    bg = "white",
    unit = "cm",
    dpi = 300
  )
}

# (6) Visualisation (environmental breadth only) ###############################

myco_types <- rev(c("AM", "EcM", "Dual", "NM"))
my_myco_colours <- c('NM' = '#d470a2', "AM" = "#E69F00", "EcM" = "#56B4E9", "Dual" = "#009E73")
tag_size <- 16
strip_size <- 10
title_size <- 8
text_size <- 7

common_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    axis.ticks = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.title = element_text(face = "bold", size = strip_size, hjust = 0.5, vjust = 0),
    strip.text = element_text(face = "bold", size = strip_size),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    aspect.ratio = 1.2
  )

# Create mycorrhizal type legend
dummy_myco_data <- data.frame(
  mycorrhizal_type = factor(myco_types, levels = myco_types),
  value = 1:4
)

# Add theme if I want to use the plot as a legend
dummy_plot_myco <- ggplot(
  dummy_myco_data, 
  aes(x = value, y = value, fill = mycorrhizal_type)
) +
  geom_point(
    colour = alpha('grey30', 0.5),
    shape = 22,
    size = 5
  ) +
  scale_fill_manual(values = my_myco_colours, name = "Mycorrhizal\ntype") +
  theme_void() +
  scale_y_continuous(limits = c(0,0))

legend_myco <- cowplot::get_legend(dummy_plot_myco)

# Calculate axis limits
get_all_x_limits <- function(data_type, response, subset_name) {
  all_x_values <- c()
  for (dataset_name in names(processed_results_all[[data_type]])) {
    density_data <- processed_results_all[[data_type]][[dataset_name]][[response]][[subset_name]]$effect$density
    all_x_values <- c(all_x_values, density_data$x)
  }
  return(c(min(all_x_values, na.rm = TRUE), max(all_x_values, na.rm = TRUE)))
}

x_limits <- list()
for (data_type in c("breadth", "position")) {
  x_limits[[data_type]] <- list()
  for (response in response_vars[[data_type]]) {
    x_limits[[data_type]][[response]] <- list()
    for (subset_name in c("all", "tropical", "nontropical")) {
      x_limits[[data_type]][[response]][[subset_name]] <- get_all_x_limits(data_type, response, subset_name)
    }
  }
}

# Effects plot function
effects_plot <- function(density_data, effect_data, intercept_data, y_title, show_plot_title = FALSE, plot_title = "", show_x_axis = TRUE, y_axis_title = FALSE, x_limits, subset_name = NULL, dataset_name = NULL) {
  intercept_annotation <- tibble(
    x = x_limits[1],
    y = "AM",
    label = paste0("β<sub>0</sub> = ", intercept_data$intercept_formatted)
  )
  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_density_ridges(
      data = density_data,
      aes(x = x, y = mycorrhizal_type, height = y, fill = mycorrhizal_type),
      stat = "identity", scale = 0.9, colour = alpha("grey30", 0.5), alpha = 0.9
    ) +
    geom_errorbarh(
      data = effect_data,
      aes(y = mycorrhizal_type, xmin = lower, xmax = upper),
      height = 0, size = 0.8
    ) +
    geom_point(
      data = effect_data,
      aes(x = mean, y = mycorrhizal_type),
      shape = 21, fill = "white", stroke = 0.5, size = 1.2
    ) +
    geom_richtext(
      data = intercept_annotation,
      aes(x = x, y = y, label = label),
      hjust = 0.05, vjust = 1, size = 1.85,
      fill = NA, label.color = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = my_myco_colours, limits = myco_types) +
    scale_y_discrete(limits = rev(myco_types), labels = NULL) +
    scale_x_continuous(
      limits = x_limits, 
      breaks = seq(-1, 1, 0.05)
      ) +
    common_theme +
    theme(
      axis.text.x = if(show_x_axis) element_text(angle = 45, hjust = 1) else element_blank(),
      axis.title.x = if(subset_name == "tropical" && dataset_name == "sensitivity_5") element_text() else element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = if(y_axis_title) element_text(face = "bold") else element_blank(),
      plot.title = if(show_plot_title) element_text() else element_blank(),
      plot.margin = if(show_plot_title) margin(3, 1, 1, 1, "pt") else margin(1, 1, 1, 1, "pt")
    ) +
    labs(
      title = if(show_plot_title) plot_title else NULL,
      y = if(y_axis_title) y_title else "",
      x = if(subset_name == "tropical" && dataset_name == "sensitivity_5") "Effect size (mean dfference)" else ""
      )
  return(p)
}

# Create environmental breadth panel with subsets as columns
create_env_breadth_panel <- function() {
  dataset_names <- c("original", paste0("sensitivity_", 1:5))
  subset_names <- c("all", "tropical", "nontropical")
  subset_titles <- c("All", "Tropical", "Non-tropical")
  plot_list <- list()
  
  for (i in 1:length(dataset_names)) {
    dataset_name <- dataset_names[i]
    
    for (j in 1:length(subset_names)) {
      subset_name <- subset_names[j]
      subset_title <- subset_titles[j]
      
      show_plot_title <- (i == 1)  # Show title only in first row
      show_x_axis <- (i == length(dataset_names))  # Show x-axis only in last row
      y_axis_title <- (j == 1)  # Show y-axis title only in first column (all)
      
      results <- processed_results_all[["breadth"]][[dataset_name]][["env_breadth"]][[subset_name]]
      
      # Y-axis title text for the "all" column
      y_title_text <- if (y_axis_title) {
        if (i == 1) "Main model" else paste0("Sensitivity ", i-1)
      } else ""
      
      plot_list[[(i-1)*length(subset_names) + j]] <- effects_plot(
        density_data = results$effect$density,
        effect_data = results$effect$summary,
        intercept_data = results$intercept,
        y_title = y_title_text,
        show_plot_title = show_plot_title,
        plot_title = if (show_plot_title) subset_title else "",
        show_x_axis = show_x_axis,
        y_axis_title = y_axis_title,
        x_limits = x_limits[["breadth"]][["env_breadth"]][[subset_name]],
        subset_name = subset_name,
        dataset_name = dataset_name
      )
    }
  }
  
  return(wrap_plots(plot_list, nrow = 6, ncol = 3, byrow = TRUE))
}


# Create and save the environmental breadth plot
env_breadth_plot <- create_env_breadth_panel()

# Add the legend
env_breadth_plot_final <- cowplot::plot_grid(
  env_breadth_plot, legend_myco, rel_widths = c(1, 0.11)
)

# Display the plot
print(env_breadth_plot_final)

ggsave(
  filename = "output/supplimentary_niche/environmental_breadth_sensitivity.png",
  plot = env_breadth_plot,
  width = 10,
  height = 21,
  bg = "white",
  unit = "cm",
  dpi = 300
)
