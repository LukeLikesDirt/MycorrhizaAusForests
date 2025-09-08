# VIF ##########################################################################
VIF = function(x) {
  x = as.data.frame(x)
  # VIF calculation
  form    = formula(paste("fooy ~ ", paste(strsplit(names(x), " "), collapse = " + ")))
  x       = data.frame(fooy = 1 + rnorm(nrow(x)) ,x)
  lm_mod  = lm(form, x)
  # End VIF calculation
  cat("\n\nVariance inflation factors\n\n")
  print(supVIF(lm_mod))
}
# VIF function dependency
supVIF = function(mod) {
  v = vcov(mod)
  assign = attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v = v[-1, -1]
    assign = assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms = labels(terms(mod))
  n.terms = length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor) = 0
    if (any(tmp_cor == 1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R = cov2cor(v)
  detR = det(R)
  result = matrix(0, n.terms, 3)
  rownames(result) = terms
  colnames(result) = c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs = which(assign == term)
    result[term, 1] = det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] = length(subs)
  }
  if (all(result[, 2] == 1)) {
    result = data.frame(GVIF=result[, 1])
  } else {
    result[, 3] = result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

# Standardise ##################################################################
std = function(x) {(x - mean(x)) / sd(x)}
unstd = function(z, original_mean, original_sd) {
  (z * original_sd) + original_mean
}

# Function to standardise values to 100 percentiles for visualization
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

# Fill missing values using k-nearest neighbors ################################

fill_missing_knn <- function(data, x, y, variable, k = 1) {
  library(RANN)
  # Input validation
  if (!all(c(x, y, variable) %in% names(data))) {
    stop("One or more specified columns not found in data")
  }
  
  # Create a copy of the data to avoid modifying the original
  data_copy <- data
  
  # Identify rows with missing values in the target variable
  missing_idx <- which(data_copy[[variable]] == "" | is.na(data_copy[[variable]]))
  
  # If no missing values, return original data
  if (length(missing_idx) == 0) return(data_copy)
  
  # Get complete cases for reference
  complete_idx <- which(data_copy[[variable]] != "" & !is.na(data_copy[[variable]]))
  
  # Create coordinate matrices
  missing_coords <- cbind(data_copy[[x]][missing_idx], 
                          data_copy[[y]][missing_idx])
  reference_coords <- cbind(data_copy[[x]][complete_idx], 
                            data_copy[[y]][complete_idx])
  
  # Find k nearest neighbors
  nn <- nn2(reference_coords, missing_coords, k = k)
  
  # For k = 1, directly assign nearest neighbor values
  if (k == 1) {
    data_copy[[variable]][missing_idx] <- data_copy[[variable]][complete_idx][nn$nn.idx]
  } else {
    # For k > 1, use majority voting for categorical variables
    # or mean for numeric variables
    for (i in seq_along(missing_idx)) {
      neighbor_values <- data_copy[[variable]][complete_idx][nn$nn.idx[i,]]
      
      if (is.numeric(data_copy[[variable]])) {
        # For numeric variables, use weighted mean based on distance
        weights <- 1 / (nn$nn.dists[i,] + .Machine$double.eps)  # Add small constant to avoid division by zero
        data_copy[[variable]][missing_idx[i]] <- weighted.mean(neighbor_values, weights)
      } else {
        # For categorical variables, use majority voting
        data_copy[[variable]][missing_idx[i]] <- names(sort(table(neighbor_values), decreasing = TRUE))[1]
      }
    }
  }
  
  return(data_copy)
}

# Example usage:
# result <- fill_missing_knn(
#   data = covariates,
#   x = "x_albers",
#   y = "y_albers",
#   variable = "ecoregion",
#   k = 3
# )

# INLA model prams #############################################################

control_compute <- list(
  config = TRUE,
  waic = TRUE,
  dic = TRUE,
  residuals = TRUE
)

# Function to calculate study area dimensions and recommended mesh parameters
estimate_mesh_params <- function(coords) {
  # Calculate bounding box
  xrange <- range(coords[,1])
  yrange <- range(coords[,2])
  
  # Calculate dimensions
  width <- diff(xrange)
  height <- diff(yrange)
  
  # Calculate diameter (using diagonal of bounding box)
  diameter <- sqrt(width^2 + height^2)
  
  # Calculate recommended max.edge values
  max_edge_1_10 <- diameter/10
  max_edge_1_20 <- diameter/20
  
  # Calculate study area center
  center_x <- mean(xrange)
  center_y <- mean(yrange)
  
  # Return results
  results <- list(
    dimensions = list(
      width = width,
      height = height,
      diameter = diameter,
      center = c(x = center_x, y = center_y)
    ),
    recommendations = list(
      max_edge_1_10 = max_edge_1_10,
      max_edge_1_20 = max_edge_1_20,
      suggested_cutoff_1_10 = max_edge_1_10/10,
      suggested_cutoff_1_20 = max_edge_1_20/10
    )
  )
  
  # Print summary
  cat("Study Area Dimensions:\n")
  cat("Width:", width, "\n")
  cat("Height:", height, "\n")
  cat("Diameter:", diameter, "\n")
  cat("\nRecommended mesh parameters:\n")
  cat("max.edge (1/10 of diameter):", max_edge_1_10, "\n")
  cat("max.edge (1/20 of diameter):", max_edge_1_20, "\n")
  cat("Suggested cutoff (1/10 rule):", max_edge_1_10/10, "to", max_edge_1_20/10, "\n")
  
  return(invisible(results))
}

# Variogram ####################################################################

predict_range <- function(data, formula, x, y, crs, plot = TRUE) {
  library(sp)
  library(gstat)
  library(ggplot2)
  library(tidyr)
  
  # Work on a copy of the data to avoid modifying the original object
  data_sp <- data
  
  # Build the coordinates formula dynamically using the specified column names
  coord_formula <- as.formula(paste("~", x, "+", y))
  
  # Convert the data frame to a SpatialPointsDataFrame
  coordinates(data_sp) <- coord_formula
  
  # Assign the coordinate reference system
  proj4string(data_sp) <- CRS(crs)
  
  # Compute the empirical variogram using the specified formula
  vgm_emp <- variogram(formula, data = data_sp)
  
  # Fit an exponential variogram model to the empirical variogram
  vgm_fit <- fit.variogram(vgm_emp, model = vgm("Exp"))
  
  # Combine empirical and fitted variogram for plotting
  fitted_variogram <- variogramLine(vgm_fit, maxdist = max(vgm_emp$dist))
  
  # Optionally, plot the empirical variogram with the fitted model overlay
  if (plot) {
    plot <- ggplot() +
      geom_point(data = as.data.frame(vgm_emp), aes(x = dist, y = gamma), color = "blue") +
      geom_line(data = as.data.frame(fitted_variogram), aes(x = dist, y = gamma), color = "red") +
      labs(
        x = "Distance",
        y = "Semivariance",
        title = "Empirical Variogram with Fitted Model"
      ) +
      theme_minimal()
    
    print(plot)
  }
  
  # Extract the range parameter from the fitted variogram model
  range_estimate <- vgm_fit$range[2]
  
  # Return the estimated range
  return(range_estimate)
}

# INLA diagnostics #############################################################

# Inla plot posterior predictive densities
plot_posterior_predictive <- function(model, response, data, n_samples = 1000) {
  
  # Extract predictive values and format them for ggplot2
  set.seed(1986)
  predictive_values <- sapply(inla.posterior.sample(n_samples, model), function(sample) {
    sample$latent[grep("Predictor", rownames(sample$latent))]
  }) %>% 
    as.data.frame(t(.)) %>%
    pivot_longer(
      everything(),
      names_to = "sample",
      values_to = "predictive_value"
    )
  
  # Create the plot
  p <- ggplot() +
    geom_density(
      data = data,
      aes_string(x = response),
      color = "blue", size = 1
    ) +
    geom_density(
      data = predictive_values,
      aes(x = exp(predictive_value), group = sample),
      color = "red",
      alpha = 0.05, size = 0.05
    ) +
    labs(x = response, y = "Density", title = "Observed vs. Posterior Predictive Densities") +
    theme_minimal() +
    theme(legend.position = c(0.85, 0.85)) +
    scale_color_manual(values = c("blue", "red"),
                       labels = c("Observed", "Posterior Prediction")) +
    guides(color = guide_legend(title = "Legend"))
  
  return(p)
}

inla_spat_diagnostic_data <- function(model_names, observed_values, stack_names, backtransform = "none", n_samples = 1000) {
  
  # Input validation
  if (!is.character(model_names) || !is.character(stack_names)) {
    stop("model_names and stack_names must be character vectors")
  }
  if (length(model_names) != length(stack_names)) {
    stop("model_names and stack_names must have the same length")
  }
  if (!backtransform %in% c("none", "exp", "plogis")) {
    stop("backtransform must be one of: 'none', 'exp', 'plogis'")
  }
  if (length(observed_values) != length(model_names)) {
    stop("observed_values must have the same length as model_names")
  }
  
  # Initialize vectors
  observed <- c()
  fitted <- c()
  upper <- c()
  lower <- c()
  residuals <- c()
  pearson_residuals <- c()
  posterior_mean <- c()
  posterior_lower <- c()
  posterior_upper <- c()
  models <- c()
  
  # Loop through each model
  for (i in seq_along(model_names)) {
    # Get model object from its name with error handling
    model <- try(get(model_names[i]), silent = TRUE)
    if (inherits(model, "try-error")) {
      stop(paste("Could not find model:", model_names[i]))
    }
    
    # Extract model number from the suffix
    model_suffix <- sub(".*?_", "", model_names[i])
    
    # Define the stack name
    stack_name <- try(get(stack_names[i]), silent = TRUE)
    if (inherits(stack_name, "try-error")) {
      stop(paste("Could not find stack:", stack_names[i]))
    }
    
    # Estimation data indexes and number of observations
    index_est <- inla.stack.index(stack_name, "est")$data
    n_observations <- length(index_est)
    
    # Extract fitted values
    fitted_values <- model$summary.fitted.values[index_est, "mean"]
    lower_values <- model$summary.fitted.values[index_est, "0.025quant"]
    upper_values <- model$summary.fitted.values[index_est, "0.975quant"]
    
    # Get the corresponding observed values for this model
    current_observed <- observed_values[[i]][index_est]
    
    # Calculate residuals
    residuals_i <- current_observed - fitted_values
    
    # Calculate pearson residuals
    pearson_residuals_i <- residuals_i / sqrt(fitted_values)
    
    # Generate posterior predictive samples
    tryCatch({
      set.seed(1986)
      posterior_samples <- inla.posterior.sample(n_samples, model)
      predictive_values <- sapply(posterior_samples, function(sample) {
        sample$latent[grep("Predictor", rownames(sample$latent))]
      })
      predictive_mean <- apply(predictive_values, 1, mean)
      predictive_quantiles <- apply(predictive_values, 1, quantile, probs = c(0.025, 0.5, 0.975))
    }, error = function(e) {
      stop(paste("Error in posterior sampling for model", model_names[i], ":", e$message))
    })
    
    # Apply backtransformation
    if (backtransform == "exp") {
      posterior_pred_df <- data.frame(
        posterior_mean = exp(predictive_mean[index_est]),
        lower = exp(predictive_quantiles[1, index_est]),
        upper = exp(predictive_quantiles[3, index_est])
      )
    } else if (backtransform == "plogis") {
      posterior_pred_df <- data.frame(
        posterior_mean = plogis(predictive_mean[index_est]),
        lower = plogis(predictive_quantiles[1, index_est]),
        upper = plogis(predictive_quantiles[3, index_est])
      )
    } else {
      posterior_pred_df <- data.frame(
        posterior_mean = predictive_mean[index_est],
        lower = predictive_quantiles[1, index_est],
        upper = predictive_quantiles[3, index_est]
      )
    }
    
    # Append values to vectors
    observed <- c(observed, current_observed)
    fitted <- c(fitted, fitted_values)
    upper <- c(upper, upper_values)
    lower <- c(lower, lower_values)
    residuals <- c(residuals, residuals_i)
    pearson_residuals <- c(pearson_residuals, pearson_residuals_i)
    posterior_mean <- c(posterior_mean, posterior_pred_df$posterior_mean)
    posterior_lower <- c(posterior_lower, posterior_pred_df$lower)
    posterior_upper <- c(posterior_upper, posterior_pred_df$upper)
    models <- c(models, rep(paste0("Model ", model_suffix), n_observations))
  }
  
  # Check for consistent lengths
  expected_length <- length(observed)
  if (!all(sapply(list(fitted, upper, lower, residuals, pearson_residuals,
                       posterior_mean, posterior_lower, posterior_upper),
                  length) == expected_length)) {
    stop("Mismatch in lengths of computed values")
  }
  
  # Create the tibble
  data_diagnostics <- tibble(
    observed = observed,
    fitted = fitted,
    upper = upper,
    lower = lower,
    residuals = residuals,
    pearson_residuals = pearson_residuals,
    posterior_mean = posterior_mean,
    posterior_lower = posterior_lower,
    posterior_upper = posterior_upper,
    model = models
  )
  
  return(data_diagnostics)
}

# Inla fixed effects
inla_fixed_effects <- function(model_names) {
  # List to store results
  results <- list()
  
  # Iterate over each model
  for (i in seq_along(model_names)) {
    # Get model object from its name
    model <- get(model_names[i])
    
    # Extract model suffix from the suffix
    model_suffix <- sub(".*?_", "", model_names[i])
    
    # Extract the fixed effects summary, convert to data frame, and add model identifier
    model_df <- model$summary.fixed %>%
      as.data.frame() %>%
      rownames_to_column("parameter") %>%
      select(parameter, `0.025quant`, `0.5quant`, `0.975quant`) %>%
      mutate(
        model = paste0("Model ", model_suffix),
        parameter = case_when(
          grepl("log\\((.*)\\)", parameter) ~ gsub("log\\((.*)\\)", "\\1", parameter),
          TRUE ~ parameter
        ))
    
    # Add the data frame to the results list
    results[[i]] <- model_df
  }
  
  # Combine all data frames into one
  combined_df <- bind_rows(results)
  
  return(combined_df)
}

inla_random_effects <- function(model_names, parameter) {
  # Input validation
  if (!is.character(model_names)) {
    stop("model_names must be a character vector")
  }
  if (missing(parameter)) {
    stop("parameter must be specified")
  }
  
  # List to store results
  results <- list()
  
  # Iterate over each model
  for (i in seq_along(model_names)) {
    # Get model object from its name with error handling
    model <- try(get(model_names[i]), silent = TRUE)
    if (inherits(model, "try-error")) {
      stop(paste("Could not find model:", model_names[i]))
    }
    
    # Check if the random effect parameter exists in the model
    if (is.null(model$summary.random[[parameter]])) {
      stop(paste("Random effect", parameter, "not found in model", model_names[i]))
    }
    
    # Extract model suffix
    model_suffix <- sub(".*?_", "", model_names[i])
    
    # Extract the random effects summary and convert to data frame
    model_df <- model$summary.random[[parameter]] %>%
      as.data.frame() %>%
      select(ID = 1, `0.025quant`, `0.5quant`, `0.975quant`) %>%
      mutate(
        model = paste0("Model ", model_suffix),
        parameter = parameter
      )
    
    # Add the data frame to the results list
    results[[i]] <- model_df
  }
  
  # Combine all results into a single data frame
  combined_results <- bind_rows(results)
  
  return(combined_results)
}

compare_inla_models <- function(model_names, observed_values, stack_names) {
  # List to store results
  results <- list()
  
  # Iterate over each model
  for (i in seq_along(model_names)) {
    # Get model object from its name
    model <- get(model_names[i])
    
    # Define the stack name
    stack_name <- get(stack_names[i])
    
    # Estimation data indexes
    index_est <- inla.stack.index(stack_name, "est")$data
    n_observations <- length(index_est)
    
    # Extract fitted values for estimation
    fitted_values <- model$summary.fitted.values[index_est, "mean"]
    
    # Compute total sum of squares (sst) using only estimation data
    sst <- sum((observed_values[index_est] - mean(observed_values[index_est]))^2)
    
    # Compute residual sum of squares (ssr)
    ssr <- sum((observed_values[index_est] - fitted_values)^2)
    
    # Calculate RMSE
    rmse <- sqrt(mean((observed_values[index_est] - fitted_values)^2))
    
    # Calculate R-squared
    r2 <- 1 - (ssr / sst)
    
    # Handle missing WAIC and DIC values
    waic_val <- ifelse(!is.null(model$waic), model$waic$waic, NA)
    dic_val <- ifelse(!is.null(model$dic), model$dic$dic, NA)
    
    # Store results
    results[[paste0("Model_", i)]] <- list(
      waic = waic_val,
      dic = dic_val,
      rmse = rmse,
      r2 = r2
    )
  }
  
  # Convert list to a data frame
  results_df <- do.call(rbind, lapply(names(results), function(model) {
    data.frame(
      model = model,
      waic = results[[model]]$waic,
      dic = results[[model]]$dic,
      rmse = results[[model]]$rmse,
      r2 = results[[model]]$r2
    )
  }))
  
  return(results_df)
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

# Example usage:
# spatial_stats <- list(kappa = 0.5)
# x <- rnorm(100)
# y <- rnorm(100)
# plot_matern_correlation(spatial_stats, x, y, range = 200)
# plot_matern_correlation(spatial_stats, x, y)  # Without range


# Example usage:
# spatial_stats <- list(kappa = 0.5)
# x <- rnorm(100)
# y <- rnorm(100)
# plot_matern_correlation(spatial_stats, x, y, range = 200)
# plot_matern_correlation(spatial_stats, x, y)  # Without range

  
# Phylogenetic signal ##########################################################

# Define a function to compute the statistics for each variable
compute_phylo_stats <- function(residuals, phylo_prox, phylo_vcv) {
  cmean <- adephylo::abouheif.moran(residuals, phylo_prox)
  lambda <- phylosignal::lambdaTest(residuals, phylo_vcv)
  
  list(
    cmean_obs = cmean$obs,
    cmean_pvalue = cmean$pvalue,
    lambda_value = lambda$Lambda,
    lambda_pvalue = lambda$pvalue
  )
}

# Gamma model validation #######################################################

# This approach is specific to gamma models
# Adapted from: https://github.com/hrue/r-inla/blob/devel/rinla/vignettes/AA-group-cv.pdf

DIC <- function(x) {
  data.frame(mean.deviance = x$dic$mean.deviance,
             p.eff = x$dic$p.eff,
             DIC = x$dic$dic,
             WAIC = x$waic$waic)
}

CV_LOOCV <- function(x, O) {
  # Log Score (LS)
  LS <- -mean(log(x$cv))
  
  # Computation of E[Y_i | y_{I_i}] using trapezoid numerical integration
  expectation = numeric(length(O))
  for (i in 1:length(O)) {
    mu = x$mean[i]
    sd = x$sd[i]
    xx = seq(mu - 6 * sd, mu + 6 * sd, length.out = 100)
    yy = xx * dnorm(xx, mean = mu, sd = sd)
    expectation[i] = pracma::trapz(xx, yy)
  }
  
  # Mean Squared Predictive Error (MSPE)
  MSPE <- mean((expectation - O)^2)
  
  data.frame(`LS.LOOCV` = LS, `MSPE.LOOCV` = MSPE)
}

CV_LGOCV <- function(x, O) {
  # Log Score (LS)
  LS <- -mean(log(x$cv))
  
  # Computation of E[Y_i | y_{I_i}] using trapezoid numerical integration
  expectation = numeric(length(O))
  for (i in 1:length(O)) {
    mu = x$mean[i]
    sd = x$sd[i]
    xx = seq(mu - 6 * sd, mu + 6 * sd, length.out = 100)
    yy = xx * dnorm(xx, mean = mu, sd = sd)
    expectation[i] = pracma::trapz(xx, yy)
  }
  
  # Mean Squared Predictive Error (MSPE)
  MSPE <- mean((expectation - O)^2)
  
  data.frame(`LS.LGOCV` = LS, `MSPE.LGOCV` = MSPE)
}

# INLA summaries ###############################################################

inla_summary <- function(inla_obj) {
  # Extract summary statistics for fixed effects
  summary_fixed <- inla_obj$summary.fixed
  
  # Subset summary statistics for mean, standard deviation, and quantiles
  summary_subset <- summary_fixed[, c('mean', 'sd', '0.025quant', '0.975quant')]
  
  # Return the subsetted summary statistics
  return(summary_subset)
}

# Posterior Marginal Plot ######################################################

# Define the custom theme
my_theme <- function() {
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_line(colour = 'black', linewidth = 0.25),
    plot.title = element_text(hjust = 0.5, vjust = 1, size = rel(0.9)),
    legend.position = c(1, 1),  # Place legend at the top-right corner
    legend.justification = c(1.1, 1.1),  # Justify legend to the top-right corner
    legend.box.just = "right",  # Align legend box to the right
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Set margin to 0
    legend.background = element_blank()
  )
}

# Support function for the labels
sigma_dist_labs <- function() {
  labs(
    x = expression(sigma),
    y = expression(paste("P(", sigma, " | Data)"))
  )
}
 
# # Define the function to handle two specific models
# posterior_marginals <- function(model_1, model_2) {
#   model_1_name <- deparse(substitute(model_1))
#   model_2_name <- deparse(substitute(model_2))
#   
#   marginals_model_1 <- inla.tmarginal(
#     fun = function(x) exp(-0.5 * x),
#     marg = model_1$internal.marginals.hyperpar[[1]]
#   ) %>%
#     as_tibble() %>%
#     mutate(
#       Parameter = model_1_name
#     )
#   
#   marginals_model_2 <- inla.tmarginal(
#     fun = function(x) exp(-0.5 * x),
#     marg = model_2$internal.marginals.hyperpar[[1]]
#   ) %>%
#     as_tibble() %>%
#     mutate(
#       Parameter = model_2_name
#     )
#   
#   marginals_df <- bind_rows(marginals_model_1, marginals_model_2)
#   
#   ggplot(marginals_df, aes(x, y, group = Parameter, color = Parameter)) +
#     geom_line() +
#     ggtitle("Posterior marginal distribution of sigma") +
#     my_theme() +
#     sigma_dist_labs()
# }

# Define the function to handle an arbitrary number of models
posterior_marginals <- function(...) {
  models <- list(...)
  
  marginals_list <- lapply(seq_along(models), function(i) {
    model <- models[[i]]
    model_name <- paste0("Model ", i)
    
    inla.tmarginal(
      fun = function(x) exp(-0.5 * x),
      marg = model$internal.marginals.hyperpar[[1]]
    ) %>%
      as_tibble() %>%
      mutate(
        Parameter = model_name
      )
  })
  
  marginals_df <- bind_rows(marginals_list)
  
  ggplot(marginals_df, aes(x, y, group = Parameter, color = Parameter)) +
    geom_line() +
    ggtitle("Posterior marginal distribution of sigma") +
    my_theme() +
    sigma_dist_labs()
}

# Calculate niche position and breadth (RCs) #############################################

# Directly from the model coefficient list

# Function to calculate tree niche
calculate_tree_niche <- function(variable) {
  position_key <- paste0(variable, "_position")
  bredth_key <- paste0(variable, "_breadth")
  
  fitted_mean <- model_coefficients[[position_key]] %>%
    slice(1)
  fitted_sd <- model_coefficients[[bredth_key]] %>%
    slice(1)
  
  # Calculate interval values
  mean_val <- fitted_mean$mean
  lower <- mean_val - fitted_sd$mean
  upper <- mean_val + fitted_sd$mean
  
  # Combine the data into a data frame
  data.frame(
    mean = mean_val,
    lower = lower,
    upper = upper
  ) %>%
    setNames(paste0(variable, "_", names(.)))
}

# Function to read and process data for each variable for mycorrhizal types
calculate_mycorrhizal_niche <- function(
    variable, 
    selected_models, 
    data, 
    fitted_values,
    x
    ) {
  position_key <- paste0(variable, "_position")
  bredth_key <- paste0(variable, "_breadth")
  
  # Read in the fitted means (i.e. niche position)
  beta_mean <- selected_models[[position_key]]$summary.fixed[1:4, c("mean", "sd", "0.025quant", "0.975quant")]
  
  fitted_mean <- data %>%
    select(mycorrhizal_type) %>%
    mutate(
      fitted_mean = x %*% beta_mean[,1],
      sd = x %*% beta_mean[,2],
      lower_CI = x %*% beta_mean[,3],
      upper_CI = x %*% beta_mean[,4]
    ) %>%
    mutate(
      mycorrhizal_type = factor(
        mycorrhizal_type, levels = c("AM", "EcM", "EcM-AM", "NM")
      )) %>%
    unique(.)
  
  # Read in the fitted sd (i.e. niche width)
  beta_sd <- selected_models[[bredth_key]]$summary.fixed[1:4, c("mean", "sd", "0.025quant", "0.975quant")]
  
  fitted_sd <- data %>%
    select(mycorrhizal_type) %>%
    mutate(
      fitted_mean = x %*% beta_sd[,1],
      sd = x %*% beta_sd[,2],
      lower_CI = x %*% beta_sd[,3],
      upper_CI = x %*% beta_sd[,4]
    ) %>%
    mutate(
      mycorrhizal_type = factor(
        mycorrhizal_type, levels = c("AM", "EcM", "EcM-AM", "NM")
      )) %>%
    unique(.)
  
  # Merge means and SDs by mycorrhizal type
  merged_data <- inner_join(
    fitted_mean %>%
      select(mycorrhizal_type, fitted_mean) %>%
      rename(mean = fitted_mean),
    fitted_sd %>%
      select(mycorrhizal_type, fitted_mean) %>%
      rename(sd = fitted_mean),
    by = "mycorrhizal_type") %>%
    group_by(mycorrhizal_type) %>%
    summarise(
      mean = mean,
      lower = mean - sd,
      upper = mean + sd
    ) %>%
    ungroup() %>%
    rename_at(vars(mean, lower, upper), ~paste0(variable, "_", .))
  
  # Store fitted values in the list
  fitted_values[[position_key]] <<- fitted_mean
  fitted_values[[bredth_key]] <<-  fitted_sd
  
  return(merged_data)
}

# Function to calculate polygons for trees
calculate_tree_polygon <- function(data) {
  data %>%
    rowwise() %>%
    mutate(
      longitude_vertices = list(c(longitude_lower, longitude_upper, longitude_upper, longitude_lower, longitude_lower)),
      latitude_vertices = list(c(latitude_lower, latitude_lower, latitude_upper, latitude_upper, latitude_lower)),
      RC1_vertices = list(c(RC1_lower, RC1_upper, RC1_upper, RC1_lower, RC1_lower)),
      RC2_vertices = list(c(RC2_lower, RC2_lower, RC2_upper, RC2_upper, RC2_lower)),
      RC3_vertices = list(c(RC3_lower, RC3_upper, RC3_upper, RC3_lower, RC3_lower)),
      RC4_vertices = list(c(RC4_lower, RC4_lower, RC4_upper, RC4_upper, RC4_lower))
    ) %>%
    unnest(c(
      longitude_vertices, latitude_vertices,
      RC1_vertices, RC2_vertices, RC3_vertices, RC4_vertices
    )) %>%
    select(
      longitude_vertices, latitude_vertices,
      RC1_vertices, RC2_vertices, RC3_vertices, RC4_vertices
    )
}

# Function to calculate polygons for mycorrhizas
calculate_mycorrhiza_polygon <- function(niche_space_mycorrhizas) {
  niche_space_mycorrhizas %>%
    group_by(mycorrhizal_type) %>%
    rowwise() %>%
    mutate(
      longitude_vertices = list(c(longitude_lower, longitude_upper, longitude_upper, longitude_lower, longitude_lower)),
      latitude_vertices = list(c(latitude_lower, latitude_lower, latitude_upper, latitude_upper, latitude_lower)),
      RC1_vertices = list(c(RC1_lower, RC1_upper, RC1_upper, RC1_lower, RC1_lower)),
      RC2_vertices = list(c(RC2_lower, RC2_lower, RC2_upper, RC2_upper, RC2_lower)),
      RC3_vertices = list(c(RC3_lower, RC3_upper, RC3_upper, RC3_lower, RC3_lower)),
      RC4_vertices = list(c(RC4_lower, RC4_lower, RC4_upper, RC4_upper, RC4_lower))
    ) %>%
    unnest(c(
      longitude_vertices, latitude_vertices,
      RC1_vertices, RC2_vertices, RC3_vertices, RC4_vertices
    )) %>%
    select(
      longitude_vertices, latitude_vertices,
      RC1_vertices, RC2_vertices, RC3_vertices, RC4_vertices
    ) %>%
    ungroup()
}


# # Calculate niche position and breadth (raw variables fitted on Gamma) #######
# 
# # # !!! With uncertany intervals !!!
# # calculate_tree_niche <- function(variable) {
# #   # Read in the fitted means (i.e. niche optimum)
# #   fitted_mean <- fread(paste0("output/supplimentary_niche/coefficient_", variable, "_position.csv")) %>%
# #     slice(1) %>%
# #     select(mean, lower_ci_mean = "0.025quant", upper_ci_mean = "0.975quant")
# #   
# #   # Read in the fitted SD (i.e. the niche width)
# #   fitted_sd <- fread(paste0("output/supplimentary_niche/coefficient_", variable, "_breadth.csv")) %>%
# #     slice(1) %>%
# #     select(sd = mean, lower_ci_sd = "0.025quant", upper_ci_sd = "0.975quant")
# #   
# #   # Calculate the exponential values
# #   exp_mean <- exp(fitted_mean$mean)
# #   exp_lower_ci_mean <- exp(fitted_mean$lower_ci_mean)
# #   exp_upper_ci_mean <- exp(fitted_mean$upper_ci_mean)
# #   exp_sd <- exp(fitted_sd$sd)
# #   exp_lower_ci_sd <- exp(fitted_sd$lower_ci_sd)
# #   exp_upper_ci_sd <- exp(fitted_sd$upper_ci_sd)
# #   
# #   # Calculate the lower and upper bounds of the niche space
# #   lower <- exp_mean - exp_sd
# #   upper <- exp_mean + exp_sd
# #   
# #   # Calculate the confidence intervals for the lower and upper bounds
# #   lower_upper_ci_sd <- exp_mean - exp_lower_ci_sd
# #   lower_lower_ci_sd <- exp_mean - exp_upper_ci_sd
# #   upper_uppwe_ci_sd <- exp_mean + exp_upper_ci_sd
# #   upper_lower_ci_sd <- exp_mean + exp_lower_ci_sd
# #   
# #   # Combine the data into a data frame
# #   data.frame(
# #     mean = exp_mean,
# #     lower_ci_mean = exp_lower_ci_mean,
# #     upper_ci_mean = exp_upper_ci_mean,
# #     lower = lower,
# #     upper = upper,
# #     lower_upper_ci_sd = lower_upper_ci_sd,
# #     lower_lower_ci_sd = lower_lower_ci_sd,
# #     upper_upper_ci_sd = upper_uppwe_ci_sd,
# #     upper_lower_ci_sd = upper_lower_ci_sd
# #   ) %>%
# #     setNames(paste0(variable, "_", names(.)))
# # }
# 
# # !!! Directly from the model coefficient list !!!
# 
# # Function to calculate tree niche
# calculate_tree_niche <- function(variable) {
#   position_key <- paste0(variable, "_position")
#   bredth_key <- paste0(variable, "_breadth")
#   
#   fitted_mean <- model_coefficients[[position_key]] %>%
#     slice(1)
#   fitted_sd <- model_coefficients[[bredth_key]] %>%
#     slice(1)
#   
#   # Calculate the exponential values
#   mean_val <- exp(fitted_mean$mean)
#   lower <- mean_val - exp(fitted_sd$mean)
#   upper <- mean_val + exp(fitted_sd$mean)
#   
#   # Combine the data into a data frame
#   data.frame(
#     mean = mean_val,
#     lower = lower,
#     upper = upper
#   ) %>%
#     setNames(paste0(variable, "_", names(.)))
# }
# 
# # Function to calculate standardised tree niche
# calculate_tree_niche_std <- function(variable) {
#   position_key <- paste0(variable, "_position")
#   bredth_key <- paste0(variable, "_breadth")
#   
#   fitted_position <- model_coefficients[[position_key]] %>%
#     slice(1) %>%
#     pull(mean) %>%
#     exp()
#   fitted_breadth <- model_coefficients[[bredth_key]] %>%
#     slice(1) %>%
#     pull(mean) %>%
#     exp()
#   
#   sd_position <- sd(data[[position_key]])
#   sd_breadth <- sd(data[[bredth_key]])
#   
#   # Calculate the exponential values
#   mean_val <- (fitted_position - fitted_position) / sd_position
#   lower <- mean_val - (fitted_breadth / sd_breadth)
#   upper <- mean_val + (fitted_breadth / sd_breadth)
#   
#   # Combine the data into a data frame
#   data.frame(
#     mean = mean_val,
#     lower = lower,
#     upper = upper
#   ) %>%
#     setNames(paste0(variable, "_", names(.)))
# }
# 
# # Function to read and process data for each variable for mycorrhizal types
# calculate_mycorrhizal_niche <- function(variable, data, x) {
#   position_key <- paste0(variable, "_position")
#   bredth_key <- paste0(variable, "_breadth")
#   
#   # Read in the fitted means (i.e. niche position)
#   beta_mean <- selected_models[[position_key]]$summary.fixed[1:4, c("mean", "sd", "0.025quant", "0.975quant")]
#   
#   fitted_mean <- data %>%
#     select(mycorrhizal_type) %>%
#     mutate(
#       fitted_mean = x %*% beta_mean[,1],
#       sd = x %*% beta_mean[,2],
#       lower_CI = x %*% beta_mean[,3],
#       upper_CI = x %*% beta_mean[,4]
#     ) %>%
#     mutate(
#       across(where(is.numeric), ~ exp(.)),
#       mycorrhizal_type = factor(
#         mycorrhizal_type, levels = c("AM", "EcM", "EcM-AM", "NM")
#     )) %>%
#     unique(.)
#   
#   # Read in the fitted sd (i.e. niche width)
#   beta_sd <- selected_models[[bredth_key]]$summary.fixed[1:4, c("mean", "sd", "0.025quant", "0.975quant")]
#   
#   fitted_sd <- data %>%
#     select(mycorrhizal_type) %>%
#     mutate(
#       fitted_mean = x %*% beta_sd[,1],
#       sd = x %*% beta_sd[,2],
#       lower_CI = x %*% beta_sd[,3],
#       upper_CI = x %*% beta_sd[,4]
#     ) %>%
#     mutate(
#       across(where(is.numeric), ~ exp(.)),
#       mycorrhizal_type = factor(
#         mycorrhizal_type, levels = c("AM", "EcM", "EcM-AM", "NM")
#       )) %>%
#     unique(.)
#   
#   # Merge means and SDs by mycorrhizal type
#   merged_data <- inner_join(
#     fitted_mean %>%
#       select(mycorrhizal_type, fitted_mean) %>%
#       rename(mean = fitted_mean),
#     fitted_sd %>%
#       select(mycorrhizal_type, fitted_mean) %>%
#       rename(sd = fitted_mean),
#     by = "mycorrhizal_type") %>%
#     group_by(mycorrhizal_type) %>%
#     summarise(
#       mean = mean,
#       lower = mean - sd,
#       upper = mean + sd
#     ) %>%
#     ungroup() %>%
#     rename_at(vars(mean, lower, upper), ~paste0(variable, "_", .))
#   
#   # Store fitted values in the list
#   fitted_values[[position_key]] <<- fitted_mean
#   fitted_values[[bredth_key]] <<-  fitted_sd
#   
#   return(merged_data)
# }
# 
# # Function to calculate polygons for trees
# calculate_tree_polygon <- function(data) {
#   data %>%
#     rowwise() %>%
#     mutate(
#       bio1_vertices = list(c(bio1_lower, bio1_upper, bio1_upper, bio1_lower, bio1_lower)),
#       bio12_vertices = list(c(bio12_lower, bio12_lower, bio12_upper, bio12_upper, bio12_lower)),
#       POC_proportion_vertices = list(c(POC_proportion_lower, POC_proportion_upper, POC_proportion_upper, POC_proportion_lower, POC_proportion_lower)),
#       MAOC_proportion_vertices = list(c(MAOC_proportion_lower, MAOC_proportion_lower, MAOC_proportion_upper, MAOC_proportion_upper, MAOC_proportion_lower)),
#       CN_ratio_vertices = list(c(CN_ratio_lower, CN_ratio_upper, CN_ratio_upper, CN_ratio_lower, CN_ratio_lower)),
#       N_total_vertices = list(c(N_total_lower, N_total_lower, N_total_upper, N_total_upper, N_total_lower)),
#       CP_ratio_vertices = list(c(CP_ratio_lower, CP_ratio_upper, CP_ratio_upper, CP_ratio_lower, CP_ratio_lower)),
#       P_total_vertices = list(c(P_total_lower, P_total_lower, P_total_upper, P_total_upper, P_total_lower))
#     ) %>%
#     unnest(c(
#       bio1_vertices, bio12_vertices, CN_ratio_vertices, CP_ratio_vertices, 
#       N_total_vertices, P_total_vertices, POC_proportion_vertices, MAOC_proportion_vertices
#     )) %>%
#     select(
#       bio1_vertices, bio12_vertices, CN_ratio_vertices, CP_ratio_vertices, 
#       N_total_vertices, P_total_vertices, POC_proportion_vertices, MAOC_proportion_vertices
#     )
# }
# 
# # Function to calculate polygons for mycorrhizas
# calculate_mycorrhiza_polygon <- function(niche_space_mycorrhizas) {
#   niche_space_mycorrhizas %>%
#     group_by(mycorrhizal_type) %>%
#     rowwise() %>%
#     mutate(
#       bio1_vertices = list(c(bio1_lower, bio1_upper, bio1_upper, bio1_lower, bio1_lower)),
#       bio12_vertices = list(c(bio12_lower, bio12_lower, bio12_upper, bio12_upper, bio12_lower)),
#       POC_proportion_vertices = list(c(POC_proportion_lower, POC_proportion_upper, POC_proportion_upper, POC_proportion_lower, POC_proportion_lower)),
#       MAOC_proportion_vertices = list(c(MAOC_proportion_lower, MAOC_proportion_lower, MAOC_proportion_upper, MAOC_proportion_upper, MAOC_proportion_lower)),
#       CN_ratio_vertices = list(c(CN_ratio_lower, CN_ratio_upper, CN_ratio_upper, CN_ratio_lower, CN_ratio_lower)),
#       N_total_vertices = list(c(N_total_lower, N_total_lower, N_total_upper, N_total_upper, N_total_lower)),
#       CP_ratio_vertices = list(c(CP_ratio_lower, CP_ratio_upper, CP_ratio_upper, CP_ratio_lower, CP_ratio_lower)),
#       P_total_vertices = list(c(P_total_lower, P_total_lower, P_total_upper, P_total_upper, P_total_lower))
#     ) %>%
#     unnest(c(
#       bio1_vertices, bio12_vertices, CN_ratio_vertices, CP_ratio_vertices, 
#       N_total_vertices, P_total_vertices, MAOC_proportion_vertices, POC_proportion_vertices
#     )) %>%
#     select(
#       bio1_vertices, bio12_vertices, CN_ratio_vertices, CP_ratio_vertices, 
#       N_total_vertices, P_total_vertices, MAOC_proportion_vertices, POC_proportion_vertices
#     ) %>%
#     ungroup()
# }
# 
# # Function to read in and center effect sizes for plotting
# centre_effect_sizes <- function(variable) {
#   # Read in the effect size data
#   effect_sizes_optima <- fread(
#     paste0("data/niche_analysis/effects_", variable, "_position.csv")
#   ) %>%
#     select(mycorrhizal_type, mean_optima = mean, lower_optima = "0.025quant", upper_optima = "0.975quant")
#   
#   # Read in the fitted SD (i.e. the niche width)
#   effect_sizes_breadth <- fread(
#     paste0("data/niche_analysis/effects_", variable, "_breadth.csv")
#   ) %>%
#     select(mycorrhizal_type, mean_breadth = mean, lower_breadth = "0.025quant", upper_breadth = "0.975quant")
#   
#   # Center the data by subtracting the mean from each effect size
#   effect_sizes <- inner_join(
#     effect_sizes_optima, effect_sizes_breadth, by = "mycorrhizal_type"
#   ) %>%
#     mutate(
#       !!paste0(variable, "_mean_optima") := mean_optima - mean(mean_optima),
#       !!paste0(variable, "_lower_optima") := lower_optima - mean(mean_optima),
#       !!paste0(variable, "_upper_optima") := upper_optima - mean(mean_optima),
#       !!paste0(variable, "_mean_breadth") := mean_breadth - mean(mean_breadth),
#       !!paste0(variable, "_lower_breadth") := lower_breadth - mean(mean_breadth),
#       !!paste0(variable, "_upper_breadth") := upper_breadth - mean(mean_breadth)
#     )
#   
#   # Select only the renamed centered columns
#   effect_sizes <- effect_sizes %>%
#     select(mycorrhizal_type, starts_with(paste0(variable, "_")))
#   
#   return(effect_sizes)
#   
# }
# 
# # Function to read in and center effect size distributions for plotting
# centre_effect_dists <- function(variable) {
#   # Read in the effect size data for optimum
#   effect_dists_position <- fread(paste0("data/niche_analysis/effects_dist_", variable, "_position.csv")) %>%
#     select(x, y, mycorrhizal_type = parameter)
#   
#   effect_size_mean_position <- fread(paste0("data/niche_analysis/effects_", variable, "_position.csv")) %>%
#     select(mycorrhizal_type, mean)
#   
#   # Center the x values by subtracting the mean of the means for the optimum data
#   effect_dists_position <- effect_dists_position %>%
#     left_join(effect_size_mean_position, by = "mycorrhizal_type") %>%
#     mutate(
#       x = x - mean(mean),
#       type = "optimum",
#       mycorrhizal_type = factor(
#         mycorrhizal_type,
#         levels = c("AM", "ECM", "ECM-AM", "NM"))
#     ) %>%
#     select(-mean)
#   
#   # Read in the effect size data for width
#   effect_dists_breadth <- fread(paste0("data/niche_analysis/effects_dist_", variable, "_breadth.csv")) %>%
#     select(x, y, mycorrhizal_type = parameter)
#   
#   effect_size_mean_breadth <- fread(paste0("data/niche_analysis/effects_", variable, "_breadth.csv")) %>%
#     select(mycorrhizal_type, mean)
#   
#   # Center the x values by subtracting the mean of the means for the width data
#   effect_dists_breadth <- effect_dists_breadth %>%
#     left_join(effect_size_mean_breadth, by = "mycorrhizal_type") %>%
#     mutate(
#       x = x - mean(mean),
#       type = "width",
#       mycorrhizal_type = factor(
#         mycorrhizal_type,
#         levels = c("AM", "ECM", "ECM-AM", "NM"))
#     ) %>%
#     select(-mean)
#   
#   # Combine the centered data for optimum and width
#   effect_dists <- bind_rows(effect_dists_position, effect_dists_breadth)
#   
#   # Return the centered dataframe
#   return(effect_dists)
# }
# 

# Create phylo_listw ##########################################################

make_phylo_listw <- function(tree, data) {
  # Ensure tip order matches data
  species_ordered <- match(tree$tip.label, data$species)
  
  if (any(is.na(species_ordered))) {
    warning("Some tree tips are not found in the data. Check species names.")
  }
  
  # Compute phylogenetic proximity (patristic distances)
  prox <- adephylo::proxTips(tree, method = "patristic")
  
  # Convert proximity matrix to listw object for MEMs
  listw <- spdep::mat2listw(as.matrix(prox), style = "B")
  
  return(listw)
}

# Prune phylo listw ##########################################################

prune_phylo_listw <- function(tree, species_df) {
  species_keep <- intersect(tree$tip.label, species_df$species)
  tree_pruned <- ape::drop.tip(tree, setdiff(tree$tip.label, species_keep))
  prox <- adephylo::proxTips(tree_pruned, method = "patristic")
  listw <- spdep::mat2listw(as.matrix(prox), style = "B")
  return(listw)
}

# Sensitivity Analysis: Mycorrhizal Type Reassignment ##########################

reassign_mycorrhizal_type <- function(data, type_probs, prop_reassign = 0.10) {
  
  # Get species eligible for reassignment
  species_pool <- data %>%
    distinct(species, .keep_all = TRUE)
  
  n_reassign <- floor(nrow(species_pool) * prop_reassign)
  
  # Randomly select species to reassign
  species_to_reassign <- sample(species_pool$species, size = n_reassign)
  
  # Create reassignment lookup table
  reassigned <- tibble(
    species = species_to_reassign,
    mycorrhizal_type = sample(
      type_probs$mycorrhizal_type,
      size = n_reassign,
      replace = TRUE,
      prob = type_probs$prob
    )
  )
  
  # Replace the selected species' mycorrhizal types
  data_new <- data %>%
    left_join(reassigned, by = "species", suffix = c("", "_new")) %>%
    mutate(mycorrhizal_type = ifelse(!is.na(mycorrhizal_type_new), mycorrhizal_type_new, mycorrhizal_type)) %>%
    select(-mycorrhizal_type_new)
  
  return(data_new)
}
