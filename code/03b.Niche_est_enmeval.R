# 03b.Niche_est_enmeval.R
# ==============================================================================
# MAIN FUNCTIONS FOR ENMevaluate SPECIES DISTRIBUTION MODELING
# ==============================================================================
# This script contains the core functions for species distribution modeling
# using ENMevaluate, including spatial cross-validation, model fitting,
# and niche breadth calculation.
#
# Adapted to extract niche position values (suitability-weighted means) for
# each RC axis and geographic coordinates (longitude/latitude).
#
# ==============================================================================

# Load required libraries
library(terra)
library(data.table)
library(fst)
library(blockCV)
library(ENMeval)
library(dsmextra)
library(ecospat)
library(maxnet)
library(lhs)  # Required for Latin Hypercube Sampling in env.breadth.maxent

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Calculate B2 breadth metric from ENMTools
#' 
#' This function computes the B2 niche breadth metric, which is a measure of niche specialization
#' based on a vector of suitability scores or probabilities. It normalizes the input vector to sum to 1
#' and calculates the inverse Simpson index adjusted for sample size.
#' 
#' @param x Numeric vector of values (e.g. suitability scores). Non-finite and NA values are removed.
#' @return B2 niche breadth value (numeric). Returns NA if input is empty after cleaning.
calc.B2 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x <- x/sum(x)
  return((1/sum(x^2) - 1) / (length(x)-1))
}

#' Calculate B2 breadth metric: Adapted from env.breadth function and ENMTools
#' (see https://github.com/danlwarren/ENMTools/blob/master/R/env.breadth.R) to
#' work with maxent models (i.e. predict with type = "cloglog") and with 
#' SpatRaster objects
#' 
#' This function estimates environmental niche breadth (B2) using Latin Hypercube Sampling (LHS)
#' to generate points in environmental space and predict suitability using a Maxent model.
#' It iteratively samples until convergence based on a tolerance threshold.
#' 
#' @param model Maxent model object from maxnet.
#' @param env SpatRaster of environmental variables (must match model predictors).
#' @param tolerance Tolerance for convergence in mean B2 (default 1e-04).
#' @param max.reps Maximum repetitions to find initial suitable conditions (default 10).
#' @param chunk.size Size of each LHS sampling chunk (default 1e+05).
#' @return List containing environmental breadth (env.B2) estimate (mean B2 after convergence).
env.breadth.maxent <- function (model, env, tolerance = 1e-04, max.reps = 10, chunk.size = 1e+05) {
  
  if (inherits(env, "SpatRaster")) {
    mins <- terra::minmax(env)[1, ]
    maxes <- terra::minmax(env)[2, ]
  }
  continue <- FALSE
  n.reps <- 0
  while (continue == FALSE & n.reps < max.reps) {
    gens <- chunk.size
    this.lhs <- lhs::randomLHS(chunk.size, length(names(env)))
    predict.table <- t(t(this.lhs) * (maxes - mins) + mins)
    colnames(predict.table) <- names(env)
    pred <- as.numeric(predict(model, newdata = data.frame(predict.table), 
                               type = "cloglog"))
    if (max(pred) == 0) {
      this.B2 <- NA
    }
    else {
      this.B2 <- calc.B2(pred)
    }
    if (!is.na(this.B2)) {
      continue <- TRUE
    }
    else {
      n.reps <- n.reps + 1
    }
  }
  if (n.reps == max.reps) {
    warning("\n\nCould not find suitable starting conditions for environmental breadth, returning NA\n\n")
    return(list(env.B2 = NA, B2.plot = NA))
  }
  else {
    delta <- 1
    while (delta > tolerance) {
      this.lhs <- lhs::randomLHS(chunk.size, length(names(env)))
      predict.table <- t(t(this.lhs) * (maxes - mins) + 
                           mins)
      colnames(predict.table) <- names(env)
      pred <- as.numeric(predict(model, newdata = data.frame(predict.table), type = "cloglog"))
      if (max(pred) == 0) {
        next
      }
      else {
        this.B2 <- c(this.B2, calc.B2(pred))
        gens <- c(gens, max(gens) + chunk.size)
      }
      delta <- abs(mean(this.B2) - mean(this.B2[-length(this.B2)]))
    }
  }
  output <- list(env.B2 = mean(this.B2)
  )
  return(output)
}

#' Calculate observed B2 values from occurrence points in 3D environmental space
#' using kernel density estimation. Designed for RC1–RC3 rasters.
#'
#' This function estimates the observed environmental niche breadth using kernel density estimation (KDE)
#' in 3D environmental space. It scales the environmental values, computes KDE, and then calculates B2
#' from the normalized density estimates.
#'
#' @param occ_vect SpatVector of observed species occurrence points.
#' @param env_rast SpatRaster of 3 environmental predictors (e.g., RC1–RC3).
#' @return A list with raw and normalized B2 niche breadth estimates.
env.breadth.observed <- function(occ_vect, env_rast) {
  
  # --- Input checks
  if (!inherits(occ_vect, "SpatVector")) stop("occ_vect must be a SpatVector")
  if (!inherits(env_rast, "SpatRaster")) stop("env_rast must be a SpatRaster")
  if (terra::nlyr(env_rast) != 3) stop("env_rast must have exactly 3 layers (RC1–RC3)")
  
  # --- Extract and clean environmental values
  env_vals <- terra::extract(env_rast, occ_vect, ID = FALSE)
  env_vals <- na.omit(env_vals)
  if (nrow(env_vals) < 10) stop("Too few points with complete data to estimate niche breadth.")
  
  # --- Standardise predictors
  env_vals_scaled <- scale(env_vals)
  
  # --- KDE in 3D
  if (!requireNamespace("ks", quietly = TRUE)) {
    stop("Please install the 'ks' package to use this function: install.packages('ks')")
  }
  kde_result <- ks::kde(x = env_vals_scaled, compute.cont = FALSE, gridsize = rep(50, 3))
  prob <- as.vector(kde_result$estimate)
  prob <- prob / sum(prob)  # Normalise to sum to 1
  
  # --- Compute and return B2
  return(list(
    env.B2 = calc.B2(prob)
  ))
}

#' Calculate continuous Boyce Index for model evaluation with robust error handling
#' 
#' This function computes the continuous Boyce Index to evaluate model predictive performance.
#' It extracts background and occurrence predictions, cleans non-finite values, and performs
#' checks for sufficient data and variation before calculating the index using ecospat.
#' Returns NA on errors or insufficient data, with warnings for diagnostics.
#' 
#' @param pred SpatRaster of model predictions over the background area.
#' @param occ_vect SpatVector of occurrence points.
#' @return Numeric Boyce Index value (-1 to 1; >0.5 indicates good predictive performance), or NA on failure.
calc_boyce <- function(pred, occ_vect) {
  tryCatch({
    # Extract all prediction values (background)
    pred_all <- terra::values(pred, na.rm = TRUE)[, 1]
    
    # Check if we have valid background predictions
    if (length(pred_all) == 0) {
      cat("Warning: No valid background predictions found\n")
      return(NA)
    }
    
    # Remove non-finite values from background
    pred_all_finite <- pred_all[is.finite(pred_all)]
    if (length(pred_all_finite) == 0) {
      cat("Warning: No finite background predictions\n")
      return(NA)
    }
    
    # Extract predictions at occurrence points
    pred_occ_raw <- terra::extract(pred, occ_vect)
    if (ncol(pred_occ_raw) < 2) {
      cat("Warning: No prediction values extracted at occurrence points\n")
      return(NA)
    }
    
    pred_occ <- pred_occ_raw[, 2]  # Column 2 contains the values
    
    # Remove non-finite values from occurrence predictions
    pred_occ_finite <- pred_occ[is.finite(pred_occ)]
    if (length(pred_occ_finite) == 0) {
      cat("Warning: No finite occurrence predictions\n")
      return(NA)
    }
    
    # Check for sufficient variation in predictions
    var_all <- var(pred_all_finite)
    var_occ <- var(pred_occ_finite)
    if (is.na(var_all) || is.na(var_occ) || var_all == 0 || var_occ == 0) {
      cat("Warning: No variation in predictions (constant values)\n")
      return(NA)
    }
    
    # Check ranges
    bg_range <- range(pred_all_finite)
    occ_range <- range(pred_occ_finite)
    
    if (diff(bg_range) == 0 || diff(occ_range) == 0) {
      cat("Warning: Prediction range is zero\n")
      return(NA)
    }
    
    # Ensure we have reasonable sample sizes
    if (length(pred_all_finite) < 100) {
      cat("Warning: Too few background predictions (", length(pred_all_finite), ")\n")
      return(NA)
    }
    
    if (length(pred_occ_finite) < 5) {
      cat("Warning: Too few occurrence predictions (", length(pred_occ_finite), ")\n")
      return(NA)
    }
    
    # Calculate Boyce Index using ecospat
    boyce_res <- ecospat::ecospat.boyce(
      fit = pred_all_finite, 
      obs = pred_occ_finite, 
      nclass = 0,  # Use continuous method
      PEplot = FALSE
    )
    
    # Check if result is valid
    if (is.null(boyce_res) || is.null(boyce_res$cor)) {
      cat("Warning: ecospat.boyce returned NULL result\n")
      return(NA)
    }
    
    boyce_value <- boyce_res$cor
    
    # Final validation of result
    if (!is.finite(boyce_value)) {
      cat("Warning: Boyce index is not finite\n")
      return(NA)
    }
    
    return(boyce_value)
    
  }, error = function(e) {
    cat("Error in calc_boyce:", e$message, "\n")
    return(NA)
  })
}

#' Estimate optimal buffer radius for background sampling using iterative Maxent modeling
#' 
#' This function determines a data-driven buffer radius for background sampling by fitting
#' simple Maxent models on expanding concentric buffers around the occurrence centroid.
#' It evaluates model performance using the Boyce Index and stops when performance stabilizes
#' or reaches a threshold. Includes robust error handling and fallbacks for stability.
#' 
#' @param occ_vect SpatVector of occurrence points.
#' @param env SpatRaster of environmental layers (RC1-RC3).
#' @param spatial_folds List containing spatial fold information from blockCV.
#' @param start_radius_km Starting buffer radius in km (default 50).
#' @param max_radius_km Maximum buffer radius in km (default 500).
#' @param step_km Step size for increasing radius in km (default 50).
#' @param boyce_threshold Minimum Boyce Index for acceptable performance (default 0.3).
#' @param max_iters Maximum number of iterations (default 10).
#' @param min_bg_ratio Minimum ratio of background points to available cells (default 0.1).
#' @param range_estimate Spatial autocorrelation range estimate in meters (required for fallback).
#' @param fallback_km Fallback radius in km if no valid radius found (required).
#' @param stabilisation_delta Delta threshold for stabilisation to exit early (default 0.01).
#' @param fit_all_models If TRUE, fits all models and selects best stabilised radius (default FALSE).
#' @param n_iterations Number of iterations to compute mean CBI for robustness (default 5).
#' @return Optimal buffer radius in meters.
# Custom function to find minimal effective buffer radius
estimate_optimal_buffer <- function(occ_vect, env, spatial_folds, 
                                    start_radius_km = 50, max_radius_km = 500, step_km = 50, 
                                    boyce_threshold = 0.5, max_iters = 10, 
                                    min_bg_ratio = 0.1, range_estimate = NULL, fallback_km,
                                    stabilisation_delta = 0.01, 
                                    fit_all_models = FALSE, n_iterations = 5) {
  
  # Validate input parameters
  if (!is.finite(start_radius_km) || start_radius_km <= 0) {
    cat("Warning: Invalid start_radius_km, using 50 km\n")
    start_radius_km <- 50
  }
  if (!is.finite(max_radius_km) || max_radius_km <= start_radius_km) {
    max_radius_km <- start_radius_km + 200
    cat("Warning: Invalid max_radius_km, using", max_radius_km, "km\n")
  }
  if (!is.finite(step_km) || step_km <= 0) {
    step_km <- 50
    cat("Warning: Invalid step_km, using 50 km\n")
  }
  if (!is.finite(min_bg_ratio) || min_bg_ratio <= 0 || min_bg_ratio > 1) {
    min_bg_ratio <- 0.1
    cat("Warning: Invalid min_bg_ratio, using 0.1\n")
  }
  if (missing(fallback_km) || !is.finite(fallback_km) || fallback_km <= 0) {
    stop("fallback_km must be provided and be a positive finite number")
  }
  if (!is.finite(stabilisation_delta) || stabilisation_delta <= 0) {
    stabilisation_delta <- 0.05
    cat("Warning: Invalid stabilisation_delta, using 0.05\n")
  }
  
  # Validate n_iterations parameter
  if (!is.finite(n_iterations) || n_iterations < 1) {
    n_iterations <- 1
    cat("Warning: Invalid n_iterations, using 1\n")
  }
  n_iterations <- as.integer(n_iterations)
  
  if (n_iterations > 1) {
    cat("Using", n_iterations, "iterations to compute mean CBI\n")
  }
  
  # Convert fallback_km to meters for internal use
  fallback_r <- fallback_km * 1000
  
  radii <- seq(start_radius_km, max_radius_km, by = step_km) * 1000  # convert to meters
  boyce_scores <- numeric(length(radii))
  prev_boyce <- -Inf
  
  occ_crds <- crds(occ_vect)
  if (nrow(occ_crds) == 0) {
    cat("No occurrence points. Using fallback:", fallback_r / 1000, "km\n")
    return(fallback_r)
  }
  
  # Extract spatial fold information
  blocks_sf <- spatial_folds$blocks
  fold_assignments <- spatial_folds$occ_folds
  
  # Define function to process a single radius
  process_radius <- function(r, radius_index) {
    cat("Testing radius:", round(r / 1000), "km\n")
    
    # Create per-point buffering instead of centroid-based buffering
    sp_buffer <- tryCatch(terra::buffer(occ_vect, width = r), error = function(e) NULL)
    if (is.null(sp_buffer)) {
      cat("Skipping radius due to buffer error\n")
      return(list(index = radius_index, boyce = NA, error = "buffer_error"))
    }
    
    # Create mask for buffered area
    buffer_rast <- tryCatch({
      terra::rasterize(sp_buffer, env[[1]], field = 1)
    }, error = function(e) NULL)
    if (is.null(buffer_rast)) {
      cat("Skipping radius due to rasterize error\n")
      return(list(index = radius_index, boyce = NA, error = "rasterize_error"))
    }
    
    # Create buffered mask
    buffered_mask <- tryCatch({
      !is.na(buffer_rast) & !is.na(env[[1]])
    }, error = function(e) NULL)
    if (is.null(buffered_mask)) {
      cat("Skipping radius due to mask creation error\n")
      return(list(index = radius_index, boyce = NA, error = "mask_error"))
    }
    
    # Apply mask to environmental layers
    env_mask <- tryCatch({
      mask(env, buffered_mask, maskvalues = FALSE)
    }, error = function(e) NULL)
    if (is.null(env_mask)) {
      cat("Skipping radius due to environmental masking error\n")
      return(list(index = radius_index, boyce = NA, error = "env_mask_error"))
    }
    
    # Count available cells for background sampling
    available_cells <- sum(!is.na(values(env_mask[[1]])))
    if (available_cells == 0) {
      cat("No available cells in buffer, skipping\n")
      return(list(index = radius_index, boyce = NA, error = "no_cells"))
    }
    
    # Adaptive background sampling
    # Sample 50% of available background points, with minimum of 100
    target_bg <- max(100, floor(0.5 * available_cells))
    
    # Use all available cells if fewer than target
    n_bg <- min(available_cells, target_bg)
    
    cat("Available cells:", available_cells, "| Target bg:", target_bg, "| Using:", n_bg, "\n")
    
    # Sample background points
    bg_points <- tryCatch({
      spatSample(env_mask[[1]], size = n_bg, method = "random", na.rm = TRUE, xy = TRUE)
    }, error = function(e) NULL)
    
    if (is.null(bg_points) || nrow(bg_points) < 50) {
      cat("Too few background points (", ifelse(is.null(bg_points), 0, nrow(bg_points)), "), skipping\n")
      return(list(index = radius_index, boyce = NA, error = "insufficient_bg"))
    }
    
    # Extract environmental data
    pres_env <- tryCatch(terra::extract(env, occ_vect)[, -1], error = function(e) NULL)
    bg_env <- tryCatch(terra::extract(env, bg_points[, c("x", "y")])[, -1], error = function(e) NULL)
    
    if (is.null(pres_env) || is.null(bg_env)) {
      cat("Failed to extract environmental data, skipping\n")
      return(list(index = radius_index, boyce = NA, error = "extraction_failed"))
    }
    
    # Clean data (remove NA values)
    pres_env <- pres_env[complete.cases(pres_env), , drop = FALSE]
    bg_env <- bg_env[complete.cases(bg_env), , drop = FALSE]
    
    if (nrow(pres_env) < 5 || nrow(bg_env) < 50) {
      cat("Not enough clean presence/background data, skipping\n")
      return(list(index = radius_index, boyce = NA, error = "insufficient_data"))
    }
    
    # Prepare data for modeling with spatial fold information
    p <- c(rep(1, nrow(pres_env)), rep(0, nrow(bg_env)))
    data_all <- rbind(pres_env, bg_env)
    
    # Add fold information (background points get fold = 0 or separate handling)
    fold_info <- c(fold_assignments, rep(0, nrow(bg_env)))
    
    # Fit MaxEnt model
    model <- tryCatch({
      maxnet(p, data_all, f = maxnet.formula(p, data_all, classes = "lq"))
    }, error = function(e) NULL)
    
    if (is.null(model)) {
      cat("Model fitting failed, skipping\n")
      return(list(index = radius_index, boyce = NA, error = "model_failed"))
    }
    
    # Predict across the study area
    cells_df <- tryCatch(as.data.frame(env_mask, na.rm = TRUE), error = function(e) NULL)
    if (is.null(cells_df) || nrow(cells_df) == 0) {
      cat("No valid cells for prediction, skipping\n")
      return(list(index = radius_index, boyce = NA, error = "no_pred_cells"))
    }
    
    pred_vals <- tryCatch(predict(model, newdata = cells_df, type = "cloglog"), error = function(e) NULL)
    if (is.null(pred_vals)) {
      cat("Prediction failed, skipping\n")
      return(list(index = radius_index, boyce = NA, error = "prediction_failed"))
    }
    
    # Create prediction raster
    pred <- env_mask[[1]] * NA
    non_na_cells <- which(!is.na(values(env_mask[[1]])))
    if (length(pred_vals) >= length(non_na_cells)) {
      pred[non_na_cells] <- pred_vals[1:length(non_na_cells)]
    } else {
      cat("Mismatch in prediction cells, skipping\n")
      return(list(index = radius_index, boyce = NA, error = "pred_mismatch"))
    }
    
    # Calculate Boyce index (with iterations if specified)
    if (n_iterations == 1) {
      # Single iteration
      boyce_score <- tryCatch({
        score <- calc_boyce(pred, occ_vect)
        if (is.na(score) || is.nan(score)) NA else score
      }, error = function(e) {
        cat("Boyce calculation failed: ", e$message, "\n")
        return(NA)
      })
    } else {
      # Multiple iterations - compute mean CBI
      boyce_scores_iter <- numeric(n_iterations)
      
      for (iter in 1:n_iterations) {
        # Re-sample background points for each iteration
        bg_points_iter <- tryCatch({
          spatSample(env_mask[[1]], size = n_bg, method = "random", na.rm = TRUE, xy = TRUE)
        }, error = function(e) NULL)
        
        if (is.null(bg_points_iter) || nrow(bg_points_iter) < 50) {
          boyce_scores_iter[iter] <- NA
          next
        }
        
        # Extract environmental data for this iteration
        bg_env_iter <- tryCatch(terra::extract(env, bg_points_iter[, c("x", "y")])[, -1], error = function(e) NULL)
        if (is.null(bg_env_iter)) {
          boyce_scores_iter[iter] <- NA
          next
        }
        
        # Clean data
        bg_env_iter <- bg_env_iter[complete.cases(bg_env_iter), , drop = FALSE]
        if (nrow(bg_env_iter) < 50) {
          boyce_scores_iter[iter] <- NA
          next
        }
        
        # Prepare data for modeling
        p_iter <- c(rep(1, nrow(pres_env)), rep(0, nrow(bg_env_iter)))
        data_all_iter <- rbind(pres_env, bg_env_iter)
        
        # Fit MaxEnt model for this iteration
        model_iter <- tryCatch({
          maxnet(p_iter, data_all_iter, f = maxnet.formula(p_iter, data_all_iter, classes = "lq"))
        }, error = function(e) NULL)
        
        if (is.null(model_iter)) {
          boyce_scores_iter[iter] <- NA
          next
        }
        
        # Predict for this iteration
        pred_vals_iter <- tryCatch(predict(model_iter, newdata = cells_df, type = "cloglog"), error = function(e) NULL)
        if (is.null(pred_vals_iter)) {
          boyce_scores_iter[iter] <- NA
          next
        }
        
        # Create prediction raster for this iteration
        pred_iter <- env_mask[[1]] * NA
        if (length(pred_vals_iter) >= length(non_na_cells)) {
          pred_iter[non_na_cells] <- pred_vals_iter[1:length(non_na_cells)]
        } else {
          boyce_scores_iter[iter] <- NA
          next
        }
        
        # Calculate Boyce for this iteration
        boyce_scores_iter[iter] <- tryCatch({
          score <- calc_boyce(pred_iter, occ_vect)
          if (is.na(score) || is.nan(score)) NA else score
        }, error = function(e) {
          return(NA)
        })
      }
      
      # Compute mean CBI across valid iterations
      valid_scores <- boyce_scores_iter[!is.na(boyce_scores_iter)]
      if (length(valid_scores) > 0) {
        boyce_score <- mean(valid_scores)
        cat("CBI across", length(valid_scores), "valid iterations: mean =", round(boyce_score, 3), 
            "| sd =", round(sd(valid_scores), 3), "\n")
      } else {
        boyce_score <- NA
        cat("All iterations failed for this radius\n")
      }
    }
    
    cat("Radius:", round(r / 1000), "km | Boyce:", round(boyce_score, 3), 
        "| BG points:", nrow(bg_env), "\n")
    
    return(list(index = radius_index, boyce = boyce_score, bg_points = nrow(bg_env)))
  }
  
  # Process all radii
  for (i in seq_along(radii)) {
    r <- radii[i]
    result <- process_radius(r, i)
    
    if (!is.null(result$boyce) && !is.na(result$boyce)) {
      boyce_scores[i] <- result$boyce
      
      # Early stopping criterion (only if not fitting all models)
      if (!fit_all_models) {
        delta <- abs(boyce_scores[i] - prev_boyce)
        if (boyce_scores[i] >= boyce_threshold && delta < stabilisation_delta) {
          cat("Stabilised at radius:", round(r / 1000), "km | Boyce:", round(boyce_scores[i], 3), 
              "| Delta:", round(delta, 4), "\n")
          return(r)
        }
        prev_boyce <- boyce_scores[i]
      }
    }
  }
  
  # Select best radius from valid Boyce scores
  valid_idx <- which(!is.na(boyce_scores) & is.finite(boyce_scores) & boyce_scores > 0)
  
  if (length(valid_idx) > 0) {
    if (fit_all_models) {
      # Find smallest radius where performance stabilizes
      # Look for first radius where BI >= threshold and improvement is minimal
      stable_radius <- NULL
      
      for (i in valid_idx) {
        current_boyce <- boyce_scores[i]
        
        # Check if current radius meets threshold
        if (current_boyce >= boyce_threshold) {
          # Look ahead to check if performance has stabilised
          stabilised <- FALSE
          
          # Check next few radii to confirm stabilisation
          look_ahead <- min(3, length(valid_idx) - which(valid_idx == i) + 1)
          remaining_idx <- valid_idx[valid_idx > i][1:look_ahead]
          remaining_idx <- remaining_idx[!is.na(remaining_idx)]
          
          if (length(remaining_idx) > 0) {
            future_scores <- boyce_scores[remaining_idx]
            # Check if all future improvements are below threshold
            improvements <- future_scores - current_boyce
            if (all(improvements <= stabilisation_delta, na.rm = TRUE)) {
              stabilised <- TRUE
            }
          } else {
            # If no future radii, consider this as stabilised if it meets threshold
            stabilised <- TRUE
          }
          
          if (stabilised) {
            stable_radius <- radii[i]
            cat("Performance stabilised at radius:", round(stable_radius / 1000, 1), 
                "km (Boyce:", round(current_boyce, 3), ")\n")
            break
          }
        }
      }
      
      # If no stable radius found, use the best performing one
      if (is.null(stable_radius)) {
        best_i <- valid_idx[which.max(boyce_scores[valid_idx])]
        stable_radius <- radii[best_i]
        cat("No stabilisation detected; using best performing radius:", round(stable_radius / 1000, 1), 
            "km (Boyce:", round(boyce_scores[best_i], 3), ")\n")
      }
      
      return(stable_radius)
      
    } else {
      # Original logic for early stopping mode
      best_i <- valid_idx[which.max(boyce_scores[valid_idx])]
      best_r <- radii[best_i]
      cat("No early stabilisation; using best:", round(best_r / 1000, 1), 
          "km (Boyce:", round(boyce_scores[best_i], 3), ")\n")
      return(best_r)
    }
  } else {
    # Use enhanced fallback logic
    if (!is.null(range_estimate) && is.finite(range_estimate)) {
      # Use range_estimate if provided, bounded and rounded to nearest 50km
      enhanced_fallback <- min(
        500000,  # max 500km
        max(
          fallback_r,  # minimum fallback_r
          round(range_estimate / 50000) * 50000  # rounded to nearest 50km
        )
      )
      cat("No valid Boyce scores > 0; using range-based fallback:", enhanced_fallback / 1000, "km\n")
      return(enhanced_fallback)
    } else {
      cat("No valid Boyce scores > 0; using provided fallback:", fallback_r / 1000, "km\n")
      return(fallback_r)
    }
  }
}

# ==============================================================================
# MAIN MODELING FUNCTION
# ==============================================================================

#' Compute species distribution model using ENMevaluate with niche breadth metrics
#' 
#' This function performs the full SDM workflow for a single species: data preparation,
#' spatial autocorrelation estimation, CV fold creation, optimal buffer estimation for
#' background sampling, model tuning with ENMeval, best model selection, predictions,
#' niche breadth calculations, extrapolation correction, and niche position estimation.
#' 
#' @param species_name Scientific name of the species.
#' @param occ_env_fst Path to FST file containing occurrence and environmental data.
#' @param env_file Path to environmental raster file.
#' @param n_cores Number of cores for parallel processing.
#' @return data.table with model results and niche metrics, or NULL on failure.
compute_SDM_enmeval <- function(species_name, occ_env_fst, env_file, n_cores) {
  
  # Error handling wrapper
  tryCatch({
    
    # Load environmental data fresh for each call
    env <- rast(env_file)
    
    # ============================================================================
    # 1. PREPARE OCCURRENCE DATA
    # ============================================================================
    cat("1. Loading occurrence data...\n")
    
    # Load occurrence data
    occ_env <- read_fst(occ_env_fst, as.data.table = TRUE)
    sp_dt <- occ_env[scientific_name == species_name, .(x_albers, y_albers)]
    
    # Check minimum records
    if (nrow(sp_dt) < 10) {
      cat("Insufficient records (", nrow(sp_dt), " < 10). Skipping species.\n")
      return(NULL)
    }
    
    cat("Found", nrow(sp_dt), "records\n")
    
    # Create spatial objects
    sp_vect <- vect(sp_dt, geom = c("x_albers", "y_albers"), crs = "EPSG:3577")
    
    # Filter for valid environmental data
    vals <- terra::extract(env, sp_vect)
    valid <- complete.cases(vals)
    
    if (sum(valid) < 10) {
      cat("Insufficient valid records (", sum(valid), " < 10). Skipping species.\n")
      return(NULL)
    }
    
    sp_vect_valid <- vect(sp_dt[valid], geom = c("x_albers", "y_albers"), crs = "EPSG:3577")
    sp_vect_valid$occ <- 1
    
    cat("Valid records:", length(sp_vect_valid), "\n")
    
    # ============================================================================
    # 2. CALCULATE OBSERVED METRICS
    # ============================================================================
    cat("2. Calculating observed metrics...\n")
    
    # Area of Occupancy (AOO)
    AOO <- nrow(unique(sp_dt[valid])) * 100  # 100 km² per grid cell
    
    # Extent of Occurrence (EOO)
    EOO <- tryCatch({
      if (nrow(unique(sp_dt[valid])) >= 3) {
        expanse(convHull(sp_vect_valid), "km")
      } else {
        NA
      }
    }, error = function(e) NA)
    
    cat("AOO:", AOO, "km²\n")
    cat("EOO:", round(EOO, 0), "km²\n")
    
    # ============================================================================
    # 3. ESTIMATE SPATIAL AUTOCORRELATION RANGE
    # ============================================================================
    cat("3. Estimating spatial autocorrelation range...\n")
    
    auto_range <- tryCatch({
      cv_spatial_autocor(
        x = sp_vect_valid,
        r = env,
        column = "occ",
        plot = FALSE
      )
    }, error = function(e) {
      cat("Range estimation failed:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(auto_range)) {
      cat("Cannot proceed without range estimation. Skipping species.\n")
      return(NULL)
    }
    
    range_estimate <- mean(auto_range$range, na.rm = TRUE)
    cat("Range estimate:", round(range_estimate/1000, 1), "km\n")
    
    # ============================================================================
    # 4. CREATE SPATIAL CROSS-VALIDATION FOLDS
    # ============================================================================
    cat("4. Creating spatial CV folds...\n")
    
    k_folds <- min(5, length(unique(auto_range$plots$data$folds)))
    
    spatial_folds <- tryCatch({
      cv_spatial(
        x = sp_vect_valid,
        r = env,
        size = range_estimate,
        k = k_folds,
        column = "occ",
        plot = FALSE
      )
    }, error = function(e) {
      cat("Spatial CV failed:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(spatial_folds)) {
      cat("Cannot proceed without spatial CV. Skipping species.\n")
      return(NULL)
    }
    
    cat("Spatial CV folds:", k_folds, "\n")
    
    # ============================================================================
    # 4a. ESTIMATE OPTIMAL BUFFER FOR BACKGROUND SAMPLING
    # ============================================================================
    cat("4a. Estimating optimal buffer with custom function...\n")
    
    buffer_distance <- estimate_optimal_buffer(
      sp_vect_valid, 
      env, 
      spatial_folds = spatial_folds,
      range_estimate = range_estimate, 
      fallback_km = 250,
      n_iterations = 5,
      stabilisation_delta = 0.02,
      fit_all_models = FALSE
    )
    
    # ============================================================================
    # 5. PREPARE SPATIALLY-BALANCED BACKGROUND POINTS
    # ============================================================================
    cat("5. Sampling background points from optimal buffered area...\n")
    
    # Create buffered sampling area around occurrences
    sp_buffer <- terra::buffer(sp_vect_valid, width = buffer_distance)
    
    # Create mask for buffered area
    buffer_rast <- terra::rasterize(sp_buffer, env[[1]], field = 1)
    buffered_mask <- !is.na(buffer_rast) & !is.na(env[[1]])
    
    # Sample background points from buffered area
    nback <- min(10000, sum(values(buffered_mask), na.rm = TRUE))
    
    bg_points <- tryCatch({
      spatSample(
        mask(env[[1]], buffered_mask), 
        size = nback, 
        method = "random", 
        as.points = TRUE, 
        na.rm = TRUE
      )
    }, error = function(e) {
      cat("Buffered sampling failed:", e$message, "\n")
      NULL
    })
    
    if (is.null(bg_points)) {
      cat("Buffered sampling failed. Skipping species.\n")
      return(NULL)
    }
    
    bg_df <- terra::as.data.frame(bg_points, geom = "xy")[c("x", "y")]
    
    # Assign CV folds to background points (using existing blocks)
    bg_vect <- vect(bg_df, geom = c("x", "y"), crs = crs(env))
    bg_block_ids <- terra::extract(block_rast, bg_vect)$folds
    
    # Handle points outside CV blocks (assign random fold)
    na_indices <- is.na(bg_block_ids)
    if (any(na_indices)) {
      bg_block_ids[na_indices] <- sample(1:k_folds, sum(na_indices), replace = TRUE)
      cat("Assigned", sum(na_indices), "background points outside CV blocks to random folds\n")
    }
    
    cat("Background points:", nrow(bg_df), "\n")
    cat("Background points per fold:", table(bg_block_ids), "\n")
    
    # ============================================================================
    # 6. PREPARE DATA FOR ENMeval
    # ============================================================================
    cat("6. Preparing data for ENMeval...\n")
    
    # Create occurrence dataframe
    occ_df <- terra::as.data.frame(sp_vect_valid, geom = "xy")[c("x", "y")]
    
    # Extract fold assignments for occurrence points
    occ_folds <- spatial_folds$folds_ids
    
    # Combine fold assignments
    all_folds <- list(
      occs.grp = occ_folds,
      bg.grp = bg_block_ids
    )
    
    # ============================================================================
    # 7. RUN ENMevaluate
    # ============================================================================
    cat("7. Running ENMevaluate...\n")
    
    e <- tryCatch({
      ENMevaluate(
        occs = occ_df,
        envs = env,
        algorithm = "maxnet",
        tune.args = list(
          fc = c("L", "LQ"),         # Linear, Linear+Quadratic
          rm = seq(0.5, 4, 0.5)      # Regularisation multipliers
        ),
        bg = bg_df,
        partitions = "user",         # Custom block design CV
        user.grp = all_folds,
        parallel = TRUE,
        numCores = n_cores,
        clamp = TRUE
      )
    }, error = function(e) {
      cat("ENMevaluate failed:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(e)) {
      cat("ENMevaluate returned NULL. Skipping species.\n")
      return(NULL)
    }
    
    cat("ENMevaluate complete. Models evaluated:", nrow(e@results), "\n")
    
    # ============================================================================
    # 8. SELECT BEST MODEL
    # ============================================================================
    cat("8. Selecting best model...\n")
    
    # Round relevant metrics
    e@results$or.10p.avg <- round(e@results$or.10p.avg, 2)
    e@results$auc.val.avg <- round(e@results$auc.val.avg, 2)
    e@results$cbi.val.avg <- round(e@results$cbi.val.avg, 2)
    e@results$delta.AICc <- round(e@results$delta.AICc)
    
    # Define valid models (updated thresholds to match main workflow)
    valid_models <- subset(e@results,
                           or.10p.avg >= 0.05 & or.10p.avg <= 0.3 & auc.val.avg >= 0.6
    )
    
    # Decide which models to select from
    if (nrow(valid_models) == 0) {
      cat("Warning: No models met omission and AUC thresholds; selecting from all models\n")
      candidate_models <- e@results
    } else {
      candidate_models <- valid_models
      cat("Selecting best model from valid subset\n")
    }
    
    # Check if we have any candidate models
    if (nrow(candidate_models) == 0) {
      cat("ERROR: No candidate models available. Skipping species.\n")
      return(NULL)
    }
    
    # Step 1: Lowest omission
    min_omission <- min(candidate_models$or.10p.avg, na.rm = TRUE)
    step1 <- candidate_models[!is.na(candidate_models$or.10p.avg) & 
                                candidate_models$or.10p.avg == min_omission, ]
    
    # Step 2: Highest AUC
    max_auc <- max(step1$auc.val.avg, na.rm = TRUE)
    step2 <- step1[!is.na(step1$auc.val.avg) & 
                     step1$auc.val.avg == max_auc, ]
    
    # Step 3: Highest CBI
    max_cbi <- max(step2$cbi.val.avg, na.rm = TRUE)
    step3 <- step2[!is.na(step2$cbi.val.avg) & 
                     step2$cbi.val.avg == max_cbi, ]
    
    # If CBI is all NA, use step2
    if (nrow(step3) == 0) {
      cat("Warning: CBI values are all NA, using AUC-selected models\n")
      step3 <- step2
    }
    
    # Step 4: First in list
    best_idx <- as.integer(rownames(step3)[1])
    best_model_params <- e@results[best_idx, ]
    
    # Additional validation of best model parameters
    if (is.na(best_model_params$auc.val.avg)) {
      cat("ERROR: Best model has NA AUC. Skipping species.\n")
      return(NULL)
    }
    
    # Output best model
    cat("Best model - Features:", as.character(best_model_params$fc), 
        "| Regularization:", best_model_params$rm, 
        "| AUC:", best_model_params$auc.val.avg,
        "| Omission:", best_model_params$or.10p.avg, "\n")
    
    # ============================================================================
    # 9. EXTRACT MODEL FEATURES AND COEFFICIENTS
    # ============================================================================
    cat("9. Extracting model features and coefficients...\n")
    
    # Get best model and extract coefficients
    best_model_obj <- e@models[[best_idx]]
    
    # Check if model object exists and has coefficients
    if (is.null(best_model_obj) || is.null(best_model_obj$betas)) {
      cat("ERROR: Best model object is NULL or has no coefficients. Skipping species.\n")
      return(NULL)
    }
    
    model_betas <- best_model_obj$betas
    feature_names <- names(model_betas)
    
    cat("Model coefficients (betas):\n")
    print(model_betas)
    
    # Generate formula
    create_maxnet_formula <- function(betas) {
      terms <- names(betas)
      vars <- c("RC1", "RC2", "RC3")
      
      # Categorize terms (excluding hinge terms)
      linear <- terms[terms %in% vars]
      quadratic <- terms[grepl("\\^2", terms)]
      product <- terms[grepl(":", terms)]
      threshold <- terms[grepl("threshold", terms)]
      
      # Combine all non-hinge terms
      all_terms <- c(linear, quadratic, product, threshold)
      
      if (length(all_terms) > 0) {
        return(paste(all_terms, collapse = " + "))
      } else {
        return("No terms")
      }
    }
    
    # Generate the formula
    best_model_formula <- create_maxnet_formula(model_betas)
    cat("Best model formula:\n")
    cat(best_model_formula, "\n")
    
    # ============================================================================
    # 10. MAKE PREDICTIONS
    # ============================================================================
    cat("10. Making predictions...\n")
    
    # Extract predictions
    predictions <- tryCatch({
      eval.predictions(e)[[best_idx]]
    }, error = function(e) {
      cat("ERROR: Could not extract predictions:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(predictions)) {
      cat("ERROR: Predictions are NULL. Skipping species.\n")
      return(NULL)
    }
    
    # Convert predictions to SpatRaster for niche position estimation
    predictions_SpatRaster <- terra::rast(predictions)
    
    # Calculate 10th percentile threshold
    train_coords <- terra::as.data.frame(sp_vect_valid, geom = "xy")[c("x", "y")]
    train_preds <- terra::extract(predictions, train_coords, ID = FALSE)
    threshold_10p <- quantile(train_preds[[1]], 0.10, na.rm = TRUE)
    
    # ============================================================================
    # 11. CALCULATE PREDICTED METRICS
    # ============================================================================
    cat("11. Calculating predicted metrics...\n")
    
    # Geographic B2 predicted
    geo_B2 <- tryCatch({
      suit_vals <- values(predictions_SpatRaster)
      calc.B2(suit_vals)
    }, error = function(e) {
      cat("Error calculating geo_B2:", e$message, "\n")
      return(NA)
    })
    
    # Environmental B2 predicted
    env_B2 <- tryCatch({
      result <- env.breadth.maxent(best_model_obj, env)
      as.numeric(result$env.B2)
    }, error = function(e) {
      cat("Error calculating env_B2:", e$message, "\n")
      return(NA)
    })
    
    cat("Geo. B2 (predicted):", round(geo_B2, 4), "\n")
    cat("Env. B2 (predicted):", round(env_B2, 4), "\n")
    
    # ============================================================================
    # 11a. EXTRAPOLATION DETECTION AND CORRECTION
    # ============================================================================
    cat("11a. Performing extrapolation detection and correction...\n")
    
    # Prepare reference (occurrences with env values)
    sp_reference <- as.data.frame(sp_vect_valid, geom = "xy")
    sp_reference <- cbind(sp_reference, terra::extract(env, sp_reference[, c("x", "y")]))
    
    # Prepare prediction grid (all env cells)
    env_df <- as.data.frame(env, xy = TRUE, na.rm = TRUE)
    
    # Compute ExDet
    sp_exdet <- tryCatch({
      compute_extrapolation(
        samples = sp_reference,
        covariate.names = c("RC1", "RC2", "RC3"),
        prediction.grid = env_df,
        coordinate.system = as.character(crs(env))
      )
    }, error = function(e) {
      cat("Error computing ExDet:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(sp_exdet)) {
      ex_dent <- NA
      predictions_masked <- predictions_SpatRaster  # Fallback, no masking
    } else {
      exdet_vec <- sp_exdet$data$all$ExDet
      ex_dent <- sum(exdet_vec >= 0, na.rm = TRUE) / sum(!is.na(exdet_vec))
      
      # Create mask raster: 1 where ExDet >=0, 0 otherwise
      mask_rast <- env[[1]] * NA
      non_na_cells <- terra::cellFromXY(env[[1]], env_df[, c("x", "y")])
      mask_vals <- ifelse(exdet_vec >= 0, 1, 0)
      mask_rast[non_na_cells] <- mask_vals
      
      # Mask suitability
      predictions_masked <- predictions_SpatRaster * mask_rast
    }
    
    cat("ExDent proportion:", round(ex_dent, 3), "\n")
    
    # Corrected geo B2 (on masked)
    geo_B2_corrected <- tryCatch({
      calc.B2(values(predictions_masked))
    }, error = function(e) {
      cat("Error calculating corrected geo_B2:", e$message, "\n")
      return(NA)
    })
    
    # Corrected env B2 (multiply by ex_dent)
    env_B2_corrected <- env_B2 * ex_dent
    
    cat("Geo. B2 (corrected):", round(geo_B2_corrected, 4), "\n")
    cat("Env. B2 (corrected):", round(env_B2_corrected, 4), "\n")
    
    # ============================================================================
    # 11b. CALCULATE NICHE POSITIONS AND OPTIMA
    # ============================================================================
    cat("11b. Calculating niche positions (weighted means) and optima...\n")
    
    # For RC axes positions (use masked suitability)
    total_suit_masked <- global(predictions_masked, "sum", na.rm = TRUE)[1, 1]
    if (total_suit_masked == 0 || is.na(total_suit_masked)) {
      cat("Error: Total masked suitability is zero or NA. Skipping positions.\n")
      rc1_position <- NA
      rc2_position <- NA
      rc3_position <- NA
    } else {
      rc1_position <- global(env[[1]] * predictions_masked, "sum", na.rm = TRUE)[1, 1] / total_suit_masked
      rc2_position <- global(env[[2]] * predictions_masked, "sum", na.rm = TRUE)[1, 1] / total_suit_masked
      rc3_position <- global(env[[3]] * predictions_masked, "sum", na.rm = TRUE)[1, 1] / total_suit_masked
    }
    
    cat("RC1 position:", round(rc1_position, 3), "\n")
    cat("RC2 position:", round(rc2_position, 3), "\n")
    cat("RC3 position:", round(rc3_position, 3), "\n")
    
    # For longitude/latitude (use masked suitability)
    suit_df <- as.data.frame(predictions_masked, xy = TRUE, na.rm = TRUE)
    names(suit_df) <- c("x_proj", "y_proj", "suit")
    if (nrow(suit_df) == 0 || sum(suit_df$suit) == 0) {
      cat("Error: No valid masked suitability values for geo positions.\n")
      lon_position <- NA
      lat_position <- NA
    } else {
      coords_vect <- terra::vect(suit_df, geom = c("x_proj", "y_proj"), crs = crs(env))
      coords_ll <- project(coords_vect, "EPSG:4326")
      suit_df$lon <- crds(coords_ll)[, 1]
      suit_df$lat <- crds(coords_ll)[, 2]
      lon_position <- sum(suit_df$lon * suit_df$suit) / sum(suit_df$suit)
      lat_position <- sum(suit_df$lat * suit_df$suit) / sum(suit_df$suit)
    }
    
    cat("Longitude position:", round(lon_position, 3), "\n")
    cat("Latitude position:", round(lat_position, 3), "\n")
    
    rc1_resp <- maxnet::response.plot(best_model_obj, v = "RC1", type = "cloglog", plot = FALSE)
    rc1_optimum <- rc1_resp$RC1[which.max(rc1_resp$pred)]
    rc1_opt_type <- if (rc1_optimum == min(rc1_resp$RC1)) "min_boundary" else if (rc1_optimum == max(rc1_resp$RC1)) "max_boundary" else "interior"
    
    rc2_resp <- maxnet::response.plot(best_model_obj, v = "RC2", type = "cloglog", plot = FALSE)
    rc2_optimum <- rc2_resp$RC2[which.max(rc2_resp$pred)]
    rc2_opt_type <- if (rc2_optimum == min(rc2_resp$RC2)) "min_boundary" else if (rc2_optimum == max(rc2_resp$RC2)) "max_boundary" else "interior"
    
    rc3_resp <- maxnet::response.plot(best_model_obj, v = "RC3", type = "cloglog", plot = FALSE)
    rc3_optimum <- rc3_resp$RC3[which.max(rc3_resp$pred)]
    rc3_opt_type <- if (rc3_optimum == min(rc3_resp$RC3)) "min_boundary" else if (rc3_optimum == max(rc3_resp$RC3)) "max_boundary" else "interior"
    
    cat("RC1 optimum:", round(rc1_optimum, 3), " (type: ", rc1_opt_type, ")\n")
    cat("RC2 optimum:", round(rc2_optimum, 3), " (type: ", rc2_opt_type, ")\n")
    cat("RC3 optimum:", round(rc3_optimum, 3), " (type: ", rc3_opt_type, ")\n")
    
    # ============================================================================
    # 12. COMPILE RESULTS
    # ============================================================================
    cat("12. Compiling results...\n")
    
    # Create results data.table
    results <- data.table(
      species = sp,
      n_occ = length(sp_vect_valid),
      auc = best_model_params$auc.val.avg,
      omission = best_model_params$or.10p.avg,
      cbi = best_model_params$cbi.val.avg,
      AOO = AOO,
      EOO = EOO,
      env_B2 = env_B2,
      geo_B2 = geo_B2,
      ex_dent = ex_dent,
      env_B2_corrected = env_B2_corrected,
      geo_B2_corrected = geo_B2_corrected,
      threshold_10p = threshold_10p,
      regularisation = best_model_params$rm,
      features = best_model_params$fc,
      formula = best_model_formula,
      n_coefs = best_model_params$ncoef,
      coefficients = paste(names(model_betas), "=", round(model_betas, 2), collapse = ";"),
      spatial_range_km = round(range_estimate/1000),
      cv_folds = k_folds,
      bg_points = nrow(bg_df),
      RC1_position = rc1_position,
      RC2_position = rc2_position,
      RC3_position = rc3_position,
      lon_position = lon_position,
      lat_position = lat_position,
      RC1_optimum = rc1_optimum,
      RC2_optimum = rc2_optimum,
      RC3_optimum = rc3_optimum,
      RC1_opt_type = rc1_opt_type,
      RC2_opt_type = rc2_opt_type,
      RC3_opt_type = rc3_opt_type
    )
    
    # Add validation criteria
    results[, pass_validation := !is.na(auc) & !is.na(omission) & 
              omission >= 0.05 & omission <= 0.30 & auc >= 0.6]
    
    cat("Species", species_name, "completed successfully\n")
    cat("Model performance: AUC =", round(results$auc, 2), 
        "| Omission =", round(results$omission, 2), 
        "| Pass validation =", results$pass_validation, "\n")
    
    return(results)
    
  }, error = function(e) {
    cat("ERROR processing species", species_name, ":", e$message, "\n")
    return(NULL)
  })
}

# ==============================================================================
# WRAPPER FUNCTION FOR BATCH PROCESSING
# ==============================================================================

#' Wrapper function for processing multiple species
#' 
#' This function processes a list of species in batch mode, calling compute_SDM_enmeval
#' for each and combining the results into a single data.table. Handles NULL results
#' by removing them from the final output.
#' 
#' @param species_list Vector of species names to process.
#' @param occ_env_fst Path to FST file containing occurrence data.
#' @param env_file Path to environmental raster file.
#' @return data.table with combined results from all successful species models.
process_species_batch <- function(species_list, occ_env_fst, env_file) {
  
  cat("Processing batch of", length(species_list), "species\n")
  
  # Process each species
  results_list <- lapply(species_list, function(sp) {
    compute_SDM_enmeval(sp, occ_env_fst, env_file, n_cores)
  })
  
  # Remove NULL results
  results_list <- results_list[!sapply(results_list, is.null)]
  
  if (length(results_list) > 0) {
    # Combine results
    final_results <- rbindlist(results_list, fill = TRUE)
    cat("Batch processing complete:", nrow(final_results), "successful models\n")
    return(final_results)
  } else {
    cat("No successful models in batch\n")
    return(NULL)
  }
}

cat("ENMevaluate modeling functions loaded successfully\n")