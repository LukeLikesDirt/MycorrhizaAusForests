# 03c.Niche_run_enmeval.R
# ==============================================================================
# MAIN WORKFLOW FOR RUNNING ENMevaluate SPECIES DISTRIBUTION MODELS
# ==============================================================================
# This script runs the complete ENMevaluate workflow for species distribution
# modeling, including data loading, model fitting, and results compilation.
# It processes species in batches to manage memory usage and enable parallel processing.
#
# Run in caffeinate from R project home directory:
# > caffeinate -dims Rscript code/03c.Niche_run_enmeval.R
#
# ==============================================================================

cat("=== ENMevaluate Species Distribution Modeling Workflow ===\n")
cat("Starting batch processing of species distribution models...\n")

# ==============================================================================
# 1. LOAD REQUIRED LIBRARIES AND FUNCTIONS
# ==============================================================================
cat("\n--- Loading libraries and functions ---\n")

# Load required libraries
library(dplyr)
library(future.apply)
library(fst)
library(terra)
library(data.table)
library(ENMeval)
library(blockCV)

# Source the main modelling functions
source("code/03b.Niche_est_enmeval.R")

# ==============================================================================
# 2. SETUP DIRECTORIES AND PATHS
# ==============================================================================
cat("\n--- Setting up directories and paths ---\n")

# Create results directory
results_dir <- "data/niche_estimates_enmeval/"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  cat("Created results directory:", results_dir, "\n")
}

# Create output directories for individual species results (optional)
species_output_dir <- "output/species_enmeval/"
if (!dir.exists(species_output_dir)) {
  dir.create(species_output_dir, recursive = TRUE)
  cat("Created species output directory:", species_output_dir, "\n")
}

cat("Results will be written to:", normalizePath(results_dir), "\n")

# ==============================================================================
# 3. CONFIGURE SYSTEM SETTINGS
# ==============================================================================
cat("\n--- Configuring system settings ---\n")

# Set terra temp directory to SSD if available
if (Sys.getenv("R_TEMP_DIR") != "") {
  terraOptions(tempdir = Sys.getenv("R_TEMP_DIR"))
  cat("Using custom temp directory:", Sys.getenv("R_TEMP_DIR"), "\n")
} else {
  # Try to use a fast temp directory
  if (dir.exists("/Users/lukelikesdirt/ssd_tmp")) {
    terraOptions(tempdir = "/Users/lukelikesdirt/ssd_tmp")
    cat("Using SSD temp directory: /Users/lukelikesdirt/ssd_tmp\n")
  } else {
    cat("Using default temp directory\n")
  }
}

# Set random seed for reproducibility
set.seed(1986)
cat("Random seed set to 1986 for reproducibility\n")

# ==============================================================================
# 4. LOAD DATA
# ==============================================================================
cat("\n--- Loading cached data ---\n")

# Define data paths
occ_fst <- "cache/occ_with_env_enmeval.fst"
env_file <- "cache/env_predictors_enmeval.tif"

# Verify cache files exist
if (!file.exists(occ_fst)) {
  stop("Occurrence cache file not found: ", occ_fst, 
       "\nPlease run 03a.Niche_cache_enmeval.R first")
}

if (!file.exists(env_file)) {
  stop("Environmental cache file not found: ", env_file, 
       "\nPlease run 03a.Niche_cache_enmeval.R first")
}

# Load environmental data to check dimensions (don't keep in memory for parallel)
cat("Checking environmental predictors...\n")
env_check <- rast(env_file)
cat("Environmental layers available:", names(env_check), "\n")
cat("Raster dimensions:", dim(env_check), "\n")
rm(env_check)  # Remove from memory

# Load occurrence data to get species list
cat("Loading occurrence data...\n")
occ_env <- read_fst(occ_fst, as.data.table = TRUE)
cat("Occurrence records loaded:", nrow(occ_env), "\n")

# Get species list with sufficient records
species_counts <- occ_env[, .N, by = scientific_name]
all_species <- species_counts[N >= 10, scientific_name]
# Random sample for testing:
#all_species <- sample(all_species, 10)
cat("Species with ≥10 records:", length(all_species), "\n")

# Clean up occurrence data from memory (will be loaded fresh in each function)
rm(occ_env)
gc()

# ==============================================================================
# 5. CONFIGURE BATCH PROCESSING
# ==============================================================================
cat("\n--- Configuring batch processing ---\n")

# Define batch size (reduced for ENMevaluate complexity and memory management)
batch_size <- 25  # Smaller batches for stability
total_batches <- ceiling(length(all_species) / batch_size)

# Create species batches
species_batches <- split(all_species, ceiling(seq_along(all_species) / batch_size))

cat("Total species to process:", length(all_species), "\n")
cat("Batch size:", batch_size, "\n")
cat("Number of batches:", total_batches, "\n")

# ==============================================================================
# 6. SETUP PARALLEL PROCESSING
# ==============================================================================
cat("\n--- Setting up parallel processing ---\n")

# Configure parallel processing
# Use sequential processing to avoid external pointer issues
# This is more stable for terra objects in ENMevaluate
cat("Using sequential processing for stability with terra objects\n")

# Set up sequential plan (no parallel processing to avoid terra issues)
plan(sequential)

# ==============================================================================
# 7. PROCESS SPECIES IN BATCHES
# ==============================================================================
cat("\n--- Processing species in batches ---\n")

# Track processing time
start_time <- Sys.time()

# Initialize species counter for progress tracking
species_counter <- 0
total_species <- length(all_species)

# Process each batch sequentially
for (i in seq_along(species_batches)) {
  
  batch_start <- Sys.time()
  cat("Processing batch", i, "of", length(species_batches), "\n")
  
  # Get species for this batch
  sp_batch <- species_batches[[i]]
  cat("Batch", i, "contains", length(sp_batch), "species\n")
  
  # Process species in batch with progress tracking
  results <- tryCatch({
    lapply(sp_batch, function(sp) {
      # Increment counter and show progress
      species_counter <<- species_counter + 1
      cat("\n=== Processing species:", sp, "(", species_counter, "of", total_species, ") ===\n")
      
      # Call the modelling function
      compute_SDM_enmeval(sp, occ_fst, env_file, n_cores = 10)
    })
  }, error = function(e) {
    cat("ERROR in batch", i, ":", e$message, "\n")
    return(list())
  })
  
  # Remove NULL results (failed species)
  results <- results[!sapply(results, is.null)]
  
  if (length(results) > 0) {
    # Combine results into data.table
    batch_dt <- rbindlist(results, fill = TRUE)
    
    # Save batch results
    batch_file <- sprintf("%s/batch_%02d.fst", results_dir, i)
    write_fst(batch_dt, batch_file)
    
    # Calculate processing time
    batch_time <- as.numeric(difftime(Sys.time(), batch_start, units = "mins"))
    
    cat("Batch", i, "completed successfully\n")
    cat("  Processed:", nrow(batch_dt), "species\n")
    cat("  Processing time:", round(batch_time, 2), "minutes\n")
    cat("  Saved to:", batch_file, "\n")
    
    # Clean up memory
    rm(results, batch_dt)
    gc()
    
  } else {
    cat("Batch", i, "failed - no successful models\n")
  }
}

# Calculate total processing time
total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("All batches completed in", round(total_time, 2), "minutes\n")

# ==============================================================================
# 8. AGGREGATE RESULTS
# ==============================================================================
cat("\n--- Aggregating batch results ---\n")

# Read all batch files
batch_files <- list.files(results_dir, pattern = "batch_.*\\.fst$", full.names = TRUE)
cat("Found", length(batch_files), "batch files\n")

if (length(batch_files) > 0) {
  # Load and combine all batch results
  all_results <- rbindlist(lapply(batch_files, function(file) {
    tryCatch({
      read_fst(file, as.data.table = TRUE)
    }, error = function(e) {
      cat("Error reading", file, ":", e$message, "\n")
      return(NULL)
    })
  }), fill = TRUE)
  
  # Remove any NULL results
  all_results <- all_results[!is.null(all_results)]
  
  if (nrow(all_results) > 0) {
    # Save aggregated results
    final_file <- paste0(results_dir, "niche_model_results_enmeval.fst")
    write_fst(all_results, final_file)
    
    cat("Aggregated results saved to:", final_file, "\n")
    cat("Total successful models:", nrow(all_results), "\n")
  } else {
    cat("No successful results to aggregate\n")
  }
} else {
  cat("No batch files found for aggregation\n")
}

# ==============================================================================
# 9. GENERATE SUMMARY STATISTICS
# ==============================================================================
if (exists("all_results") && nrow(all_results) > 0) {
  cat("\n--- Generating summary statistics ---\n")
  
  # Model performance summary
  cat("Model Performance Summary:\n")
  cat("  Mean AUC:", round(mean(all_results$auc, na.rm = TRUE), 3), "\n")
  cat("  Median AUC:", round(median(all_results$auc, na.rm = TRUE), 3), "\n")
  cat("  SD AUC:", round(sd(all_results$auc, na.rm = TRUE), 3), "\n")
  cat("  Mean Omission:", round(mean(all_results$omission, na.rm = TRUE), 3), "\n")
  cat("  Median Omission:", round(median(all_results$omission, na.rm = TRUE), 3), "\n")
  cat("  SD Omission:", round(sd(all_results$omission, na.rm = TRUE), 3), "\n")
  cat("  Mean CBI:", round(mean(all_results$cbi, na.rm = TRUE), 3), "\n")
  cat("  Median CBI:", round(median(all_results$cbi, na.rm = TRUE), 3), "\n")
  cat("  SD CBI:", round(sd(all_results$cbi, na.rm = TRUE), 3), "\n")
  cat("  Models passing validation:", sum(all_results$pass_validation, na.rm = TRUE), 
      "(", round(sum(all_results$pass_validation, na.rm = TRUE)/nrow(all_results)*100, 1), "%)\n")
  
  # Feature class usage
  cat("\nFeature Class Usage:\n")
  fc_table <- table(all_results$features)
  for (fc in names(fc_table)) {
    cat("  ", fc, ":", fc_table[fc], "(", round(fc_table[fc]/nrow(all_results)*100, 1), "%)\n")
  }

  # Niche breadth summary
  cat("\nNiche Breadth Summary:\n")
  cat("  Mean Geo. B2:", round(mean(all_results$geo_B2, na.rm = TRUE), 4), "\n")
  cat("  Median Geo. B2:", round(median(all_results$geo_B2, na.rm = TRUE), 4), "\n")
  cat("  SD Geo. B2:", round(sd(all_results$geo_B2, na.rm = TRUE), 4), "\n\n")
  cat("  Mean Geo. B2 (corrected):", round(mean(all_results$geo_B2_corrected, na.rm = TRUE), 4), "\n")
  cat("  Median Geo. B2 (corrected):", round(median(all_results$geo_B2_corrected, na.rm = TRUE), 4), "\n")
  cat("  SD Geo. B2 (corrected):", round(sd(all_results$geo_B2_corrected, na.rm = TRUE), 4), "\n\n")
  cat("  Mean Env. B2:", round(mean(all_results$env_B2, na.rm = TRUE), 4), "\n")
  cat("  Median Env. B2:", round(median(all_results$env_B2, na.rm = TRUE), 4), "\n")
  cat("  SD Env. B2:", round(sd(all_results$env_B2, na.rm = TRUE), 4), "\n\n")
  cat("  Mean Env. B2 (corrected):", round(mean(all_results$env_B2_corrected, na.rm = TRUE), 4), "\n")
  cat("  Median Env. B2 (corrected):", round(median(all_results$env_B2_corrected, na.rm = TRUE), 4), "\n")
  cat("  SD Env. B2 (corrected):", round(sd(all_results$env_B2_corrected, na.rm = TRUE), 4), "\n\n")
  cat("  Mean RC1 position:", round(mean(all_results$RC1_position, na.rm = TRUE), 4), "\n")
  cat("  Median RC1 position:", round(median(all_results$RC1_position, na.rm = TRUE), 4), "\n")
  cat("  SD RC1 position:", round(sd(all_results$RC1_position, na.rm = TRUE), 4), "\n\n")
  cat("  Mean RC2 position:", round(mean(all_results$RC2_position, na.rm = TRUE), 4), "\n")
  cat("  Median  RC2 position:", round(median(all_results$RC2_position, na.rm = TRUE), 4), "\n\n")
  cat("  SD RC2 position:", round(sd(all_results$RC2_position, na.rm = TRUE), 4), "\n")
  cat("  Mean RC3 position:", round(mean(all_results$RC3_position, na.rm = TRUE), 4), "\n")
  cat("  Median RC3 position:", round(median(all_results$RC3_position, na.rm = TRUE), 4), "\n")
  cat("  SD RC3 position:", round(sd(all_results$RC3_position, na.rm = TRUE), 4), "\n\n")
  
  # Save summary statistics
  summary_stats <- list(
    total_species_processed = nrow(all_results),
    processing_time_minutes = total_time,
    mean_auc = mean(all_results$auc, na.rm = TRUE),
    median_auc = median(all_results$auc, na.rm = TRUE),
    sd_auc = sd(all_results$auc, na.rm = TRUE),
    mean_omission = mean(all_results$omission, na.rm = TRUE),
    median_omission = median(all_results$omission, na.rm = TRUE),
    sd_omission = sd(all_results$omission, na.rm = TRUE),
    mean_cbi = mean(all_results$cbi, na.rm = TRUE),
    median_cbi = median(all_results$cbi, na.rm = TRUE),
    sd_cbi = sd(all_results$cbi, na.rm = TRUE),
    models_passing_validation = sum(all_results$pass_validation, na.rm = TRUE),
    completion_time = Sys.time()
  )
  
  saveRDS(summary_stats, paste0(results_dir, "workflow_summary.rds"))
  cat("Summary statistics saved to:", paste0(results_dir, "workflow_summary.rds"), "\n")
}

# ==============================================================================
# 10. CLEANUP AND FINAL MESSAGES
# ==============================================================================
cat("\n--- Workflow cleanup ---\n")

# Clean memory
gc()

# Final summary
cat("\n=== ENMevaluate Workflow Complete ===\n")
cat("Processing started:", format(start_time), "\n")
cat("Processing completed:", format(Sys.time()), "\n")
cat("Total processing time:", round(total_time, 2), "minutes\n")

if (exists("all_results") && nrow(all_results) > 0) {
  cat("Successfully processed:", nrow(all_results), "species\n")
  cat("Models passing validation:", sum(all_results$pass_validation, na.rm = TRUE), "\n")
  cat("Results saved to:", results_dir, "\n")
} else {
  cat("No successful models generated\n")
}

cat("ENMevaluate workflow completed!\n")