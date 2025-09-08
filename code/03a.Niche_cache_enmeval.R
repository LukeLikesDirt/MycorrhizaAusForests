# 03a.Niche_cache_enmeval.R
# ==============================================================================
# CACHE OCCURRENCE AND ENVIRONMENTAL DATA FOR ENMevaluate WORKFLOW
# ==============================================================================
# This script prepares and caches occurrence data with environmental predictor
# values for species distribution modeling using ENMevaluate. It's adapted from
# the original ENMTools caching workflow but optimized for ENMevaluate requirements.
#
# Run in caffeinate from R project home directory:
# > caffeinate -dims Rscript code/03a.Niche_cache_enmeval.R
#
# ==============================================================================

# Load required libraries
library(terra)
library(data.table)
library(fst)

cat("=== ENMevaluate Data Caching Workflow ===\n")
cat("Starting data preparation and caching...\n")

# ==============================================================================
# 1. LOAD AND FILTER OCCURRENCE DATA
# ==============================================================================
cat("\n--- Loading occurrence data ---\n")

# Load presence points with at least 10 records per species
# This ensures we have enough data for robust model fitting and cross-validation
occ_data <- fread("data/presence/trees_10.txt")
cat("Initial records loaded:", nrow(occ_data), "\n")

# Filter to keep only species with at least 10 records
# This is critical for ENMevaluate as it requires sufficient data for k-fold CV
occ_data_filtered <- occ_data[, .SD[.N >= 10], by = scientific_name]
cat("Records after filtering (>=10 per species):", nrow(occ_data_filtered), "\n")
cat("Number of species:", length(unique(occ_data_filtered$scientific_name)), "\n")

# Create spatial vector object for coordinate reference system
# Using Australian Albers projection (EPSG:3577) for accurate area calculations
vect_pts <- vect(occ_data_filtered, 
                 geom = c("x_albers", "y_albers"), 
                 crs = "EPSG:3577")

cat("Spatial vector created with", length(vect_pts), "points\n")

# ==============================================================================
# 2. LOAD ENVIRONMENTAL PREDICTORS
# ==============================================================================
cat("\n--- Loading environmental predictors ---\n")

# Load environmental predictors optimized for forest species
# Using reduced complexity model with 3 principal components for computational efficiency
env <- rast("data/aus_forests_23/predictors_10.tif")
cat("Environmental raster loaded with", nlyr(env), "layers\n")
cat("Raster dimensions:", dim(env), "\n")
cat("Raster resolution:", res(env), "m\n")

# Select only the first 3 principal components for ENMevaluate
# This reduces model complexity while maintaining most environmental variation
env_subset <- env[[c("RC1", "RC2", "RC3")]]
cat("Using environmental layers:", names(env_subset), "\n")

# Set temporary directory for terra operations (use SSD if available)
# This is important for large raster operations
if (Sys.getenv("R_TEMP_DIR") != "") {
  terraOptions(tempdir = Sys.getenv("R_TEMP_DIR"))
  cat("Using custom temp directory:", Sys.getenv("R_TEMP_DIR"), "\n")
} else {
  cat("Using default temp directory\n")
}

# ==============================================================================
# 3. EXTRACT ENVIRONMENTAL VALUES
# ==============================================================================
cat("\n--- Extracting environmental values ---\n")

# Extract environmental values at all occurrence points
# This is the most computationally intensive step
cat("Extracting environmental values for", length(vect_pts), "points...\n")
vals <- extract(env_subset, vect_pts)
cat("Environmental extraction complete\n")

# Check for missing values
missing_vals <- sum(!complete.cases(vals))
cat("Records with missing environmental data:", missing_vals, "\n")
cat("Percentage missing:", round(missing_vals/nrow(vals)*100, 2), "%\n")

# Combine occurrence data with environmental values
occ_env <- cbind(occ_data_filtered, vals)
cat("Combined dataset dimensions:", dim(occ_env), "\n")

# Set data.table key for efficient species-based operations
setkey(occ_env, scientific_name)

# ==============================================================================
# 4. DATA QUALITY CHECKS
# ==============================================================================
cat("\n--- Data quality checks ---\n")

# Check species with complete environmental data
species_complete <- occ_env[complete.cases(occ_env), .N, by = scientific_name]
species_incomplete <- occ_env[!complete.cases(occ_env), .N, by = scientific_name]

cat("Species with complete environmental data:", nrow(species_complete), "\n")
cat("Records with complete environmental data:", sum(species_complete$N), "\n")

if (nrow(species_incomplete) > 0) {
  cat("Species with some missing environmental data:", nrow(species_incomplete), "\n")
  cat("Records with missing environmental data:", sum(species_incomplete$N), "\n")
}

# Check for duplicate coordinates within species
duplicates <- occ_env[, .N, by = .(scientific_name, x_albers, y_albers)][N > 1]
if (nrow(duplicates) > 0) {
  cat("Warning: Found", nrow(duplicates), "duplicate coordinate pairs\n")
} else {
  cat("No duplicate coordinates detected\n")
}

# Summary statistics for environmental variables
cat("\n--- Environmental variable summary ---\n")
for (var in c("RC1", "RC2", "RC3")) {
  var_summary <- summary(occ_env[[var]])
  cat("Variable", var, ":\n")
  cat("  Range:", round(var_summary[1], 3), "to", round(var_summary[6], 3), "\n")
  cat("  Mean:", round(var_summary[4], 3), "±", round(sd(occ_env[[var]], na.rm = TRUE), 3), "\n")
  cat("  Missing:", sum(is.na(occ_env[[var]])), "\n")
}

# ==============================================================================
# 5. CACHE PREPARED DATA
# ==============================================================================
cat("\n--- Caching prepared data ---\n")

# Create cache directory if it doesn't exist
if (!dir.exists("cache")) {
  dir.create("cache")
  cat("Created cache directory\n")
}

# Save to FST format with maximum compression
# FST is faster than RDS for large datasets and preserves data.table structure
cache_file <- "cache/occ_with_env_enmeval.fst"
write_fst(occ_env, cache_file, compress = 100)

# Verify cache file was created
if (file.exists(cache_file)) {
  file_size <- file.info(cache_file)$size / 1024^2  # Convert to MB
  cat("Cache file created successfully:", cache_file, "\n")
  cat("File size:", round(file_size, 2), "MB\n")
} else {
  stop("Failed to create cache file")
}

# ==============================================================================
# 6. CACHE ENVIRONMENTAL RASTER SUBSET
# ==============================================================================
cat("\n--- Caching environmental raster subset ---\n")

# Save the subset environmental raster for consistent use across the workflow
env_cache_file <- "cache/env_predictors_enmeval.tif"
writeRaster(env_subset, env_cache_file, overwrite = TRUE)

if (file.exists(env_cache_file)) {
  file_size <- file.info(env_cache_file)$size / 1024^2  # Convert to MB
  cat("Environmental raster cached:", env_cache_file, "\n")
  cat("File size:", round(file_size, 2), "MB\n")
} else {
  stop("Failed to cache environmental raster")
}

# ==============================================================================
# 7. CREATE METADATA SUMMARY
# ==============================================================================
cat("\n--- Creating metadata summary ---\n")

# Create summary metadata for the cached dataset
metadata <- list(
  cache_date = Sys.time(),
  total_records = nrow(occ_env),
  total_species = length(unique(occ_env$scientific_name)),
  complete_records = sum(complete.cases(occ_env)),
  environmental_layers = names(env_subset),
  coordinate_system = "EPSG:3577",
  min_records_per_species = 10,
  cache_files = c(cache_file, env_cache_file)
)

# Save metadata
metadata_file <- "cache/enmeval_cache_metadata.rds"
saveRDS(metadata, metadata_file)

# Print final summary
cat("\n=== Caching Complete ===\n")
cat("Dataset summary:\n")
cat("  Total records:", metadata$total_records, "\n")
cat("  Total species:", metadata$total_species, "\n")
cat("  Complete records:", metadata$complete_records, "\n")
cat("  Environmental layers:", paste(metadata$environmental_layers, collapse = ", "), "\n")
cat("  Cache files created:\n")
for (file in metadata$cache_files) {
  cat("    -", file, "\n")
}
cat("  Metadata file:", metadata_file, "\n")

cat("\nData caching workflow completed successfully!\n")
cat("Cached data is ready for ENMevaluate species distribution modeling.\n")