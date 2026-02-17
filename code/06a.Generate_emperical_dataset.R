# ============================================================================== #
# MYCORRHIZAL TYPE CLASSIFICATION FOR AUSTRALIAN NATIVE SPECIES
# ============================================================================== #
# This script classifies plant species into mycorrhizal types based on empirical
# fungal root observations. Classifications follow these rules:
#
# AM:   Species with AM type, even if found NM elsewhere (facultative AM)
# EcM:  Species with EcM type, even if found NM elsewhere (facultative EcM)
# Dual: Species with BOTH AM and EcM types. AM must be observed in established
#       individuals (many EcM plants form AM associations as seedlings only)
# NM:   Species with ONLY NM type across ≥2 observations (to avoid misclassifying
#       facultative mycorrhizal species). Exception: 1 observation is sufficient
#       if cluster roots are present.
# ============================================================================== #

# Required packages
require(data.table)
require(tidyverse)

# ============================================================================== #
# (1) LOAD AND PREPARE DATA ----------------------------------------------------
# ============================================================================== #

# Load fungal root data and filter to native species only
fungal_root_emp <- fread("data/mycorrhizal_type_emperical/fungal_root_curated.txt") %>%
  filter(native_status == "native") %>%
  select(
    family, genus, species, mycorrhizal_type, host_age_group, 
    am_evaluated, ecm_evaluated, am_structures, ecm_structures, nm_structures
  )

# Load environmental breadth estimates
niche_breadth <- fread("data/niche_estimates_enmeval/niche_estimates.txt")

# ============================================================================== #
# (2) IDENTIFY POTENTIAL DUAL MYCORRHIZAL GENERA -------------------------------
# ============================================================================== #
# The dual mycorrhizal trait is generally conserved within genera. We need to
# identify these genera to apply stricter evaluation criteria (i.e. require
# AM in established individuals).

dual_genera <- fread("data/FungalRoot/harmonised_fungal_root_aus_trees.txt") %>%
  filter(
    # Only include genera present in our empirical dataset
    genus %in% fungal_root_emp$genus,
    # Filter to genera known to have the dual mycorrhizal trait
    mycorrhizal_type == "EcM-AM"
  ) %>%
  pull(genus) %>%
  unique() %>%
  # Add eucalypt genera (also known to have both AM and EcM species)
  c("Angophora", "Corymbia", "Eucalyptus")

# ============================================================================== #
# (3) VECTORS EVALUATED AM & ECM TYPES -----------------------------------------
# ============================================================================== #
# For dual genera, we need to know which species were actually examined for
# mycorrhizal structures to avoid false negatives.

# Species evaluated for AM structures (in ESTABLISHED individuals from dual genera)
# Note: We only track AM evaluation in established plants because seedlings of
# EcM species often have temporally transient AM colonisation

# Species with AM observed in ESTABLISHED individuals (for dual genera checks)
am_evaluated_species <- fungal_root_emp %>%
  filter(
    am_evaluated == "yes",
    host_age_group == "established"
  ) %>%
  select(genus, species, mycorrhizal_type) %>%
  unique()

# Species evaluated for EcM structures (any age)
ecm_evaluated_species <- fungal_root_emp %>%
  filter(ecm_evaluated == "yes") %>%
  select(genus, species, mycorrhizal_type) %>%
  unique()

# Species that were evaluated for AM and found to HAVE AM structures
am_eval_am_species <- am_evaluated_species %>%
  filter(mycorrhizal_type == "AM" | mycorrhizal_type == "EcM-AM") %>%
  pull(species)

# Species that were evaluated for EcM and found to HAVE EcM structures
ecm_eval_ecm_species <- ecm_evaluated_species %>%
  filter(mycorrhizal_type == "EcM" | mycorrhizal_type == "EcM-AM") %>%
  pull(species)

# ============================================================================== #
# (4) SPECIES VECTORS BY MYCORRHIZAL TYPE --------------------------------------
# ============================================================================== #

# All species observed with AM (ANY age - for general AM presence)
am_species_any_age <- fungal_root_emp %>%
  filter(mycorrhizal_type == "AM") %>%
  pull(species) %>%
  unique()

# Species with AM observed in ESTABLISHED individuals only
# (used for dual genera checks and EcM exclusion)
am_species_established <- fungal_root_emp %>%
  filter(
    mycorrhizal_type == "AM",
    host_age_group == "established"
  ) %>%
  pull(species) %>%
  unique()

# All species observed with EcM (regardless of other types)
ecm_species <- fungal_root_emp %>%
  # If EcM-AM is NOT in an established individual, mutate to EcM
  # (AM must be observed in established individuals, but EcM at any age)
  mutate(
    mycorrhizal_type = ifelse(
      mycorrhizal_type == "EcM-AM" & host_age_group != "established",
      "EcM",
      mycorrhizal_type
    )
  ) %>%
  filter(mycorrhizal_type == "EcM") %>%
  pull(species) %>%
  unique()

# Species classified as DUAL mycorrhizal (have both AM and EcM)
# Important: AM must be observed in ESTABLISHED individuals
dual_species <- fungal_root_emp %>%
  filter(
    # Condition 1: Explicitly coded as dual type IN ESTABLISHED individuals
    (mycorrhizal_type == "EcM-AM" & host_age_group == "established") |
      # Condition 2: Found with both AM (in established, when evaluated) AND 
      # EcM (any age, when evaluated) across different observations
      (species %in% am_eval_am_species & species %in% ecm_eval_ecm_species)
  ) %>%
  pull(species) %>%
  unique()

# All species observed as NM (regardless of other types)
nm_species <- fungal_root_emp %>%
  filter(mycorrhizal_type == "NM") %>%
  pull(species) %>%
  unique()

# ============================================================================== #
# (5) CLASSIFY SPECIES MYCORRHIZAL TYPES ---------------------------------------
# ============================================================================== #

#### AM ####
# Species are AM-only if:
# 1. They have mycorrhizal_type == "AM" (at least once)
# 2. They have NOT been observed with EcM
# 3. They are NOT dual species
# 4. For non-dual genera: just need one AM observation
# 5. For dual genera: must have AM in established plants, been evaluated for 
#    EcM, AND have ≥2 AM observations (to avoid false positives)

# Count AM observations per species (for dual genera check)
am_obs_counts <- fungal_root_emp %>%
  filter(mycorrhizal_type == "AM") %>%
  group_by(species) %>%
  summarise(am_obs_count = n()) %>%
  ungroup()

# Create FungalRoot dataset for AM species
fungal_root_am <- fungal_root_emp %>%
  # Keep only species observed with AM but not with EcM or dual
  filter(
    species %in% am_species_any_age &  # ← Use ANY age for general AM presence
      !species %in% ecm_species & 
      !species %in% dual_species
  ) %>%
  # Add observation counts
  left_join(am_obs_counts, by = "species") %>%
  # For species in dual genera, require:
  # - AM was observed in established individuals (for that species)
  # - Species was evaluated for EcM
  # - At least 2 AM observations
  filter(
    # If NOT in dual genus: keep all
    !genus %in% dual_genera |
      # If IN dual genus: require established AM, EcM evaluation, AND ≥2 observations
      (species %in% am_species_established & 
         species %in% ecm_evaluated_species$species &
         am_obs_count >= 2)
  ) %>%
  # Keep only AM observations
  filter(mycorrhizal_type == "AM") %>%
  select(family, genus, species, mycorrhizal_type) %>%
  unique()

#### EcM ####
# Species are EcM-only if:
# 1. They have mycorrhizal_type == "EcM" (at least once)
# 2. They have NOT been observed with AM (in established individuals)
# 3. They are NOT dual species
# 4. For non-dual genera: just need one EcM observation
# 5. For dual genera: must have been evaluated for AM in established plants
#    AND have ≥2 EcM observations (to avoid false positives)

# Count EcM observations per species (for dual genera check)
ecm_obs_counts <- fungal_root_emp %>%
  filter(mycorrhizal_type == "EcM" | 
           (mycorrhizal_type == "EcM-AM" & host_age_group != "established")) %>%
  group_by(species) %>%
  summarise(ecm_obs_count = n()) %>%
  ungroup()

fungal_root_ecm <- fungal_root_emp %>%
  # Keep only species observed with EcM but not with AM or dual
  filter(
    species %in% ecm_species & 
      !species %in% am_species_established &  # ← Use ESTABLISHED for EcM exclusion
      !species %in% dual_species
  ) %>%
  # Add observation counts
  left_join(ecm_obs_counts, by = "species") %>%
  # For species in dual genera, require AM evaluation in established individuals
  # AND at least 2 EcM observations
  filter(
    # NOT in a dual genus (keep all) OR evaluated for AM in established + ≥2 obs
    !genus %in% dual_genera | 
      (species %in% am_evaluated_species$species & ecm_obs_count >= 2)
  ) %>%
  # Keep only EcM observations
  filter(mycorrhizal_type == "EcM") %>%
  select(family, genus, species, mycorrhizal_type) %>%
  unique()

#### Dual ####
# Species already identified in dual_species list
fungal_root_dual <- fungal_root_emp %>%
  filter(species %in% dual_species) %>%
  select(family, genus, species) %>%
  # Assign all as "EcM-AM" type
  mutate(mycorrhizal_type = "EcM-AM") %>%
  unique()

#### NM ####
# Species are NM-only if:
# 1. They have been observed as NM
# 2. They have NOT been observed with AM, EcM, or as dual
# 3. Either: ≥2 observations of NM type, OR 1 observation with cluster roots,
#    OR in Proteaceae (cluster root family)

fungal_root_nm <- fungal_root_emp %>%
  # Keep only species observed as NM but not with AM, EcM, or dual
  filter(
    species %in% nm_species & 
      !species %in% am_species_any_age & 
      !species %in% ecm_species & 
      !species %in% dual_species
  ) %>%
  # Group by species to count NM observations
  group_by(species) %>%
  mutate(nm_count = n()) %>%
  ungroup() %>%
  # Keep species with ≥2 NM observations OR cluster roots OR Proteaceae
  filter(
    nm_count >= 2 | nm_structures == "cr" | family == "Proteaceae"
  ) %>%
  select(family, genus, species, mycorrhizal_type) %>%
  unique()

#### COMBINE DATA ####
# Base dataset
fungal_root_curated <- bind_rows(
  fungal_root_am,
  fungal_root_ecm,
  fungal_root_dual,
  fungal_root_nm
) %>%
  arrange(family, genus, species) %>%
  unique()

# Environmental breadth dataset
fungal_root_curated_env <- fungal_root_curated %>%
  # Add niche breadth estimates
  left_join(
    niche_breadth %>% select(-mycorrhizal_type), 
    by = c("family", "genus", "species")
  ) %>%
  # Remove species without niche breadth estimates
  filter(!is.na(env_B2_corrected))

# ============================================================================== #
# (6) SUMMARY STATISTICS -------------------------------------------------------
# ============================================================================== #

cat("\n=== MYCORRHIZAL TYPE CLASSIFICATION SUMMARY ===\n\n")
cat("Total species classified:", nrow(fungal_root_curated), "\n")
cat("  AM (identified):  ", sum(fungal_root_curated$mycorrhizal_type == "AM"), "\n")
cat("  AM (with breadth est):  ", sum(fungal_root_curated_env$mycorrhizal_type == "AM"), "\n")
cat("  EcM (identified): ", sum(fungal_root_curated$mycorrhizal_type == "EcM"), "\n")
cat("  EcM (with breadth est): ", sum(fungal_root_curated_env$mycorrhizal_type == "EcM"), "\n")
cat("  Dual (identified): ", sum(fungal_root_curated$mycorrhizal_type == "EcM-AM"), "\n")
cat("  Dual (with breadth est): ", sum(fungal_root_curated_env$mycorrhizal_type == "EcM-AM"), "\n")
cat("  NM (identified):  ", sum(fungal_root_curated$mycorrhizal_type == "NM"), "\n")
cat("  NM (with breadth est):  ", sum(fungal_root_curated_env$mycorrhizal_type == "NM"), "\n")

# Check for potential dual genera species that weren't evaluated
cat("=== QUALITY CHECK: Dual Genera Evaluation Status ===\n\n")
for (gen in dual_genera) {
  species_in_genus <- fungal_root_emp %>% 
    filter(genus == gen) %>% 
    pull(species) %>% 
    unique()
  
  classified_species <- fungal_root_curated %>% 
    filter(genus == gen) %>% 
    pull(species) %>% 
    unique()
  
  missing_species <- setdiff(species_in_genus, classified_species)
  
  if (length(missing_species) > 0) {
    cat("Genus:", gen, "\n")
    cat("  Species in data but NOT classified:", length(missing_species), "\n")
    cat("  These may need evaluation for both AM and EcM:\n")
    cat("    ", paste(missing_species, collapse = "\n     "), "\n\n")
  }
}

# Summerise mycorrhizal types by family and genus
cat("=== MYCORRHIZAL TYPE SUMMARY BY FAMILY ===\n\n")
fungal_root_curated %>%
  group_by(family, mycorrhizal_type) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = mycorrhizal_type, values_from = count, values_fill = 0) %>%
  print(n = Inf)
cat("\n=== MYCORRHIZAL TYPE SUMMARY BY GENUS ===\n\n")
fungal_root_curated %>%
  group_by(genus, mycorrhizal_type) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = mycorrhizal_type, values_from = count, values_fill = 0) %>%
  print(n = Inf)

# Summerise niche breadth by mycorrhizal type
cat("=== ENVIRONMENTAL BREADTH SUMMARY BY MYCORRHIZAL TYPE ===\n\n")
fungal_root_curated_env %>%
  group_by(mycorrhizal_type) %>%
  summarise(
    count = n(),
    mean_breadth = mean(env_B2_corrected, na.rm = TRUE),
    median_breadth = median(env_B2_corrected, na.rm = TRUE),
    sd_breadth = sd(env_B2_corrected, na.rm = TRUE)
  ) %>%
  print()

# Linear models of niche breadth by mycorrhizal type
lm_breadth <- lm(env_B2_corrected ~ mycorrhizal_type, data = fungal_root_curated_env)
parameters::parameters(lm_breadth) %>%
  print()

# Save the data
fwrite(
  fungal_root_curated, 
  "data/mycorrhizal_type_emperical/fungal_root_emp_classified.txt",
  sep = "\t"
)
fwrite(
  fungal_root_curated_env, 
  "data/mycorrhizal_type_emperical/niche_estimates_emperical.txt",
  sep = "\t"
)