
# This subset data to extract sources for mycorrhizal types in FungalRoot.
# Unfortunately this is messy and I have taken the time to make a clean, 
# logical, and well annotated script.

# Required packages
require(data.table)
require(tidyverse)

# (1) Organise the global data -------------------------------------------------

### Read in the data ###

# Get species names
species_names <- fread("data/FungalRoot/harmonised_tree_occurrences.txt") %>%
  select(ID, family, genus, species = scientific_name)

# Get species measurements
species_measurements <- inner_join(
  species_names,
  fread("data/FungalRoot/measurements.csv") %>%
    select(ID = `Core ID`, measurementType, measurementValue) %>%
    pivot_wider(names_from = measurementType, values_from = measurementValue) %>%
    rename(
      source = Source,
      name = Name,
      mycorrhizal_type = `Mycorrhiza type`,
      AM_intensity = `AM intensity`,
      AM_method = `AM method`,
      elevation_min = `Elevation min`,
      ECM_method = `EcM method`,
      host_age = `Host age`,
      ECM_intensity = `EcM intensity`,
      remark = `Remark: type`,
      ECM_frequency = `EcM frequency`,
      AM_frequency = `AM frequency`,
      ERM_intensity = `ErM intensity`,
      ERM_method = `ErM method`,
      ORM_method = `OrM method`,
      ORM_intensity = `OrM intensity`,
      elevation_max = `Elevation max`
    ) %>%
    mutate(
      mycorrhizal_type = case_when(
        mycorrhizal_type == "EcM, AM undetermined" ~ "EcM",
        mycorrhizal_type == "non-ectomycorrhizal (AM undetermined)" ~ "uncertain",
        mycorrhizal_type == "EcM, no AM colonization" ~ "EcM",
        mycorrhizal_type == "non-mycorrhizal" ~ "NM",
        mycorrhizal_type == "Other" ~ "uncertain",
        mycorrhizal_type == "EcM,AM" ~ "EcM-AM",
        mycorrhizal_type == "AM-like (non-vascular plants)" ~ "AM-like_non-vascular_plant",
        mycorrhizal_type == "ErM,EcM" ~ "ErM",
        mycorrhizal_type == "ErM,AM" ~ "ErM",
        TRUE ~ mycorrhizal_type
      )),
  by = "ID"
) %>%
  as_tibble()

# How many unique sources are there?
n_distinct(species_measurements$source)

### Compute mycorrhizal type estimates ###

# Potential dual species: Retain all species with either EcM-AM or one record of 
# each for AM and EcM type
potential_dual <- species_measurements %>%
  filter(
    mycorrhizal_type %in% c("EcM-AM", "AM", "EcM")
  ) %>%
  # Group by genus and name to identify species with both types
  group_by(species) %>%
  mutate(
    n_types = n_distinct(mycorrhizal_type)
  ) %>%
  filter(
    n_types == 2 | mycorrhizal_type == "EcM-AM"
  ) %>%
  ungroup() %>%
  select(ID, family, genus, species, mycorrhizal_type, source) %>%
  distinct() %>%
  arrange(family, genus, species, mycorrhizal_type) %>%
  print(n = Inf)

# Facilitative dual species: Retain all "potential_dual" species with NM records
faculative_dual <- species_measurements %>%
  filter(
    mycorrhizal_type == "NM"
  ) %>%
  # Join and filter faculative duals
  inner_join(
    potential_dual %>%
      select(species) %>%
      distinct(),
    by = "species"
  ) %>%
  select(ID, family, genus, species, mycorrhizal_type, source) %>%
  arrange(family, genus, species, mycorrhizal_type) %>%
  print(n = Inf)

# Facilitative AM species: Retain all species with AM and NM records across 
# different measurements
faculative_AM <- species_measurements %>%
  filter(
    !species %in% potential_dual$species,    # Remove the potentially dual species
    mycorrhizal_type %in% c("AM", "NM")      # Keep only AM and NM types
  ) %>%
  # Group by species and count types
  group_by(species) %>%
  mutate(
    n_types = n_distinct(mycorrhizal_type)
  ) %>%
  # Filter "faculative AM" species with both AM and NM records
  filter(
    n_types == 2
  ) %>%
  ungroup() %>%
  select(ID, family, genus, species, mycorrhizal_type, source) %>%
  distinct() %>%
  arrange(family, genus, species, mycorrhizal_type) %>%
  print(n = Inf)

# Facilitative EcM species: Retain all species with EcM and NM records across
# different measurements
faculative_EcM <- species_measurements %>%
  filter(
    !species %in% potential_dual$species,    # Remove the potentially dual species
    mycorrhizal_type %in% c("EcM", "NM")     # Keep only EcM and NM types
  ) %>%
  # Group by species and count types
  group_by(species) %>%
  mutate(
    n_types = n_distinct(mycorrhizal_type)
  ) %>%
  # Filter "faculative EcM" species with both EcM and NM records
  filter(
    n_types == 2
  ) %>%
  ungroup() %>%
  select(ID, family, genus, species, mycorrhizal_type, source) %>%
  distinct() %>%
  arrange(family, genus, species, mycorrhizal_type) %>%
  print(n = Inf)

#### Compute global estimates ####

global_estimates <- species_measurements %>%
  # Remove uncertain mycorrhizal types
  filter(
    mycorrhizal_type != "uncertain"
  ) %>%
  # Count the number of observations per species
  group_by(species) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  mutate(
    mycorrhizal_type = case_when(
      # Mutate potential duals to EcM-AM
      species %in% potential_dual$species ~ "EcM-AM",
      # Mutate species with "AM and "NM" to "AM"
      species %in% faculative_AM$species & mycorrhizal_type == "NM" ~ "AM",
      # Mutate species with "EcM" and "NM" to "EcM"
      species %in% faculative_EcM$species & mycorrhizal_type == "NM" ~ "EcM",
      TRUE ~ mycorrhizal_type
    )
  ) %>%
  select(family, genus, species, mycorrhizal_type, n_obs) %>%
  unique() %>%
  arrange(family, genus, species) %>%
  print(n = Inf)

# Compute the proportion of each mycorrhizal type in global trees
global_prop <- global_estimates %>%
  mutate(
    total_species = n_distinct(species)
  ) %>%
  group_by(mycorrhizal_type) %>%
  summarise(
    n_species = n_distinct(species),
    proportion = n_species / first(total_species)
  )

# Compute the propotion of each mycorrhizal type in each family
global_family_prop <- global_estimates %>%
  group_by(family) %>%
  mutate(
    total_family_species = n_distinct(species)
  ) %>%
  group_by(family, mycorrhizal_type) %>%
  summarise(
    n_species = n_distinct(species),
    proportion = n_species / first(total_family_species)
  ) %>%
  arrange(family, mycorrhizal_type) %>%
  print(n = Inf)

#### Dual mycorrhizal sources with host age ####

# Fungal root measurements by age of host plant
fungal_root_host_age <- species_measurements %>%
  filter(
    # Remove uncertain mycorrhizal types
    mycorrhizal_type != "uncertain",
    # Retain potential duals
    species %in% potential_dual$species
  ) %>%
  # Mutate host age to a manageable format
  mutate(host_age = case_when(
    host_age == ">12 months" | host_age == "12 months" ~ ">1 yr",
    host_age == "<1 month" | host_age == "1 month" | host_age == "2 months" |
      host_age == "3 months" | host_age == "4 months" | host_age == "5 months" |
      host_age == "8 months" | host_age == "9 months" | host_age == "10 months" |
      host_age == "<12 months" ~ "<1 yr",
    is.na(host_age) ~ "uncertain",
    TRUE ~ as.character(host_age)
  )) %>%
  select(ID, family, genus, species, host_age, mycorrhizal_type, source) %>%
  as_tibble() %>%
  print(n = Inf)

# How many unique "potentially dual" species are there?
n_distinct(fungal_root_host_age$species)

# (2) Mycorrhizal type estimates for Australian trees --------------------------

# Species data for Australian trees
australian_estimates <- species_measurements %>%
  filter(species %in% (fread(
    "generated_data/global_tree_mycorrhizal_types.txt"
  ) %>%
    filter(native_status == "native") %>%
    pull(scientific_name))
  ) %>%
  # Remove uncertain mycorrhizal types
  filter(
    mycorrhizal_type != "uncertain"
  ) %>%
  # Count the number of observations per species
  group_by(species) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  mutate(
    mycorrhizal_type = case_when(
      # Mutate potential duals to EcM-AM
      species %in% potential_dual$species ~ "EcM-AM",
      # Mutate species with "AM and "NM" to "AM"
      species %in% faculative_AM$species & mycorrhizal_type == "NM" ~ "AM",
      # Mutate species with "EcM" and "NM" to "EcM"
      species %in% faculative_EcM$species & mycorrhizal_type == "NM" ~ "EcM",
      TRUE ~ mycorrhizal_type
    )
  ) %>%
  select(family, genus, species, mycorrhizal_type, n_obs) %>%
  unique() %>%
  arrange(family, genus, species) %>%
  print(n = Inf)
# Compute the proportion of each mycorrhizal type in Australian trees
australian_estimates %>%
  mutate(
    total_species = n_distinct(species)
  ) %>%
  group_by(mycorrhizal_type) %>%
  summarise(
    n_species = n_distinct(species),
    proportion = n_species / first(total_species)
  )

# Get niche estimate info
niche_estimates <- fread("generated_data/niche_estimates.txt") %>%
  select(
    family, genus, species, env_breadth = env_B2_corrected, biome, climate_zone
  ) %>%
  # Join australian estimates
  inner_join(
    australian_estimates,
    by = c("family", "genus", "species")
  ) %>%
  # Mutate "EcM-AM" to "Dual"
  mutate(
    mycorrhizal_type = case_when(
      mycorrhizal_type == "EcM-AM" ~ "Dual",
      TRUE ~ mycorrhizal_type
    )
  )

# (3) Sources for mycorrhizal types in Australian trees ------------------------

# NOTE: Here I will evaluate non-native Australia trees when they are found in
# sources along with Australian trees in case I decide to expand the data set to
# global trees later

# Save the Australian estimate data with sources
australian_sources <- species_measurements %>%
  filter(species %in% australian_estimates$species) %>%
  # Remove uncertain mycorrhizal types
  filter(
    mycorrhizal_type != "uncertain"
  ) %>%
  # Compute n_sources and n_obs per species
  group_by(species) %>%
  mutate(
    # Add a native column
    native_status = "native"
  ) %>%
  select(ID, family, genus, species, mycorrhizal_type, native_status, host_age, remark, source) %>%
  arrange(source, family, genus, species) %>%
  print(n = Inf)

# Generate a list of non-native Australian tree species from the sources
global_sources <- species_measurements %>%
  filter(
    # Filter to sources used for Australian tree species
    source %in% australian_sources$source,
    # Remove Australian tree species
    !species %in% australian_estimates$species
  ) %>%
  # Remove uncertain mycorrhizal types
  filter(
    mycorrhizal_type != "uncertain"
  ) %>%
  # Compute n_sources and n_obs per species
  group_by(species) %>%
  mutate(
    # Add a native column
    native_status = "non-native"
  ) %>%
  select(ID, family, genus, species, mycorrhizal_type, native_status, host_age, remark, source) %>%
  arrange(source, family, genus, species) %>%
  print(n = Inf)

# Join the two data sets
australian_sources <- bind_rows(
  australian_sources,
  global_sources
) %>%
  arrange(source, family, genus, species) %>%
  print(n = Inf)

# How many unique sources are there for Australian tree species?
n_distinct(australian_sources$source)

# Check the number of sources for the species in the niche estimate data
niche_estimate_sources <- australian_sources %>%
  filter(
    species %in% niche_estimates$species
  ) %>%
  pull(source) %>%
  unique() %>%
  glimpse()

# Generate a data set with sources for niche estimates and additional sources
australian_sources_niche <- australian_sources %>%
  # Filter to species used for niche estimates
  filter(source %in% niche_estimate_sources) %>%
  print(n = Inf)
australian_sources_additional <- australian_sources %>%
  # Get sources not used for niche estimates
  filter(!source %in% niche_estimate_sources) %>%
  print(n = Inf)

# How many unique sources for niche estimates and additional sources?
n_distinct(australian_sources_niche$source)
n_distinct(australian_sources_additional$source)

# Create a named list of sources
list_of_sources <- list(
  all_sources = australian_sources,
  niche_estimate_sources = australian_sources_niche,
  additional_sources = australian_sources_additional
)

# Save as xlsx
writexl::write_xlsx(
  list_of_sources,
  "data/mycorrhizal_type_emperical/australian_tree_sources.xlsx"
)

# (4) Combine curated and uncurated data ---------------------------------------

# Australian tree species list
australian_trees <- fread("generated_data/global_tree_mycorrhizal_types.txt") %>%
  filter(native_status == "native") %>%
  pull(scientific_name)

# Fungal root curated: Data that has been curated for mycorrrhizal structures
# and host life stage
fungal_root_curated <- fread("data/mycorrhizal_type_emperical/fungal_root_curated.txt") %>%
  select(-native_status) %>%
  filter(species %in% australian_trees) %>%
  arrange(source, family, species)

# Fungal root unassessed: Fungal root data that could not be accessed and curated
# due to language barriers or paywalls
fungal_root_uncurated <- fread("data/mycorrhizal_type_emperical/fungal_root_unassessed.txt") %>%
  select(-native_status) %>%
  mutate(replicates = as.numeric(replicates)) %>%
  filter(species %in% australian_trees) %>%
  arrange(source, family, species)

# Fungal root uncurated: Fungal root data that has not been curated due to lack of time
fungal_root_uncurated_no_time <- readxl::read_excel(
  "data/mycorrhizal_type_emperical/australian_tree_sources.xlsx",
  sheet = "all_sources"
) %>%
  select(fungal_root_id = ID, family, genus, species, mycorrhizal_type, source) %>%
  filter(
    # Retain native Australian tree species
    species %in% australian_trees,
    # Remove curated and uncurated data
    !fungal_root_id %in% fungal_root_curated$fungal_root_id,
    !fungal_root_id %in% fungal_root_uncurated$fungal_root_id
  ) %>%
  mutate(
    species_verbatim = "",
    replicates = as.numeric(""),
    habitat = "",
    host_age = "",
    host_age_group = "",
    am_evaluated = case_when(
      mycorrhizal_type %in% c("AM", "EcM-AM") ~ "yes",
      TRUE ~ "-"
    ),
    am_criteria = "",
    am_structures = "",
    ecm_evaluated = case_when(
      mycorrhizal_type %in% c("EcM", "EcM-AM") ~ "yes",
      TRUE ~ "-"
    ),
    ecm_structures = "",
    erm_structure = "",
    nm_structures = "",
    latitude = "",
    longitude = "",
    comment = "uncurated: no time to curate",
    link = ""
  ) %>%
  select(
    fungal_root_id, family, genus, species, species_verbatim, mycorrhizal_type,
    replicates, habitat, host_age, host_age_group, am_evaluated, am_criteria,
    am_structures, ecm_evaluated, ecm_structures, erm_structure, nm_structures,
    latitude, longitude, comment, source, link
  ) %>%
  arrange(source, family, species)

# Check if the order of the names maches across the three data sets
names(fungal_root_curated)==names(fungal_root_uncurated)
names(fungal_root_curated)==names(fungal_root_uncurated_no_time)

# Combine the three data sets
fungal_root_combined <- bind_rows(
  fungal_root_curated,
  fungal_root_uncurated,
  fungal_root_uncurated_no_time
) %>%
  fwrite(
    "generated_data/emperical_measurments.txt",
    sep = "\t"
  )

