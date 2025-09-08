# Required packages
require(data.table)
require(tidyverse)

# (1) Organise the data ############################################################

# Global Tree Search data
global_trees <- fread("data/GlobalTreeSearch/harmonised_global_tree.txt")

# Compute unique genus-family pairs
global_genus_families <- global_trees %>%
  distinct(genus, family)

# Compute unique species-genus-family pairs
global_species_genus_families <- global_trees %>%
  distinct(scientific_name, genus, family)

# Fungal root estimates for genera (curated by experts)
fungal_root_genus <- fread("data/FungalRoot/harmonised_fungal_root.txt") %>%
  # Filter to tree genera
  filter(
    genus %in% global_trees$genus
  ) %>%
  # Update family based on Global Tree Search
  select(
    -family
  ) %>%
  # Add family information based on Global Tree Search
  left_join(global_genus_families, by = "genus") %>%
  select(
    family, genus, mycorrhizal_type
  ) %>%
  # Simplify mycorrhizal types
  mutate(
    mycorrhizal_type = case_when(
      mycorrhizal_type == "NM-AM, rarely EcM" ~ "NM-AM",
      str_detect(mycorrhizal_type, "species-specific") ~ "uncertain",
      TRUE ~ mycorrhizal_type
    )
  ) %>%
  # Exclude uncertain mycorrhizal types - the annotation we will resolve
  filter(
    mycorrhizal_type != "uncertain"
  ) %>%
  # Resolve conflicts by keeping the most frequent mycorrhizal type for each genera
  group_by(genus, mycorrhizal_type) %>%
  mutate(
    count = n()
  ) %>%
  group_by(genus) %>%
  filter(
    count == max(count)
  ) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    -count
  )

# Check unique mycorrhizal types
unique(fungal_root_genus$mycorrhizal_type)

# Check for missing families
fungal_root_genus %>%
  filter(is.na(family)) %>%
  print(n = Inf)

# Fungal root measurements for tree species (some errors in mycorrhizal assinments)
fungal_root_species <- fread("data/FungalRoot/harmonised_measurements.txt") %>%
  # Filter to tree species
  filter(
    scientific_name %in% global_trees$scientific_name
  ) %>%
  # Updated family and genus information based on Global Tree Search
  select(
    -family, -genus
  ) %>%
  left_join(global_species_genus_families, by = "scientific_name") %>%
  select(family, genus, scientific_name, mycorrhizal_type) %>%
  # Simplify mycorrhizal types: This was actually done in the previous step 
  # "01a.Harmonise_plant_databases.R", but included it here for clarity of 
  # assumptions.
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
    )
  ) %>%
  # Exclude uncertain mycorrhizal types - the annotation we want to resolve
  filter(
    mycorrhizal_type != "uncertain" 
  ) %>%
  # Resolve conflicts by keeping the most frequent mycorrhizal type for each species
  group_by(scientific_name, mycorrhizal_type) %>%
  mutate(
    count = n()
  ) %>%
  group_by(scientific_name) %>%
  filter(
    count == max(count)
  ) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    -count
  )

# Check unique mycorrhizal types
unique(fungal_root_species$mycorrhizal_type)

# Check for missing families
fungal_root_species %>%
  filter(is.na(family)) %>%
  print(n = Inf)

# (2) Family-wise estimates of mycorrhizal type ################################

# For genera in GlobalTreeSearch with uncertain or unknown mycorrhizal types,
# we resolve theseusing a 90% consensus approach at the family-level.
# We compute family-wise estimates for both genus and species level data, giving 
# preference to genus-level data (expert curated, taking into account such as 
# observations for only a single mycorrhizal type) over species-level data
# in case of conflicts.

# Family-wise estimates for species-level data
fungal_root_family_1 <- fungal_root_genus %>%
  # Count the total number of genera per family
  group_by(family) %>%
  mutate(
    total = n()
  ) %>%
  # Compute the proportion of genera for each mycorrhizal type
  group_by(family, mycorrhizal_type) %>%
  summarise(
    count = n(),
    proportion = n() / first(total),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# Family-wise estimates for species-level data
fungal_root_family_2 <- fungal_root_species %>%
  # Count the total number of species per family
  group_by(family) %>%
  mutate(
    total = n()
  ) %>%
  # Compute the proportion of species for each mycorrhizal type
  group_by(family, mycorrhizal_type) %>%
  summarise(
    count = n(),
    proportion = n() / first(total),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# Take family-wise estimates at 90% concensus
fungal_root_family_consensus <- bind_rows(
  fungal_root_family_1 %>% mutate(dataset = "genus"),
  fungal_root_family_2 %>% mutate(dataset = "species")
  ) %>%
  # Filter to families with at least 90% consensus
  filter(
    proportion >= 0.9
  ) %>%
  # Dereplicate
  distinct(family, mycorrhizal_type, .keep_all = TRUE) %>%
  # Take genus (expert curated) estimates for residual duplicates/conflicts
  group_by(family) %>%
  mutate(n = n()) %>%
  filter(
    n == 1 | (n > 1 & dataset == "genus")
  ) %>%
  select(-n) %>%
  arrange(family) %>%
  print(n = Inf)

# Check for duplicates in family-wise estimates
fungal_root_family_consensus %>%
  group_by(family) %>%
  filter(
    n() > 1
  ) %>%
  arrange(family, mycorrhizal_type) %>%
  print(n = Inf)

# Family-wise unresolved mycorrhizal types
fungal_root_family_unresolved <- bind_rows(
  fungal_root_family_1 %>% mutate(dataset = "genus"),
  fungal_root_family_2 %>% mutate(dataset = "species")
) %>%
  # Filter to families with unresolved mycorrhizal types
  filter(
    !family %in% fungal_root_family_consensus$family
  ) %>%
  arrange(family) %>%
  print(n = Inf)

# How many families have resolved and unresolved mycorrhizal types at the 
# family-level?
n_distinct(fungal_root_family_consensus$family)
n_distinct(fungal_root_family_unresolved$family)

# Join data with single "uncertain" mycorrhizal type for unresolved families
fungal_root_family_estimates <- bind_rows(
  fungal_root_family_consensus,
  fungal_root_family_unresolved %>%
    mutate(mycorrhizal_type = "uncertain"
    )
  ) %>%
  select(
    family, mycorrhizal_type_family = mycorrhizal_type
  ) %>%
  distinct() %>%
  print(n = Inf)

# (3) Assign mycorrhizal type to GlobalTreeSearch ##############################

global_mycorrhizal_types <- global_trees %>%
  # Add fungal root estimates for genera
  left_join(fungal_root_genus, by = c("genus", "family")) %>%
  # Add family estimates for unresolved mycorrhizal types (NAs)
  left_join(fungal_root_family_estimates, by = "family") %>%
  mutate(
    estimate = case_when(
      is.na(mycorrhizal_type) ~ "family",
      !is.na(mycorrhizal_type) ~ "genus",
      TRUE ~ NA_character_
    ),
    mycorrhizal_type = case_when(
      !is.na(mycorrhizal_type) ~ mycorrhizal_type,
      is.na(mycorrhizal_type) ~ mycorrhizal_type_family,
      TRUE ~ NA_character_
    ),
    # Residual NAs in estimate represent families not in fungal root
    estimate = case_when(
      is.na(mycorrhizal_type) ~ NA_character_,
      TRUE ~ estimate
    ),
    mycorrhizal_type = case_when(
      is.na(mycorrhizal_type) ~ "uncertain",
      TRUE ~ mycorrhizal_type
    )
  ) %>%
  select(
    scientific_name, genus, family, mycorrhizal_type, estimate
  )

# Check for residual unknowns
unknown_mycorrhizal_types <- global_mycorrhizal_types %>% 
  filter(is.na(mycorrhizal_type) | is.na(estimate)) %>%
  as_tibble() %>%
  print(n = Inf)

# Double check these are genera not in fungal root
fungal_root_genus %>%
  filter(family %in% unknown_mycorrhizal_types$family)
fungal_root_genus %>%
  filter(genus %in% unknown_mycorrhizal_types$genus)
fungal_root_species %>%
  filter(family %in% unknown_mycorrhizal_types$family)
fungal_root_species %>%
  filter(genus %in% unknown_mycorrhizal_types$genus)

# Check distinct mycorrhizal types
unique(global_mycorrhizal_types$mycorrhizal_type)

# (4) Amend errors in mycorrhizal types assignments ############################

# Check if any uncertain species are in the fungal_root_species dataset:
# Note Pisonia and Daviesia are species-specific:
species_specific_estimates <- global_mycorrhizal_types %>%
  filter(mycorrhizal_type == "uncertain") %>%
  select(scientific_name) %>%
  inner_join(fungal_root_species, by = "scientific_name") %>%
  select(scientific_name, mycorrhizal_type) %>%
  # Pisonia umbellifera is likely AM: https://doi.org/10.1016/j.funeco.2014.09.001
  bind_rows(
    tibble(
      scientific_name = "Pisonia umbellifera",
      mycorrhizal_type = "AM"
    )
  ) %>%
  as_tibble() %>%
  print(n = Inf)

# Check ambiguous assignments for OM trees
global_mycorrhizal_types %>%
  filter(mycorrhizal_type == "OM")

# Lobelia OM status appears to be a mistake
fungal_root_species %>%
  filter(genus == "Lobelia")
fungal_root_genus %>%
  filter(genus == "Lobelia")

# Check eucalypt assignments given eucalypt dual status is generally considered 
# temporally dependent, transitioning from AM or dual-mycorrhizal as seedlings 
# to ECM as mature trees:
# Brundrett & Tedersoo. 2018. Evolutionary history of mycorrhizal symbioses and global host plant diversity. New Phytologist 220: 1108–1115.
# Teste et al. 2020. Dual-mycorrhizal plants: their ecology and relevance. New Phytologist 225: 1835–1851.
# Note: Eucalypts are cloesly related and are often considered a single 
# high-level clade with 13 subgenera (including Eucalyptusm, Corymbia and 
# Angophora), hence lumping for mycorrhizal estimates:
# Brooker. 2000. A new classification of the genus Eucalyptus L'Her. (Myrtaceae). Australian systematic botany 13(1): 79-148.

# Check for eucalypt assignments
fungal_root_genus %>%
  filter(genus %in% c("Eucalyptus", "Corymbia", "Angophora")) %>%
  select(family, genus, mycorrhizal_type) %>%
  print(n = Inf)

# Check the temporal mycorrhizal type assumption on induviduals less than 1 year
# old against species older than 1 year old
fungal_root_eucalypts <- inner_join(
  fread("data/FungalRoot/harmonised_occurrences.txt"),
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
  filter(
    genus %in% c("Eucalyptus", "Corymbia", "Angophora"),
    !is.na(host_age)
  ) %>%
  # Mutate host age to a manageable format
  mutate(host_age = case_when(
    host_age == ">12 months" |  host_age == "12 months" ~ ">12",
    host_age == "<1 month" | host_age == "1 month" | host_age == "2 months" |
      host_age == "3 months" | host_age == "4 months" | host_age == "5 months" |
      host_age == "8 months" | host_age == "9 months" | host_age == "10 months" ~ "<12",
    TRUE ~ as.character(host_age)
  )) %>%
  select(family, genus, host_age, mycorrhizal_type) %>%
  as_tibble() %>%
  print(n = Inf)

# Compute stats to confirm the temporal mycorrhizal type assumption
bind_rows(
  # Seedlings
  fungal_root_eucalypts %>%
    filter(host_age == "<12") %>%
    mutate(total = n()) %>%
    # Calculate the relative portions of findings
    group_by(mycorrhizal_type) %>%
    summarise(
      n_obs = n(),
      rel_obs = round(first(n_obs)/first(total), digits = 3)
    ) %>%
    mutate(age = "<12"),
  # Older trees
  fungal_root_eucalypts %>%
    filter(host_age == ">12") %>%
    mutate(total = n()) %>%
    # Calculate the relative portions of findings
    group_by(mycorrhizal_type) %>%
    summarise(
      n_obs = n(),
      rel_obs = round(first(n_obs)/first(total), digits = 3)
    ) %>%
    mutate(age = ">12")
  ) %>%
  print(n = Inf)

# Generate the final dataset
global_mycorrhizal_types_final <- global_mycorrhizal_types %>%
  mutate(
    mycorrhizal_type = case_when(
      # Fix the ambiguous "OM" assignment for Lobelia
      genus == "Lobelia" ~ "AM",
      # Fix the Pisonia species-specific assignments
      scientific_name %in% species_specific_estimates$scientific_name ~ 
        species_specific_estimates$mycorrhizal_type[
          match(scientific_name, species_specific_estimates$scientific_name)
        ],
      # Fix the eucalypt  assignments
      genus %in% c("Eucalyptus", "Corymbia", "Angophora") ~ "EcM",
      TRUE ~ mycorrhizal_type
    ),
    # Update to species-specific estimates for Pisonia
    estimate = case_when(
      scientific_name %in% species_specific_estimates$scientific_name ~ "species",
      TRUE ~ estimate
    )
  )

# (5) Save Global Tree Search with mycorrhizal types #############################

# Compute stats for species assignments
global_mycorrhizal_types_final %>%
  group_by(mycorrhizal_type) %>%
  summarise(
    n_obs = n(),
    rel_obs = round(n_obs / nrow(global_mycorrhizal_types_final), digits = 3)
  ) %>%
  print(n = Inf)

# Compute stats for genus assignments
global_mycorrhizal_types_final %>%
  select(genus, mycorrhizal_type) %>%
  distinct() %>%
  mutate(total = n()) %>%
  group_by(mycorrhizal_type) %>%
  summarise(
    n_obs = n(),
    rel_obs = round(n_obs / first(total), digits = 3)
  ) %>%
  print(n = Inf)

# Count destinct species and genera
n_distinct(global_mycorrhizal_types_final$scientific_name)
n_distinct(global_mycorrhizal_types_final$genus)

# Save the final dataset
global_mycorrhizal_types_final %>%
  fwrite(
    "output/generated_data/global_tree_mycorrhizal_types.txt", sep = "\t"
  )

# (6) Assign mycorrhizal types to Australian native trees ######################

# Read in Australian tree
australian_trees <- fread("data/APC/harmonised_apc_tree_list.txt") %>%
  inner_join(
    global_mycorrhizal_types_final %>%
      select(scientific_name, mycorrhizal_type),
    by = "scientific_name"
  )

# Check missing and uncertain taxa
australian_trees %>%
  filter(is.na(mycorrhizal_type) | mycorrhizal_type == "uncertain") %>%
  as_tibble() %>%
  print(n = Inf)

# Address Myrtaceae on a species-by-species basis due to uncertainty
australian_trees_final <- australian_trees %>%
  mutate(
    mycorrhizal_type = case_when(
      # The Myrtaceae species "Lindsayomyrtus racemoides" & "Osbornia octodonta"
      # are likely AM because of thier mangrove and rainforest habitat
      scientific_name == "Lindsayomyrtus racemoides" ~ "AM",
      scientific_name == "Osbornia octodonta" ~ "AM",
      scientific_name == "Pisonia brunoniana" ~ "AM",
      # Australian Proteaceae are functionally NM bust can be colonised by AM
      family == "Proteaceae" ~ "NM",
      TRUE ~ mycorrhizal_type
    )) %>%
  as_tibble() %>%
  print(n = Inf)

# Save the Australian native trees with mycorrhizal types
australian_trees_final %>%
  unique() %>%
  fwrite("output/generated_data/australian_tree_mycorrhizal_types.txt", sep = "\t")
