
# Build the phylogenetic tree of Australian native tree species, used for
# phylogenetic eigenvector predictors in 05a.Species_niche_breadth_differences.R,
# 05b.Species_niche_position_differences.R, 07d.Niche_breadth_differences_emperical.R,
# and to draw the tree in Figure_1.R.
# Run this after 01.Harmonise_tree_data.R (it reads that script's output) and
# before any of the four scripts above.

# Required packages
require(V.PhyloMaker2)
require(ape)
require(data.table)
require(tidyverse)

# Read in the harmonised tree list and restrict to native species with a
# resolved mycorrhizal type, matching the species population used throughout
# the rest of this pipeline (see Figure_1.R)
data <- fread("generated_data/global_tree_mycorrhizal_types.txt") %>%
  filter(
    native_status == "native",
    mycorrhizal_type != "uncertain"
  ) %>%
  mutate(species = gsub(" ", "_", scientific_name))

# Create a species list for phylo.maker
species_list <- data %>%
  select(species, genus, family) %>%
  unique() %>%
  # Change family names to suit the phylo.maker backbone dataset
  mutate(
    family = case_when(
      family == "Viburnaceae" ~ "Adoxaceae",
      family == "Orchidaceae" ~ "Sapotaceae",
      TRUE ~ family
    )
  ) %>%
  as.data.frame()

# How many species are there in total?
n_distinct(species_list$species)

# Create and save the phylogenetic tree
phylo_tree <- phylo.maker(species_list, scenarios = "S3")$scenario.3
write.tree(phylo_tree, file = "generated_data/phylo_tree_mycorrhizal_types.tre")
