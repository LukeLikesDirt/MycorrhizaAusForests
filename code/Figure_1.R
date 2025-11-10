
# Load libraries
require(data.table)
library(ggtext)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ape)
library(cowplot)
require(tidyverse)

# Constants and helper functions ###############################################

# Plot constants
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 9
myco_colors <- c(
  'ErM' = "grey40",
  'NM' = '#d470a2',
  "AM" = "#E69F00",
  "NM-AM" = "#D55E00",
  "EcM" = "#56B4E9",
  "Dual" = "#009E73"
)

# Function to get node range for a family in the phylogenetic tree
get_family_range <- function(tree_data, family_name) {
  family_tips <- data_matched$species[data_matched$family == family_name]
  family_tips <- family_tips[family_tips %in% tree_data$tip.label]
  
  if(length(family_tips) < 2) return(NULL)
  
  tip_numbers <- match(family_tips, tree_data$tip.label)
  return(range(tip_numbers))
}

# Helper function to add family strip with mycorrhizal type species distribution
add_family_strip <- function(plot_obj, tree, family_name, hjust_val, offset_text_val) {
  family_range <- get_family_range(tree, family_name)
  if(is.null(family_range) || length(family_range) != 2) return(plot_obj)
  
  # Get family data for mycorrhizal type breakdown
  family_data <- data_matched[data_matched$family == family_name, ]
  myco_counts <- table(family_data$mycorrhizal_type)
  
  # Add main family label
  plot_obj <- plot_obj +
    geom_strip(family_range[1], family_range[2], 
               label = family_name, 
               color = 'black',  # Changed to black for abundance
               hjust = hjust_val, 
               offset.text = offset_text_val, 
               fontsize = 4, 
               offset = 35, 
               barsize = 0.6)
  
  # Add mycorrhizal type strips
  myco_types <- names(myco_counts)[myco_counts > 0]
  current_start <- family_range[1]
  
  for(myco_type in myco_types) {
    type_count <- myco_counts[myco_type]
    type_end <- current_start + type_count - 1
    
    if(type_end <= family_range[2]) {
      plot_obj <- plot_obj +
        geom_strip(current_start, type_end, 
                   color = myco_colors[myco_type], 
                   hjust = 0.5, 
                   fontsize = 3, 
                   barsize = 3, 
                   offset = 70)
    }
    current_start <- type_end + 1
  }
  
  return(plot_obj)
}

# Mycorrhizal distributions ####################################################

#### * Global tree data * ####

global_trees <- fread("output/generated_data/global_tree_mycorrhizal_types.txt") %>% 
  mutate(
    # Mutate EcM-AM to Dual
    mycorrhizal_type = ifelse(mycorrhizal_type == "EcM-AM", "Dual", mycorrhizal_type)
  )

# Summarise global tree data:
global_tree_summary <- global_trees %>%
  # Remove uncertain mycorrhizal types
  filter(mycorrhizal_type != "uncertain") %>%
  mutate(
    n_species = n_distinct(scientific_name)
  ) %>%
  group_by(mycorrhizal_type) %>%
  summarise(
    abs_richness = n_distinct(scientific_name),
    rel_richness = n_distinct(scientific_name) / first(n_species),
    dataset = "Global Trees"
  ) %>%
  ungroup() %>%
  mutate(
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "NM-AM", "EcM", "Dual", "ErM", "NM")
    )
  ) %>%
  print()

# How many species are there in total?
n_distinct(global_trees$scientific_name)

# How many genera are there in total?
n_distinct(global_trees$genus)

#### * Australian tree data * ####

australian_trees <- fread("output/generated_data/australian_tree_mycorrhizal_types.txt") %>% 
  mutate(
    # Mutate EcM-AM to Dual
    mycorrhizal_type = ifelse(mycorrhizal_type == "EcM-AM", "Dual", mycorrhizal_type),
    # Create species column with underscores for tree matching
    species = gsub(" ", "_", scientific_name)
  )

# Summarise Australian tree data:
australian_tree_summary <- australian_trees %>%
  # Calculate relative richness
  mutate(
    n_species = n_distinct(scientific_name)
  ) %>%
  group_by(mycorrhizal_type) %>%
  summarise(
    abs_richness = n_distinct(scientific_name),
    rel_richness = n_distinct(scientific_name) / first(n_species),
    dataset = "Aust. Trees"
  ) %>%
  ungroup() %>%
  mutate(
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "NM-AM", "EcM", "Dual", "ErM", "NM")
    )
  ) %>%
  print()

# How many species are there in total?
n_distinct(australian_trees$scientific_name)

# How many genera are there in total?
n_distinct(australian_trees$genus)

#### * Summarise diversity of mycorrhizal types * ####

# Summarise richness by mycorrhizal type
global_summary <- global_tree_summary %>%
  group_by(mycorrhizal_type) %>%
  summarise(global_richness = sum(abs_richness, na.rm = TRUE), .groups = "drop")

australian_summary <- australian_tree_summary %>%
  group_by(mycorrhizal_type) %>%
  summarise(australian_richness = sum(abs_richness, na.rm = TRUE), .groups = "drop")

# Join summaries and calculate proportion
summary_table <- full_join(global_summary, australian_summary, by = "mycorrhizal_type") %>%
  mutate(
    australian_richness = coalesce(australian_richness, 0),
    proportion_in_aus = australian_richness / global_richness
  )

# Add combined "AM + NM-AM" row
am_nmam_global <- sum(global_tree_summary$abs_richness[global_tree_summary$mycorrhizal_type %in% c("AM", "NM-AM")], na.rm = TRUE)
am_nmam_aus <- sum(australian_tree_summary$abs_richness[australian_tree_summary$mycorrhizal_type %in% c("AM", "NM-AM")], na.rm = TRUE)
am_nmam_proportion <- am_nmam_aus / am_nmam_global

am_nmam_row <- tibble(
  mycorrhizal_type = "AM + NM-AM",
  global_richness = am_nmam_global,
  australian_richness = am_nmam_aus,
  proportion_in_aus = am_nmam_proportion
)

# Combine into final table
final_table <- bind_rows(summary_table, am_nmam_row) %>%
  arrange(mycorrhizal_type)

# Print table
print(final_table)

#### * Plot Australian trees vs global trees * ####

mycorrhizal_diversity_plot <- ggplot(
  tibble(
    Dataset = c(australian_tree_summary$dataset, global_tree_summary$dataset),
    relative_richness = c(australian_tree_summary$rel_richness, global_tree_summary$rel_richness),
    mycorrhizal_type = c(australian_tree_summary$mycorrhizal_type, global_tree_summary$mycorrhizal_type)
  ),
  aes(Dataset, relative_richness, fill = mycorrhizal_type)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = TRUE)) +
  theme_void() +
  theme(
    plot.tag.location = "panel",
    plot.tag.position = c(-0.1, 1.1),
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    legend.title = element_blank(),
    legend.text = element_text(size = text_size),
    legend.justification = 'center',
    legend.margin = margin(t = 10, b = 5),  # add space above legend and below
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = title_size, margin = margin(r = 3)),  # space to the right of y-axis text
    axis.text.x = element_text(size = text_size, margin = margin(t = 3)),  # space below x-axis text
    axis.title.x = element_text(size = title_size, margin = margin(t = 5)),  # space above x-axis title
    plot.tag = element_markdown(size = tag_size),
    aspect.ratio = 0.25,
    plot.margin = margin(t = 0, r = 10, b = 0, l = 5),
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    panel.grid.minor.x = element_line(colour = "grey90", linewidth = 0.35),
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = FALSE)) +
  scale_x_discrete(labels = c("Australia", "Global")) +
  scale_fill_manual(
    values = myco_colors,
    labels = c(
      'AM' = 'Arbuscular mycorrhizal (AM)',
      'NM-AM' = 'Weakly arbuscular\nmycorrhizal (NM–AM)',
      'EcM' = 'Ectomycorrhizal (EcM)',
      'Dual' = 'Dual-mycorrhizal',
      'ErM' = 'Ericoid mycorrhizal (ErM)',
      'NM' = 'Non-mycorrhizal (NM)'
    )
  ) +
  coord_flip() +
  labs(x = NULL, y = 'Relative richness', tag = "(**b**)") + 
  guides(fill = guide_legend(nrow = 2)) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent,
    expand = c(0.01, 0.01)  # shrinks padding on both ends
  )

# Display the plot
print(mycorrhizal_diversity_plot)

# Phylogenetic tree for Australian tree species ################################
# Read in the phylogenetic tree
phylo_tree <- read.tree("output/generated_data/phylo_tree_mycorrhizal_types.tre")

# Your color palette
myco_colors <- c(
  'ErM' = "grey40",
  'NM' = '#d470a2',
  "AM" = "#E69F00",
  "NM-AM" = "#D55E00",
  "EcM" = "#56B4E9",
  "Dual" = "#009E73"
)

# Top 10 families
top_families <- australian_trees %>%
  group_by(family) %>%
  summarise(n = n_distinct(scientific_name), .groups = "drop") %>%
  arrange(desc(n)) %>%
  slice_head(n = 10) %>%
  pull(family)

# Prepare species list matching tree tips
species_list <- australian_trees %>%
  select(species, genus, family, mycorrhizal_type) %>%
  distinct(species, .keep_all = TRUE)

# Match species to tree tips and drop unmatched tips
spmatch <- match(phylo_tree$tip.label, species_list$species)
pruned_tree <- drop.tip(phylo_tree, phylo_tree$tip.label[is.na(spmatch)])

# Reorder data to match pruned tree
spmatch_clean <- match(pruned_tree$tip.label, species_list$species)
data_matched <- species_list[spmatch_clean, ]

# Verify matching
all(pruned_tree$tip.label == data_matched$species)

# Create a dummy column for abundance (like the 'mete' column in original)
data_matched$abundance <- 1

# Create the base circular tree plot
tree_plot <- ggtree(
  pruned_tree, 
  ladderize = FALSE, 
  layout = "circular", 
  color = "grey50", 
  size = 0.12
  )

# Add abundance bars (in black) and mycorrhizal type coloring
tree_plot1 <- tree_plot +
  geom_fruit(
    data = data_matched,
    geom = geom_bar,
    mapping = aes(x = abundance, y = species, fill = mycorrhizal_type),
    orientation = "y",
    stat = "identity",
    colour = NA,
    pwidth = 0.1,
    offset = 0
  ) +
  scale_fill_manual(
    values = myco_colors,
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 1),
    name = "Mycorrhizal type"
  )


# Manual positioning values with family names for easy adjustment
family_positions <- list(
  "Myrtaceae" = list(hjust = 1, offset_text = 60),
  "Fabaceae" = list(hjust = 0.5, offset_text = 75),
  "Proteaceae" = list(hjust = 0, offset_text = 55),
  "Lauraceae" = list(hjust = 0, offset_text = 50),
  "Sapindaceae" = list(hjust = 1, offset_text = 55),
  "Rubiaceae" = list(hjust = 0, offset_text = 50),
  "Rutaceae" = list(hjust = 1, offset_text = 55),
  "Malvaceae" = list(hjust = 1, offset_text = 55),
  "Euphorbiaceae" = list(hjust = 1, offset_text = 60),
  "Moraceae" = list(hjust = 1, offset_text = 60)
)

# Add family strips for top families
tree_plot2 <- tree_plot1

# Add strips for each top family using named positioning
for(family in top_families) {
  # Get positioning for this family, with defaults if not specified
  if(family %in% names(family_positions)) {
    hjust_val <- family_positions[[family]]$hjust
    offset_val <- family_positions[[family]]$offset_text
  }
  
  tree_plot2 <- add_family_strip(tree_plot2, pruned_tree, family, hjust_val, offset_val)
}

# Final plot adjustments
tree_plot_final <- tree_plot2 +
  labs(
    tag = "(**a**)"
  ) +
  theme(
    plot.tag = element_markdown(size = 14),
    plot.tag.location = "panel",
    plot.tag.position = c(-0.01, 0.88),
    plot.margin = margin(t = -25, r = -15, b = -25, l = 0, "pt"),
    legend.position = "none",
    aspect.ratio = 1
  )

# Display the plot
print(tree_plot_final)

# Join and save the plots ######################################################

# Final layout: top row stacked over diversity plot
final_plot <- plot_grid(
  tree_plot_final,
  mycorrhizal_diversity_plot,
  rel_heights = c(2.2, 1),
  ncol = 1
)

# Save png
ggsave(
  "output/figure1.png",
  final_plot,
  width = 16, 
  height = 20, 
  bg = "white",
  units = "cm", 
  dpi = 300
)

# Save tiff
ggsave(
  "output/figure1.tiff",
  final_plot,
  width = 16, 
  height = 20, 
  bg = "white",
  units = "cm"
)

