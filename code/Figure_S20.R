require(data.table)
require(dplyr)

fread("output/generated_data/australian_tree_mycorrhizal_types.txt") %>%
  filter(mycorrhizal_type == "NM") %>%
  group_by(family) %>%
  summarise(n = n())
fread("data/presence/trees_10.txt") %>%
  select(mycorrhizal_type, family, genus, species = scientific_name) %>%
  distinct() %>%
  filter(mycorrhizal_type == "NM") %>%
  group_by(family) %>%
  summarise(n = n())
    
# Trees in our data ############################################################

trees<-fread("output/generated_data/australian_tree_mycorrhizal_types.txt") %>%
  filter(family == "Proteaceae") %>%
  filter(
    family == "Proteaceae",
    genus != "Protea",
    genus != "Knightia",
    genus != "Brabejum",
    genus != "Embothrium",
    genus != "Serruria",
    genus != "Kermadecia",
    genus != "Finschia",
    genus != "Roupala"
  ) %>%
  select(family, genus, species = scientific_name) %>%
  distinct() %>%
  # Mutate synonyms to accepted names
  mutate(
    genus = case_when(
      genus == "Conchium" ~ "Hakea",
      genus == "Linkia" ~ "Persoonia",
      genus == "Hylogyne" ~ "Telopea",
      genus == "Pentadactylon" ~ "Persoonia",
      genus == "Pycnonia" ~ "Persoonia",
      genus == "Anadenia" ~ "Grevillea",
      genus == "Sirmuellera" ~ "Banksia",
      genus == "Tricondylus" ~ "Lomatia",
      genus == "Cyanocarpus" ~ "Helicia",
      genus == "Lysanthe" ~ "Grevillea",
      genus == "Molloya" ~ "Strangea",
      TRUE ~ genus
    )
  ) %>%
  mutate(
    vegetation_type = case_when(
      genus %in% c(
        "Grevillea", "Hakea", "Persoonia", "Dryandra", "Banksia", "Synaphea", 
        "Petrophile", "Conospermum", "Adenanthos", "Isopogon", "Lomatia", 
        "Lambertia", "Stirlingia", "Xylomelum", "Telopea", "Strangea", 
        "Franklandia", "Symphionema", "Acidonia", "Bellendena"
      ) ~ "sclerophyll",
      genus %in% c(
        "Stenocarpus", "Helicia", "Orites", "Macadamia", "Hollandaea", 
        "Lasjia", "Alloxylon", "Hicksbeachia", "Bleasdalea", "Buckinghamia", 
        "Austromuellera", "Musgravea", "Darlingia", "Eidothea", "Cardwellia",
        "Athertonia", "Catalepidia", "Opisthiolepis", "Floydia", "Neorites",
        "Megahertzia","Carnarvonia","Sphalmium", "Cenarrhenes", "Agastachys",
        "Placospermum", "Nothorites", "Triunia"
      ) ~ "rainforest",
    )
  ) %>%
  as_tibble() %>% 
  print(n = Inf)
trees %>%
  group_by(vegetation_type) %>%
  summarise(n = n()) %>%
  print(n = Inf)
trees %>%
  select(genus, vegetation_type) %>%
  distinct() %>%
  group_by(vegetation_type) %>%
  summarise(n = n()) %>%
  print(n = Inf)

# Trees in the Australian flora ################################################

trees <- fread("output/generated_data/plants_10.txt")

# Plants in the Australian flora ###############################################

plants <- fread("data/APC/harmonised_apc_flora_list.txt") %>%
  filter(
    family == "Proteaceae",
    genus != "Protea",
    genus != "Knightia",
    genus != "Brabejum",
    genus != "Embothrium",
    genus != "Serruria",
    genus != "Kermadecia",
    genus != "Finschia",
    genus != "Roupala"
    ) %>%
  select(family, genus, species = scientific_name) %>%
  distinct() %>%
  # Mutate synonyms to accepted names
  mutate(
    genus = case_when(
      genus == "Conchium" ~ "Hakea",
      genus == "Linkia" ~ "Persoonia",
      genus == "Hylogyne" ~ "Telopea",
      genus == "Pentadactylon" ~ "Persoonia",
      genus == "Pycnonia" ~ "Persoonia",
      genus == "Anadenia" ~ "Grevillea",
      genus == "Sirmuellera" ~ "Banksia",
      genus == "Tricondylus" ~ "Lomatia",
      genus == "Cyanocarpus" ~ "Helicia",
      genus == "Lysanthe" ~ "Grevillea",
      genus == "Molloya" ~ "Strangea",
      TRUE ~ genus
    )
  ) %>%
  mutate(
    vegetation_type = case_when(
      genus %in% c(
        "Grevillea", "Hakea", "Persoonia", "Dryandra", "Banksia", "Synaphea", 
        "Petrophile", "Conospermum", "Adenanthos", "Isopogon", "Lomatia", 
        "Lambertia", "Stirlingia", "Xylomelum", "Telopea", "Strangea", 
        "Franklandia", "Symphionema", "Acidonia", "Bellendena"
      ) ~ "sclerophyll",
      genus %in% c(
        "Stenocarpus", "Helicia", "Orites", "Macadamia", "Hollandaea", 
        "Lasjia", "Alloxylon", "Hicksbeachia", "Bleasdalea", "Buckinghamia", 
        "Austromuellera", "Musgravea", "Darlingia", "Eidothea", "Cardwellia",
        "Athertonia", "Catalepidia", "Opisthiolepis", "Floydia", "Neorites",
        "Megahertzia","Carnarvonia","Sphalmium", "Cenarrhenes", "Agastachys",
        "Placospermum", "Nothorites", "Triunia"
      ) ~ "rainforest",
      TRUE ~ "unknown"
    )
  ) %>%
  as_tibble() %>% 
  print(n = Inf)

# Shrubs in the Australian flora ###############################################

# Remove trees from the plants data to get shrubs
shrubs <- anti_join(
  plants,
  trees,
  by = c("family", "genus", "species")
) %>%
  print(n = Inf)

# Plot by growth form ###############################################################

plants <- plants %>%
  # Define tree and 
  mutate(
    growth_form = case_when(
      paste(family, genus, species) %in% paste(trees$family, trees$genus, trees$species) ~ "tree",
      TRUE ~ "shrub"
    )
  ) 

plants %>%
  group_by(genus) %>%
  summarise(
    genus = first(genus),
    n_species = n(),
    n_shrub_species = sum(growth_form == "shrub"),
    n_tree_species = sum(growth_form == "tree"),
    vegetation_type = first(vegetation_type),
  ) %>%
  print(n = Inf)

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtext)

# Calculate proportions for each chart
# Chart 1: Overall growth form proportions
growth_props <- plants %>%
  count(growth_form) %>%
  mutate(proportion = n / sum(n))

# Chart 2: Vegetation type within shrubs
shrub_veg <- plants %>%
  filter(growth_form == "shrub") %>%
  count(vegetation_type) %>%
  mutate(proportion = n / sum(n))

# Chart 3: Vegetation type within trees
tree_veg <- plants %>%
  filter(growth_form == "tree") %>%
  count(vegetation_type) %>%
  mutate(proportion = n / sum(n))

# Create the three charts
p1 <- ggplot(
  growth_props,
  aes(x = "", y = proportion, fill = growth_form)) +
  geom_col(width = 1, color = "white", size = 1) +
  scale_fill_manual(
    values = c("tree" = "#6B8E23", "shrub" = "#CC7722"),
    labels = c("tree" = "Trees", "shrub" = "Shrubs")
    ) +
  labs(
    y = NULL,
    x = NULL,
    tag = ("(**a**)"),
    fill = "Growth\nForm"
    ) +
  coord_flip() +
  # Reverse the axis
  scale_y_reverse() +
  theme_void() +
  theme(
    plot.tag = element_markdown(size = 14),
    aspect.ratio = 0.2,
    plot.margin = margin(1, 1, 1, 1)
    )

p2 <- ggplot(
  shrub_veg, 
  aes(x = "", y = proportion, fill = vegetation_type)) +
  geom_col(width = 1, color = "white", size = 1) +
  scale_fill_manual(
    values = c("rainforest" = "forestgreen", "sclerophyll" = "#cc4c02"),
    labels = c("rainforest" = "Rainforest", "sclerophyll" = "Sclerophyll"),
    breaks = c("sclerophyll", "rainforest") 
  
    ) +
  labs(
    y = NULL, 
    x = NULL, 
    tag = ("(**b**)"),
    fill = "Vegetation\nType"
    ) +
  coord_flip() +
  theme_void() +
  theme(
    plot.tag = element_markdown(size = 14),
    aspect.ratio = 0.2,
    plot.margin = margin(1, 1, 1, 1)
  )

p3 <- ggplot(
  tree_veg,
  aes(x = "", y = proportion, fill = vegetation_type)) +
  geom_col(width = 1, color = "white", size = 1) +
  scale_fill_manual(
    values = c("rainforest" = "forestgreen", "sclerophyll" = "#cc4c02"),
    labels = c("rainforest" = "Rainforest", "sclerophyll" = "Sclerophyll")
  ) +
  labs(y = NULL, x = NULL, tag = ("(**c**)")) +
  coord_flip() +
  # Scale axis to percentages
  scale_y_continuous(
    labels = scales::percent,
    breaks = seq(0, 1, by = 0.2)
    ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.tag = element_markdown(size = 14),
    axis.text.x = element_text(size = 10),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    aspect.ratio = 0.2,
    plot.margin = margin(1, 1, 1, 1)
  )

# Combine the three plots
combined_plot <- p1 / p2 / p3 

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave(
  filename = "output/suplimentary_proteaceae/figure.png",
  plot = combined_plot,
  width = 16,
  height = 10,
  units = "cm",
  dpi = 300
)

# Optional: Print summary statistics
cat("\nSummary Statistics:\n")
cat("==================\n\n")
cat("Growth Form Distribution:\n")
print(growth_props)
cat("\nVegetation Types in Shrubs:\n")
print(shrub_veg)
cat("\nVegetation Types in Trees:\n")
print(tree_veg)
