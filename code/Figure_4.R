
# Required packages
require(ggtext)
require(tidyverse)

# Load the data
load('output/generated_data/figure_4.RData')

#### Theme ####

# Myco types
myco_types <- c("AM", "EcM", "Dual", "NM")

# Colours
myco_colours <- c(
  'NM' = '#d470a2',
  "AM" = "#E69F00",
  "EcM" = "#56B4E9",
  "Dual" = "#009E73"
)

# Plot constants
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 9

# Define common plot theme for consistent styling
common_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    axis.ticks = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.title = element_text(face = "bold", size = strip_size, hjust = 0.5),
    plot.tag = element_markdown(size = tag_size),
    strip.text = element_text(face = "bold", size = strip_size),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    aspect.ratio = 1
  )

# Create mycorrhizal type legend
dummy_myco_data <- data.frame(
  mycorrhizal_type = factor(myco_types, levels = myco_types),
  value = 1:4
)

# Add theme if I want to use the plot as a legend
dummy_plot_myco <- ggplot(
  dummy_myco_data, 
  aes(x = value, y = value, colour = mycorrhizal_type)
) +
  geom_point(
    fill = alpha('lightgrey', 0.33),
    shape = 22,
    size = 7,
    stroke = 1.25
  ) +
  scale_colour_manual(values = myco_colours, name = "Mycorrhizal\ntype") +
  theme_void() +
  scale_y_continuous(limits = c(0,0))

legend_myco <- cowplot::get_legend(dummy_plot_myco)

# # Load and prepare niche data
# data_position <- data.table::fread("data/data_position_enmeval/niche_estimates.txt") %>%
#   filter(mycorrhizal_type != "ErM") %>%
#   mutate(
#     # Rename EcM-AM to Dual
#     mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
#     mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM")),
#     species = str_replace_all(species, " ", "_"),
#     # Scale the RC variables for the full dataset
#     RC1 = as.numeric(scale(RC1_position)),
#     RC2 = as.numeric(scale(RC2_position)),
#     RC3 = as.numeric(scale(RC3_position)),
#     biome = ifelse(dominant_biome == "Tropical", "tropical", "nontropical"),
#     biome = factor(
#       biome,
#       levels = c("tropical", "nontropical")
#     )
#   ) %>%
#   select(
#     family, genus, species, mycorrhizal_type, biome,
#     RC1, RC2, RC3
#   )


# Figure a: Niche overlap kernel density estimates #############################

# Define plot parameters
rc_vars <- c("RC1", "RC2", "RC3")
rc_labels <- c("RC1 (temperature & decomposition)", "RC2 (soil moisture)", "RC3 (soil phosphorus)")
biomes <- c("all", "tropical", "nontropical")
biome_titles <- c("All", "Tropical", "Temperate")

# Breaks and limit density peaks:
# !!! These need to be adjusted based visual inspection of kernel density estimates !!!
y_limits <- list(c(0, 1.42), c(0, 0.8), c(0, 0.675))
y_breaks <- list(seq(0, 2, by = 0.4), seq(0, 2, by = 0.2), seq(0, 2, by = 0.2))

# Function to create density plots
create_density_plot <- function(rc_var, biome, rc_idx, biome_idx) {
  
  # Filter data based on biome
  data = data_position
  plot_data <- if(biome == "all") data else data %>% filter(biome == !!biome)
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = !!sym(rc_var), colour = mycorrhizal_type)) +
    geom_density(aes(fill = mycorrhizal_type), alpha = 0.33, fill = "lightgrey", linewidth = 1.25) +
    scale_colour_manual(values = myco_colours) +
    scale_y_continuous(
      limits = y_limits[[rc_idx]],
      breaks = y_breaks[[rc_idx]]
    ) +
    scale_x_continuous(
      limits = c(min(data[[rc_var]]), max(data[[rc_var]])),
      breaks = scales::pretty_breaks()
    ) +
    common_theme
  
  # Set labels based on position
  x_label <- if(rc_idx == 1 && biome == "tropical") rc_labels[rc_idx] else NULL
  x_label <- if(rc_idx == 2 && biome == "tropical") rc_labels[rc_idx] else x_label
  x_label <- if(rc_idx == 3 && biome == "tropical") rc_labels[rc_idx] else x_label
  
  y_label <- if(biome == "all" && rc_idx == 2) "Density" else NULL
  
  title <- if(rc_idx == 1) biome_titles[biome_idx] else NULL
  
  tag <- if(rc_var == "RC1" && biome == "all") "(**a**)" else NULL
  
  # Apply labels
  p <- p + labs(x = x_label, y = y_label, title = title, tag = tag)
  
  # Apply theme modifications
  if(biome != "all") {
    p <- p + theme(axis.text.y = element_blank())
  }
  
  if(rc_idx == 1 && biome == "all") {
    p <- p + theme(plot.tag.position = c(0.05, 0.95))
  }
  
  if(biome == "tropical" && !is.null(x_label)) {
    p <- p + theme(plot.margin = margin(1, 1, 3, 1, "pt"))
  }
  
  return(p)
}

# Generate all plots using nested loops
density_plots <- list()
for(i in seq_along(rc_vars)) {
  for(j in seq_along(biomes)) {
    plot_name <- paste0(rc_vars[i], "_", biomes[j], "_overlap")
    density_plots[[plot_name]] <- create_density_plot(rc_vars[i], biomes[j], i, j)
  }
}

# Figure b: Niche overlap Schoener's D #########################################

# Split into environmental and geographic groups
env_list <- split(
  data_position %>%
    select("RC1", "RC2", "RC3"),
  data_position %>%
    pull(mycorrhizal_type)
)
env_list_tropical <- split(
  data_position %>%
    filter(biome == "tropical") %>%
    select("RC1", "RC2", "RC3"),
  data_position %>%
    filter(biome == "tropical") %>%
    pull(mycorrhizal_type)
)
env_list_nontropical <- split(
  data_position %>%
    filter(biome == "nontropical") %>%
    select("RC1", "RC2", "RC3"),
  data_position %>%
    filter(biome == "nontropical") %>%
    pull(mycorrhizal_type)
)

# Function to compute 3D Schoener's D using kernel density estimation
compute_schoeners_d_3d <- function(data_list, space_type = "Environmental") {
  n <- length(data_list)
  groups <- names(data_list)
  result <- list()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      group1_data <- data_list[[i]]
      group2_data <- data_list[[j]]
      
      schoeners_d <- tryCatch({
        # Convert to matrices
        group1_matrix <- as.matrix(group1_data)
        group2_matrix <- as.matrix(group2_data)
        combined_matrix <- rbind(group1_matrix, group2_matrix)
        
        # Define grid boundaries
        x_range <- range(combined_matrix[,1])
        y_range <- range(combined_matrix[,2])
        z_range <- range(combined_matrix[,3])
        
        # Create 3D grid
        grid_res <- 20 # Grid resolution
        x_seq <- seq(x_range[1], x_range[2], length.out = grid_res)
        y_seq <- seq(y_range[1], y_range[2], length.out = grid_res)
        z_seq <- seq(z_range[1], z_range[2], length.out = grid_res)
        
        grid_points <- as.matrix(expand.grid(x = x_seq, y = y_seq, z = z_seq))
        
        # Calculate kernel density for each group
        n1 <- nrow(group1_matrix)
        n2 <- nrow(group2_matrix)
        d <- ncol(group1_matrix)  # dimensions
        
        # Silverman's bandwidth for multivariate case
        h1 <- (4/(d+2))^(1/(d+4)) * n1^(-1/(d+4)) * apply(group1_matrix, 2, sd)
        h2 <- (4/(d+2))^(1/(d+4)) * n2^(-1/(d+4)) * apply(group2_matrix, 2, sd)
        
        # Calculate densities at grid points
        density1 <- rep(0, nrow(grid_points))
        density2 <- rep(0, nrow(grid_points))
        
        for (k in 1:nrow(grid_points)) {
          # Density for group 1
          for (l in 1:nrow(group1_matrix)) {
            diff1 <- (grid_points[k,] - group1_matrix[l,]) / h1
            density1[k] <- density1[k] + exp(-0.5 * sum(diff1^2)) / (prod(h1) * (2*pi)^(d/2))
          }
          density1[k] <- density1[k] / n1
          
          # Density for group 2
          for (l in 1:nrow(group2_matrix)) {
            diff2 <- (grid_points[k,] - group2_matrix[l,]) / h2
            density2[k] <- density2[k] + exp(-0.5 * sum(diff2^2)) / (prod(h2) * (2*pi)^(d/2))
          }
          density2[k] <- density2[k] / n2
        }
        
        # Normalise densities
        density1 <- density1 / sum(density1)
        density2 <- density2 / sum(density2)
        
        # Calculate Schoener's D
        schoeners_d <- 1 - 0.5 * sum(abs(density1 - density2))
        
      }, error = function(e) {
        warning(paste("Error computing Schoener's D for", groups[i], "vs", groups[j], ":", e$message))
        NA
      })
      
      result[[length(result) + 1]] <- tibble(
        group1 = groups[i],
        group2 = groups[j],
        schoeners_d = as.numeric(schoeners_d),
        type = space_type
      )
    }
  }
  
  bind_rows(result)
}

# Compute Schoener's D for environmental data (3D)
schoeners_env <- compute_schoeners_d_3d(env_list, "Environmental")
schoeners_env_tropical <- compute_schoeners_d_3d(env_list_tropical, "Environmental (Tropical)")
schoeners_env_nontropical <- compute_schoeners_d_3d(env_list_nontropical, "Environmental (Nontropical)")

# Priority for ordering pairs
priority <- c("AM" = 1, "EcM" = 2, "Dual" = 3, "NM" = 4)
order_pair <- function(g1, g2) {
  if (priority[g1] < priority[g2]) c(g1, g2)
  else if (priority[g1] > priority[g2]) c(g2, g1)
  else sort(c(g1, g2))
}

# Prepare data for plotting
schoeners_df <- schoeners_env %>%
  rowwise() %>%
  mutate(
    ordered_pair = paste(order_pair(group1, group2), collapse = " ∩ ")
  ) %>%
  ungroup() %>%
  distinct(ordered_pair, type, .keep_all = TRUE) %>%
  mutate(pair_label = ordered_pair) %>%
  # Levels for pair_label
  mutate(pair_label = factor(pair_label, levels = unique(ordered_pair)))
schoeners_df_tropical <- schoeners_env_tropical %>%
  rowwise() %>%
  mutate(
    ordered_pair = paste(order_pair(group1, group2), collapse = " ∩ ")
  ) %>%
  ungroup() %>%
  distinct(ordered_pair, type, .keep_all = TRUE) %>%
  mutate(pair_label = ordered_pair) %>%
  # Levels for pair_label
  mutate(pair_label = factor(pair_label, levels = unique(ordered_pair)))
schoeners_df_nontropical <- schoeners_env_nontropical %>%
  rowwise() %>%
  mutate(
    ordered_pair = paste(order_pair(group1, group2), collapse = " ∩ ")
  ) %>%
  ungroup() %>%
  distinct(ordered_pair, type, .keep_all = TRUE) %>%
  mutate(pair_label = ordered_pair) %>%
  # Levels for pair_label
  mutate(pair_label = factor(pair_label, levels = unique(ordered_pair)))

# Check the levels of pair_label
print(levels(schoeners_df$pair_label))

# Define thresholds
low_thresh <- 0.2
moderate_thresh <- 0.4
high_thresh <- 0.6
very_high_thresh <- 0.8

# Schoener's D plot for all trees
schoener_all <- ggplot(schoeners_df, aes(x = pair_label, y = schoeners_d, fill = type, colour = type)) +
  geom_col(linewidth = 0.25) +
  scale_fill_manual(values = c("Environmental" = "grey")) +
  scale_colour_manual(values = c("Environmental" = "black")) +
  geom_hline(yintercept = low_thresh, linetype = "dotted", color = "#de2d26") +
  geom_hline(yintercept = moderate_thresh, linetype = "dotted", color = "#de2d26") +
  geom_hline(yintercept = high_thresh, linetype = "dotted", color = "#de2d26") +
  geom_hline(yintercept = very_high_thresh, linetype = "dotted", color = "#de2d26") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "Schoener's D",
    tag = "(**b**)"
  ) +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# Schoener's D plot for tropical trees
schoener_tropical <- ggplot(schoeners_df_tropical, aes(x = pair_label, y = schoeners_d, fill = type, colour = type)) +
  geom_col(linewidth = 0.25) +
  scale_fill_manual(values = c("Environmental (Tropical)" = "grey")) +
  scale_colour_manual(values = c("Environmental (Tropical)" = "black")) +
  geom_hline(yintercept = low_thresh, linetype = "dotted", color = "#de2d26") +
  geom_hline(yintercept = moderate_thresh, linetype = "dotted", color = "#de2d26") +
  geom_hline(yintercept = high_thresh, linetype = "dotted", color = "#de2d26") +
  geom_hline(yintercept = very_high_thresh, linetype = "dotted", color = "#de2d26") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )

# Schoener's D plot for temperate trees
schoener_nontropical <- ggplot(schoeners_df_nontropical, aes(x = pair_label, y = schoeners_d, fill = type, colour = type)) +
  geom_col(linewidth = 0.25) +
  scale_fill_manual(values = c("Environmental (Nontropical)" = "grey")) +
  scale_colour_manual(values = c("Environmental (Nontropical)" = "black")) +
  geom_hline(yintercept = low_thresh, linetype = "dotted", color = "#de2d26") +
  geom_hline(yintercept = moderate_thresh, linetype = "dotted", color = "#de2d26") +
  geom_hline(yintercept = high_thresh, linetype = "dotted", color = "#de2d26") +
  geom_hline(yintercept = very_high_thresh, linetype = "dotted", color = "#de2d26") +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )

# Join and save ################################################################

figure4 <- patchwork::wrap_plots(
  density_plots[["RC1_all_overlap"]], density_plots[["RC1_tropical_overlap"]], density_plots[["RC1_nontropical_overlap"]],
  density_plots[["RC2_all_overlap"]], density_plots[["RC2_tropical_overlap"]], density_plots[["RC2_nontropical_overlap"]],
  density_plots[["RC3_all_overlap"]], density_plots[["RC3_tropical_overlap"]], density_plots[["RC3_nontropical_overlap"]],
  schoener_all, schoener_tropical, schoener_nontropical,
  nrow = 4
)

figure4_final <- cowplot::plot_grid(
  figure4, legend_myco, rel_widths = c(1, 0.16)
)

# Save the plot
ggsave(
  "output/figure4.png",
  plot = figure4_final,
  width = 15.5,
  height = 20.75, 
  bg = "white",
  units = "cm",
  dpi = 300
)

# Save the plot
ggsave(
  "output/figure4.tif",
  plot = figure4_final,
  width = 15.5,
  height = 20.75, 
  bg = "white",
  units = "cm"
)

# Figure b: Heat map option ####################################################

# Define the colours
my_colors <- rev(paletteer::paletteer_c("grDevices::Spectral", 100))

# Dummy data for heatmap legend
dummy_data_schoener <- tibble(x = 1, y = 1, z = 0:1)

# Create dummy plot for Schoener's D legend
plot_schoener <- ggplot(dummy_data_schoener, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = my_colors,
    name = "\nSchoener's D",
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    #labels = c("Very low", "Low", "Moderate", "High", "Very high")
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0)
  )


# Extract legends and convert to gtables
gt_myco <- cowplot::get_legend(dummy_plot_myco)
gt_schoener <- ggplotGrob(plot_schoener)$grobs[[
  which(sapply(ggplotGrob(plot_schoener)$grobs, function(x) x$name) == "guide-box")
]]

# Combine legends with controlled spacing
legend_combined <- gridExtra::gtable_rbind(gt_myco, gt_schoener, size = "max")
legend_combined$heights[2] <- unit(0, "pt")  # Adjust this value to control gap size

# Schoener's D heatmap for all trees
schoener_heatmap_all <- ggplot(schoeners_df, aes(x = group1, y = group2, fill = schoeners_d)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", schoeners_d)), 
            size = 3, color = "black", fontface = "bold") +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    name = "Schoener's D"
  ) +
  scale_x_discrete(
    limits = c("AM", "EcM", "Dual")
  ) +
  scale_y_discrete(
    limits = rev(c("EcM", "Dual", "NM"))
  ) +
  labs(
    x = NULL,
    y = NULL,
    tag = "(**b**)"
  ) +
  common_theme +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Schoener's D heatmap for tropical trees
schoener_heatmap_tropical <- ggplot(schoeners_df_tropical, aes(x = group1, y = group2, fill = schoeners_d)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", schoeners_d)), 
            size = 3, color = "black", fontface = "bold") +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    name = "Schoener's D"
  ) +
  scale_x_discrete(
    limits = c("AM", "EcM", "Dual")
  ) +
  scale_y_discrete(
    limits = rev(c("EcM", "Dual", "NM"))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  common_theme +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank()
  )

# Schoener's D heatmap for nontropical trees
schoener_heatmap_nontropical <- ggplot(schoeners_df_nontropical, aes(x = group1, y = group2, fill = schoeners_d)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", schoeners_d)), 
            size = 3, color = "black", fontface = "bold") +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    name = "Schoener's D"
  ) +
  scale_x_discrete(
    limits = c("AM", "EcM", "Dual")
  ) +
  scale_y_discrete(
    limits = rev(c("EcM", "Dual", "NM"))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  common_theme +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank()
  )

figure4 <- patchwork::wrap_plots(
  density_plots[["RC1_all_overlap"]], density_plots[["RC1_tropical_overlap"]], density_plots[["RC1_nontropical_overlap"]],
  density_plots[["RC2_all_overlap"]], density_plots[["RC2_tropical_overlap"]], density_plots[["RC2_nontropical_overlap"]],
  density_plots[["RC3_all_overlap"]], density_plots[["RC3_tropical_overlap"]], density_plots[["RC3_nontropical_overlap"]],
  schoener_heatmap_all, schoener_heatmap_tropical, schoener_heatmap_nontropical,
  nrow = 4
)

figure4_final <- cowplot::plot_grid(
  figure4, legend_combined, rel_widths = c(1, 0.18)
)

# Save the plot
ggsave(
  "output/figure4_heat.png",
  plot = figure4_final,
  width = 15.5,
  height = 19,   
  bg = "white",
  units = "cm",
  dpi = 300
)

# Save the plot
ggsave(
  "output/figure4_heat.tif",
  plot = figure4_final,
  width = 15.5,
  height = 19,   
  bg = "white",
  units = "cm"
)