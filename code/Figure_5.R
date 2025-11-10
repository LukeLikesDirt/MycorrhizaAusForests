
# Required packages
require(ggtext)
require(tidyverse)

# Load the data
load('output/generated_data/figure_5.RData')

# Access effect results
density_effect_all <- env_breadth_results[["all"]]$effect$density
effect_all <- env_breadth_results[["all"]]$effect$summary
density_effect_tropical <- env_breadth_results[["tropical"]]$effect$density
effect_tropical <- env_breadth_results[["tropical"]]$effect$summary
density_effect_nontropical <- env_breadth_results[["nontropical"]]$effect$density
effect_nontropical <- env_breadth_results[["nontropical"]]$effect$summary

# Access intercept results
intercept_all <- env_breadth_results[["all"]]$intercept
intercept_tropical <- env_breadth_results[["tropical"]]$intercept
intercept_nontropical <- env_breadth_results[["nontropical"]]$intercept

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
    panel.grid.major = element_line(linewidth = 0.3),
    panel.grid.minor = element_line(linewidth = 0.15),
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
  aes(x = value, y = value, fill = mycorrhizal_type)
) +
  geom_point(
    colour = alpha('grey30', 0.5),
    shape = 22,
    size = 5
  ) +
  scale_fill_manual(values = myco_colours, name = "Mycorrhizal\ntype") +
  theme_void() +
  scale_y_continuous(limits = c(0,0))

legend_myco <- cowplot::get_legend(dummy_plot_myco)

# Axis limits for effects plot
x_limits <- c(
  # Get all x values
  density_effect_all %>% pull(x),
  density_effect_tropical %>% pull(x),
  density_effect_nontropical %>% pull(x)
)

# Define min and max limits, with slight padding for visual clarity
x_min <- min(x_limits)
x_max <- max(x_limits)

# Function to create individual effects plots
effects_plot <- function(density_data, effect_data, intercept_data, title, breadth_type = "", x_lab = "", plot_tag = "") {
  
  # Determine y-axis text visibility
  show_y_axis <- title == "All"
  
  # Create intercept annotation data
  intercept_annotation <- tibble(
    x = -Inf,
    y = "AM",  # Position at the top mycorrhizal type
    label = paste0("&beta;<sub>0</sub> = ", intercept_data$intercept_formatted)
  )
  
  ggplot() +
    # Baseline (the grand mean)
    geom_vline(
      xintercept = 0,
      linetype = "dotted"
    ) +
    # Density ridges
    ggridges::geom_density_ridges(
      data = density_data,
      aes(x = x, y = mycorrhizal_type, height = y, fill = mycorrhizal_type),
      stat = "identity", scale = 0.85, colour = alpha("grey30", 0.5), alpha = 0.9
    ) +
    # Error bars
    geom_errorbarh(
      data = effect_data,
      aes(y = mycorrhizal_type, xmin = lower, xmax = upper),
      height = 0, size = 1
    ) +
    # Mean points
    geom_point(
      data = effect_data,
      aes(x = mean, y = mycorrhizal_type),
      shape = 21, fill = "white", stroke = 0.5, size = 1.5
    ) +
    # Intercept annotation
    geom_richtext(
      data = intercept_annotation,
      aes(y = y, x = x, label = label),
      hjust = 0.15, vjust = 0.1, size = 2.5,
      fill = NA, label.color = NA,
      inherit.aes = FALSE
    ) +
    # Styling
    scale_fill_manual(values = myco_colours, limits = myco_types) +
    scale_y_discrete(limits = myco_types) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = seq(-1, 1, 0.1)
    ) +
    coord_flip() +
    common_theme +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = if(show_y_axis) element_text() else element_blank(),
      plot.title = element_blank(),
      axis.title.y = if(show_y_axis) element_text(margin = margin(r = -2.5)) else element_blank()
    ) +
    labs(x = if(show_y_axis) x_lab else "", tag = plot_tag)
}

# Function to plots raw data
raw_data_plot <- function(raw_data, title, breadth_type = "", x_lab = "", plot_tag = "") {
  # Prepare raw data
  raw_data <- raw_data %>%
    mutate(
      mycorrhizal_type = factor(mycorrhizal_type, levels = myco_types),
      mycorrhizal_type_nudged = as.numeric(mycorrhizal_type) - 0.3
    ) 
  
  # Means
  mean_data <- raw_data %>%
    group_by(mycorrhizal_type) %>%
    mutate(
      mean = mean(.data[[breadth_type]]),
      sd = sd(.data[[breadth_type]])
    ) %>%
    ungroup()
  
  # Determine y-axis text visibility
  show_y_axis <- title == "All"
  
  ggplot() +
    # Raw points (left)
    geom_point(
      data = raw_data,
      aes(x = .data[[breadth_type]], y = mycorrhizal_type, color = mycorrhizal_type),
      position = position_jitter(width = 0, height = 0.2),
      alpha = 0.2, size = 1
    ) +
    # Error bars
    geom_errorbarh(
      data = mean_data,
      aes(y = mycorrhizal_type, xmin = mean - sd, xmax = mean + sd),
      height = 0, size = 1
    ) +
    # Mean points
    geom_point(
      data = mean_data,
      aes(x = mean, y = mycorrhizal_type),
      shape = 21, fill = "white", stroke = 0.5, size = 1.5
    ) +
    # Styling
    scale_color_manual(values = myco_colours) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    coord_flip() +
    common_theme +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = if(show_y_axis) element_text() else element_blank(),
      axis.title.y = if(show_y_axis) element_text() else element_blank(),
      plot.title = if(breadth_type == "env_breadth") element_text() else element_blank(),
      plot.tag.position = c(0.045, 0.92), # Adjust for title
      panel.grid.minor = element_blank()
    ) +
    labs(x = if(show_y_axis) x_lab else "", 
         title = if(breadth_type == "env_breadth") title else "",
         tag = plot_tag)
}

# Create raw data plots
raw_plot_all <- raw_data_plot(data_breadth, "All", "env_breadth", "Environmental breadth", "(**a**)")
raw_plot_trop <- raw_data_plot(data_breadth_tropical, "Tropical", "env_breadth")
raw_plot_nontrop <- raw_data_plot(data_breadth_nontropical, "Temperate", "env_breadth")

# Create effects panels
effect_plot_all <- effects_plot(density_effect_all, effect_all, intercept_all, "All", x_lab = "Effect size", plot_tag = "(**b**)")
effect_plot_trop <- effects_plot(density_effect_tropical, effect_tropical, intercept_tropical, "Tropical")
effect_plot_nontrop <- effects_plot(density_effect_nontropical, effect_nontropical, intercept_nontropical, "Temperate")

# Wrap the plots
figure_5 <- patchwork::wrap_plots(
  raw_plot_all, raw_plot_trop, raw_plot_nontrop,
  effect_plot_all, effect_plot_trop, effect_plot_nontrop,
  nrow = 2
)

figure_5_final <- cowplot::plot_grid(
  figure_5, legend_myco, rel_widths = c(1, 0.15)
)

# Save the plot
ggsave(
  "output/figure5.png",
  plot = figure_5_final,
  width = 16,
  height = 9.25, 
  bg = "white",
  units = "cm",
  dpi = 300
)
ggsave(
  "output/figure5.tif",
  plot = figure_5_final,
  width = 16,
  height = 9.25, 
  bg = "white",
  units = "cm"
)
