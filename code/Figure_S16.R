
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

# Effect size plots ############################################################

# Function to calculate axis limits for an RC
calc_limits <- function(results) {
  x_vals <- c(
    results$density_effect_all %>% pull(x),
    results$density_effect_tropical %>% pull(x),
    results$density_effect_nontropical %>% pull(x)
  )
  c(min(x_vals), max(x_vals))
}

# Function to create effect plot
create_effect_plot <- function(results, climate, rc_num, limits, plot_type = "middle") {
  # Get the appropriate data based on climate
  density_data <- switch(climate,
                         "all" = results$density_effect_all,
                         "tropical" = results$density_effect_tropical,
                         "nontropical" = results$density_effect_nontropical
  )
  
  effect_data <- switch(climate,
                        "all" = results$effect_all,
                        "tropical" = results$effect_tropical,
                        "nontropical" = results$effect_nontropical
  )
  
  intercept_data <- switch(climate,
                           "all" = results$intercept_all,
                           "tropical" = results$intercept_tropical,
                           "nontropical" = results$intercept_nontropical
  )
  
  # Create intercept annotation
  intercept_annotation <- tibble(
    x = Inf,
    y = "AM",
    label = paste0("&beta;<sub>0</sub> = ", intercept_data$intercept_formatted)
  )
  
  # Base plot
  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ggridges::geom_density_ridges(
      data = density_data,
      aes(x = x, y = mycorrhizal_type, height = y, fill = mycorrhizal_type),
      stat = "identity", scale = 0.9, colour = alpha("grey30", 0.5), alpha = 0.9
    ) +
    geom_errorbarh(
      data = effect_data,
      aes(y = mycorrhizal_type, xmin = lower, xmax = upper),
      height = 0, size = 1
    ) +
    geom_point(
      data = effect_data,
      aes(x = mean, y = mycorrhizal_type),
      shape = 21, fill = "white", stroke = 0.5, size = 1.5
    ) +
    geom_richtext(
      data = intercept_annotation,
      aes(y = y, x = x, label = label),
      hjust = 0.15, vjust = 1, size = 2.4,
      fill = NA, label.color = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = myco_colours, limits = myco_types) +
    scale_y_discrete(limits = myco_types) +
    coord_flip() +
    common_theme
  
  # Apply plot-specific styling and axis limits
  if (plot_type == "left") {
    # Left column plots (all climate)
    p <- p +
      scale_x_continuous(
        limits = limits,
        breaks = scales::pretty_breaks()
      ) +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      )
    
    # Add special labels for specific plots
    if (rc_num == 1 && climate == "all") {
      p <- p + labs(x = "", tag = "(**b**)")
    } else if (rc_num == 2 && climate == "all") {
      p <- p + labs(x = "Effect size")
    } else if (rc_num == 3 && climate == "all") {
      p <- p + 
        labs(x = "") +
        theme(
          axis.text.y = element_text(),
          plot.title = element_blank(),
          axis.title.y = element_text(margin = margin(r = -2.5))
        )
    } else {
      p <- p + labs(x = "")
    }
    
  } else if (plot_type == "middle") {
    # Middle column plots (tropical)
    p <- p +
      scale_x_continuous(
        limits = limits,
        breaks = scales::pretty_breaks()
      ) +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank()
      ) +
      labs(x = "")
    
  } else if (plot_type == "right") {
    # Right column plots (nontropical)
    p <- p +
      scale_x_continuous(
        limits = limits,
        breaks = scales::pretty_breaks(),
        sec.axis = sec_axis(~ ., name = paste0("RC", rc_num))
      ) +
      theme(
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = title_size, vjust = 0)
      ) +
      labs(x = if (rc_num == 1) "RC1" else "")
  }
  
  return(p)
}

# Calculate limits for all RCs
rc_results <- list(RC1_results, RC2_results, RC3_results)
rc_limits <- map(rc_results, calc_limits)

# Create all plots using nested loops
climates <- c("all", "tropical", "nontropical")
plot_types <- c("left", "middle", "right")
effect_plots <- list()

for (i in 1:3) {  # RC1, RC2, RC3
  for (j in 1:3) {  # all, tropical, nontropical
    plot_name <- paste0("effect_plot_RC", i, "_", climates[j])
    effect_plots[[plot_name]] <- create_effect_plot(
      results = rc_results[[i]],
      climate = climates[j],
      rc_num = i,
      limits = c(-1, 1.38),  # Using a fixed limit for cross-comparison
      #limits = rc_limits[[i]],
      plot_type = plot_types[j]
    )
  }
}

# Raw data plots ###############################################################

# Define axis limits for each RC
x_limits <- list()
for (pos in c("RC1", "RC2", "RC3")) {
  raw_min <- min(data_position[[pos]], na.rm = TRUE)
  raw_max <- max(data_position[[pos]], na.rm = TRUE)
  x_limits[[pos]] <- c(raw_min, raw_max)
}

raw_data_plot <- function(raw_data, title, position_type = "", x_lab = "", plot_tag = "", x_limits) {
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
      mean = mean(.data[[position_type]]),
      sd = sd(.data[[position_type]])
    ) %>%
    ungroup()
  
  show_y_axis <- title == "All"
  axis_limits <- x_limits[[position_type]]
  
  ggplot() +
    geom_point(
      data = raw_data,
      aes(x = .data[[position_type]], y = mycorrhizal_type, color = mycorrhizal_type),
      position = position_jitter(width = 0, height = 0.2),
      alpha = 0.2, size = 1
    ) +
    geom_errorbarh(
      data = mean_data,
      aes(y = mycorrhizal_type, xmin = mean - sd, xmax = mean + sd),
      height = 0, size = 1
    ) +
    geom_point(
      data = mean_data,
      aes(x = mean, y = mycorrhizal_type),
      shape = 21, fill = "white", stroke = 0.5, size = 1.5
    ) +
    scale_color_manual(values = myco_colours) +
    scale_x_continuous(
      limits = axis_limits,
      breaks = seq(floor(axis_limits[1]), ceiling(axis_limits[2]), by = 1)
    ) +
    coord_flip() +
    common_theme +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = if (show_y_axis) element_text() else element_blank(),
      axis.title.y = if (show_y_axis) element_text() else element_blank(),
      plot.title = if (position_type == "RC1") element_text() else element_blank(),
      plot.tag.position = c(0.05, 0.95),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = if (show_y_axis) x_lab else "", 
      title = if (position_type == "RC1") title else "",
      tag = plot_tag
    )
}

# Combine the plots ############################################################

# Figure 5
figureS16 <- patchwork::wrap_plots(
  raw_data_plot(data_position, "All", "RC1", "RC1", "(**a**)", x_limits = x_limits),
  raw_data_plot(data_position_tropical, "Tropical", "RC1", x_limits = x_limits),
  raw_data_plot(data_position_nontropical, "Temperate", "RC1", x_limits = x_limits),
  raw_data_plot(data_position, "All", "RC2", "RC2", x_limits = x_limits),
  raw_data_plot(data_position_tropical, "Tropical", "RC2", x_limits = x_limits),
  raw_data_plot(data_position_nontropical, "Temperate", "RC2", x_limits = x_limits),
  raw_data_plot(data_position, "All", "RC3", "RC3", x_limits = x_limits),
  raw_data_plot(data_position_tropical, "Tropical", "RC3", x_limits = x_limits),
  raw_plot_RC3_nontropical <- raw_data_plot(data_position_nontropical, "Temperate", "RC3", x_limits = x_limits),
  effect_plots[["effect_plot_RC1_all"]], effect_plots[["effect_plot_RC1_tropical"]], effect_plots[["effect_plot_RC1_nontropical"]],
  effect_plots[["effect_plot_RC2_all"]], effect_plots[["effect_plot_RC2_tropical"]], effect_plots[["effect_plot_RC2_nontropical"]],
  effect_plots[["effect_plot_RC3_all"]], effect_plots[["effect_plot_RC3_tropical"]], effect_plots[["effect_plot_RC3_nontropical"]],
  nrow = 6
)

# Add the legend to the alternate figure
figureS16_final <- cowplot::plot_grid(
  figureS16, legend_myco, rel_widths = c(1, 0.16)
)

# Save the plot
ggsave(
  "output/figureS16.png",
  plot = figureS16_final,
  width = 14.6,
  height = 24, 
  bg = "white",
  units = "cm",
  dpi = 300
)
ggsave(
  "output/figureS16.tif",
  plot = figureS16_final,
  width = 14.6,
  height = 24, 
  bg = "white",
  units = "cm",
)
