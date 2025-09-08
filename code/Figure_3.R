
# Required packages
library(ggtext)
library(tidyverse)

# Plot constants
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 8
common_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.grid.minor = element_line(linewidth = 0.15),
    axis.ticks = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
    plot.tag = element_markdown(size = tag_size),
    strip.text = element_text(face = "bold", size = strip_size),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    aspect.ratio = 1
  )

# Figure 3 #####################################################################

# Load the data:
load("output/generated_data/figure_3.RData")

# Take a glimpse of the data
glimpse(relative_richness_data)
glimpse(relative_richness_marginal_effects)
glimpse(relative_richness_coeficients)

# Loop through RC1 to RC4 to create marginal effects plots
marginal_effects_plots <- list()
for(i in 1:3) {
  rc_param <- paste0("RC", i)
  
  # Filter data for the current RC
  df_current <- relative_richness_data %>%
    select(mycorrhizal_type, response, all_of(rc_param)) %>%
    rename(RC = rc_param) %>%
    mutate(
      mycorrhizal_type = if_else(mycorrhizal_type == "EcM-AM", "Dual", mycorrhizal_type),
      mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
      )
  
  # For the marginal effects filtering, use the actual RC value (RC1, RC2, etc.)
  filtered_effects <- relative_richness_marginal_effects %>%
    filter(variable == rc_param) %>%
    rename(RC = predictor) %>%
    mutate(
      mycorrhizal_type = if_else(mycorrhizal_type == "EcM-AM", "Dual", mycorrhizal_type),
      mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
    )
  
  # For the coefficients filtering, use the actual RC value (RC1, RC2, etc.)
  filtered_coefficients <- relative_richness_coeficients %>%
    filter(parameter == rc_param) %>%
    rename(RC = parameter) %>%
    mutate(
      mycorrhizal_type = if_else(mycorrhizal_type == "EcM-AM", "Dual", mycorrhizal_type),
      mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
    )
  
  # Create panel label
  panel_label <- if (i == 1) {
    "(**a**)"
  } else if (i == 2) {
    "(**b**)"
  } else if (i == 3) {
    "(**c**)"
  }
  
  p <- ggplot(data = df_current, aes(x = RC, y = response)) +
    # geom_hex(bins = 40, aes(alpha = after_stat(ndensity)), fill = "black") +
    # scale_alpha(range = c(0.1, 1)) +
    geom_hex(bins = 40, aes(
      alpha = after_stat(ndensity),
      fill = after_stat(ndensity)
    )) +
    scale_alpha(range = c(0.8, 1)) +
    scale_fill_gradient(
      low = "#EEEEEE",
      high = "#212121"
    ) +
    # stat_smooth(
    #   method = "loess",
    #   colour = "#009E73",
    #   linewidth = 0.25,
    # ) +
    geom_hline(yintercept = 0, linetype = 'dotted', colour = '#de2d26', linewidth = 0.5) +
    geom_ribbon(
      data = filtered_effects,
      aes(ymin = lower, ymax = upper), fill = "#3366FF", alpha = 0.2,
    ) +
    geom_line(
      data = filtered_effects,
      aes(y = response), colour = "#3366FF", linewidth = 0.5) +
    # Add the coefficient annotations
    geom_text(
      data = filtered_coefficients,
      aes(x = x_pos, y = y_pos, label = annotation),
      size = 2.25, 
      hjust = 1.05,
      vjust = 2.5
    ) +
    facet_wrap(
      ~mycorrhizal_type, nrow = 1
    ) +
    scale_y_continuous(
      limits = c(-2.25, 2),
      breaks = c(-2, -1, 0, 1, 2),
      # Labels on the original scale
      labels = c("1:100", "1:10", "1:1", "10:1", "100:1")
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
      x = case_when(
        rc_param == "RC1" ~ "RC1 (temperature & decomposition)",
        rc_param == "RC2" ~ "RC2 (soil moisture)",
        rc_param == "RC3" ~ "RC3 (soil phosphorus)",
        TRUE ~ rc_param
      ),
      y = NULL,
      # # Add tag for annotation
      # tag = panel_label
    ) +
    theme(
      # Position the tag in the top left
      # plot.tag = element_text(size = tag_size, face = "bold", hjust = 0),
      # theme(plot.tag.position = "plot")
    )
  
  # Conditionally adjust the strip text formatting:
  if(i == 1) {
    p <- p + common_theme + theme(
      strip.text = element_text(face = "bold", size = strip_size)
    )
  } else {
    p <- p + common_theme + theme(
      strip.text = element_blank()
    )
  }
  
  # Save to the correctly named list
  marginal_effects_plots[[i]] <- p
}

# Combine the plots into a single plot
combined_plot <- patchwork::wrap_plots(
  marginal_effects_plots[[1]],
  marginal_effects_plots[[2]] + ylab("Relative richness"),
  marginal_effects_plots[[3]],
  nrow = 3
)

# # Use the tag label as a y-axis label
# combined_marginal_effects_plot <- patchwork::wrap_elements(combined_plot) +
#   labs(tag = expression("Relative richness")) +
#   theme(
#     plot.margin = margin(t = 0, r = 0, b = 0, l = 10, "pt"),
#     plot.tag = element_text(size = title_size, angle = 90),
#     plot.tag.position = "left"
#   )

# Save the plot
ggsave(
  filename = "output/figure3.png",
  plot = combined_plot,
  width = 16,
  height = 15.25,
  bg = "white",
  units = "cm",
  dpi = 300
)

# Save the plot
ggsave(
  filename = "output/figure3.tif",
  plot = combined_plot,
  width = 16,
  height = 15.25,
  bg = "white",
  units = "cm"
)

