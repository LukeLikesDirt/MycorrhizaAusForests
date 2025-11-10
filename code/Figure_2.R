
# Required packages
library(terra)
library(ggtext)
library(tidyverse)
source("code/map_australia.R")

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


# Figure 2 #####################################################################

# Load the data:
# Richness (log-ratio) is standardised to percentiles in figure 2a, 2c and 2d
# Richness ratio is on the log10-scale in figure 2b
data_figure_2a <- rast("output/generated_data/figure_2a.tif")
load("output/generated_data/figure_2b.RData")
load("output/generated_data/figure_2c.RData")
load("output/generated_data/figure_2c.RData")

#### Legend ####

# Define the colours
my_colors_figure_2 <- rev(paletteer::paletteer_c("grDevices::Spectral", 100))

# Create a dummy tibble
dummy_data_figure_2 <- tibble(x = 1, y = 1, z = 1:100)

# Create a dummy plot to generate the legend
dummy_plot_figure_2 <- ggplot(dummy_data_figure_2, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = my_colors_figure_2,
    name = "Relative\nrichness",
    breaks = c(5, 95),
    labels = c("Low", "High")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.ticks = element_blank(),
    legend.title = element_text(size = title_size),
    legend.text = element_text(size = text_size)
  )

# Extract the legend from the dummy plot
legend_figure_2 <- cowplot::get_legend(dummy_plot_figure_2)

#### 2a ####

# Richnees of mycorrhizal types
figure_2a <- aus_map +
  tidyterra::geom_spatraster(
    data = terra::rast(list(
      AM = data_figure_2a[["AM"]],
      EcM = data_figure_2a[["EcM"]],
      `Dual` = data_figure_2a[["EcM-AM"]],
      NM = data_figure_2a[["NM"]]
    )),
    na.rm = TRUE
  ) +
  scale_fill_stepsn(
    colors = rev(paletteer::paletteer_c("grDevices::Spectral", 20)),
    breaks = 1:100,
    labels = c("Low", rep("", 98), "High"),
    na.value = "transparent"
  ) +
  facet_wrap(~ lyr, nrow = 1) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = strip_size),
    strip.background = element_blank(),
    plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
    axis.text = element_text(size = text_size),
    plot.tag = element_markdown(size = tag_size),
    plot.tag.location = "plot",
    # Adjust the position of the tag to account for the strip
    plot.tag.position = c(0.025, 0.875)
  ) +
  labs(x = "Longitude", y = "Latitude", tag = "(**a**)")

# Display the plot
print(figure_2a)

#### 2b ####
 
# Then use this in your plotting code
figure_2b <- ggplot(
  latitude_gradient_data,
  aes(x = latitude, y = relative_richness)
  ) +
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
  geom_hline(
    yintercept = 0, linetype = 'dotted', colour = '#de2d26'
  ) +
  geom_ribbon(
    data = latitude_gradient_marginal_effects,
    aes(ymin = lower, ymax = upper), fill = "#3366FF", alpha = 0.2
  ) +
  geom_line(
    data = latitude_gradient_marginal_effects,
    aes(y = relative_richness), colour = "#3366FF", linewidth = 0.5) +
  scale_y_continuous(
    limits = c(-2.25, 2),
    breaks = c(-2, -1, 0, 1, 2),
    # Back-transformation of the log-ratio to the original scale
    labels = c("1:100", "1:10", "1:1", "10:1", "100:1")
  ) +
  # Add the coefficient annotations
  geom_text(
    data = latitude_gradient_coeficients,
    aes(x = x_pos, y = y_pos, label = annotation),
    size = 2.25, 
    hjust = 1.05,
    vjust = 2.5
  ) +
  scale_x_continuous(
    labels = function(x) paste0(abs(x), "°S")
  ) +
  common_theme +
  theme(strip.text = element_blank()) +
  labs(x = "Latitude", y = "Relative richness", tag = "(**b**)") +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# Display the plot
print(figure_2b)

#### 2c ####

# RC1-RC2 environmental space plot
figure_2c <- data_figure_2c %>%
  mutate(
    # Back transform bins to the original scale
    RC1_original = data_figure_2c_limits$RC1_min + (RC1_bin - 0.5) * (data_figure_2c_limits$RC1_max - data_figure_2c_limits$RC1_min) / 20,
    RC2_original = data_figure_2c_limits$RC2_min + (RC2_bin - 0.5) * (data_figure_2c_limits$RC2_max - data_figure_2c_limits$RC2_min) / 20
  ) %>%
  pivot_longer(
    cols = -c(RC1_bin, RC2_bin, RC1_original, RC2_original),
    names_to = "mycorrhizal_type",
    values_to = "predicted"
  ) %>%
  ggplot(aes(x = RC1_original, y = RC2_original, fill = predicted)) +
  geom_tile() +
  # Add 0 lines for the RC1 and RC2 axes
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  scale_fill_stepsn(
    colors = rev(paletteer::paletteer_c("grDevices::Spectral", 100)),
    breaks = 1:100,
    labels = c("Low", rep("", 98), "High"),
    na.value = "transparent"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  common_theme +
  theme(strip.text = element_blank()) +
  labs(
    x = "RC1 (temperature & decomposition)",
    y = "RC2 (soil moisture)", 
    tag = "(**c**)"
  ) +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# Display the plot
print(figure_2c)

#### Save figure 2 ####

# Combine the plots
figure_2 <- patchwork::wrap_plots(
  figure_2a,
  figure_2b,
  figure_2c,
  ncol = 1
)

figure_2_final <- cowplot::plot_grid(
  figure_2, legend_figure_2, rel_widths = c(1, 0.11)
)

# Display the final plot
print(figure_2_final)

# Save the final plot
ggsave(
  filename = "output/figure2.png",
  plot = figure_2_final,
  width = 16,
  height = 13.5,
  bg = "white",
  units = "cm",
  dpi = 300
)
ggsave(
  filename = "output/figure2.tif",
  plot = figure_2_final,
  width = 16,
  height = 13.5,
  bg = "white",
  units = "cm"
)

