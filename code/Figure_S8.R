
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

# Figure 2 #####################################################################

# Load the data:
data_figure_2a <- rast("output/supplimentary_absolute_richness/main_figure_a.tif")
load("output/supplimentary_absolute_richness/main_figure_b.RData")
load("output/supplimentary_absolute_richness/main_figure_c.RData")

#### Legend ####

# Define the colours
my_colours <- rev(paletteer::paletteer_c("grDevices::Spectral", 100))

# Create a dummy tibble
dummy_data <- tibble(x = 1, y = 1, z = 1:100)

# Create a dummy plot to generate the legend
dummy_plot <- ggplot(dummy_data, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = my_colours,
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
legend <- cowplot::get_legend(dummy_plot)

#### Figure a ####

# Richnees of mycorrhizal types
figure_a <- aus_map +
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
    #panel.grid = element_line(colour = "white", linewidth = 0.125),
    axis.text = element_text(size = text_size),
    plot.tag = element_markdown(size = tag_size),
    plot.tag.location = "plot",
    # Adjust the position of the tag to account for the strip
    plot.tag.position = c(0.025, 0.875)
  ) +
  labs(x = "Longitude", y = "Latitude", tag = "(**a**)")

# Display the plot
print(figure_a)

#### Figure b ####

# Then use this in your plotting code
figure_b <- ggplot(
  latitude_gradient_data_abs %>%
    # Add 1 to avoid log(0)
    mutate(richness = richness + 1),
  aes(x = latitude, y = richness)
) +
  # geom_hex(bins = 30, aes(alpha = after_stat(ndensity)), fill = "black") +
  # scale_alpha(range = c(0.1, 1)) +
  geom_hex(bins = 40, aes(
    alpha = after_stat(ndensity),
    fill = after_stat(ndensity)
  )) +
  scale_alpha(range = c(0.6, 1)) +
  scale_fill_gradient(
    low = "#EEEEEE",
    high = "#212121"
  ) +
  geom_ribbon(
    data = latitude_gradient_marginal_effects_abs,
    aes(ymin = lower, ymax = upper), fill = "#3366FF", alpha = 0.2
  ) +
  geom_line(
    data = latitude_gradient_marginal_effects_abs,
    aes(y = richness), colour = "#3366FF", linewidth = 0.5) +
  # Add the coefficient annotations
  geom_text(
    data = latitude_gradient_coeficients_abs,
    aes(x = x_pos, y = y_pos, label = annotation),
    size = 2.25, 
    hjust = 1.05,
    vjust = 2.5
  ) +
  scale_x_continuous(
    labels = function(x) paste0(abs(x), "°S")
  ) +
  scale_y_log10(
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    aspect.ratio = 1,
    strip.text = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.tag = element_markdown(size = tag_size),
    plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
  ) +
  labs(x = "Latitude", y = "Absolute richness", tag = "(**b**)") +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# Display the plot
print(figure_b)

#### Figure c ####

# RC1-RC2 environmental space plot
figure_c <- data_figure_c %>%
  mutate(
    # Back transform bins to the original scale
    RC1_original = data_figure_c_limits$RC1_min + (RC1_bin - 0.5) * (data_figure_c_limits$RC1_max - data_figure_c_limits$RC1_min) / 20,
    RC2_original = data_figure_c_limits$RC2_min + (RC2_bin - 0.5) * (data_figure_c_limits$RC2_max - data_figure_c_limits$RC2_min) / 20
  ) %>%
  pivot_longer(
    cols = -c(RC1_bin, RC2_bin, RC1_original, RC2_original),
    names_to = "mycorrhizal_type",
    values_to = "predicted"
  ) %>%
  ggplot(aes(x = RC1_original, y = RC2_original, fill = predicted)) +
  geom_tile() +
  scale_fill_stepsn(
    colors = rev(paletteer::paletteer_c("grDevices::Spectral", 100)),
    breaks = 1:100,
    labels = c("Low", rep("", 98), "High"),
    na.value = "transparent"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "grey90", fill = NA),
    legend.position = "none",
    strip.text = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.tag = element_markdown(size = tag_size),
    aspect.ratio = 1,
    plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
  ) +
  labs(
    x = "RC1 (temperature & decomposition)",
    y = "RC2 (soil moisture)", 
    tag = "(**c**)"
  ) +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# Display the plot
print(figure_c)

#### Save figure ####

# Combine the plots
figure <- patchwork::wrap_plots(
  figure_a,
  figure_b,
  figure_c,
  ncol = 1
)

figure_final <- cowplot::plot_grid(
  figure, legend, rel_widths = c(1, 0.11)
)

# Display the final plot
print(figure_final)

# Save the final plot
ggsave(
  filename = "output/figureS8.png",
  plot = figure_final,
  width = 16,
  height = 13.5,
  bg = "white",
  units = "cm",
  dpi = 300
)
