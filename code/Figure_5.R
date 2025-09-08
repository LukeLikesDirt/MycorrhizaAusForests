
# Required packages
library(terra)
library(cowplot)
library(gridExtra)
library(gtable)
library(grid)
library(ggtext)
library(tidyverse)
source("code/map_australia.R")

# Load the data:
# Env breadth is standardised to percentiles in figure 5a and 5c
data_figure_5a <- rast("output/generated_data/figure_5a.tif")
load("output/generated_data/figure_5b.RData")
load("output/generated_data/figure_5c.RData")

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

# --- Legend -------------------------------------------------------------------

# Colors
my_colors_figure_5 <- rev(paletteer::paletteer_c("grDevices::Spectral", 100))

# Dummy data
dummy_data_figure_5 <- tibble(x = 1, y = 1, z = 1:100) %>%
  mutate(group = factor(
    ifelse(z <= 50, "Tropical", "Temperate"),
    levels = c("Tropical", "Temperate"))
    )

# 1. Relative richness legend (fill gradient)
plot_breadth <- ggplot(dummy_data_figure_5, aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = my_colors_figure_5,
    name = "\nNiche\nbreadth", # <-- Here is where I create the legend spacing
    breaks = c(5, 95),
    labels = c("Low", "High")
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.ticks = element_blank()
  )

# 2. Climate zone legend (colour discrete)
plot_climate <- ggplot(dummy_data_figure_5, aes(x, y, colour = group)) +
  geom_point(shape = 15, size = 4) +
  scale_colour_manual(
    values = c("Tropical" = "#de2d26", "Temperate" = "#3366FF"),
    name = "Climate\nzone",
    guide = guide_legend(
      label.hjust = 0,
      override.aes = list(size = 4),
      keywidth = unit(1, "lines"),        # Width of legend key
      keyheight = unit(1, "lines")        # Height of legend key
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.width = unit(1, "lines"),         # Overall key area width
    legend.key.height = unit(1, "lines"),        # Overall key area height  
    legend.spacing.y = unit(0, "pt"),            # Vertical spacing between items
    legend.text = element_text(margin = margin(l = 0, unit = "pt")),  # Move text closer
    legend.margin = margin(0, 0, 0, 0)
  )

# Extract legends and convert to gtables
gt_climate <- ggplotGrob(plot_climate)$grobs[[
  which(sapply(ggplotGrob(plot_climate)$grobs, function(x) x$name) == "guide-box")
]]
gt_breadth <- ggplotGrob(plot_breadth)$grobs[[
  which(sapply(ggplotGrob(plot_breadth)$grobs, function(x) x$name) == "guide-box")
]]

# Combine legends with controlled spacing
legend_figure_5 <- gtable_rbind(gt_climate, gt_breadth, size = "max")
legend_figure_5$heights[2] <- unit(0, "pt")  # Adjust this value to control gap size

# --- Figure 5a ----------------------------------------------------------------

# Richnees of mycorrhizal types
figure_5a <- aus_map +
  tidyterra::geom_spatraster(
    data = terra::rast(list(
      AM = data_figure_5a[["AM"]],
      EcM = data_figure_5a[["EcM"]],
      Dual = data_figure_5a[["Dual"]],
      NM = data_figure_5a[["NM"]]
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
  # Add a h line to separate tropical and temperate
  geom_hline(yintercept = -23.45, linetype = "dotted", linewidth = 0.3) +
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
print(figure_5a)

# --- Figure 5b ----------------------------------------------------------------

# Remove trailing spaces from annotations
coefficients <- coefficients %>%
  mutate(annotation = str_trim(annotation))

# Then use this in your plotting code
figure_5b <- ggplot(raw_data, aes(x = predictor, y = response)) +
  # geom_hex(bins = 40, aes(alpha = after_stat(ndensity)), fill = "grey30") +
  # scale_alpha(range = c(0.075, 0.95)) +
  geom_hex(bins = 30, aes(
    alpha = after_stat(ndensity),
    fill = after_stat(ndensity)
  )) +
  scale_alpha(range = c(0.8, 1)) +
  scale_fill_gradient(
    low = "#EEEEEE",
    high = "#212121"
  ) +
  geom_ribbon(
    data = marginal_effects %>% filter(model %in% c("AM", "EcM", "Dual", "NM")),
    aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.2
  ) +
  geom_line(
    data = marginal_effects %>% filter(model %in% c("AM", "EcM", "Dual", "NM")),
    aes(y = response), colour = "black", linewidth = 0.5
  ) +
  geom_ribbon(
    data = marginal_effects %>% filter(grepl("_tropical", model)),
    aes(ymin = lower, ymax = upper), fill = "#de2d26", alpha = 0.2
  ) +
  geom_line(
    data = marginal_effects %>% filter(grepl("_tropical", model)),
    aes(y = response), colour = "#de2d26", linewidth = 0.5
  ) +
  geom_richtext(
    data = coefficients %>% filter(grepl("_tropical", model)),
    aes(x = x_pos, y = y_pos, label = annotation),
    colour = "#de2d26",
    size = 2.25, hjust = 1, vjust = 1,
    fill = NA, label.colour = NA
  ) +
  geom_ribbon(
    data = marginal_effects %>% filter(grepl("_nontropical", model)),
    aes(ymin = lower, ymax = upper), fill = "#3366FF", alpha = 0.2
  ) +
  geom_line(
    data = marginal_effects %>% filter(grepl("_nontropical", model)),
    aes(y = response), colour = "#3366FF", linewidth = 0.5
  ) +
  geom_richtext(
    data = coefficients %>% filter(grepl("_nontropical", model)),
    aes(x = x_pos, y = y_pos, label = annotation),
    colour = "#3366FF",
    size = 2.25, hjust = 0.98, vjust = 1.65,
    fill = NA, label.colour = NA
  ) +
  common_theme +
  theme(
    strip.text = element_blank()
    ) +
  ylim(0, 0.8) +
  scale_x_continuous(
    labels = function(x) paste0(abs(x), "°S")
  ) +
  labs(x = "Latitude", y = "Niche breadth", tag = "(**b**)") +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# Display the plot
print(figure_5b)

# --- Figure 5c ----------------------------------------------------------------

# bio4-bio31 environmental space plot
figure_5c <- data_figure_5c %>%
  mutate(
    # Back transform bins to the original scale
    bio4_original = data_figure_5c_limits$bio4_std_min + (bio4_bin - 0.5) * (data_figure_5c_limits$bio4_std_max - data_figure_5c_limits$bio4_std_min) / 20,
    bio31_original = data_figure_5c_limits$bio31_std_min + (bio31_bin - 0.5) * (data_figure_5c_limits$bio31_std_max - data_figure_5c_limits$bio31_std_min) / 20
  ) %>%
  pivot_longer(
    cols = -c(bio4_bin, bio31_bin, bio4_original, bio31_original),
    names_to = "mycorrhizal_type",
    values_to = "predicted"
  ) %>%
  ggplot(aes(x = bio4_original, y = bio31_original, fill = predicted)) +
  geom_tile() +
  # Add 0 lines for the axes
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
  theme(
    strip.text = element_blank()
  ) +
  labs(
    x = "Temperature seasonality",
    y = "Moisture seasonality", 
    tag = "(**c**)"
  ) +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# Display the plot
print(figure_5c)

# --- Join & save --------------------------------------------------------------

# Combine the plots
figure_5 <- patchwork::wrap_plots(
  figure_5a,
  figure_5b,
  figure_5c,
  ncol = 1
)

# Combine with main figure
figure_5_final <- plot_grid(
  figure_5, 
  legend_figure_5, 
  rel_widths = c(1, 0.14)
) +
  theme(
    plot.margin = margin(0, 0, 0, -6.5, "pt")  # Negative left margin trims left side
  )

# Display the final plot
print(figure_5_final)

# Save the final plot
ggsave(
  filename = "output/figure5.png",
  plot = figure_5_final,
  width = 16,
  height = 13.75,
  bg = "white",
  units = "cm",
  dpi = 300
)
ggsave(
  filename = "output/figure5.tif",
  plot = figure_5_final,
  width = 16,
  height = 13.5,
  bg = "white",
  units = "cm"
)

