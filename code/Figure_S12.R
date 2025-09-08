# Required packages
require(ggtext)
require(tidyverse)
source("code/functions.R")

# Set theme
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 9
common_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    axis.ticks = element_blank(),
    axis.text = element_text(size = text_size),
    axis.title = element_text(size = title_size),
    plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
    plot.tag = element_markdown(size = tag_size),
    strip.text = element_text(face = "bold", size = strip_size),
    plot.margin = margin(1, 1, 1, 1, "pt"),
    aspect.ratio = 1
  )

# Read the estimation data
data_est <- data.table::fread(
  "data/presence/sites_relative_richness_est_10.txt",
  stringsAsFactors = TRUE
) %>%
  # Compute the log ratio of richness
  mutate(
    x_km = x_albers / 1000, y_km = y_albers / 1000,
    
    # Compute the relative richness
    richness_AM_prop = (richness_AM + 0.5) / (richness + 1),
    richness_EcM_prop = (richness_EcM + 0.5) / (richness + 1),
    richness_EcM_AM_prop = (richness_EcM_AM + 0.5) / (richness + 1),
    richness_NM_prop = (richness_NM + 0.5) / (richness + 1),
    
    # Apply the log-ratio transformation
    richness_AM_log_ratio = log(richness_AM_prop / (1 - richness_AM_prop)),
    richness_EcM_log_ratio = log(richness_EcM_prop / (1 - richness_EcM_prop)),
    richness_EcM_AM_log_ratio = log(richness_EcM_AM_prop / (1 - richness_EcM_AM_prop)),
    richness_NM_log_ratio = log(richness_NM_prop / (1 - richness_NM_prop)),

  )

# Does the log-ratio transformation account for sampling effort?

# Number of sites per grid cell
sample_effort_plot_log_ratio <- data_est %>%
  select(
    n_obs, AM = richness_AM_log_ratio, EcM = richness_EcM_log_ratio, 
    Dual = richness_EcM_AM_log_ratio, NM = richness_NM_log_ratio) %>%
  pivot_longer(
    cols = -n_obs,
    names_to = "richness_type",
    values_to = "richness_log_ratio"
  ) %>%
  ggplot(aes(n_obs, richness_log_ratio)) +
  geom_point(alpha = 0.15, size = 0.3) +
  geom_smooth(method = "loess") +
  stat_smooth(method = "lm", colour = "red") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label)), colour = "red", size = 3) +
  labs(
    x = "Number of Sites (log-scale)",
    y = "Relative richness (log-ratio)",
  ) +
  scale_x_log10() +
  scale_y_continuous(
    breaks = c(-4, -2, 0, 2, 4),
    limits = c(-4.595120, 4.595120)
    ) +
  labs(
    tag = "(**a**)"
  ) +
  theme_minimal() +
  common_theme +
  theme(
    strip.text = element_text(face = "bold", size = strip_size),
    plot.tag = element_markdown(size = tag_size),
    # Adjust the tag position to account for the strip
    plot.tag.position = c(0.02, 0.93),
  ) +
  facet_wrap(~ richness_type, nrow = 1)

# How about absolute richness?

sample_effort_plot_absolute <- data_est %>%
  select(
    n_obs, AM = richness_AM, EcM = richness_EcM, 
    Dual = richness_EcM_AM, NM = richness_NM) %>%
  pivot_longer(
    cols = -n_obs,
    names_to = "richness_type",
    values_to = "richness"
  ) %>%
  # Add a pseudo-count to avoid log(0)
  mutate(richness = richness + 1) %>%
  ggplot(aes(n_obs, richness)) +
  geom_point(alpha = 0.15, size = 0.3) +
  geom_smooth(method = "loess") +
  stat_smooth(method = "lm", colour = "red") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label)), colour = "red", size = 3) +
  labs(
    x = NULL,
    y = "Richness (log-scale)",
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    tag = "(**b**)"
  ) +
  common_theme +
  theme(
    strip.text = element_blank(),
    plot.tag = element_markdown(size = tag_size),
    axis.text.x = element_blank()
  ) +
  facet_wrap(~ richness_type, nrow = 1)

# Join the plots
sample_effort_plot <- sample_effort_plot_log_ratio / sample_effort_plot_absolute

# Save the plot
ggsave(
  "output/figureS12png",
  sample_effort_plot, height = 12, width = 20, units = "cm", dpi = 300
)
