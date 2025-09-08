# Required packages
library(ggtext)
library(tidyverse)

# Plot constants
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 8

# Load the data
load("output/generated_data/figure_2b.RData")
load("output/supplimentary_absolute_richness/main_figure_b.RData")

# Compute max log10(upper richness) for log-transformed scaling for absolute richness
absolute_max <- max(latitude_gradient_marginal_effects_abs$upper)

log_upper_max <- latitude_gradient_marginal_effects_abs %>%
  pull(upper) %>%
  log10() %>%
  max(na.rm = TRUE)

log_upper_min <- latitude_gradient_marginal_effects_abs %>%
  pull(lower) %>%
  log10() %>%
  min(na.rm = TRUE)

# Compute min and max for relative richness
relative_richness_min <- latitude_gradient_marginal_effects %>%
  pull(relative_richness) %>%
  min(na.rm = TRUE)
relative_richness_max <- latitude_gradient_marginal_effects %>%
  pull(relative_richness) %>%
  max(na.rm = TRUE)

# Define min-max scaling function
min_max_scale <- function(
    x, old_min, old_max,
    new_min = relative_richness_min, new_max = relative_richness_max) {
  if (old_max == old_min) {
    return(rep(new_min, length(x)))
  }
  new_min + (x - old_min) * (new_max - new_min) / (old_max - old_min)
}

# Prepare and combine datasets
fitted_combined <- latitude_gradient_marginal_effects_abs %>%
  select(
    mycorrhizal_type, latitude,
    fitted_absolute = richness, lower_absolute = lower, upper_absolute = upper
  ) %>%
  mutate(
    ID = row_number(),
    lower_absolute = lower_absolute,
    # Recode EcM-AM to Dual
    mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
    # Level the mycorrhizal_type
    mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
  ) %>%
  inner_join(
    latitude_gradient_marginal_effects %>%
      select(
        mycorrhizal_type, latitude,
        fitted_relative = relative_richness,
        lower_relative = lower,
        upper_relative = upper
      ) %>%
      mutate(
        ID = row_number(),
        # Recode EcM-AM to Dual
        mycorrhizal_type = recode(mycorrhizal_type, "EcM-AM" = "Dual"),
        # Level the mycorrhizal_type
        mycorrhizal_type = factor(mycorrhizal_type, levels = c("AM", "EcM", "Dual", "NM"))
      ),
    by = c("mycorrhizal_type", "ID")
  ) %>%
  mutate(
    latitude = (latitude.x + latitude.y) / 2
  ) %>%
  select(-ID, -latitude.x, -latitude.y) %>%
  mutate(
    # Regular min-max scaling for absolute richness
    fitted_absolute_scaled = min_max_scale(fitted_absolute, old_min = 0, old_max = absolute_max),
    lower_absolute_scaled = min_max_scale(lower_absolute, old_min = 0, old_max = absolute_max),
    upper_absolute_scaled = min_max_scale(upper_absolute, old_min = 0, old_max = absolute_max),
    
    # Min-max scale log10(absolute richness)
    fitted_absolute_log_scaled = min_max_scale(log10(fitted_absolute), old_min = log_upper_min, old_max = log_upper_max),
    lower_absolute_log_scaled = min_max_scale(log10(lower_absolute), old_min = log_upper_min, old_max = log_upper_max),
    upper_absolute_log_scaled = min_max_scale(log10(upper_absolute), old_min = log_upper_min, old_max = log_upper_max)
  )

# ==== Plot 1: Linear Absolute Richness (already working) ====
plot_linear <- ggplot(fitted_combined, aes(x = latitude)) +
  geom_ribbon(aes(ymin = upper_relative, ymax = lower_relative), fill = "#de2d26", alpha = 0.2) +
  geom_line(aes(y = fitted_relative), colour = "#de2d26", linewidth = 1) +
  geom_ribbon(aes(ymin = lower_absolute_scaled, ymax = upper_absolute_scaled),
              fill = "#3366FF", alpha = 0.2) +
  geom_line(aes(y = fitted_absolute_scaled),
            colour = "#3366FF", linewidth = 1) +
  scale_y_continuous(
    name = "Relative richness",
    limits = c(relative_richness_min, 0.5),
    breaks = c(-1, -0.477, 0, 0.477),
    labels = c("1:10", "1:3", "0", "3:1"),
    sec.axis = sec_axis(
      trans = ~ ((. - (relative_richness_min)) / (relative_richness_max - (relative_richness_min))) * (max(fitted_combined$upper_absolute) - 0) + 0,
      name = "Absolute richness",
      labels = function(x) format(x, big.mark = " ", scientific = FALSE),
      breaks = c(0, 10, 20, 30)
    )
  ) +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°S")) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    axis.ticks.x = element_line(colour = "grey90", size = 0.5),
    axis.ticks.y.left = element_line(color = alpha("#de2d26", 0.6), size = 0.5),
    axis.ticks.y.right = element_line(color = alpha("#3366FF", 0.6), size = 0.5),
    axis.ticks.length = unit(-2, "pt"),
    aspect.ratio = 1,
    axis.text.y.left = element_text(color = "#de2d26"),
    axis.text.y.right = element_text(color = "#3366FF"),
    strip.text = element_text(size = strip_size, face = "bold")
  ) +
  labs(x = "Latitude") +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# ==== Plot 2: Log10 Absolute Richness Min-Max Scaled ====
plot_log <- ggplot(fitted_combined, aes(x = latitude)) +
  geom_ribbon(aes(ymin = upper_relative, ymax = lower_relative), fill = "#de2d26", alpha = 0.2) +
  geom_line(aes(y = fitted_relative), colour = "#de2d26", linewidth = 1) +
  geom_ribbon(aes(ymin = lower_absolute_log_scaled, ymax = upper_absolute_log_scaled),
              fill = "#3366FF", alpha = 0.2) +
  geom_line(aes(y = fitted_absolute_log_scaled),
            colour = "#3366FF", linewidth = 1) +
  scale_y_continuous(
    name = "Relative richness",
    limits = c(relative_richness_min - 0.2, 0.55),
    breaks = c(-1, -0.477, 0, 0.477),
    labels = c("1:10", "1:3", "0", "3:1"),
    sec.axis = sec_axis(
      trans = ~ 10^(((. - (relative_richness_min)) / (relative_richness_max - (relative_richness_min))) * (log_upper_max - log_upper_min) + log_upper_min),
      name = "Absolute richness",
      breaks = c(1, 3, 10, 30),
      labels = function(x) format(x, big.mark = " ", scientific = FALSE)
    )
  ) +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°S")) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    axis.ticks.x = element_line(colour = "grey90", size = 0.5),
    axis.ticks.y.left = element_line(color = alpha("#de2d26", 0.6), size = 0.5),
    axis.ticks.y.right = element_line(color = alpha("#3366FF", 0.6), size = 0.5),
    axis.ticks.length = unit(-2, "pt"),
    aspect.ratio = 1,
    axis.text.y.left = element_text(color = "#de2d26"),
    axis.text.y.right = element_text(color = "#3366FF"),
    strip.text = element_text(size = strip_size, face = "bold")
  ) +
  labs(x = "Latitude") +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# Display the plot
print(plot_linear)
print(plot_log)

# Save the plots
ggsave(
  "output/figureS7.png",
  plot = plot_log,
  bg = "white",
  width = 16, height = 5.25, units = "cm", dpi = 300
)
