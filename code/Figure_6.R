
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

# Helper function to visualise marginal effects
create_marginal_effects_plot <- function(raw_data, marginal_effects,
                                         coefficients, var) {
  
  ggplot(
    raw_data %>% filter(variable == var),
    aes(x = predictor, y = response)) +
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
      data = marginal_effects %>% filter(variable == var),
      aes(ymin = lower, ymax = upper), fill = "#3366FF", alpha = 0.2
    ) +
    geom_line(
      data = marginal_effects %>% filter(variable == var),
      aes(y = response), colour = "#3366FF", linewidth = 0.5
    ) +
    geom_richtext(
      data = coefficients %>% filter(parameter == var),
      aes(x = x_pos, y = y_pos, label = annotation),
      size = 2.5, hjust = 1, vjust = 1,
      fill = NA, label.colour = NA
    ) +
    common_theme +
    scale_x_continuous(
      breaks = if(var == "richness") seq(-4, 4, 2) else seq(-4, 4, 1),
    ) +
    ylim(0, 0.7) +
    labs(
      x = NULL, y = NULL
    ) + 
    facet_wrap(
      ~mycorrhizal_type, nrow = 1,
    )
}

# Figure 3 #####################################################################

# Load the data:
load("output/generated_data/figure_6.RData")

# Organise the raw data
raw_data <- raw_data %>%
  filter(
    # Use the standardised variables
    str_detect(variable, "_std")
  ) %>%
  mutate(
    # Remove the "_std" suffix from the variable names
    variable = str_remove(variable, "_std"),
    # Level mycorrhizal type
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "EcM", "Dual", "NM")
  ))

# Organise the marginal effects data
marginal_effects_data <- marginal_effects_data %>%
  # Use the standardised variables
  select(
    mycorrhizal_type, variable, response,
    lower, upper,
    predictor = predictor_std
    ) %>%
  mutate(
    # Level the mycorrhizal type
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "EcM", "Dual", "NM")
  ))

# Level the mycorrhizal type in the coefficients
coeficients_data <- coeficients_data %>%
  mutate(
    mycorrhizal_type = factor(
      mycorrhizal_type,
      levels = c("AM", "EcM", "Dual", "NM")
  ))

# Create the marginal effects plots
figure_6 <- patchwork::wrap_plots(
  create_marginal_effects_plot(
    raw_data,
    marginal_effects_data,
    coeficients_data,
    "bio4"
  ) +
    xlab("Temperature seasonality"),
  
  create_marginal_effects_plot(
    raw_data,
    marginal_effects_data,
    coeficients_data,
    "bio31"
  ) +
    theme(strip.text = element_blank()) +
    xlab("Moisture index seasonality") +
    ylab("Environmental breadth")
    ,
  
  create_marginal_effects_plot(
    raw_data,
    marginal_effects_data,
    coeficients_data,
    "richness"
  ) +
    theme(strip.text = element_blank()) +
    xlab("Tree richness"),
  
  ncol = 1
)

# Display the plot
print(figure_6)

# Save the plot
ggsave(
  filename = "output/figure6.png",
  plot = figure_6,
  width = 16,
  height = 15.25,
  bg = "white",
  units = "cm",
  dpi = 300
)

# Save the plot
ggsave(
  filename = "output/figure6.tif",
  plot = figure_6,
  width = 16,
  height = 15.25,
  bg = "white",
  units = "cm"
)
