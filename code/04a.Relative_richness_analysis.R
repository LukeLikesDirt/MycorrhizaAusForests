# Required packages
require(INLA) # inla.list.models()
require(fmesher)
require(terra)
require(ggtext)
require(tidyverse)
source("code/functions.R")
source("code/map_australia.R")

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

# (1) Read data ################################################################

# Forest raster
forest_raster <- rast("data/aus_forests_23/aus_for23_masked_10.tif")

# Read the estimation data
data_est <- data.table::fread(
  "data/presence/sites_relative_richness_est_10.txt",
  stringsAsFactors = TRUE
) %>%
  # Compute the log ratio of richness
  mutate(
    x_km = x_albers / 1000, y_km = y_albers / 1000,
    richness_AM_log_ratio = log((richness_AM + 1) / (richness - richness_AM + 1)),
    richness_EcM_log_ratio = log((richness_EcM + 1) / (richness - richness_EcM + 1)),
    richness_EcM_AM_log_ratio = log((richness_EcM_AM + 1) / (richness - richness_EcM_AM + 1)),
    richness_NM_log_ratio = log((richness_NM + 1) / (richness - richness_NM + 1)),
  )

# Count the number of cells for estimation
n_est <- nrow(data_est) %>%
  print()

# Read in the covariates for prediction
data_pred <- data.table::fread(
  "data/presence/sites_relative_richness_pred_10.txt",
  stringsAsFactors = TRUE
) %>%
  mutate(
    x_km = x_albers / 1000, y_km = y_albers / 1000,
  )

# Count the number of cells for prediction
n_pred <- nrow(data_pred) %>%
  print()

# Standardise the covariates (rotated components are off centre)
RC1_mean <- mean(c(data_est$RC1, data_pred$RC1), na.rm = TRUE)
RC1_sd <- sd(c(data_est$RC1, data_pred$RC1), na.rm = TRUE)
RC2_mean <- mean(c(data_est$RC2, data_pred$RC2), na.rm = TRUE)
RC2_sd <- sd(c(data_est$RC2, data_pred$RC2), na.rm = TRUE)
RC3_mean <- mean(c(data_est$RC3, data_pred$RC3), na.rm = TRUE)
RC3_sd <- sd(c(data_est$RC3, data_pred$RC3), na.rm = TRUE)

# Apply standardization to both datasets using combined parameters
data_est <- data_est %>%
  mutate(
    RC1 = (RC1 - RC1_mean) / RC1_sd,
    RC2 = (RC2 - RC2_mean) / RC2_sd,
    RC3 = (RC3 - RC3_mean) / RC3_sd
  )

data_pred <- data_pred %>%
  mutate(
    RC1 = (RC1 - RC1_mean) / RC1_sd,
    RC2 = (RC2 - RC2_mean) / RC2_sd,
    RC3 = (RC3 - RC3_mean) / RC3_sd
  )

# (2) Model formulas ###########################################################

# Plot richness log ratio against the covariates to check for non-linearities
data_est %>%
  select(RC1, RC2, RC3, richness_AM_log_ratio) %>%
  pivot_longer(cols = -richness_AM_log_ratio) %>%
  ggplot(aes(value, richness_AM_log_ratio)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ name, scales = "free")
data_est %>%
  select(RC1, RC2, RC3, richness_EcM_log_ratio) %>%
  pivot_longer(cols = -richness_EcM_log_ratio) %>%
  ggplot(aes(value, richness_EcM_log_ratio)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ name, scales = "free")
data_est %>%
  select(RC1, RC2, RC3, richness_EcM_AM_log_ratio) %>%
  pivot_longer(cols = -richness_EcM_AM_log_ratio) %>%
  ggplot(aes(value, richness_EcM_AM_log_ratio)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ name, scales = "free")
data_est %>%
  select(RC1, RC2, RC3, richness_NM_log_ratio) %>%
  pivot_longer(cols = -richness_NM_log_ratio) %>%
  ggplot(aes(value, richness_NM_log_ratio)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ name, scales = "free")

# Does the log-ratio transformation account for sampling effort?

# Number of sites per grid cell
sample_effort_plot <- data_est %>%
  select(
    n_obs, AM = richness_AM_log_ratio, EcM = richness_EcM_log_ratio, 
    Dual = richness_EcM_AM_log_ratio, NM = richness_NM_log_ratio) %>%
  pivot_longer(
    cols = -n_obs,
    names_to = "richness_type",
    values_to = "richness_log_ratio"
  ) %>%
  ggplot(aes(n_obs, richness_log_ratio)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  #stat_smooth(method = "lm", colour = "red") +
  #ggpubr::stat_cor(aes(label = after_stat(rr.label)), colour = "red", size = 3) +
  labs(
    x = "Number of Sites (log-scale)",
    y = "Relative richness (log-ratio)",
  ) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    aspect.ratio = 1,
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "pt")
  ) +
  facet_wrap(~ richness_type, ncol = 2)

# Save the plot
ggsave(
  "output/supplimentary_richness_log_ratio/sample_effort_plot.png",
  sample_effort_plot, height = 15, width = 15, units = "cm", dpi = 300
)

# Notes:
#   - For random walk models, set the recommended pc.prior and scale the random 
#     walk for a smoother fit that is less prone to over-fitting
#   - 'rw2' is less responsive to abrupt changes compared to the 'rw1'

# Set prior on precision:
#   - rw1 = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
#   - rw2 = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))
#pcprec <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))

# Define the formulas
formula_AM <- as.formula(paste0(
  "y_AM ~ -1 + intercept + ",
  "RC1 + RC2 + RC3 +",
  "f(spatial, model = spde)"
))
formula_EcM <- as.formula(paste0(
  "y_EcM ~ -1 + intercept + ",
  "RC1 + RC2 + RC3 +",
  "f(spatial, model = spde)"
))

formula_EcM_AM <- as.formula(paste0(
  "y_EcM_AM ~ -1 + intercept + ",
  "RC1 + RC2 + RC3 +",
  "f(spatial, model = spde)"
))
formula_NM <- as.formula(paste0(
  "y_NM ~ -1 + intercept + ",
  "RC1 + RC2 + RC3 +",
  "f(spatial, model = spde)"
))

# (3) Prepare the mesh #########################################################

# Main tools to control the mesh:
#   (1) 'max.edge'    The largest allowed triangle length.
#                     The lower the value for 'max.edge' the higher the resolution. 
#   (2) 'quantileoff' Sites with a distance value smaller than the quantileoff 
#                     are modelled by a single vertex in the mesh.
#   (3) 'boundary'    Useful for islands and fjords.
#   (4) 'offset'      Determines how far the inner and outer boundaries extend.

# Site locations for the estimation and prediction data
loc_est <- data_est %>%
  select(x_km, y_km) %>%
  as.matrix()
loc_pred <- data_pred %>%
  select(x_km, y_km) %>%
  as.matrix()

# Distances between estimation site.
hist(dist(loc_est), main = "", xlab = "Distance between sites (km)")
# Small scale looks to 2000-3000 km at the continental scale and 0 - 1000 km 
# at the regional scale.
small_scale <- 300

# Maximun edge for the triangulation (about 1/3 of the spatial correlation range)
max_edge <- small_scale / 5

# Prepare the mesh
mesh <- fm_mesh_2d_inla(
  loc = loc_est,
  max.edge = c(1, 5) * max_edge,
  cutoff = max_edge / 5
)

# Define the mesh CSR
fm_crs(mesh) <- crs(aus_map_albers$data)

# Visualise the mesh using ggplot
ggplot() +
  theme_minimal() +
  fmesher::geom_fm(data = mesh) +
  geom_point(
    data = data_est,
    aes(x_km, y_km), size = 0.2, alpha = 0.5
  ) +
  geom_sf(
    data = aus_map_albers$data,
    fill = "transparent",
    col = "red"
  ) +
  labs(
    x = NULL,
    y = NULL
  )

# Define projector matrix for the estimation and prediction data
A_est <- inla.spde.make.A(
  mesh = mesh,
  loc = loc_est
)
dim(A_est)
A_pred <- inla.spde.make.A(
  mesh = mesh,
  loc = loc_pred
)
dim(A_pred)

# (4) Define the SPDE ##########################################################

# Range and sigma estimations for the priors
range_est <- 80
sigma_est <- 1

# Define the SPDE for the estimation
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  # Range of the spatial correlation is unlikely to be less range_est at p=0.01
  prior.range = c(range_est, 0.01),
  # Prior for the standard deviation is unlikely to be larger than sigma_est at p=0.01
  prior.sigma = c(sigma_est, 0.01)
)

# Define the SPDE index
index_spatial <- inla.spde.make.index(
  name = "spatial",
  n.spde = spde$n.spde
)

# The size of the mesh and n.spde should match
mesh$n
spde$n.spde

# (5) Stack the data ###########################################################

# Covariates for the estimation and prediction (spatial and marginal effects)
covariates <- bind_rows(
  tibble(
    RC1 = c(data_est$RC1, data_pred$RC1),
    RC2 = c(data_est$RC2, data_pred$RC2),
    RC3 = c(data_est$RC3, data_pred$RC3),
    ecoregion = c(data_est$ecoregion, data_pred$ecoregion),
    tag = c(rep("est", n_est), rep("pred", n_pred))
  ),
  tibble(
    RC1 = seq(
      min(data_est$RC1), max(data_est$RC1), length.out = 100
    ),
    RC2 = median(data_est$RC2),
    RC3 = median(data_est$RC3),
    tag = "pred_RC1"
  ),
  tibble(
    RC1 = median(data_est$RC1),
    RC2 = seq(
      min(data_est$RC2), max(data_est$RC2), length.out = 100
    ),
    RC3 = median(data_est$RC3),
    tag = "pred_RC2"
  ),
  tibble(
    RC1 = median(data_est$RC1),
    RC2 = median(data_est$RC2),
    RC3 = seq(
      min(data_est$RC3), max(data_est$RC3), length.out = 100
    ),
    tag = "pred_RC3"
  )
)

# Stack the estimation data
stack_est <- inla.stack(
  tag = "est",
  data = list(
    y_AM = data_est$richness_AM_log_ratio,
    y_EcM = data_est$richness_EcM_log_ratio,
    y_EcM_AM = data_est$richness_EcM_AM_log_ratio,
    y_NM = data_est$richness_NM_log_ratio
  ),
  # The order needs to matches the order in effects
  A = list(1, 1, A_est),
  effects = list(
    intercept = rep(1, n_est),
    X = covariates %>% filter(tag == "est") %>% select(-tag),
    spatial = index_spatial
  )
)

# Stack the prediction data
stack_pred <- inla.stack(
  tag = "pred",
  data = list(
    y_AM = rep(NA_real_, n_pred),
    y_EcM = rep(NA_real_, n_pred),
    y_EcM_AM = rep(NA_real_, n_pred),
    y_NM = rep(NA_real_, n_pred)
  ),
  # The order needs to matches the order in effects
  A = list(1, 1, A_pred),
  effects = list(
    intercept = rep(1, n_pred),
    X = covariates %>% filter(tag == "pred") %>% select(-tag),
    spatial = index_spatial
  )
)

# Stack the data for the marginal effects
stack_RC1 <- inla.stack(
  tag = "pred_RC1",
  data = list(
    y_AM = rep(NA_real_, 100),
    y_EcM = rep(NA_real_, 100),
    y_EcM_AM = rep(NA_real_, 100),
    y_NM = rep(NA_real_, 100)
  ),
  # The order needs to matches the order in effects
  A = list(1, 1),
  effects = list(
    intercept = rep(1, 100),
    X = covariates %>% filter(tag == "pred_RC1") %>% select(-tag)
  )
)
stack_RC2 <- inla.stack(
  tag = "pred_RC2",
  data = list(
    y_AM = rep(NA_real_, 100),
    y_EcM = rep(NA_real_, 100),
    y_EcM_AM = rep(NA_real_, 100),
    y_NM = rep(NA_real_, 100)
  ),
  # The order needs to matches the order in effects
  A = list(1, 1),
  effects = list(
    intercept = rep(1, 100),
    X = covariates %>% filter(tag == "pred_RC2") %>% select(-tag)
  )
)
stack_RC3 <- inla.stack(
  tag = "pred_RC3",
  data = list(
    y_AM = rep(NA_real_, 100),
    y_EcM = rep(NA_real_, 100),
    y_EcM_AM = rep(NA_real_, 100),
    y_NM = rep(NA_real_, 100)
  ),
  # The order needs to matches the order in effects
  A = list(1, 1),
  effects = list(
    intercept = rep(1, 100),
    X = covariates %>% filter(tag == "pred_RC3") %>% select(-tag)
  )
)

# Join the stacks
stack <- inla.stack.join(
  stack_est,
  stack_pred,
  stack_RC1,
  stack_RC2,
  stack_RC3
)

# Retrieve the stack indexes
index_est <- inla.stack.index(stack, "est")$data
index_pred <- inla.stack.index(stack, "pred")$data
index_RC1 <- inla.stack.index(stack, "pred_RC1")$data
index_RC2 <- inla.stack.index(stack, "pred_RC2")$data
index_RC3 <- inla.stack.index(stack, "pred_RC3")$data

# (6) Fit the models ###########################################################

# AM
model_AM <- inla(
  formula_AM,
  family = 'gaussian',
  data = inla.stack.data(stack),
  control.predictor = list(
    compute = TRUE, link = 1,
    A = inla.stack.A(stack)
  ),
  control.compute = control_compute
)

# EcM
model_EcM <- inla(
  formula_EcM,
  family = 'gaussian',
  data = inla.stack.data(stack),
  control.predictor = list(
    compute = TRUE, link = 1,
    A = inla.stack.A(stack)
  ),
  control.compute = control_compute
)

# EcM-AM
model_EcM_AM <- inla(
  formula_EcM_AM,
  family = 'gaussian',
  data = inla.stack.data(stack),
  control.predictor = list(
    compute = TRUE, link = 1,
    A = inla.stack.A(stack)
  ),
  control.compute = control_compute
)

# NM
model_NM <- inla(
  formula_NM,
  family = 'gaussian',
  data = inla.stack.data(stack),
  control.predictor = list(
    compute = TRUE, link = 1,
    A = inla.stack.A(stack)
  ),
  control.compute = control_compute
)

# (7) Model diagnostics ######################################################

# Compare sigma
posterior_marginals(model_AM, model_EcM, model_EcM_AM, model_NM)

#### * Compute diagnostic data * ####

# Diagnostic data: Compute fitted values, residuals, and posterior predictive 
# checks
data_diagnostics <- inla_spat_diagnostic_data(
  model_names = c("model_AM", "model_EcM", "model_EcM_AM", "model_NM"),
  stack_names = c("stack_est", "stack_est", "stack_est", "stack_est"),
  observed_values = list(
    data_est$richness_AM_log_ratio,
    data_est$richness_EcM_log_ratio,
    data_est$richness_EcM_AM_log_ratio,
    data_est$richness_NM_log_ratio
  )
) %>%
  mutate(
    model = if_else(model == "Model EcM_AM", "Model EcM-AM", model)
  )

#### * Normality of residuals * ####

normality_plot <- ggplot(data_diagnostics, aes(x = residuals)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey",
    colour = "black"
  ) +
  geom_density(
    colour = "red",
    linewidth = 0.5
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    aspect.ratio = 1,
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "pt"),
    strip.background = element_rect(fill = "grey90", colour = "grey90")
  ) +
  facet_wrap(~ model, nrow = 1) +
  labs(tag = "a", x = "Residuals", y = "Density")
print(normality_plot)

#### * Homogeneity of variance * ####

homogeneity_plot <- ggplot(data_diagnostics, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.3, shape = 1) +
  geom_smooth(method = "loess", colour = "red") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  facet_wrap(~ model, nrow = 1, scales = "free_x") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    aspect.ratio = 1,
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "pt"),
    strip.text = element_blank()
  ) +
  labs(tag = "b", x = "Fitted values", y = "Residuals")

# Display the plot
print(homogeneity_plot)

#### * Observed vs predicted * ####
limits_list <- list(
  "Model AM" = c(
    min(data_diagnostics %>% filter(model == "Model AM") %>% pull(observed), 
        data_diagnostics %>% filter(model == "Model AM") %>% pull(fitted)),
    max(data_diagnostics %>% filter(model == "Model AM") %>% pull(observed), 
        data_diagnostics %>% filter(model == "Model AM") %>% pull(fitted))
  ),
  "Model EcM" = c(
    min(data_diagnostics %>% filter(model == "Model EcM") %>% pull(observed), 
        data_diagnostics %>% filter(model == "Model EcM") %>% pull(fitted)),
    max(data_diagnostics %>% filter(model == "Model EcM") %>% pull(observed), 
        data_diagnostics %>% filter(model == "Model EcM") %>% pull(fitted))
  ),
  "Model EcM_AM" = c(
    min(data_diagnostics %>% filter(model == "Model EcM-AM") %>% pull(observed),
        data_diagnostics %>% filter(model == "Model EcM-AM") %>% pull(fitted)),
    max(data_diagnostics %>% filter(model == "Model EcM-AM") %>% pull(observed),
        data_diagnostics %>% filter(model == "Model EcM-AM") %>% pull(fitted))
  ),
  "Model NM" = c(
    min(data_diagnostics %>% filter(model == "Model NM") %>% pull(observed), 
        data_diagnostics %>% filter(model == "Model NM") %>% pull(fitted)),
    max(data_diagnostics %>% filter(model == "Model NM") %>% pull(observed), 
        data_diagnostics %>% filter(model == "Model NM") %>% pull(fitted))
  )
)

# Observed vs predicted plot
obs_fit_plot <- ggplot(data_diagnostics, aes(x = posterior_mean, y = observed)) +
  geom_point(alpha = 0.3, shape = 1, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour = "red") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label)), colour = "red", size = 3) +
  facet_wrap(~ model, nrow = 1, scales = "free") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    aspect.ratio = 1,
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "pt"),
    strip.text = element_blank()
  ) +
  labs(tag = "c", x = "Posterior mean predictions", y = "Observed values") +
  ggh4x::facetted_pos_scales(
    x = list(
      "Model AM" = scale_x_continuous(limits = limits_list[["Model AM"]]),
      "Model EcM" = scale_x_continuous(limits = limits_list[["Model EcM"]]),
      "Model EcM-AM" = scale_x_continuous(limits = limits_list[["Model EcM_AM"]]),
      "Model NM" = scale_x_continuous(limits = limits_list[["Model NM"]])
    ),
    y = list(
      "Model AM" = scale_y_continuous(limits = limits_list[["Model AM"]]),
      "Model EcM" = scale_y_continuous(limits = limits_list[["Model EcM"]]),
      "Model EcM-AM" = scale_y_continuous(limits = limits_list[["Model EcM_AM"]]),
      "Model NM" = scale_y_continuous(limits = limits_list[["Model NM"]])
    )
  )

# Display the plot
print(obs_fit_plot)

#### * Residuals vs predictors * ####

tibble(
  residuals = data_diagnostics %>% filter(model == "Model NM") %>% pull("residuals"),
  RC1 = data_est[["RC1"]],
  RC2 = data_est[["RC2"]],
  RC3 = data_est[["RC3"]]
) %>%
  pivot_longer(
    -residuals
  ) %>%
  ggplot(aes(x = value, y = residuals)) +
  geom_point() +
  stat_smooth(method = "loess") +
  facet_wrap(~ name, scales = "free_x") +
  theme_minimal() +
  theme(aspect.ratio = 1)

#### * Save the diagnostic plots * ####

# Join the plots
diagnostic_plot <- patchwork::wrap_plots(
  normality_plot,
  homogeneity_plot,
  obs_fit_plot,
  ncol = 1
)

# Save the plots
ggsave(
  "output/supplimentary_richness_log_ratio/validation_plot.png",
  diagnostic_plot, 
  height = 14.75, width = 15, units = "cm",
  dpi = 300
)

# (8) Summarise model outputs #######################################################

#### * Fixed effects: RCs * ####

# Extract model coefficients
fixed_effects_data <- inla_fixed_effects(
  model_names = c("model_AM", "model_EcM", "model_EcM_AM", "model_NM")
) %>%
  rename(median = `0.5quant`, lower = `0.025quant`, upper = `0.975quant`) %>%
  mutate(
    model = case_when(
      model == "Model AM"     ~ "AM",
      model == "Model EcM"    ~ "EcM",
      model == "Model EcM_AM" ~ "EcM-AM",
      model == "Model NM"     ~ "NM",
      TRUE ~ model
    ),
    colour = case_when(
      lower > 0 & upper > 0 ~ "#3366FF",
      lower < 0 & upper < 0 ~ "#ca0020",
      TRUE ~ "grey"
    ),
    parameter = factor(parameter, levels = rev(levels(factor(parameter))))
  )

#### * Predicted values * ####

# Define the model names and get observed values
model_info <- list(
  AM = list(
    name = "model_AM",
    observed = data_est$richness_AM_log_ratio
  ),
  EcM = list(
    name = "model_EcM",
    observed = data_est$richness_EcM_log_ratio
  ),
  EcM_AM = list(
    name = "model_EcM_AM",
    observed = data_est$richness_EcM_AM_log_ratio
  ),
  NM = list(
    name = "model_NM",
    observed = data_est$richness_NM_log_ratio
  )
)

# Create an empty list to store the predicted values for each model
fitted_values_list <- list()

# Loop through each model
for (model_type in names(model_info)) {
  # Get the current model object
  current_model <- get(model_info[[model_type]]$name)
  
  # Create initial predicted values tibble
  current_fitted <- tibble(
    model = model_type,
    predicted = c(
      current_model$summary.fitted$mean[index_est],
      current_model$summary.fitted$mean[index_pred]
    ),
    upper = c(
      current_model$summary.fitted$`0.975quant`[index_est],
      current_model$summary.fitted$`0.975quant`[index_pred]
    ),
    lower = c(
      current_model$summary.fitted$`0.025quant`[index_est],
      current_model$summary.fitted$`0.025quant`[index_pred]
    ),
    sd = c(
      current_model$summary.fitted$sd[index_est],
      current_model$summary.fitted$sd[index_pred]
    ),
    observed = c(
      model_info[[model_type]]$observed,
      rep(NA_integer_, n_pred)
    ),
    latitude = c(
      data_est$latitude,
      data_pred$latitude
    ),
    longitude = c(
      data_est$longitude,
      data_pred$longitude
    ),
    x_albers = c(
      data_est$x_albers,
      data_pred$x_albers
    ),
    y_albers = c(
      data_est$y_albers,
      data_pred$y_albers
    ),
  ) %>%
    mutate(
      uncertainty = upper - lower,
      residuals = observed - predicted,
      pearson_residuals = residuals / sd,
      x_albers = as.numeric(x_albers),
      y_albers = as.numeric(y_albers)
    )
  
  # Create projection matrix
  A_posterior <- fmesher::fm_basis(
    mesh,
    loc = matrix(
      c(current_fitted$x_albers / 1000, current_fitted$y_albers / 1000),
      ncol = 2
    )
  )
  
  # Extract spatial field for current model
  spatial_field <- tibble(
    spatial_field = as.vector(A_posterior %*% current_model$summary.random$spatial[,"mean"]),
    x_albers = current_fitted$x_albers,
    y_albers = current_fitted$y_albers
  )
  
  # Join spatial field to predicted values
  fitted_values_list[[model_type]] <- current_fitted %>%
    left_join(spatial_field, by = c("x_albers", "y_albers"))
}

# Combine all models into one dataframe
fitted_values <- bind_rows(fitted_values_list)

#### * Rasterise predicted values * ####

# AM raster
fitted_rast_AM <- rasterize(
  vect(
    fitted_values %>% filter(model == "AM"),
    geom = c("x_albers", "y_albers"),
    crs = crs(forest_raster)
  ),
  forest_raster,
  field = names(fitted_values %>% select(-latitude, -longitude, -x_albers, -y_albers)),
  background = NA
) %>%
  project(crs(aus_map$data))

# EcM raster
fitted_rast_EcM <- rasterize(
  vect(
    fitted_values %>% filter(model == "EcM"),
    geom = c("x_albers", "y_albers"),
    crs = crs(forest_raster)
  ),
  forest_raster,
  field = names(fitted_values %>% select(-latitude, -longitude, -x_albers, -y_albers)),
  background = NA
) %>%
  project(crs(aus_map$data))

# EcM-AM raster
fitted_rast_EcM_AM <- rasterize(
  vect(
    fitted_values %>% filter(model == "EcM_AM"),
    geom = c("x_albers", "y_albers"),
    crs = crs(forest_raster)
  ),
  forest_raster,
  field = names(fitted_values %>% select(-latitude, -longitude, -x_albers, -y_albers)),
  background = NA
) %>%
  project(crs(aus_map$data))

# NM raster
fitted_rast_NM <- rasterize(
  vect(
    fitted_values %>% filter(model == "NM"),
    geom = c("x_albers", "y_albers"),
    crs = crs(forest_raster)
  ),
  forest_raster,
  field = names(fitted_values %>% select(-latitude, -longitude, -x_albers, -y_albers)),
  background = NA
) %>%
  project(crs(aus_map$data))

#### * Marginal effects * ####

# Create an empty list to store marginal effects for each model
marginal_effects_list <- list()

# Loop through each model
for (model_type in names(model_info)) {
  # Get the current model object
  current_model <- get(model_info[[model_type]]$name)
  
  # Compute marginal effects for current model
  marginal_effects_list[[model_type]] <- bind_rows(
    tibble(
      response = current_model[["summary.linear.predictor"]]$mean[index_RC1],
      predictor = covariates %>% filter(tag == "pred_RC1") %>% pull(RC1),
      upper = current_model[["summary.linear.predictor"]]$`0.975quant`[index_RC1],
      lower = current_model[["summary.linear.predictor"]]$`0.025quant`[index_RC1],
      variable = "RC1",
      model = model_type
    ),
    tibble(
      response = current_model[["summary.linear.predictor"]]$mean[index_RC2],
      predictor = covariates %>% filter(tag == "pred_RC2") %>% pull(RC2),
      upper = current_model[["summary.linear.predictor"]]$`0.975quant`[index_RC2],
      lower = current_model[["summary.linear.predictor"]]$`0.025quant`[index_RC2],
      variable = "RC2",
      model = model_type
    ),
    tibble(
      response = current_model[["summary.linear.predictor"]]$mean[index_RC3],
      predictor = covariates %>% filter(tag == "pred_RC3") %>% pull(RC3),
      upper = current_model[["summary.linear.predictor"]]$`0.975quant`[index_RC3],
      lower = current_model[["summary.linear.predictor"]]$`0.025quant`[index_RC3],
      variable = "RC3",
      model = model_type
    )
  )
}

# Combine all marginal effects into one dataframe
marginal_effects <- bind_rows(marginal_effects_list)

# (9) Effect plots #############################################################

#### * Fixed effects: RCs * ####

# Plot the fixed effects
fixed_effects_plot <- fixed_effects_data  %>%
  filter(parameter != "intercept") %>%
  ggplot(aes(x = parameter, y = median, colour = colour)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_pointrange(
    aes(ymin = lower, ymax = upper),
    shape = 20,
    linewidth = 0.5,
    size = 0.5
  ) +
  scale_colour_identity() +
  facet_wrap(~ model, nrow = 1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  coord_flip() +
  common_theme +
  theme(plot.tag.location = "plot") +
  labs(
    x = NULL,
    y = 'Posterior standardised effect',
    tag = "(**a**)"
  )

# Display the plot
print(fixed_effects_plot)

#### * Marginal effects * ####

# Loop through RC1 to RC3 to create marginal effects plots
marginal_effects_plots <- list()
for(i in 1:3) {
  rc_param <- paste0("RC", i)
  
  # Filter data for the current RC
  df_current <- data_est %>%
    select(ends_with("_log_ratio"), all_of(rc_param)) %>%
    pivot_longer(
      -all_of(rc_param),
      names_to = "model",
      values_to = "response"
    ) %>%
    rename(RC = rc_param) %>%
    mutate(
      model = case_when(
        model == "richness_AM_log_ratio" ~ "AM",
        model == "richness_EcM_log_ratio" ~ "EcM",
        model == "richness_EcM_AM_log_ratio" ~ "EcM-AM",
        model == "richness_NM_log_ratio" ~ "NM",
        TRUE ~ model
      )
    )
  
  # For the marginal effects filtering, use the actual RC value (RC1, RC2, etc.)
  filtered_effects <- marginal_effects %>%
    filter(variable == rc_param) %>%
    rename(RC = predictor) %>%
    mutate(
      model = ifelse(model == "EcM_AM", "EcM-AM", model)
    )
  
  # Create panel label - only "b" for RC1, empty for others
  panel_label <- ifelse(i == 1, "(**b**)", "")
  
  p <- ggplot(data = df_current, aes(x = RC, y = response)) +
    geom_hex(bins = 40, aes(alpha = after_stat(ndensity)), fill = "black") +
    scale_alpha(range = c(0.1, 1)) +
    geom_hline(yintercept = 0, linetype = 'dotted', colour = '#de2d26', linewidth = 0.5) +
    geom_ribbon(
      data = filtered_effects,
      aes(ymin = lower, ymax = upper), fill = "#3366FF", alpha = 0.2,
    ) +
    geom_line(
      data = filtered_effects,
      aes(y = response), colour = "#3366FF", linewidth = 0.5) +
    facet_wrap(
      ~model, nrow = 1
    ) +
    scale_y_continuous(
      limits = c(-5, 4), 
      breaks = c(-4, -2,  0, 2, 4)
      # limits = c(-2.2, 2), 
      # breaks = c(-2, -1, 0, 1, 2)
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
      x = case_when(
        rc_param == "RC1" ~ "RC1 (Temperature & decomposition)",
        rc_param == "RC2" ~ "RC2 (Soil moisture)",
        rc_param == "RC3" ~ "RC3 (Soil phosphorus)",
        TRUE ~ rc_param
      ),
      y = if(i == 0) 'Richness (log-scaled)' else NULL,
      # Add tag for annotation
      tag = panel_label
    ) +
    theme(
      # Position the tag in the top left
      plot.tag = element_text(size = 14, face = "bold"),
      theme(plot.tag.location = "plot")
      #plot.tag.position = c(-0.02, 1)
    )
  
  # Conditionally adjust the strip text formatting:
  if(i == 0) {
    p <- p + common_theme + theme(
      plot.margin = margin(t = 5, r = 1, b = 5, l = 1, "pt"),
      strip.text = element_text(face = "bold", size = strip_size)
    )
  } else {
    p <- p + common_theme + theme(
      plot.margin = margin(t = 5, r = 1, b = 5, l = 1, "pt"),
      strip.text = element_blank()
    )
  }
  
  # Save to the correctly named list
  marginal_effects_plots[[i]] <- p
}

# (10) Richness plots ##########################################################

# Standardised richness to percentiles
standardised_rast_AM <- standardise_raster_percentiles(fitted_rast_AM[["predicted"]])
standardised_rast_EcM <- standardise_raster_percentiles(fitted_rast_EcM[["predicted"]])
standardised_rast_EcM_AM <- standardise_raster_percentiles(fitted_rast_EcM_AM[["predicted"]])
standardised_rast_NM <- standardise_raster_percentiles(fitted_rast_NM[["predicted"]])

richness_main <- aus_map_2 +
  tidyterra::geom_spatraster(
    data = terra::rast(list(
      AM = standardised_rast_AM,
      EcM = standardised_rast_EcM,
      `EcM-AM` = standardised_rast_EcM_AM,
      NM = standardised_rast_NM
    )),
    na.rm = TRUE
  ) +
  scale_fill_stepsn(
    colors = rev(paletteer::paletteer_c("grDevices::Spectral", 100)),
    breaks = 1:100,
    labels = c("Low", rep("", 98), "High"),
    na.value = "transparent"
  ) +
  facet_wrap(~ lyr, nrow = 1) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
    strip.text = element_text(face = "bold", size = strip_size),
    strip.background = element_blank(),
    panel.grid = element_line(colour = "white", linewidth = 0.125),
    axis.text = element_text(size = text_size),
    plot.tag = element_markdown(size = tag_size),
    plot.tag.location = "plot",
    plot.tag.position = c(0.025, 0.85)
  ) +
  labs(x = NULL, y = NULL, tag = "(**a**)")

# Display the plot
print(richness_main)

# (11) Env space plots #########################################################

# Load the tree data
data_env_space <- tibble(
  predicted_AM = fitted_values %>% filter(model == "AM") %>% pull(predicted),
  predicted_EcM = fitted_values %>% filter(model == "EcM") %>% pull(predicted),
  predicted_EcM_AM = fitted_values %>% filter(model == "EcM_AM") %>% pull(predicted),
  predicted_NM = fitted_values %>% filter(model == "NM") %>% pull(predicted),
  RC1 = c(data_est[["RC1"]], data_pred[["RC1"]]),
  RC2 = c(data_est[["RC2"]], data_pred[["RC2"]]),
  RC3 = c(data_est[["RC3"]], data_pred[["RC3"]])
)

#### * Create env bins * ###

# Find min and max for axis pairs
min_max_limits  <- data_env_space %>%
  summarise(
    RC1_min = min(RC1, na.rm = TRUE),
    RC1_max = max(RC1, na.rm = TRUE),
    RC2_min = min(RC2, na.rm = TRUE),
    RC2_max = max(RC2, na.rm = TRUE),
    RC3_min = min(RC3, na.rm = TRUE),
    RC3_max = max(RC3, na.rm = TRUE)
  )

# Create the break sequences outside the mutate function
RC1_breaks <- seq(min_max_limits $RC1_min, min_max_limits $RC1_max, length.out = 21)
RC2_breaks <- seq(min_max_limits $RC2_min, min_max_limits $RC2_max, length.out = 21)
RC3_breaks <- seq(min_max_limits $RC3_min, min_max_limits $RC3_max, length.out = 21)

# Create 2D bins for RC1 and RC2
RC1_RC2_bin_lookup <- expand.grid(
  RC1_bin = 1:20,
  RC2_bin = 1:20
)

# Assign unique bin numbers
RC1_RC2_bin_lookup[["bin_2d"]] <- 1:nrow(RC1_RC2_bin_lookup)

# Add bin numbers to the data
RC1_RC2_2d_bins <- data_env_space %>%
  mutate(
    # Get the bin numbers again (1-40)
    RC1_bin = cut(RC1, breaks = RC1_breaks, labels = FALSE, include.Low = TRUE),
    RC2_bin = cut(RC2, breaks = RC2_breaks, labels = FALSE, include.Low = TRUE)
  ) %>%
  left_join(RC1_RC2_bin_lookup, by = c("RC1_bin", "RC2_bin"))

#### * Summarise the data * ####

# RC1-RC2
data_RC1_RC2 <- RC1_RC2_2d_bins %>%
  group_by(RC1_bin, RC2_bin) %>%
  summarise(
    predicted_AM = mean(predicted_AM, na.rm = TRUE),
    predicted_EcM = mean(predicted_EcM, na.rm = TRUE),
    predicted_EcM_AM = mean(predicted_EcM_AM, na.rm = TRUE),
    predicted_NM = mean(predicted_NM, na.rm = TRUE)
  )

# Function to standardise richness values to percentiles for colour scale
data_RC1_RC2_quantised <- quantise_tibble_percentiles(
  data_RC1_RC2, 
  c("predicted_AM", "predicted_EcM", "predicted_EcM_AM", "predicted_NM")
)

#### * Create the plots * ####

# RC1-RC2 environmental space plot
rc12_plot <- data_RC1_RC2_quantised %>%
  mutate(
    RC1_original = min_max_limits$RC1_min + (RC1_bin - 0.5) * (min_max_limits$RC1_max - min_max_limits$RC1_min) / 20,
    RC2_original = min_max_limits$RC2_min + (RC2_bin - 0.5) * (min_max_limits$RC2_max - min_max_limits$RC2_min) / 20
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
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.tag = element_markdown(size = 16),
    aspect.ratio = 1,
    plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
  ) +
  labs(
    x = "RC1 (Temperature & decomposition)",
    y = "RC2 (Soil moisture)", 
    tag = "(**c**)"
  ) +
  facet_wrap(~mycorrhizal_type, nrow = 1)

# Display the plots
print(rc12_plot)

# (12) Supplementary plots #####################################################

#### * Predicted values: Standardised to deciles * ####

richness_sup_1 <- patchwork::wrap_plots(
  # AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = standardised_rast_AM,
      na.rm = TRUE
    ) +
    scale_fill_stepsn(
      colors = rev(paletteer::paletteer_c("grDevices::Spectral", 100)),
      breaks = 1:100,
      labels = c(rep("", 10), "Low", rep("", 78), "High", rep("", 10)),
      na.value = "transparent"
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = "Predicted richness\n(standardised to percentiles)", fill = NULL) +
    ggtitle("AM"),
  # EcM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = standardised_rast_EcM,
      na.rm = TRUE
    ) +
    scale_fill_stepsn(
      colors = rev(paletteer::paletteer_c("grDevices::Spectral", 100)),
      breaks = 1:100,
      labels = c(rep("", 10), "Low", rep("", 78), "High", rep("", 10)),
      na.value = "transparent"
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL) +
    ggtitle("EcM"),
  # EcM-AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = standardised_rast_EcM_AM,
      na.rm = TRUE
    ) +
    scale_fill_stepsn(
      colors = rev(paletteer::paletteer_c("grDevices::Spectral", 100)),
      breaks = 1:100,
      labels = c(rep("", 10), "Low", rep("", 78), "High", rep("", 10)),
      na.value = "transparent"
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL) +
    ggtitle("Dual"),
  # NM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = standardised_rast_NM,
      na.rm = TRUE
    ) +
    scale_fill_stepsn(
      colors = rev(paletteer::paletteer_c("grDevices::Spectral", 100)),
      breaks = 1:100,
      labels = c(rep("", 10), "Low", rep("", 78), "High", rep("", 10)),
      na.value = "transparent"
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL) +
    ggtitle("NM"),
  # Number of rows
  nrow = 1
)

# Display the map
print(richness_sup_1)

#### * Predicted values: Regular scale * ####

richness_sup_2 <- patchwork::wrap_plots(
  # AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_AM,
      aes(fill = predicted),
      na.rm = TRUE
    ) +
    
    scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = rel(0.9), face = "bold"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.title = element_markdown()
    ) +
    labs(x = NULL, y = "Predicted richness<br>(log-ratio)", fill = NULL),
  # EcM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM,
      aes(fill = predicted),
      na.rm = TRUE
    ) +
    scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = rel(0.9), face = "bold"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL),
  # EcM-AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM_AM,
      aes(fill = predicted),
      na.rm = TRUE
    ) +
    scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = rel(0.9), face = "bold"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL),
  # NM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_NM,
      aes(fill = predicted),
      na.rm = TRUE
    ) +
    scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = rel(0.9), face = "bold"),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL),
  # Number of rows
  nrow = 1
)

# Display the map
print(richness_sup_2)

#### * Observed richness * ####

# Plot the observed values
richness_observed <- patchwork::wrap_plots(
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_AM,
      aes(fill = observed),
      na.rm = TRUE
    ) +
    scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      axis.title.y = element_markdown()
    ) +
    labs(x = NULL, y = "Observed richness<br>(log-ratio)", fill = NULL),
  # EcM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM,
      aes(fill = observed),
      na.rm = TRUE
    ) +
    scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL),
  # EcM-AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM_AM,
      aes(fill = observed),
      na.rm = TRUE
    ) +
    scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL),
  # NM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_NM,
      aes(fill = observed),
      na.rm = TRUE
    ) +
    scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL),
  # Number of rows
  nrow = 1
)

# Display the map
print(richness_observed)

#### * Uncertainty * ####

# Plot the uncertainty
richness_uncertainty <- patchwork::wrap_plots(
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_AM,
      aes(fill = uncertainty),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      option = "B",
      na.value = "transparent"
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = "Prediction uncertainty\n(95% credible interval)", fill = NULL),
  # EcM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM,
      aes(fill = uncertainty),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      option = "B",
      na.value = "transparent",
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL),
  # EcM-AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM_AM,
      aes(fill = uncertainty),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      option = "B",
      na.value = "transparent",
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL),
  # NM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_NM,
      aes(fill = uncertainty),
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(
      option = "B",
      na.value = "transparent",
      breaks = scales::pretty_breaks()
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.text = element_text(hjust = 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(7, "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL),
  # Number of rows
  nrow = 1
)

# Display the map
print(richness_uncertainty)

#### * Observed vs predicted  * ####
obs_fit_plot_2 <- ggplot(data_diagnostics, aes(x = posterior_mean, y = observed)) +
  geom_point(alpha = 0.3, shape = 1, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, colour = "red") +
  ggpubr::stat_cor(aes(label = after_stat(rr.label)), colour = "red", size = 4) +
  facet_wrap(~ model, nrow = 1, scales = "free") +
  theme_minimal()  +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    aspect.ratio = 1,
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, "pt"),
    strip.text = element_blank()
  ) +
  labs(x = "Posterior mean predictions", y = "Observed values") +
  ggh4x::facetted_pos_scales(
    x = list(
      "Model AM" = scale_x_continuous(limits = limits_list[["Model AM"]]),
      "Model EcM" = scale_x_continuous(limits = limits_list[["Model EcM"]]),
      "Model EcM-AM" = scale_x_continuous(limits = limits_list[["Model EcM_AM"]]),
      "Model NM" = scale_x_continuous(limits = limits_list[["Model NM"]])
    ),
    y = list(
      "Model AM" = scale_y_continuous(limits = limits_list[["Model AM"]]),
      "Model EcM" = scale_y_continuous(limits = limits_list[["Model EcM"]]),
      "Model EcM-AM" = scale_y_continuous(limits = limits_list[["Model EcM_AM"]]),
      "Model NM" = scale_y_continuous(limits = limits_list[["Model NM"]])
    )
  )

# Display the plot
print(obs_fit_plot_2)

#### * Join and save the supplementary plots * ####

# Join the supplementary plots
supplementary_plots <- patchwork::wrap_plots(
  richness_sup_1,
  richness_sup_2,
  richness_observed,
  richness_uncertainty,
  obs_fit_plot_2,
  ncol = 1
) + 
  theme(
    plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
  )

# Save the supplementary plots
ggsave(
  filename = "output/supplimentary_richness_log_ratio/model_performance_plot.png",
  plot = supplementary_plots,
  width = 22,
  height = 26,
  units = "cm",
  dpi = 300
)

# (13) Spatial diagnostics #####################################################

#### * Matérn correlation function * ####

# Compute statistics
spatial_stats_AM <- spatial_params(model_AM, spde, name = "spatial")
spatial_stats_EcM <- spatial_params(model_EcM, spde, name = "spatial")
spatial_stats_EcM_AM <- spatial_params(model_EcM_AM, spde, name = "spatial")
spatial_stats_NM <- spatial_params(model_NM, spde, name = "spatial")

# Compare prior and posterior sigma:
# Re-fit models if posterior sigma is larger than prior sigma
sigma_est
spatial_stats_AM["sigma_spatial"]
spatial_stats_EcM["sigma_spatial"]
spatial_stats_EcM_AM["sigma_spatial"]
spatial_stats_NM["sigma_spatial"]

# Compare prior and posterior range:
# Re-fit models if posterior range is smaller than prior range
range_est
spatial_stats_AM["range"]
spatial_stats_EcM["range"]
spatial_stats_EcM_AM["range"]
spatial_stats_NM["range"]

# Matérn correlation function plot
matern_plot <- patchwork::wrap_plots(
  # AM plot
  plot_matern_correlation(
    spatial_stats_AM,
    x = data_est$x_km,
    y = data_est$y_km,
    range = spatial_stats_AM["range"],
    y_limits = c(0, 0.71),
    x_limits = c(0, 400), 
    x_breaks = seq(0, 400, by = 100)
  ) +
    theme(plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")),
  # EcM plot
  plot_matern_correlation(
    spatial_stats_EcM,
    x = data_est$x_km,
    y = data_est$y_km,
    range = spatial_stats_EcM["range"],
    y_limits = c(0, 0.71),
    x_limits = c(0, 400), 
    x_breaks = seq(0, 400, by = 100)
  ) + 
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ),
  # EcM-AM plot
  plot_matern_correlation(
    spatial_stats_EcM_AM,
    x = data_est$x_km,
    y = data_est$y_km,
    range = spatial_stats_EcM_AM["range"],
    y_limits = c(0, 0.71),
    x_limits = c(0, 400), 
    x_breaks = seq(0, 400, by = 100)
  ) + 
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ),
  # NM plot
  plot_matern_correlation(
    spatial_stats_NM,
    x = data_est$x_km,
    y = data_est$y_km,
    range = spatial_stats_NM["range"],
    y_limits = c(0, 0.71),
    x_limits = c(0, 400), 
    x_breaks = seq(0, 400, by = 100)
  ) + 
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ),
  # Number of rows
  nrow = 1
)

# Display the plot
print(matern_plot)

#### * Compute Morans I * ####

# Compute spatial weights
coords <- coordinates(data_est %>% select(x_albers, y_albers))
tri <- deldir::deldir(coords[,1], coords[,2])
nb <- spdep::dnearneigh(coords, d1 = 0, d2 = quantile(dist(coords), 0.1))
lw <- spdep::nb2listw(nb, style = "W")

# List to store results
moran_results <- list()
unique(data_diagnostics$model)
# Loop over each model and compute Moran’s I
for (model_name in c("AM", "EcM", "EcM-AM", "NM")) {
  # Extract residuals for the current model
  residuals <- data_diagnostics %>%
    filter(model == paste0("Model ", model_name)) %>%
    pull(residuals)
  
  # Compute Moran’s I
  moran_test <- spdep::moran.test(residuals, lw)
  
  # Store results in a list
  moran_results[[model_name]] <- list(
    Moran_I = moran_test$estimate[1],
    p_value = moran_test$p.value
  )
}

# Print results
print(moran_results)

#### * Spatial residuals * ####

# Spatial residuals plot
spatial_residuals <- patchwork::wrap_plots(
  # AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_AM,
      aes(fill = residuals),
      na.rm = TRUE
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = "Residuals", fill = NULL) +
    geom_richtext(
      aes(
        x = 114,
        y = -Inf,
        label = paste("*p*-value = ", round(moran_results[["AM"]][["p_value"]], 2))
      ),
      size = 2.5,
      hjust = 0.1,
      vjust = 0,
      label.color = NA,
      fill = NA
    ) +
    geom_richtext(
      aes(
        x = 114,
        y = -Inf,
        label = paste0("Moran's I = ", round(moran_results[["AM"]][["Moran_I"]], 3))
      ),
      size = 2.5,
      hjust = 0.1,
      vjust = -0.75,
      label.color = NA,
      fill = NA
    ),
  # EcM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM,
      aes(fill = residuals),
      na.rm = TRUE
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL) +
    geom_richtext(
      aes(
        x = 114,
        y = -Inf,
        label = paste("*p*-value = ", round(moran_results[["EcM"]][["p_value"]], 2))
      ),
      size = 2.5,
      hjust = 0.1,
      vjust = 0,
      label.color = NA,
      fill = NA
    ) +
    geom_richtext(
      aes(
        x = 114,
        y = -Inf,
        label = paste0("Moran's I = ", round(moran_results[["EcM"]][["Moran_I"]], 3))
      ),
      size = 2.5,
      hjust = 0.1,
      vjust = -0.75,
      label.color = NA,
      fill = NA
    ),
  # EcM-AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM_AM,
      aes(fill = residuals),
      na.rm = TRUE
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL) +
    geom_richtext(
      aes(
        x = 114,
        y = -Inf,
        label = paste("*p*-value = ", round(moran_results[["EcM-AM"]][["p_value"]], 2))
      ),
      size = 2.5,
      hjust = 0.1,
      vjust = 0,
      label.color = NA,
      fill = NA
    ) +
    geom_richtext(
      aes(
        x = 114,
        y = -Inf,
        label = paste0("Moran's I = ", round(moran_results[["EcM-AM"]][["Moran_I"]], 3))
      ),
      size = 2.5,
      hjust = 0.1,
      vjust = -0.75,
      label.color = NA,
      fill = NA
    ),
  # NM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_NM,
      aes(fill = residuals),
      na.rm = TRUE
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      limits = c(-1, 1),
      breaks = c(-1, 0, 1)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
    ) +
    labs(x = NULL, y = NULL, fill = NULL) +
    geom_richtext(
      aes(
        x = 114,
        y = -Inf,
        label = paste("*p*-value = ", round(moran_results[["NM"]][["p_value"]], 2))
      ),
      size = 2.5,
      hjust = 0.1,
      vjust = 0,
      label.color = NA,
      fill = NA
    ) +
    geom_richtext(
      aes(
        x = 114,
        y = -Inf,
        label = paste0("Moran's I = ", round(moran_results[["NM"]][["Moran_I"]], 3))
      ),
      size = 2.5,
      hjust = 0.1,
      vjust = -0.75,
      label.color = NA,
      fill = NA
    ),
  # Number of rows
  nrow = 1
)

# Display the plot
print(spatial_residuals)

#### * Posterior spatial effect * ####

plot_posterior_spatial <- patchwork::wrap_plots(
  # AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_AM,
      aes(fill = spatial_field),
      na.rm = TRUE
    )  +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold")
    ) +
    labs(x = NULL, y = "Posterior spatial effect", fill = NULL, title = "AM"),
  # EcM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM,
      aes(fill = spatial_field),
      na.rm = TRUE
    )  +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold")
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = "EcM"),
  # EcM-AM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_EcM_AM,
      aes(fill = spatial_field),
      na.rm = TRUE
    )  +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold")
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = "Dual"),
  # NM plot
  aus_map +
    tidyterra::geom_spatraster(
      data = fitted_rast_NM,
      aes(fill = spatial_field),
      na.rm = TRUE
    )  +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      breaks = scales::pretty_breaks(n = 3)
    ) +
    theme(
      axis.text = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0.115, 0.925),
      legend.background = element_blank(),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(hjust = 1, size = 7),
      plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = strip_size, face = "bold")
    ) +
    labs(x = NULL, y = NULL, fill = NULL, title = "NM"),
  # Number of rows
  nrow = 1
)

# Display the map
print(plot_posterior_spatial)

#### * Join and save the spatial diagnostics * ####

# Join the spatial diagnostics
spatial_diagnostics <- patchwork::wrap_plots(
  plot_posterior_spatial,
  spatial_residuals,
  matern_plot,
  ncol = 1
) + 
  theme(
    plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1), "pt")
  )

# Save the spatial diagnostics
ggsave(
  filename = "output/supplimentary_richness_log_ratio/spatial_diagnostics.png",
  plot = spatial_diagnostics,
  width = 16,
  height = 12.75,
  units = "cm",
  dpi = 300
)

# (12) Save the data ###########################################################

# Save raster data fro figure 2a
writeRaster(
  rast(c(
    AM = standardised_rast_AM,
    EcM = standardised_rast_EcM,
    `EcM-AM` = standardised_rast_EcM_AM,
    NM = standardised_rast_NM)),
  overwrite = TRUE,
  filename = "output/generated_data/figure_2a.tif"
)

# Save data for figure 2c
data_figure_2c <- data_RC1_RC2_quantised
data_figure_2c_limits <- min_max_limits
save(
  data_figure_2c,
  data_figure_2c_limits,
  file = "output/generated_data/figure_2c.RData"
)

# Organise and save the data for figure 3

# Observed data
relative_richness_data <- data_est %>%
  mutate(
    AM = log10((richness_AM + 1) / (richness - richness_AM + 1)),
    EcM = log10((richness_EcM + 1) / (richness - richness_EcM + 1)),
    `EcM-AM` = log10((richness_EcM_AM + 1) / (richness - richness_EcM_AM + 1)),
    NM = log10((richness_NM + 1) / (richness - richness_NM + 1))
  ) %>%
  select(
    RC1, RC2, RC3, AM, EcM, `EcM-AM`, NM
  ) %>%
  pivot_longer(
    cols = -c(RC1, RC2, RC3),
    names_to = "mycorrhizal_type",
    values_to = "response"
  )

# Marginal effects data
relative_richness_marginal_effects <- marginal_effects %>%
  mutate(
    response = log10(exp(response)),
    lower = log10(exp(lower)),
    upper = log10(exp(upper)),
    mycorrhizal_type = model,
    mycorrhizal_type = ifelse(
      mycorrhizal_type == "EcM_AM", "EcM-AM", mycorrhizal_type
    )
  ) %>%
  select(variable, response, predictor, lower, upper, mycorrhizal_type)

# Model coefficients
relative_richness_coeficients <- fixed_effects_data %>%
  filter(
    parameter != "intercept"
  ) %>%
  mutate(
    # Format the text with adjusted slope and CI
    annotation = sprintf("β = %.3f [%.3f, %.3f]", median, lower, upper),
    # Use Inf for positioning
    x_pos = Inf,
    y_pos = Inf,
    mycorrhizal_type = model
  ) %>%
  select(parameter, x_pos, y_pos, annotation, mycorrhizal_type)

# Save the data for plotting
save(
  relative_richness_data,
  relative_richness_marginal_effects,
  relative_richness_coeficients,
  file = "output/generated_data/figure_3.RData"
)

# Save the model coefficients
model_names <- c("model_AM", "model_EcM", "model_EcM_AM", "model_NM")
map_dfr(model_names, function(model_name) {
  model <- get(model_name)
  
  bind_rows(
    model[["summary.fixed"]] %>%
      rownames_to_column("parameter") %>%
      as_tibble() %>%
      select(
        parameter, median = `0.5quant`, lower = `0.025quant`, upper = `0.975quant`
      ),
    model[["summary.hyperpar"]] %>%
      rownames_to_column("parameter") %>%
      as_tibble() %>%
      select(
        parameter, median = `0.5quant`, lower = `0.025quant`, upper = `0.975quant`
      )
  ) %>%
    mutate(
      quantiles = sprintf("[%.3f, %.3f]", lower, upper),
      parameter = recode(parameter, 
                         "intercept" = "Intercept",
                         "Precision for the Gaussian observations" = "Observation Precision",
                         "Range for spatial" = "Spatial Range",
                         "Stdev for spatial" = "Spatial SD"
      ),
      mycorrhizal_type = str_remove(model_name, "model_"),
      mycorrhizal_type = ifelse(mycorrhizal_type == "EcM_AM", "Dual", mycorrhizal_type)
    ) %>%
    select(parameter, median, quantiles, mycorrhizal_type)
}) %>%
  data.table::fwrite(
    "output/supplimentary_richness_log_ratio/model_coefficients.txt",
    sep = "\t"
  )
