# Required packages
require(INLA) # inla.list.models()
require(fmesher)
require(terra)
require(ggtext)
require(tidyverse)
source("code/functions.R")
source("code/map_australia.R")
set.seed(1986)

# (1) Read data ################################################################

# Read the estimation data
data_est <- data.table::fread(
  "data/presence/sites_relative_richness_est_10.txt",
  stringsAsFactors = TRUE
) %>%
  # Compute the log ratio of richness
  mutate(
    x_km = x_albers / 1000, y_km = y_albers / 1000,
    latitude = abs(latitude),
    latitude_std = std(latitude),
    richness_AM_log_ratio = log((richness_AM + 1) / (richness - richness_AM + 1)),
    richness_EcM_log_ratio = log((richness_EcM + 1) / (richness - richness_EcM + 1)),
    richness_EcM_AM_log_ratio = log((richness_EcM_AM + 1) / (richness - richness_EcM_AM + 1)),
    richness_NM_log_ratio = log((richness_NM + 1) / (richness - richness_NM + 1))
  )

# Count the number of cells for estimation
n_est <- nrow(data_est) %>%
  print()

# (2) Model formulas ###########################################################

# Plot richness log ratio against latitude
data_est %>%
  select(
    latitude, richness_AM_log_ratio, richness_EcM_log_ratio,
    richness_EcM_AM_log_ratio, richness_NM_log_ratio
    ) %>%
  pivot_longer(cols = -latitude) %>%
  ggplot(aes(latitude, value)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ name, scales = "free")

# Define the formulas
formula_AM <- as.formula(paste0(
  "y_AM ~ -1 + intercept + ",
  "latitude_std +",
  "f(spatial, model = spde)"
))
formula_EcM <- as.formula(paste0(
  "y_EcM ~ -1 + intercept + ",
  "latitude_std +",
  "f(spatial, model = spde)"
))

formula_EcM_AM <- as.formula(paste0(
  "y_EcM_AM ~ -1 + intercept + ",
  "latitude_std +",
  "f(spatial, model = spde)"
))
formula_NM <- as.formula(paste0(
  "y_NM ~ -1 + intercept + ",
  "latitude_std +",
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

# Distances between estimation site.
hist(dist(loc_est), main = "", xlab = "Distance between sites (km)")
# Small scale looks to 2000-3000 km at the continental scale and 0 - 1000 km 
# at the regional scale. I'll use the 500km peak.
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

# Covariates for the estimation (spatial and marginal effects)
covariates <- bind_rows(
  tibble(
    latitude_std = data_est$latitude_std,
    ecoregion = data_est$ecoregion,
    tag = rep("est", n_est)
  ),
  tibble(
    latitude_std = seq(
      min(data_est$latitude_std), max(data_est$latitude_std), length.out = 100
    ),
    tag = "est_latitude"
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

# Stack the data for the marginal effects
stack_est_latitude <- inla.stack(
  tag = "est_latitude",
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
    X = covariates %>% filter(tag == "est_latitude") %>% select(-tag)
  )
)

# Join the stacks
stack <- inla.stack.join(
  stack_est,
  stack_est_latitude
)

# Retrieve the stack indexes
index_est <- inla.stack.index(stack, "est")$data
index_est_latitude <- inla.stack.index(stack, "est_latitude")$data

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

#### * Marginal effects * ####

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

# Create an empty list to store marginal effects for each model
marginal_effects_list <- list()

# Loop through each model
for (model_type in names(model_info)) {
  # Get the current model object
  current_model <- get(model_info[[model_type]]$name)
  
  # Compute marginal effects for current model
  marginal_effects_list[[model_type]] <- tibble(
      response = current_model[["summary.linear.predictor"]]$mean[index_est_latitude],
      predictor = covariates %>% filter(tag == "est_latitude") %>% pull(latitude_std),
      upper = current_model[["summary.linear.predictor"]]$`0.975quant`[index_est_latitude],
      lower = current_model[["summary.linear.predictor"]]$`0.025quant`[index_est_latitude],
      variable = "latitude",
      model = model_type
      ) %>%
    mutate(
      model = ifelse(model == "EcM_AM", "EcM-AM", model),
      # Back-standardise latitude
      predictor = case_when(
        variable == "latitude" ~ unstd(
          predictor, mean(data_est$latitude), sd(data_est$latitude)
          ))
    )
}

# Combine all marginal effects into one dataframe
marginal_effects <- bind_rows(marginal_effects_list)

#### * Model coefficients * ####

# Extract model coefficients
coeficient_data <- inla_fixed_effects(
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
    parameter = factor(parameter, levels = rev(levels(factor(parameter))))
  )

# (9) Organise the plot data ###################################################

# I will plot the data on the log 10 scale for intuitive interpretation

# Observed data
latitude_gradient_data <- data_est %>%
  mutate(
    latitude = abs(latitude),
    AM = log10((richness_AM + 1) / (richness - richness_AM + 1)),
    EcM = log10((richness_EcM + 1) / (richness - richness_EcM + 1)),
    `EcM-AM` = log10((richness_EcM_AM + 1) / (richness - richness_EcM_AM + 1)),
    NM = log10((richness_NM + 1) / (richness - richness_NM + 1))
  ) %>%
  select(
    latitude, AM, EcM, `EcM-AM`, NM
  ) %>%
  pivot_longer(
    cols = -latitude,
    names_to = "mycorrhizal_type",
    values_to = "relative_richness"
  )

# Marginal effects data
latitude_gradient_marginal_effects <- marginal_effects %>%
  mutate(
    latitude = abs(predictor),
    relative_richness = log10(exp(response)),
    lower = log10(exp(lower)),
    upper = log10(exp(upper)),
    mycorrhizal_type = model,
  ) %>%
  select(latitude, relative_richness, lower, upper, mycorrhizal_type)

# Model coefficients
latitude_gradient_coeficients <- coeficient_data %>%
  filter(
    parameter == "latitude_std"
  ) %>%
  mutate(
    # Format the text with adjusted slope and CI
    annotation = sprintf("β = %.3f [%.3f, %.3f]", median, lower, upper),
    # Use Inf for positioning
    x_pos = Inf,
    y_pos = Inf,
    mycorrhizal_type = model
  ) %>%
  select(x_pos, y_pos, annotation, mycorrhizal_type)

# (10) Plot the data ############################################################

ggplot(
  data = latitude_gradient_data,
  aes(x = latitude, y = relative_richness)) +
  geom_hex(bins = 50, aes(alpha = after_stat(ndensity)), fill = "grey10") +
  scale_alpha(range = c(0.1, 1)) +
  geom_hline(
    yintercept = 0, linetype = 'dotted', colour = '#de2d26', linewidth = 0.5) +
  geom_ribbon(
    data = latitude_gradient_marginal_effects,
    aes(ymin = lower, ymax = upper), fill = "#3366FF", alpha = 0.2,
  ) +
  geom_line(
    data = latitude_gradient_marginal_effects,
    aes(y = relative_richness), colour = "#3366FF", linewidth = 0.5
    ) +
  # Add the annotations
  geom_text(
    data = latitude_gradient_coeficients,
    aes(x = x_pos, y = y_pos, label = annotation),
    size = 2.5, 
    hjust = 1.05,
    vjust = 2.5
  ) +
  facet_wrap(
    ~mycorrhizal_type, nrow = 1
  ) +
  theme_minimal() +
  scale_x_continuous(
    labels = function(x) paste0(abs(x), "°S")
  ) +
  scale_y_continuous(
    limits = c(-2.25, 2),
    breaks = c(-2, -1, 0, 1, 2),
    # Labels on the original scale
    labels = c("1:100", "1:10", "1:1", "10:1", "100:1")
  ) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    plot.margin = margin(t = 5, r = 0, b = 5, l = 0, "pt"),
    axis.title.x = element_text(size = rel(0.9)),
    aspect.ratio = 1
  )

# (11) Save the data ###########################################################

save(
  latitude_gradient_data,
  latitude_gradient_marginal_effects,
  latitude_gradient_coeficients,
  file = "generated_data/figure_2b.RData"
)
