# Aus map
aus_map <- rnaturalearth::ne_countries(
  scale = 'medium', type = 'map_units', returnclass = 'sf'
  ) %>%
  filter(name == 'Australia') %>%
  ggplot() +
  geom_sf(color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(
    limits = c(113, 154),
    breaks = c(120, 130, 140, 150),
    labels = c("120°E", "", "140°E", "")
  ) +
  scale_y_continuous(
    limits = c(-43.5, -10),
    breaks = seq(-40, -10, 10)
  ) +
  theme_grey() +
  theme(
    plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  )

aus_map_2 <- rnaturalearth::ne_countries(
  scale = 'medium', type = 'map_units', returnclass = 'sf'
) %>%
  filter(name == 'Australia') %>%
  ggplot() +
  geom_sf(color = "black", fill = NA) +
  xlab(NULL) +
  ylab(NULL) +
  scale_x_continuous(
    limits = c(113, 154),
    breaks = c(120, 130, 140, 150),
    labels = c("120°E", "", "140°E", "")
    ) +
  scale_y_continuous(
    limits = c(-43.5, -10),
    breaks = seq(-40, -10, 10)
  ) +
  theme_grey() +
  theme(
    plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
    panel.border = element_rect(colour = "grey90", fill = NA, linewidth = 0.5),
    panel.grid = element_line(colour = "white", linewidth = 0.25),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  )


aus_map_3 <- rnaturalearth::ne_countries(
  scale = 'medium', type = 'map_units', returnclass = 'sf'
) %>%
  filter(name == 'Australia')
 

# Map asp ratio
map_aspect_ratio <-  ((-10) - (-43.5)) / (154 - 113)

aus_map_albers <- sf::st_transform(
  x = aus_map$data,
  crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +datum=WGS84 +units=km"
) %>%
  ggplot() +
  geom_sf() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  ) +
  scale_x_continuous(limits = c(-2000, 2250)) +
  scale_y_continuous(
    limits = c(-1100, -4800)
  ) 

