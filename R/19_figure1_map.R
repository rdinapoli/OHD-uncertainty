# 19_figure1_map.R - Rapa Nui map with study areas (Figure 1)
# Purpose: Creates map showing three study areas from Stevenson et al. (2015)
# Inputs: data/stevenson_2015.csv (for study area locations)
#         output/Figure1/Fig 1b.tif (island DEM raster, not included in repo)
#         output/Figure1/Fig1 1a globe.tif (globe inset raster, not included in repo)
# Outputs: output/figures/figure1_map.png
# Runtime: ~1 minute
#
# NOTE: This script requires two TIF raster files that are not distributed
# with the repository due to their size. The pre-generated figure1_map.png
# in output/figures/ is provided for review. To regenerate, obtain the TIF
# files and place them in output/Figure1/.

library(terra)
library(ggplot2)
library(patchwork)
library(here)

base_dir <- here::here()
output_dir <- file.path(base_dir, "output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load and prepare island DEM raster

cat("Loading island DEM raster...\n")
island <- rast(file.path(base_dir, "output/Figure1/Fig 1b.tif"))

# Downsample 4x for ggplot performance (~1800x1350 px)
island_small <- aggregate(island, fact = 4)

# Convert to normalized RGB array for annotation_raster()
island_arr <- as.array(island_small) / 255
island_arr[is.nan(island_arr)] <- 1.0

# The TIF was exported from ArcGIS with south-to-north row ordering.
# annotation_raster() places row 1 at ymax (north), so flip vertically.
island_arr <- island_arr[dim(island_arr)[1]:1, , ]

# Geographic extent calibrated from two control points:
#   Rano Kau crater center: pixel (row=1083, col=1387) -> (-27.187, -109.434)
#   Terevaka summit:        pixel (row=3500, col=2860) -> (-27.076, -109.371)
# Linear fit: lat = 4.59e-5 * row - 27.237; lon = 4.01e-5 * col - 109.490
island_ext <- list(
  xmin = -109.490, xmax = -109.201,
  ymin = -27.237,  ymax = -26.989
)

# Load globe inset raster

cat("Loading globe inset...\n")
globe <- rast(file.path(base_dir, "output/Figure1/Fig1 1a globe.tif"))
globe_small <- aggregate(globe, fact = 6)
globe_arr <- as.array(globe_small) / 255
globe_arr[is.nan(globe_arr)] <- 1.0
globe_arr <- globe_arr[dim(globe_arr)[1]:1, , ]

# Study areas

study_areas <- data.frame(
  label = c("SA1: Te Niu\n(n = 127, 805 mm/yr)",
            "SA2: Maunga O'Koro\n(n = 50, 1,690 mm/yr)",
            "SA3: Anakena West\n(n = 65, 1,460 mm/yr)"),
  short = c("SA1", "SA2", "SA3"),
  lon = c(-109.395, -109.355, -109.322),
  lat = c(-27.085, -27.140, -27.075),
  # Offset label positions: avoid globe inset (upper-left) and TIF legend (lower-right)
  label_lon = c(-109.38, -109.44, -109.27),
  label_lat = c(-27.048, -27.168, -27.048),
  stringsAsFactors = FALSE
)

# Build main island map

cat("Building island map...\n")

p_island <- ggplot() +
  # DEM base layer (includes original scale bar, north arrow, elevation legend)
  annotation_raster(island_arr, interpolate = TRUE,
                    xmin = island_ext$xmin, xmax = island_ext$xmax,
                    ymin = island_ext$ymin, ymax = island_ext$ymax) +
  # Study area markers (red circles with white stroke)
  geom_point(data = study_areas, aes(x = lon, y = lat),
             shape = 21, size = 5, fill = "firebrick",
             color = "white", stroke = 1.2) +
  geom_text(data = study_areas, aes(x = lon, y = lat, label = short),
            size = 2, color = "white", fontface = "bold") +
  # Leader lines from labels to markers
  geom_segment(data = study_areas,
               aes(x = label_lon, y = label_lat, xend = lon, yend = lat),
               color = "grey30", linewidth = 0.3, lineend = "round") +
  # Study area labels (white box, black outline, black text)
  geom_label(data = study_areas,
             aes(x = label_lon, y = label_lat, label = label),
             fill = "white", color = "black", linewidth = 0.3,
             size = 2.3, lineheight = 0.85, family = "sans",
             label.padding = unit(0.15, "lines")) +
  # Coordinate system -- ratio corrects for latitude distortion
  coord_fixed(ratio = 1 / cos(27.12 * pi / 180),
              xlim = c(-109.475, -109.225),
              ylim = c(-27.210, -27.040),
              expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

# Globe inset

p_globe <- ggplot() +
  annotation_raster(globe_arr, interpolate = TRUE,
                    xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  coord_fixed(expand = FALSE) +
  theme_void()

# Composite: island map with globe inset in upper-left

cat("Compositing figure...\n")
p_final <- p_island +
  inset_element(p_globe,
                left = -0.01, right = 0.25,
                bottom = 0.68, top = 1.01)

# Save

ggsave(file.path(output_dir, "figure1_map.pdf"), p_final,
       width = 7, height = 5.5, device = cairo_pdf)
ggsave(file.path(output_dir, "figure1_map.png"), p_final,
       width = 7, height = 5.5, dpi = 300)

cat("Figure 1 saved to:", output_dir, "\n")
