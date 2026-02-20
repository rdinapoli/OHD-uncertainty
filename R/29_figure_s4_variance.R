# 29_figure_s4_variance.R - Incremental variance decomposition (Figure S4)
# Purpose: Visualizes spatial vs. temporal contributions to prediction variance
# Inputs: output/tables/variance_decomposition.csv
# Outputs: output/figures/figure_s4_variance.png
# Runtime: ~1 minute

library(tidyverse)
library(here)

# Paths
input_file <- here("output/tables/variance_decomposition.csv")
output_dir <- here("output/figures")

# =============================================================================
# Load Data
# =============================================================================

variance_data <- read_csv(input_file, show_col_types = FALSE)

# Handle any negative temporal variances (numerical artifacts) by setting floor at 0
variance_data <- variance_data %>%
  mutate(var_temporal = pmax(var_temporal, 0))

# =============================================================================
# Calculate Summary Statistics
# =============================================================================

# Calculate median percentages for annotation
total_var <- variance_data$var_spatial + variance_data$var_temporal
pct_spatial <- median(variance_data$var_spatial / total_var) * 100
pct_temporal <- median(variance_data$var_temporal / total_var) * 100

cat(sprintf("Variance decomposition: Spatial %.2f%%, Temporal %.4f%%\n",
            pct_spatial, pct_temporal))

# =============================================================================
# Create Figure
# =============================================================================

p_variance <- ggplot(variance_data, aes(x = var_spatial, y = var_temporal)) +
  # Density shading (contours filled)
  stat_density_2d(
    geom = "polygon",
    aes(fill = after_stat(level)),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_gradient(low = "#cce5ff", high = "#004085") +
  # Points on top
  geom_point(alpha = 0.4, size = 1.2, color = "gray30") +
  # Reference line (1:1 - if temporal = spatial)
  geom_abline(slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  # Labels
  labs(
    x = "Spatial Variance",
    y = "Temporal Variance"
  ) +
  # Theme with sans-serif font
  theme_bw(base_size = 10) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# Save Figure
# =============================================================================

# PNG output
ggsave(
  file.path(output_dir, "figure_s4_variance.png"),
  p_variance,
  width = 5,
  height = 4,
  dpi = 300
)

cat("Figure S4 saved to:", file.path(output_dir, "figure_s4_variance.png"), "\n")
