# 23_figure5_bottleneck.R - Dual bottleneck heatmap (Figure 5)
# Purpose: Shows dating uncertainty as function of measurement precision and parameter uncertainty
# Inputs: output/tables/measurement_error_sweep_results.csv
# Outputs: output/figures/figure5_bottleneck.png
# Runtime: ~1 minute

library(tidyverse)
library(here)

base_dir <- here::here()
output_dir <- file.path(base_dir, "output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Shared style

theme_fig <- theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

col_actual <- "firebrick"

# Load cross-sensitivity data (6 data points: 3 sigma_meas x 2 tiers)

cat("Loading cross-sensitivity data...\n")

cross_sens <- read_csv(file.path(base_dir, "output/tables/measurement_error_sweep_results.csv"),
                       show_col_types = FALSE)

# Data: sigma_meas = 0.03, 0.01, 0.001
# tier0_uncertainty = parameters known exactly
# tier3_uncertainty = realistic parameter uncertainty

cat("Data points:\n")
print(cross_sens %>% select(sigma_meas, tier0_uncertainty, tier3_uncertainty, gap))

# Analytical interpolation model
#
# Model: total = sqrt(meas_component^2 + param_component^2)
#
# meas_component: interpolated from tier0 (power law in sigma_meas)
# param_component: param_frac * parameter_floor
#
# The parameter floor (+/-184 yr from manuscript) is best estimated from the
# sigma_meas=0.001 point where measurement is negligible.

# Fit measurement component from tier0 (known parameters)
# Use log-log linear fit: log(unc) = a + b * log(sigma)
fit_meas <- lm(log(tier0_uncertainty) ~ log(sigma_meas), data = cross_sens)
cat("\nMeasurement component fit:\n")
cat("  log(unc) =", round(coef(fit_meas)[1], 3), "+",
    round(coef(fit_meas)[2], 3), "* log(sigma)\n")

# Parameter floor: use the value at sigma_meas=0.001 where meas noise is minimal
# From the cross-sensitivity CSV, gap at 0.001 = 162 yr
# But the actual tier3 uncertainty is 178 yr, and tier0 is 15 yr
# So parameter_floor = sqrt(178^2 - 15^2) = 177 yr
param_floor <- sqrt(cross_sens$tier3_uncertainty[3]^2 - cross_sens$tier0_uncertainty[3]^2)
cat("  Parameter floor:", round(param_floor, 1), "yr\n")

# Prediction function
predict_uncertainty <- function(sigma_meas_val, param_frac) {
  log_meas <- coef(fit_meas)[1] + coef(fit_meas)[2] * log(sigma_meas_val)
  meas_component <- exp(log_meas)
  param_component <- param_frac * param_floor
  sqrt(meas_component^2 + param_component^2)
}

# Verify against data
cat("\nVerification:\n")
for (i in seq_len(nrow(cross_sens))) {
  pred_t0 <- predict_uncertainty(cross_sens$sigma_meas[i], 0)
  pred_t3 <- predict_uncertainty(cross_sens$sigma_meas[i], 1)
  cat(sprintf("  sigma=%.3f: tier0 actual=%.0f pred=%.0f | tier3 actual=%.0f pred=%.0f\n",
              cross_sens$sigma_meas[i],
              cross_sens$tier0_uncertainty[i], pred_t0,
              cross_sens$tier3_uncertainty[i], pred_t3))
}

# Generate heatmap grid

sigma_grid <- 10^seq(log10(0.0005), log10(0.04), length.out = 200)
param_grid <- seq(0, 1, length.out = 200)

grid_df <- expand.grid(sigma_meas = sigma_grid, param_frac = param_grid) %>%
  as_tibble()

# Vectorized prediction
grid_df$uncertainty <- mapply(predict_uncertainty,
                              grid_df$sigma_meas, grid_df$param_frac)

grid_df <- grid_df %>%
  mutate(uncertainty_capped = pmin(uncertainty, 350))

# Key operating points

# Use actual CSV values for labels (analytical model is for heatmap surface only)
stevenson_unc <- cross_sens$tier0_uncertainty[cross_sens$sigma_meas == 0.03]
realistic_unc <- cross_sens$tier3_uncertainty[cross_sens$sigma_meas == 0.03]
perfect_unc   <- cross_sens$tier3_uncertainty[cross_sens$sigma_meas == 0.001]

cat(sprintf("\nOperating points:\n  Stevenson (known params): \u00b1%.0f yr\n  Realistic: \u00b1%.0f yr\n  Perfect measurement: \u00b1%.0f yr\n",
            stevenson_unc, realistic_unc, perfect_unc))

# Plot

cat("Creating heatmap...\n")

fig5 <- ggplot(grid_df, aes(x = sigma_meas, y = param_frac, fill = uncertainty_capped)) +
  geom_raster(interpolate = TRUE) +
  # Color scale
  scale_fill_gradientn(
    colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE",
               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
    values = scales::rescale(c(0, 30, 60, 100, 150, 200, 250, 300, 350)),
    limits = c(0, 350),
    breaks = c(30, 100, 200, 300),
    labels = c("\u00b130", "\u00b1100", "\u00b1200", "\u00b1300"),
    name = "Dating\nUncertainty\n(\u00b1 years)"
  ) +
  # Contour lines at key thresholds
  geom_contour(aes(z = uncertainty), breaks = 30,
               color = "white", linewidth = 0.8, linetype = "dashed") +
  geom_contour(aes(z = uncertainty), breaks = round(param_floor),
               color = "grey80", linewidth = 0.6, linetype = "dashed") +
  geom_contour(aes(z = uncertainty), breaks = 299,
               color = "grey90", linewidth = 0.6, linetype = "dotted") +
  # Contour labels
  annotate("text", x = 0.0008, y = 0.03, label = "\u00b130 yr",
           color = "white", size = 2.8, fontface = "bold", hjust = 0) +
  # Parameters-known operating point (black circle, lower right)
  annotate("point", x = 0.03, y = 0,
           shape = 21, size = 4, fill = "black", color = "white", stroke = 1) +
  annotate("text", x = 0.022, y = 0.06,
           label = paste0("Parameters\nknown\n(\u00b1", round(stevenson_unc), " yr)"),
           color = "black", size = 2.8, hjust = 1, fontface = "bold") +
  # Realistic operating point (black square, upper right)
  annotate("point", x = 0.03, y = 1,
           shape = 22, size = 4, fill = "black", color = "white", stroke = 1) +
  annotate("text", x = 0.022, y = 0.94,
           label = paste0("Realistic\n(\u00b1", round(realistic_unc), " yr)"),
           color = "white", size = 2.8, hjust = 1, fontface = "bold") +
  # Perfect measurement point (black triangle, upper left)
  annotate("point", x = 0.001, y = 1,
           shape = 24, size = 4, fill = "black", color = "white", stroke = 1) +
  annotate("text", x = 0.0015, y = 0.94,
           label = paste0("Perfect\nmeasurement\n(\u00b1", round(perfect_unc), " yr)"),
           color = "black", size = 2.8, hjust = 0, fontface = "bold") +
  # Axis formatting
  scale_x_log10(
    breaks = c(0.001, 0.003, 0.01, 0.03),
    labels = c("0.001", "0.003", "0.01", "0.03")
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = c("0\n(Known)", "0.25", "0.5", "0.75", "1.0\n(Realistic)")
  ) +
  labs(
    x = expression("IR-PAS Measurement Precision, " * sigma[meas] * " (wt%)"),
    y = "Parameter Uncertainty (fraction of realistic)"
  ) +
  theme_fig +
  theme(
    legend.position = "right",
    legend.key.height = unit(1.5, "cm")
  )

# Save

ggsave(file.path(output_dir, "figure5_bottleneck.pdf"), fig5,
       width = 7, height = 5, device = cairo_pdf)
ggsave(file.path(output_dir, "figure5_bottleneck.png"), fig5,
       width = 7, height = 5, dpi = 300)

cat("Figure 5 saved to:", output_dir, "\n")
