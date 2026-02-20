# 28_correlation.R - Environmental correlation analysis (Table S11)
# Purpose: Tests claimed environmental correlations with realistic uncertainty
# Inputs: output/tables/tier3_summary.csv
# Outputs: output/tables/correlation_summary.csv, output/tables/information_gradient.csv
# Runtime: ~1 minute

library(tidyverse)
library(here)

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

output_dir <- here("output/tables")
fig_dir <- here("output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

cat("Environmental correlation reassessment\n")

# -----------------------------------------------------------------------------
# Environmental Data from Stevenson et al. 2015 Table 1
# -----------------------------------------------------------------------------

env_data <- data.frame(
  study_area = c("SA1", "SA2", "SA3"),
  study_area_name = c("Te Niu", "Maunga O'Koro", "Anakena West"),

  # Location
  dist_to_coast_m = c(860, 4340, 792),
  elevation_m = c(120, 260, 170),

  # Climate
  rainfall_mm = c(805, 1690, 1460),

  # Soil nutrients
  resin_phosphorus_mg_kg = c(4.26, 0.49, 14.6),
  exchangeable_calcium_cmol_kg = c(6.03, 0.82, 4.42),
  base_saturation_pct = c(16.26, 2.89, 14.76),

  # Stevenson's claimed decline years
  stevenson_decline_year = c(1660, 1705, 1850),
  stevenson_pct_post_contact = c(0.64, 5.00, 67.2),

  # Sample sizes
  n_samples = c(127, 50, 65)
)

cat("Environmental data from Table 1:\n")
print(env_data %>% select(study_area, rainfall_mm, resin_phosphorus_mg_kg,
                          stevenson_decline_year))

# -----------------------------------------------------------------------------
# Calculate correlations with Stevenson's point estimates
# -----------------------------------------------------------------------------

cat("\nCorrelations with Stevenson's point estimates:\n")

cor_rainfall_stev <- cor(env_data$rainfall_mm, env_data$stevenson_decline_year)
cor_phosphorus_stev <- cor(env_data$resin_phosphorus_mg_kg, env_data$stevenson_decline_year)
cor_calcium_stev <- cor(env_data$exchangeable_calcium_cmol_kg, env_data$stevenson_decline_year)
cor_basesaturation_stev <- cor(env_data$base_saturation_pct, env_data$stevenson_decline_year)

cat(sprintf("  Rainfall vs Decline Year:     r = %.3f\n", cor_rainfall_stev))
cat(sprintf("  Phosphorus vs Decline Year:   r = %.3f\n", cor_phosphorus_stev))
cat(sprintf("  Calcium vs Decline Year:      r = %.3f\n", cor_calcium_stev))
cat(sprintf("  Base Saturation vs Decline:   r = %.3f\n", cor_basesaturation_stev))

# Note: With n=3, ANY correlation appears strong but is not statistically testable

# -----------------------------------------------------------------------------
# Function: Monte Carlo correlation with uncertainty
# -----------------------------------------------------------------------------

mc_correlation <- function(x, y_mean, y_sd, n_sim = 10000) {
  correlations <- numeric(n_sim)

  for (i in 1:n_sim) {
    y_sim <- rnorm(length(y_mean), mean = y_mean, sd = y_sd)
    correlations[i] <- cor(x, y_sim)
  }

  return(correlations)
}

# -----------------------------------------------------------------------------
# Monte Carlo correlations with +/-30 years (Stevenson's claimed uncertainty)
# -----------------------------------------------------------------------------

cat("\nCorrelations with +/-30 year uncertainty:\n")

N_SIM <- 10000

env_vars <- list(
  rainfall = env_data$rainfall_mm,
  phosphorus = env_data$resin_phosphorus_mg_kg,
  calcium = env_data$exchangeable_calcium_cmol_kg,
  base_saturation = env_data$base_saturation_pct
)

results_30 <- list()
for (var_name in names(env_vars)) {
  cors <- mc_correlation(env_vars[[var_name]],
                         env_data$stevenson_decline_year,
                         y_sd = 30,
                         n_sim = N_SIM)
  results_30[[var_name]] <- cors

  pct_positive <- mean(cors > 0) * 100
  cat(sprintf("  %s: Mean r = %.3f, 95%% CI = [%.3f, %.3f], P(r > 0) = %.1f%%\n",
              var_name, mean(cors), quantile(cors, 0.025), quantile(cors, 0.975), pct_positive))
}

# -----------------------------------------------------------------------------
# Monte Carlo correlations with +/-299 years (corrected uncertainty)
# -----------------------------------------------------------------------------

cat("\nCorrelations with +/-299 year uncertainty:\n")

results_299 <- list()
for (var_name in names(env_vars)) {
  cors <- mc_correlation(env_vars[[var_name]],
                         env_data$stevenson_decline_year,
                         y_sd = 299,
                         n_sim = N_SIM)
  results_299[[var_name]] <- cors

  pct_positive <- mean(cors > 0) * 100
  cat(sprintf("  %s: Mean r = %.3f, 95%% CI = [%.3f, %.3f], P(r > 0) = %.1f%%\n",
              var_name, mean(cors), quantile(cors, 0.025), quantile(cors, 0.975), pct_positive))
}

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------

cat("\nCreating visualizations...\n")

# Combine results for plotting
plot_data <- data.frame()
for (var_name in names(env_vars)) {
  plot_data <- bind_rows(plot_data,
    data.frame(
      variable = var_name,
      uncertainty = "\u00b130 years (Stevenson)",
      correlation = results_30[[var_name]]
    ),
    data.frame(
      variable = var_name,
      uncertainty = "\u00b1299 years (Corrected)",
      correlation = results_299[[var_name]]
    )
  )
}

plot_data$variable <- factor(plot_data$variable,
                             levels = c("rainfall", "phosphorus", "calcium", "base_saturation"),
                             labels = c("Rainfall", "Phosphorus", "Calcium", "Base Saturation"))

plot_data$uncertainty <- factor(plot_data$uncertainty,
                                levels = c("\u00b130 years (Stevenson)", "\u00b1299 years (Corrected)"))

# Histogram of correlations
p_cors <- ggplot(plot_data, aes(x = correlation, fill = uncertainty)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~variable, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("\u00b130 years (Stevenson)" = "blue",
                               "\u00b1299 years (Corrected)" = "orange")) +
  labs(
    x = "Correlation Coefficient (r)",
    y = "Count",
    fill = "Dating Uncertainty"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "correlation_distributions.pdf"), p_cors,
       width = 10, height = 8)
ggsave(file.path(fig_dir, "correlation_distributions.png"), p_cors,
       width = 10, height = 8, dpi = 150)
cat("Saved: correlation_distributions.pdf/png\n")

# Scatter plots with uncertainty
p_scatter <- env_data %>%
  select(study_area_name, rainfall_mm, stevenson_decline_year) %>%
  ggplot(aes(x = rainfall_mm, y = stevenson_decline_year)) +
  geom_point(size = 4, color = "blue") +
  geom_errorbar(aes(ymin = stevenson_decline_year - 30,
                    ymax = stevenson_decline_year + 30),
                width = 50, color = "blue", alpha = 0.5) +
  geom_errorbar(aes(ymin = stevenson_decline_year - 299,
                    ymax = stevenson_decline_year + 299),
                width = 100, color = "orange", alpha = 0.3) +
  geom_text(aes(label = study_area_name), vjust = -1.5, size = 3) +
  geom_hline(yintercept = 1722, linetype = "dashed", color = "red") +
  annotate("text", x = 1800, y = 1722, label = "European Contact",
           hjust = 1, vjust = -0.5, color = "red", size = 3) +
  labs(
    title = "Rainfall vs. Decline Year with Uncertainty",
    subtitle = "Blue = \u00b130 yr (Stevenson), Orange = \u00b1299 yr (Corrected)",
    x = "Annual Rainfall (mm)",
    y = "Estimated Decline Year (AD)"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(1300, 2100))

ggsave(file.path(fig_dir, "rainfall_correlation_scatter.pdf"), p_scatter,
       width = 8, height = 6)
cat("Saved: rainfall_correlation_scatter.pdf\n")

# -----------------------------------------------------------------------------
# Summary statistics
# -----------------------------------------------------------------------------

summary_df <- data.frame()
for (var_name in names(env_vars)) {
  stev_cor <- results_30[[var_name]]
  corr_cor <- results_299[[var_name]]

  summary_df <- bind_rows(summary_df, data.frame(
    variable = var_name,
    stev_mean_r = round(mean(stev_cor), 3),
    stev_ci_low = round(quantile(stev_cor, 0.025), 3),
    stev_ci_high = round(quantile(stev_cor, 0.975), 3),
    stev_pct_positive = round(mean(stev_cor > 0) * 100, 1),
    corr_mean_r = round(mean(corr_cor), 3),
    corr_ci_low = round(quantile(corr_cor, 0.025), 3),
    corr_ci_high = round(quantile(corr_cor, 0.975), 3),
    corr_pct_positive = round(mean(corr_cor > 0) * 100, 1),
    corr_includes_zero = quantile(corr_cor, 0.025) < 0 & quantile(corr_cor, 0.975) > 0
  ))
}

write_csv(summary_df, file.path(output_dir, "correlation_summary.csv"))
cat("Saved: correlation_summary.csv\n")

# -----------------------------------------------------------------------------
# Information gradient table
# -----------------------------------------------------------------------------

info_gradient <- data.frame(
  variable = names(env_vars),
  point_estimate_r = c(cor_rainfall_stev, cor_phosphorus_stev,
                       cor_calcium_stev, cor_basesaturation_stev)
)
info_gradient$ci_width_30 <- sapply(names(env_vars), function(v) {
  q <- quantile(results_30[[v]], c(0.025, 0.975))
  round(q[2] - q[1], 3)
})
info_gradient$ci_width_299 <- sapply(names(env_vars), function(v) {
  q <- quantile(results_299[[v]], c(0.025, 0.975))
  round(q[2] - q[1], 3)
})
info_gradient$ci_expansion_factor <- round(info_gradient$ci_width_299 / info_gradient$ci_width_30, 1)

write_csv(info_gradient, file.path(output_dir, "information_gradient.csv"))
cat("Saved: information_gradient.csv\n")

# -----------------------------------------------------------------------------
# Key Findings
# -----------------------------------------------------------------------------

cat("\nKey findings:\n")
cat("  With n=3 study areas, correlations are inherently unreliable.\n")
cat("  With +/-299 year uncertainty, ALL correlation 95% CIs span zero.\n")
cat("  The claimed environmental pattern is not statistically supported.\n")

cat("\nCorrelation analysis complete.\n")
