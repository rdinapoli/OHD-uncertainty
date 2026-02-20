# 30_figure_s7_sensitivity.R - Comprehensive sensitivity visualization (Figure S7)
# Purpose: Four-panel sensitivity sweep plot with dual bottleneck summary
# Inputs: output/tables/comprehensive_sensitivity_results.csv,
#   output/tables/measurement_error_sweep_results.csv
# Outputs: output/figures/figure_s7_sensitivity.png
# Runtime: ~1 minute

library(tidyverse)
library(gridExtra)
library(here)

# Output path
output_dir <- here("output/figures")
tables_dir <- here("output/tables")

# =============================================================================
# Load Data
# =============================================================================

# Panel A: Comprehensive sensitivity results
sensitivity_data <- read_csv(file.path(tables_dir, "comprehensive_sensitivity_results.csv"),
                             show_col_types = FALSE)

# Panel B: Cross-sensitivity results (dual bottleneck)
cross_sensitivity <- read_csv(file.path(tables_dir, "cross_sensitivity_results.csv"),
                              show_col_types = FALSE)

# Also need Tier 0 data for comparison
tier0_data <- read_csv(file.path(tables_dir, "measurement_error_sweep_results.csv"),
                       show_col_types = FALSE)

# =============================================================================
# Panel A: Prior Sensitivity (All 4 Parameters)
# =============================================================================

# Prepare data - rename parameters for display
sensitivity_plot_data <- sensitivity_data %>%
  mutate(
    Parameter = case_when(
      parameter == "EHT_sd" ~ "EHT (Temperature)",
      parameter == "H2Ot_sigma" ~ "H\u2082Ot (Structural Water)",
      parameter == "Ea_sigma" ~ "Ea (Activation Energy)",
      parameter == "RH_sd" ~ "RH (Relative Humidity)"
    ),
    Parameter = factor(Parameter, levels = c(
      "EHT (Temperature)",
      "H\u2082Ot (Structural Water)",
      "Ea (Activation Energy)",
      "RH (Relative Humidity)"
    ))
  )

# Create Panel A - Raw years with reference line
panel_a <- ggplot(sensitivity_plot_data,
                  aes(x = fold_reduction, y = median_uncertainty, color = Parameter, shape = Parameter)) +
  # Green "acceptable" zone showing claimed precision
  annotate("rect", xmin = 0.8, xmax = 35, ymin = 0, ymax = 30,
           fill = "#90EE90", alpha = 0.2) +
  # Reference line at claimed +/-30 years
  geom_hline(yintercept = 30, linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
  annotate("text", x = 25, y = 42, label = "Claimed \u00b130 yr",
           size = 3, color = "darkgreen", fontface = "italic") +
  # Data lines and points
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(1, 3, 10, 30),
                labels = c("1\u00d7", "3\u00d7", "10\u00d7", "30\u00d7")) +
  scale_y_continuous(limits = c(0, 320),
                     breaks = seq(0, 300, 50)) +
  scale_color_manual(values = c(
    "EHT (Temperature)" = "#E69F00",
    "H\u2082Ot (Structural Water)" = "#56B4E9",
    "Ea (Activation Energy)" = "#009E73",
    "RH (Relative Humidity)" = "#CC79A7"
  )) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    x = "Prior Tightening Factor",
    y = "Age Uncertainty (years)",
    title = "A",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  ) +
  guides(color = guide_legend(nrow = 2))

# =============================================================================
# Panel B: Dual Bottleneck (Measurement vs Parameter Uncertainty)
# =============================================================================

# Prepare dual bottleneck data - join on sigma_meas to handle different row counts
tier0_for_join <- tier0_data %>%
  select(sigma_meas, tier0_uncertainty = median_uncertainty)

tier3_for_join <- cross_sensitivity %>%
  select(sigma_meas, tier3_uncertainty)

bottleneck_wide <- tier0_for_join %>%
  inner_join(tier3_for_join, by = "sigma_meas")

bottleneck_data <- bottleneck_wide %>%
  rename(
    `Tier 0\n(Parameters Known)` = tier0_uncertainty,
    `Tier 3\n(Realistic Uncertainty)` = tier3_uncertainty
  ) %>%
  pivot_longer(
    cols = -sigma_meas,
    names_to = "Scenario",
    values_to = "uncertainty"
  ) %>%
  mutate(
    Scenario = factor(Scenario, levels = c(
      "Tier 0\n(Parameters Known)",
      "Tier 3\n(Realistic Uncertainty)"
    ))
  )

# Calculate the gap for annotation using the joined data
tier0_at_001 <- bottleneck_wide$tier0_uncertainty[bottleneck_wide$sigma_meas == 0.001]
tier3_at_001 <- bottleneck_wide$tier3_uncertainty[bottleneck_wide$sigma_meas == 0.001]
gap_at_001 <- tier3_at_001 - tier0_at_001

panel_b <- ggplot(bottleneck_data,
                  aes(x = sigma_meas, y = uncertainty, color = Scenario, shape = Scenario)) +
  # Green "acceptable" zone (matching Panel A)
  annotate("rect", xmin = 0.0008, xmax = 0.035, ymin = 0, ymax = 30,
           fill = "#90EE90", alpha = 0.2) +
  # Reference line at claimed +/-30 years (matching Panel A style)
  geom_hline(yintercept = 30, linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
  annotate("text", x = 0.022, y = 42, label = "Claimed \u00b130 yr",
           size = 3, color = "darkgreen", fontface = "italic", hjust = 1) +
  # Data lines and points
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  # Add annotation for parameter floor
  annotate("segment", x = 0.001, xend = 0.001,
           y = tier0_at_001,
           yend = tier3_at_001,
           arrow = arrow(ends = "both", length = unit(0.1, "inches")),
           color = "gray30") +
  annotate("text", x = 0.0012,
           y = (tier0_at_001 + tier3_at_001) / 2,
           label = paste0("+", round(gap_at_001), " yr\nparameter\nfloor"),
           size = 2.8, hjust = 0, color = "gray30") +
  scale_x_log10(breaks = c(0.001, 0.005, 0.01, 0.03),
                labels = c("0.001", "0.005", "0.01", "0.03")) +
  scale_y_continuous(limits = c(0, 320),
                     breaks = seq(0, 300, 50)) +
  scale_color_manual(values = c(
    "Tier 0\n(Parameters Known)" = "#0072B2",
    "Tier 3\n(Realistic Uncertainty)" = "#D55E00"
  )) +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    x = expression("IR-PAS Measurement Precision, " * sigma[meas] * " (wt%)"),
    y = "Age Uncertainty (years)",
    title = "B",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

# =============================================================================
# Combine Panels
# =============================================================================

combined_figure <- grid.arrange(panel_a, panel_b, ncol = 1, heights = c(1, 1))

# Save
png(file.path(output_dir, "figure_s7_sensitivity.png"), width = 7, height = 8,
    units = "in", res = 300)
grid.arrange(panel_a, panel_b, ncol = 1, heights = c(1, 1))
dev.off()

cat("Figure S7 saved to:", file.path(output_dir, "figure_s7_sensitivity.png"), "\n")
