# 27_discrimination.R - Scenario discrimination analysis (Table S10)
# Purpose: Tests whether OHD can distinguish between alternative chronological scenarios
# Inputs: output/tables/tier3_summary.csv
# Outputs: output/tables/discrimination_comparisons.csv, output/tables/contact_classification.csv,
#   output/tables/individual_date_postcontact_prob.csv,
#   output/tables/classification_accuracy_results.csv, output/tables/distance_required_results.csv
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

STEVENSON_SE <- 30
OUR_SE <- 299
CONTACT_YEAR <- 1722

set.seed(42)

cat("Scenario discrimination analysis\n")

# -----------------------------------------------------------------------------
# Function: Generate synthetic OHD dates from a population model
# -----------------------------------------------------------------------------

generate_synthetic_dates <- function(n, decline_year, growth_rate = 0.005,
                                     decline_rate = 0.02, start_year = 1200,
                                     end_year = 1900) {
  years <- start_year:end_year
  intensity <- numeric(length(years))

  for (i in seq_along(years)) {
    y <- years[i]
    if (y < decline_year) {
      intensity[i] <- exp(growth_rate * (y - start_year))
    } else {
      peak <- exp(growth_rate * (decline_year - start_year))
      intensity[i] <- peak * exp(-decline_rate * (y - decline_year))
    }
  }

  prob <- intensity / sum(intensity)
  dates <- sample(years, size = n, replace = TRUE, prob = prob)

  return(dates)
}

# -----------------------------------------------------------------------------
# Function: Generate SPD from dates
# -----------------------------------------------------------------------------

generate_spd <- function(dates_ad, se, year_range = c(1000, 2000), resolution = 1) {
  years <- seq(year_range[1], year_range[2], by = resolution)
  spd <- rep(0, length(years))

  if (length(se) == 1) se <- rep(se, length(dates_ad))

  for (i in seq_along(dates_ad)) {
    if (!is.na(dates_ad[i]) && !is.na(se[i])) {
      spd <- spd + dnorm(years, mean = dates_ad[i], sd = se[i])
    }
  }

  spd <- spd / sum(spd * resolution)
  return(data.frame(year = years, density = spd))
}

# -----------------------------------------------------------------------------
# Function: Calculate overlap between two SPDs (Bhattacharyya coefficient)
# -----------------------------------------------------------------------------

bhattacharyya_overlap <- function(spd1, spd2) {
  common_years <- intersect(spd1$year, spd2$year)

  p1 <- spd1$density[spd1$year %in% common_years]
  p2 <- spd2$density[spd2$year %in% common_years]

  p1 <- p1 / sum(p1)
  p2 <- p2 / sum(p2)

  bc <- sum(sqrt(p1 * p2))

  return(bc)
}

# -----------------------------------------------------------------------------
# Function: Statistical discrimination test via simulation
# -----------------------------------------------------------------------------

discrimination_test <- function(n_samples, se, decline_years, n_sim = 100) {
  n_scenarios <- length(decline_years)
  overlap_matrix <- matrix(NA, n_scenarios, n_scenarios)
  rownames(overlap_matrix) <- paste0("AD ", decline_years)
  colnames(overlap_matrix) <- paste0("AD ", decline_years)

  for (i in 1:n_scenarios) {
    for (j in 1:n_scenarios) {
      overlaps <- numeric(n_sim)

      for (sim in 1:n_sim) {
        dates_i <- generate_synthetic_dates(n_samples, decline_years[i])
        dates_j <- generate_synthetic_dates(n_samples, decline_years[j])

        spd_i <- generate_spd(dates_i, se)
        spd_j <- generate_spd(dates_j, se)

        overlaps[sim] <- bhattacharyya_overlap(spd_i, spd_j)
      }

      overlap_matrix[i, j] <- mean(overlaps)
    }
  }

  return(overlap_matrix)
}

# -----------------------------------------------------------------------------
# Define test scenarios
# -----------------------------------------------------------------------------

decline_scenarios <- c(
  1660,  # Stevenson's SA1 claim
  1705,  # Stevenson's SA2 claim
  1722,  # European contact
  1780,  # Post-contact
  1850   # Stevenson's SA3 claim
)

scenario_labels <- c(
  "AD 1660 (SA1 claim)",
  "AD 1705 (SA2 claim)",
  "AD 1722 (Contact)",
  "AD 1780 (Post-contact)",
  "AD 1850 (SA3 claim)"
)

N_SAMPLES <- 127
N_SIM <- 50

# -----------------------------------------------------------------------------
# Run discrimination tests
# -----------------------------------------------------------------------------

cat("Testing discrimination with Stevenson's +/-30 year SE...\n")
overlap_30 <- discrimination_test(N_SAMPLES, STEVENSON_SE, decline_scenarios, N_SIM)

cat("Testing discrimination with corrected +/-299 year SE...\n")
overlap_299 <- discrimination_test(N_SAMPLES, OUR_SE, decline_scenarios, N_SIM)

# -----------------------------------------------------------------------------
# Results
# -----------------------------------------------------------------------------

cat("\nOverlap coefficients (1.0 = identical, 0.0 = completely distinct)\n")

cat("\nWith +/-30 year SE (Stevenson's claim):\n")
print(round(overlap_30, 3))

cat("\nWith +/-299 year SE (corrected):\n")
print(round(overlap_299, 3))

# -----------------------------------------------------------------------------
# Key comparisons
# -----------------------------------------------------------------------------

critical_pairs <- list(
  c(1, 3),  # 1660 vs 1722 (SA1 claim vs contact)
  c(2, 3),  # 1705 vs 1722 (SA2 claim vs contact)
  c(3, 4),  # 1722 vs 1780 (contact vs post-contact)
  c(1, 5)   # 1660 vs 1850 (SA1 vs SA3)
)

comparison_names <- c(
  "SA1 claim (1660) vs Contact (1722)",
  "SA2 claim (1705) vs Contact (1722)",
  "Contact (1722) vs Post-contact (1780)",
  "SA1 claim (1660) vs SA3 claim (1850)"
)

comparison_df <- data.frame()

for (k in seq_along(critical_pairs)) {
  i <- critical_pairs[[k]][1]
  j <- critical_pairs[[k]][2]

  overlap_stev <- overlap_30[i, j]
  overlap_corr <- overlap_299[i, j]

  cat(sprintf("\n%s:\n", comparison_names[k]))
  cat(sprintf("  +/-30 yr overlap:  %.1f%%\n", overlap_stev * 100))
  cat(sprintf("  +/-299 yr overlap: %.1f%%\n", overlap_corr * 100))

  if (overlap_corr > 0.9) {
    cat("  -> Cannot distinguish with realistic uncertainty\n")
  } else if (overlap_corr > 0.7) {
    cat("  -> Marginally distinguishable (weak)\n")
  } else {
    cat("  -> Can distinguish\n")
  }

  comparison_df <- bind_rows(comparison_df, data.frame(
    comparison = comparison_names[k],
    year_diff = abs(decline_scenarios[i] - decline_scenarios[j]),
    overlap_30yr = round(overlap_stev, 3),
    overlap_299yr = round(overlap_corr, 3),
    distinguishable_30 = overlap_stev < 0.5,
    distinguishable_299 = overlap_corr < 0.5
  ))
}

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------

cat("\nCreating visualizations...\n")

# Generate example SPDs for visualization
example_dates <- list()
for (dy in decline_scenarios) {
  example_dates[[paste0("AD_", dy)]] <- generate_synthetic_dates(N_SAMPLES, dy)
}

# SPDs at +/-30 years
spd_30_all <- data.frame()
for (i in seq_along(decline_scenarios)) {
  dy <- decline_scenarios[i]
  spd <- generate_spd(example_dates[[paste0("AD_", dy)]], STEVENSON_SE)
  spd$scenario <- scenario_labels[i]
  spd_30_all <- bind_rows(spd_30_all, spd)
}
spd_30_all$uncertainty <- "\u00b130 years (Stevenson)"

# SPDs at +/-299 years
spd_299_all <- data.frame()
for (i in seq_along(decline_scenarios)) {
  dy <- decline_scenarios[i]
  spd <- generate_spd(example_dates[[paste0("AD_", dy)]], OUR_SE)
  spd$scenario <- scenario_labels[i]
  spd_299_all <- bind_rows(spd_299_all, spd)
}
spd_299_all$uncertainty <- "\u00b1299 years (Corrected)"

# Combine
spd_all <- bind_rows(spd_30_all, spd_299_all)
spd_all$uncertainty <- factor(spd_all$uncertainty,
                              levels = c("\u00b130 years (Stevenson)", "\u00b1299 years (Corrected)"))

# Plot
p_discrimination <- ggplot(spd_all, aes(x = year, y = density, color = scenario)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_vline(xintercept = CONTACT_YEAR, linetype = "dashed", color = "red") +
  facet_wrap(~uncertainty, ncol = 1, scales = "free_y") +
  scale_color_viridis_d(option = "turbo") +
  labs(
    x = "Year AD",
    y = "Probability Density",
    color = "True Decline Year"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11)
  ) +
  coord_cartesian(xlim = c(1400, 1950)) +
  guides(color = guide_legend(nrow = 2))

ggsave(file.path(fig_dir, "discrimination_test.pdf"), p_discrimination,
       width = 12, height = 8)
ggsave(file.path(fig_dir, "discrimination_test.png"), p_discrimination,
       width = 12, height = 8, dpi = 150)
cat("Saved: discrimination_test.pdf/png\n")

# Heatmap of overlap matrices
overlap_long <- data.frame()
for (i in 1:length(decline_scenarios)) {
  for (j in 1:length(decline_scenarios)) {
    overlap_long <- bind_rows(overlap_long, data.frame(
      year_i = decline_scenarios[i],
      year_j = decline_scenarios[j],
      overlap_30 = overlap_30[i, j],
      overlap_299 = overlap_299[i, j]
    ))
  }
}

overlap_long$label_i <- factor(overlap_long$year_i, labels = scenario_labels)
overlap_long$label_j <- factor(overlap_long$year_j, labels = scenario_labels)

p_heatmap <- ggplot(overlap_long, aes(x = label_i, y = label_j)) +
  geom_tile(aes(fill = overlap_299), color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", overlap_299 * 100)), size = 3) +
  scale_fill_gradient2(low = "darkgreen", mid = "yellow", high = "red",
                       midpoint = 0.5, limits = c(0, 1),
                       labels = scales::percent) +
  labs(
    title = "SPD Overlap Matrix (\u00b1299 years)",
    subtitle = "% overlap between SPDs from different true decline years",
    x = "Scenario 1",
    y = "Scenario 2",
    fill = "Overlap"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )

ggsave(file.path(fig_dir, "overlap_heatmap.pdf"), p_heatmap,
       width = 10, height = 8)
cat("Saved: overlap_heatmap.pdf\n")

# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------

write_csv(comparison_df, file.path(output_dir, "discrimination_comparisons.csv"))
cat("Saved: discrimination_comparisons.csv\n")

# Save overlap matrices
saveRDS(list(overlap_30 = overlap_30, overlap_299 = overlap_299),
        file.path(output_dir, "overlap_matrices.rds"))

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

critical_overlap <- overlap_299[1, 3]
cat(sprintf("\nCritical test: AD 1660 vs AD 1722 (62-year difference)\n"))
cat(sprintf("  Overlap at +/-30 yr:  %.1f%%\n", overlap_30[1, 3] * 100))
cat(sprintf("  Overlap at +/-299 yr: %.1f%%\n", critical_overlap * 100))

cat("\nDiscrimination analysis complete.\n")
