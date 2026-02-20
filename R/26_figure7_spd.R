# 26_figure7_spd.R - SPD comparison figure (Figure 7)
# Purpose: Side-by-side SPD comparison at +/-30 vs +/-300 year uncertainty
# Inputs: output/tables/spd_envelopes.csv
# Outputs: output/figures/figure7_spd.png
# Runtime: ~1 minute

library(tidyverse)
library(patchwork)
library(here)

# Output directory
output_dir <- here("output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
data <- read_csv(here("data/stevenson_2015.csv"), show_col_types = FALSE)

cat("Generating bootstrap SPDs for Figure 7...\n")

# Parameters
CONTACT_YEAR <- 1722
STEVENSON_SE <- 30
CORRECTED_SE <- 299
N_BOOTSTRAP <- 500
YEAR_RANGE <- c(500, 2500)

# Study area info
sa_info <- list(
  SA1 = list(name = "SA1 (Te Niu)", claim = "AD 1660", color_stev = "#2166AC", color_corr = "#B2182B"),
  SA2 = list(name = "SA2 (Maunga O'Koro)", claim = "AD 1705", color_stev = "#2166AC", color_corr = "#B2182B"),
  SA3 = list(name = "SA3 (Anakena West)", claim = "AD 1850", color_stev = "#2166AC", color_corr = "#B2182B")
)

# Bootstrap SPD function with confidence envelopes
bootstrap_spd_with_ci <- function(dates_ad, se, n_boot = 500, year_range = c(1100, 1950)) {
  years <- seq(year_range[1], year_range[2], by = 1)
  n_dates <- length(dates_ad[!is.na(dates_ad)])
  dates_clean <- dates_ad[!is.na(dates_ad)]

  boot_matrix <- matrix(NA, nrow = n_boot, ncol = length(years))

  for (b in 1:n_boot) {
    # Resample dates with replacement
    boot_dates <- sample(dates_clean, n_dates, replace = TRUE)
    # Add measurement error
    noisy_dates <- boot_dates + rnorm(n_dates, mean = 0, sd = se)

    # Calculate SPD
    spd <- rep(0, length(years))
    for (d in noisy_dates) {
      spd <- spd + dnorm(years, mean = d, sd = se)
    }
    spd <- spd / sum(spd)  # Normalize
    boot_matrix[b, ] <- spd
  }

  # Calculate percentiles
  spd_mean <- colMeans(boot_matrix)
  spd_lo_95 <- apply(boot_matrix, 2, quantile, probs = 0.025)
  spd_hi_95 <- apply(boot_matrix, 2, quantile, probs = 0.975)
  spd_lo_50 <- apply(boot_matrix, 2, quantile, probs = 0.25)
  spd_hi_50 <- apply(boot_matrix, 2, quantile, probs = 0.75)

  tibble(
    year = years,
    spd_mean = spd_mean,
    spd_lo_95 = spd_lo_95,
    spd_hi_95 = spd_hi_95,
    spd_lo_50 = spd_lo_50,
    spd_hi_50 = spd_hi_50
  )
}

# Generate SPDs for all study areas and both uncertainty levels
spd_results <- list()

for (sa in c("SA1", "SA2", "SA3")) {
  cat(sprintf("  Processing %s...\n", sa))

  dates <- data %>%
    filter(study_area == sa) %>%
    pull(date_ad)

  # Stevenson uncertainty (+/-30 yr)
  cat(sprintf("    - %s with +/-30 years...\n", sa))
  spd_stev <- bootstrap_spd_with_ci(dates, STEVENSON_SE, N_BOOTSTRAP, YEAR_RANGE) %>%
    mutate(study_area = sa, uncertainty = "\u00b130 years (Stevenson)")

  # Corrected uncertainty (+/-299 yr)
  cat(sprintf("    - %s with +/-299 years...\n", sa))
  spd_corr <- bootstrap_spd_with_ci(dates, CORRECTED_SE, N_BOOTSTRAP, YEAR_RANGE) %>%
    mutate(study_area = sa, uncertainty = "\u00b1299 years (Corrected)")

  spd_results[[paste0(sa, "_stev")]] <- spd_stev
  spd_results[[paste0(sa, "_corr")]] <- spd_corr
}

# Combine all results
all_spd <- bind_rows(spd_results)

cat("Creating figure panels...\n")

# Create individual panels
create_panel <- function(sa, unc_label, color_fill, color_line, title_suffix = "") {

  plot_data <- all_spd %>%
    filter(study_area == sa, uncertainty == unc_label)

  sa_name <- sa_info[[sa]]$name

  ggplot(plot_data, aes(x = year)) +
    # 95% CI ribbon (lighter)
    geom_ribbon(aes(ymin = spd_lo_95, ymax = spd_hi_95),
                fill = color_fill, alpha = 0.2) +
    # 50% CI ribbon (darker)
    geom_ribbon(aes(ymin = spd_lo_50, ymax = spd_hi_50),
                fill = color_fill, alpha = 0.4) +
    # Mean line
    geom_line(aes(y = spd_mean), color = color_line, linewidth = 1) +
    # Contact line
    geom_vline(xintercept = CONTACT_YEAR, linetype = "dashed", color = "red", linewidth = 0.8) +
    # Labels
    labs(
      title = paste0(sa_name, title_suffix),
      x = NULL,
      y = "Probability Density"
    ) +
    scale_x_continuous(limits = YEAR_RANGE, breaks = seq(500, 2500, 500)) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

# Find global y max for matched axes
y_max_stev <- all_spd %>%
  filter(grepl("30", uncertainty)) %>%
  pull(spd_hi_95) %>%
  max()

y_max_corr <- all_spd %>%
  filter(grepl("299", uncertainty)) %>%
  pull(spd_hi_95) %>%
  max()

# Create panels with fixed y-axis
create_panel_fixed <- function(sa, unc_label, color_fill, color_line, y_max) {

  plot_data <- all_spd %>%
    filter(study_area == sa, uncertainty == unc_label)

  sa_name <- sa_info[[sa]]$name

  ggplot(plot_data, aes(x = year)) +
    geom_ribbon(aes(ymin = spd_lo_95, ymax = spd_hi_95),
                fill = color_fill, alpha = 0.2) +
    geom_ribbon(aes(ymin = spd_lo_50, ymax = spd_hi_50),
                fill = color_fill, alpha = 0.4) +
    geom_line(aes(y = spd_mean), color = color_line, linewidth = 1) +
    geom_vline(xintercept = CONTACT_YEAR, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      title = sa_name,
      x = NULL,
      y = "Probability Density"
    ) +
    scale_x_continuous(limits = YEAR_RANGE, breaks = seq(500, 2500, 500)) +
    scale_y_continuous(limits = c(0, y_max * 1.05)) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

# Use SAME y_max for both rows to show the flattening effect
y_max_global <- y_max_stev

p_sa1_stev_f <- create_panel_fixed("SA1", "\u00b130 years (Stevenson)", "#2166AC", "#2166AC", y_max_global)
p_sa2_stev_f <- create_panel_fixed("SA2", "\u00b130 years (Stevenson)", "#2166AC", "#2166AC", y_max_global)
p_sa3_stev_f <- create_panel_fixed("SA3", "\u00b130 years (Stevenson)", "#2166AC", "#2166AC", y_max_global)

p_sa1_corr_f <- create_panel_fixed("SA1", "\u00b1299 years (Corrected)", "#B2182B", "#B2182B", y_max_global)
p_sa2_corr_f <- create_panel_fixed("SA2", "\u00b1299 years (Corrected)", "#B2182B", "#B2182B", y_max_global)
p_sa3_corr_f <- create_panel_fixed("SA3", "\u00b1299 years (Corrected)", "#B2182B", "#B2182B", y_max_global)

combined_plot_matched <- (p_sa1_stev_f | p_sa2_stev_f | p_sa3_stev_f) /
                          (p_sa1_corr_f | p_sa2_corr_f | p_sa3_corr_f) +
  plot_annotation(
    theme = theme(plot.margin = margin(t = 5))
  )

ggsave(file.path(output_dir, "figure7_spd.pdf"), combined_plot_matched,
       width = 14, height = 8, dpi = 300)
ggsave(file.path(output_dir, "figure7_spd.png"), combined_plot_matched,
       width = 14, height = 8, dpi = 300)

cat("Saved: figure7_spd.pdf/png\n")

# Summary statistics
summary_stats <- all_spd %>%
  group_by(study_area, uncertainty) %>%
  summarize(
    peak_year = year[which.max(spd_mean)],
    ci_95_width = sum(spd_hi_95 > 0.0001),  # Approximate width
    mean_density_at_contact = spd_mean[year == CONTACT_YEAR],
    .groups = "drop"
  )

cat("\nSPD summary statistics:\n")
print(summary_stats)

cat("\nFigure 7 complete.\n")
