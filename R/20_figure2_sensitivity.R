# 20_figure2_sensitivity.R - EHT sensitivity and uncertainty propagation (Figure 2)
# Purpose: Shows how calculated date varies with assumed EHT, plus fixed vs. Bayesian comparison
# Inputs: data/stevenson_2015.csv
# Outputs: output/figures/figure2_sensitivity.png
# Runtime: ~1 minute

library(tidyverse)
library(patchwork)
library(here)

base_dir <- here::here()
output_dir <- file.path(base_dir, "output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Rate equation (from sensitivity_variance_decomposition.R lines 48-71)

a_eq <- 23.8042281295
b_eq <- 0.3782846434
c_eq <- 0.0002071707
d_eq <- -6313.2075128583
R_gas <- 8.314
e_RH_default <- 2.15

compute_age <- function(H2Ome, H2Ot, Ea, EHT, RH, e_RH = e_RH_default) {
  T_K <- EHT + 273.15
  log_rate <- a_eq + b_eq * (-Ea / (R_gas * T_K)) +
              c_eq * H2Ot + d_eq * (1 / T_K) +
              e_RH * (RH - 0.98)
  log_rate <- pmax(pmin(log_rate, 50), -50)
  rate <- exp(log_rate)
  age <- H2Ome^2 / rate
  return(age)
}

# Load data and compute representative samples

dat <- read_csv(file.path(base_dir, "data/stevenson_2015.csv"),
                show_col_types = FALSE)

# Study area parameters
sa_info <- tibble(
  study_area = c("SA1", "SA2", "SA3"),
  eht_assumed = c(22.4, 21.8, 23.3),
  label = c("SA1 (Anakena)", "SA2 (Hanga Ho'onu)", "SA3 (Ahu Vinap\u00fa)")
)

# Representative sample = median composition values per SA
rep_samples <- dat %>%
  group_by(study_area) %>%
  summarize(
    H2Ome = median(H2Ome_pct),
    H2Ot  = median(H2Ot_pct),
    Ea    = median(Ea_J_per_mol),
    .groups = "drop"
  ) %>%
  left_join(sa_info, by = "study_area")

# Verify: Stevenson's dates at assumed EHT
rep_samples <- rep_samples %>%
  mutate(
    age_stevenson = compute_age(H2Ome, H2Ot, Ea, eht_assumed, 0.98),
    date_ad_stevenson = 2015 - age_stevenson
  )

cat("Representative samples (median composition per SA):\n")
for (i in 1:nrow(rep_samples)) {
  r <- rep_samples[i, ]
  cat(sprintf("  %s: H2Ome=%.4f, H2Ot=%.4f, Ea=%.0f, EHT=%.1f\u00b0C -> age=%.0f BP (AD %.0f)\n",
              r$study_area, r$H2Ome, r$H2Ot, r$Ea, r$eht_assumed, r$age_stevenson, r$date_ad_stevenson))
}

# Shared theme

theme_pub <- theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "lines"),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(size = 10, face = "bold")
  )

# Panel A: EHT sensitivity

cat("Generating Panel A: Date shifts from EHT scenarios\n")

# Use SA1 median sample: same measurement, different temperature assumptions
sa1 <- rep_samples %>% filter(study_area == "SA1")
eht_assumed <- sa1$eht_assumed  # 22.4 deg C

# Compute point estimate at each plausible EHT
eht_scenarios <- seq(18, 27, by = 1)
pts_df <- tibble(
  EHT = eht_scenarios,
  age_bp = compute_age(sa1$H2Ome, sa1$H2Ot, sa1$Ea, eht_scenarios, 0.98),
  date_ad = 2015 - age_bp,
  is_assumed = abs(EHT - eht_assumed) < 0.5
)

# Stevenson's exact assumed point
stev_date <- 2015 - compute_age(sa1$H2Ome, sa1$H2Ot, sa1$Ea, eht_assumed, 0.98)

cat("Point estimates by assumed EHT (same sample):\n")
for (i in 1:nrow(pts_df)) {
  p <- pts_df[i, ]
  shift <- p$date_ad - stev_date
  cat(sprintf("  EHT = %d\u00b0C -> AD %.0f (%+.0f yr from Stevenson)\n",
              p$EHT, p$date_ad, shift))
}

contact_year <- 1722

p_a <- ggplot(pts_df, aes(x = date_ad, y = EHT)) +
  # Light connecting line through all points
  geom_line(color = "grey60", linewidth = 0.5) +
  # Horizontal reference line at Stevenson's assumed EHT
  geom_hline(yintercept = eht_assumed, linetype = "dotted",
             color = "grey60", linewidth = 0.4) +
  # All scenario points
  geom_point(
    data = filter(pts_df, !is_assumed),
    size = 3, shape = 16, color = "grey40"
  ) +
  # Stevenson's assumed point (distinct)
  geom_point(
    aes(x = stev_date, y = eht_assumed),
    size = 4, shape = 21, fill = "black", color = "black", stroke = 1
  ) +
  # +/-30-year error bars on Stevenson's point only
  geom_errorbar(
    aes(y = eht_assumed, xmin = stev_date - 30, xmax = stev_date + 30),
    width = 0.3, linewidth = 0.6, color = "black", orientation = "y"
  ) +
  # EHT labels on left side of each point
  geom_text(
    data = filter(pts_df, !is_assumed),
    aes(label = sprintf("%d\u00b0C", EHT)),
    hjust = 1.3, vjust = 0.4, size = 2.5, color = "grey40"
  ) +
  # European contact
  geom_vline(xintercept = contact_year, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  annotate("text", x = contact_year - 8, y = 18.5,
           label = "European contact\n(AD 1722)",
           size = 2.2, hjust = 1, color = "grey50") +
  scale_x_continuous(breaks = seq(1200, 1800, by = 100)) +
  coord_cartesian(xlim = c(1200, 1800)) +
  labs(
    x = "Calculated date (AD)",
    y = "Assumed EHT (\u00b0C)"
  ) +
  theme_pub

# Panel B: Two-column contrast (Stevenson vs. Bayesian)

cat("Generating Panel B: Two-column contrast schematic\n")

# Helper: generate density curve tibble positioned at (cx, cy)
make_density <- function(x_vals, d_vals, cx, cy, width, height) {
  x_scaled <- cx - width/2 + (x_vals - min(x_vals)) / (max(x_vals) - min(x_vals)) * width
  d_scaled <- cy + d_vals / max(d_vals) * height
  tibble(x = x_scaled, y = d_scaled, ybase = cy)
}

# Left column (Stevenson): fixed point values -> narrow +/-30 yr (their claim)
# Right column (Bayesian): distributions -> "?" (the question this paper answers)

# --- Input parameter distributions ---
eht_x <- seq(15, 30, length.out = 100)
eht_d <- dnorm(eht_x, 22.4, 2.5)

rh_x <- seq(0.45, 1.25, length.out = 100)
rh_d <- dnorm(rh_x, 0.85, 0.10)

h2ot_x <- seq(-0.14, 0.50, length.out = 100)
h2ot_d <- dnorm(h2ot_x, 0.18, 0.08)

ea_x <- seq(78, 90, length.out = 100)
ea_d <- dnorm(ea_x, 84, 2)

# --- Layout parameters ---
param_y <- c(5.0, 4.2, 3.4, 2.6)
param_names <- c("EHT", "RH", "H\u2082Ot", "Ea")
stev_values <- c("22.4\u00b0C", "0.98", "0.10 wt%", "84.4 kJ/mol")
bayes_priors <- c("\u00b12.5\u00b0C", "\u00b10.10", "\u00b10.08 wt%", "\u00b12 kJ/mol")

col_left  <- 1.8   # Stevenson
col_right <- 5.2   # Bayesian

# Stevenson side: delta functions
stev_delta <- tibble(
  x = col_left, y = param_y,
  label = stev_values, name = param_names
)

# Bayesian side: mini density curves
bayes_dens <- bind_rows(
  make_density(eht_x, eht_d, col_right, param_y[1] - 0.25, 1.4, 0.4) %>% mutate(param = "EHT"),
  make_density(rh_x, rh_d, col_right, param_y[2] - 0.25, 1.4, 0.4) %>% mutate(param = "RH"),
  make_density(h2ot_x, h2ot_d, col_right, param_y[3] - 0.25, 1.4, 0.4) %>% mutate(param = "H2Ot"),
  make_density(ea_x, ea_d, col_right, param_y[4] - 0.25, 1.4, 0.4) %>% mutate(param = "Ea")
)

# Rate equation boxes
rate_box_left  <- tibble(xmin = 0.8, xmax = 2.8, ymin = 1.6, ymax = 2.1)
rate_box_right <- tibble(xmin = 4.2, xmax = 6.2, ymin = 1.6, ymax = 2.1)

# Output: Stevenson narrow spike (their published claim, not our result)
out_x_ib <- seq(1200, 2000, length.out = 200)
out_stev <- dnorm(out_x_ib, 1597, 30)
out_stev_df <- make_density(out_x_ib, out_stev, col_left, 0.3, 1.6, 0.9)

# Output: Bayesian side gets a dashed uncertain shape + "?"
# Use a wide generic bell to suggest "unknown width", drawn with dashed line
out_unknown <- dnorm(out_x_ib, 1597, 120)
out_unknown_df <- make_density(out_x_ib, out_unknown, col_right, 0.3, 1.6, 0.9)

# Arrows
arrow_left_out  <- tibble(x = col_left, xend = col_left, y = 1.6, yend = 1.3)
arrow_right_out <- tibble(x = col_right, xend = col_right, y = 1.6, yend = 1.3)

p_b <- ggplot() +
  # Column headers
  annotate("text", x = col_left, y = 5.85,
           label = "Fixed-parameter\napproach",
           size = 2.8, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = col_left, y = 5.45,
           label = "Single assumed values",
           size = 2.0, color = "grey45") +
  annotate("text", x = col_right, y = 5.85,
           label = "Bayesian\napproach",
           size = 2.8, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = col_right, y = 5.45,
           label = "Parameter distributions",
           size = 2.0, color = "grey45") +
  # Vertical divider
  annotate("segment", x = 3.5, xend = 3.5, y = 0.1, yend = 5.9,
           linetype = "dotted", color = "grey70",
           linewidth = 0.4) +
  # Parameter names (shared, on the divider)
  annotate("text", x = 3.5, y = param_y + 0.05,
           label = param_names,
           size = 2.4, fontface = "bold", color = "grey40") +
  # Fixed-parameter side: point markers with values
  geom_point(data = stev_delta,
             aes(x = x, y = y - 0.18),
             shape = 16, size = 3, color = "grey20") +
  geom_text(data = stev_delta,
            aes(x = x, y = y - 0.40, label = label),
            size = 2.2, color = "grey40") +
  # Bayesian side: density curves
  geom_ribbon(data = bayes_dens,
              aes(x = x, ymin = ybase, ymax = y, group = param),
              fill = "grey60", alpha = 0.5, color = "grey30", linewidth = 0.3) +
  geom_text(data = tibble(x = col_right, y = param_y - 0.42, label = bayes_priors),
            aes(x = x, y = y, label = label),
            size = 2.2, color = "grey40") +
  # Rate equation boxes
  geom_rect(data = rate_box_left,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", color = "grey30", linewidth = 0.5) +
  annotate("text", x = col_left, y = 1.85, label = "Rate equation",
           size = 2.5, fontface = "bold") +
  geom_rect(data = rate_box_right,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", color = "grey30", linewidth = 0.5) +
  annotate("text", x = col_right, y = 1.85, label = "Rate equation",
           size = 2.5, fontface = "bold") +
  # Arrows: rate box to output
  geom_segment(data = arrow_left_out,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               color = "grey40", linewidth = 0.3) +
  geom_segment(data = arrow_right_out,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               color = "grey40", linewidth = 0.3) +
  # Output: Stevenson narrow spike (their published claim)
  geom_ribbon(data = out_stev_df,
              aes(x = x, ymin = ybase, ymax = y),
              fill = "grey30", alpha = 0.7, color = "grey15", linewidth = 0.3) +
  annotate("text", x = col_left, y = 0.10, label = "\u00b130 yr (claimed)",
           size = 2.0, fontface = "bold") +
  # Output: Bayesian -- dashed uncertain shape with "?"
  geom_ribbon(data = out_unknown_df,
              aes(x = x, ymin = ybase, ymax = y),
              fill = "grey85", alpha = 0.4, color = "grey45",
              linewidth = 0.4, linetype = "dashed") +
  annotate("text", x = col_right, y = 0.75, label = "?",
           size = 10, fontface = "bold", color = "grey35") +
  # Output axis label
  annotate("text", x = 3.5, y = -0.10, label = "Date estimate (AD)",
           size = 2.0, color = "grey50") +
  coord_cartesian(xlim = c(0.0, 7.0), ylim = c(-0.25, 6.15)) +
  theme_void() +
  theme(
    text = element_text(family = "sans"),
    plot.margin = margin(5, 5, 5, 5)
  )

# Compose Figure 2: Panel A (sensitivity) + Panel B (schematic)

cat("Composing Figure 2 composite\n")

fig2 <- (p_a + ggtitle("A")) /
  (p_b + ggtitle("B")) +
  plot_layout(heights = c(5, 5.5)) &
  theme(plot.title = element_text(
    size = 14, face = "bold", hjust = 0
  ))

ggsave(file.path(output_dir, "figure2_sensitivity.png"),
       fig2, width = 4, height = 10.5, dpi = 300)
ggsave(file.path(output_dir, "figure2_sensitivity.pdf"),
       fig2, width = 4, height = 10.5,
       device = cairo_pdf)

cat("Figure 2 saved to:", output_dir, "\n")
