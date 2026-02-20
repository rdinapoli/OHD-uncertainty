# 14_prior_predictive.R - Prior predictive checks
# Purpose: Samples from prior distributions to verify reasonable age ranges
# Inputs: data/stevenson_2015.csv
# Outputs: output/figures/figure_s3_prior_predictive.png
# Runtime: ~2 minutes

library(tidyverse)
library(truncnorm)
library(patchwork)
library(here)

set.seed(42)

# --- Rate equation coefficients (R^2=1.000 fit to Stevenson 2015) ---
a_eq  <- 23.8042281295
b_eq  <- 0.3782846434
c_eq  <- 0.0002071707
d_eq  <- -6313.2075128583
R_gas <- 8.314

# --- Observed data (for overlay) ---
empirical_data <- read_csv(
  here("data/stevenson_2015.csv"),
  show_col_types = FALSE
)
H2Ome_obs <- empirical_data$H2Ome_pct

# Stevenson's reported age range (years AD)
stevenson_age_min_ad <- 1100  # approximate earliest date
stevenson_age_max_ad <- 1850  # approximate latest date (some post-contact)

# --- Simulation parameters ---
N_sim <- 10000
meas_sigma <- 0.03  # measurement noise (wt%)

# --- Forward simulate for each tier ---
simulate_prior_predictive <- function(tier, N = N_sim) {
  cat(sprintf("Simulating Tier %d prior predictive (N=%d)...\n", tier, N))

  # === Composition priors (all tiers) ===
  # Hyperpriors
  H2Ot_mu    <- rnorm(N, mean = 0.18, sd = 0.08)
  H2Ot_sigma <- rexp(N, rate = 10)  # mean 0.1
  Ea_mu      <- rnorm(N, mean = 84000, sd = 2000)
  Ea_sigma   <- rexp(N, rate = 0.0005)  # mean 2000

  # Sample-level (truncated normal from hyperpriors)
  H2Ot <- rtruncnorm(N, a = 0.05, b = 0.43,
                      mean = H2Ot_mu, sd = H2Ot_sigma)
  Ea   <- rtruncnorm(N, a = 75000, b = 92000,
                      mean = Ea_mu, sd = Ea_sigma)

  # === Age prior (all tiers) ===
  # age ~ Normal(500, 200) truncated to [100, 900] BP
  age_bp <- rtruncnorm(N, a = 100, b = 900, mean = 500, sd = 200)

  # === Environment priors (tier-dependent) ===
  if (tier == 1) {
    # Fixed at Stevenson's mean values
    EHT <- rep(22.5, N)  # approximate mean of SA values
    RH  <- rep(0.98, N)
    e_RH_val <- rep(0, N)  # no RH effect (fixed at 0.98, term cancels)
    temporal_climate <- rep(0, N)

  } else if (tier == 2) {
    # Hierarchical EHT
    EHT_mean <- rnorm(N, mean = 22.5, sd = 2.5)
    EHT_sd   <- rexp(N, rate = 0.33)  # mean 3
    EHT      <- rtruncnorm(N, a = 18, b = 28,
                            mean = EHT_mean, sd = EHT_sd)
    # Hierarchical RH
    RH_mean  <- rnorm(N, mean = 0.85, sd = 0.1)
    RH_sd    <- rexp(N, rate = 5)  # mean 0.2
    RH       <- rtruncnorm(N, a = 0.7, b = 1.0,
                            mean = RH_mean, sd = RH_sd)
    # e_RH coefficient
    e_RH_val <- rnorm(N, mean = 2.15, sd = 0.5)
    temporal_climate <- rep(0, N)

  } else {
    # Tier 3: same as Tier 2 + temporal climate
    EHT_mean <- rnorm(N, mean = 22.5, sd = 2.5)
    EHT_sd   <- rexp(N, rate = 0.33)
    EHT      <- rtruncnorm(N, a = 18, b = 28,
                            mean = EHT_mean, sd = EHT_sd)
    RH_mean  <- rnorm(N, mean = 0.85, sd = 0.1)
    RH_sd    <- rexp(N, rate = 5)
    RH       <- rtruncnorm(N, a = 0.7, b = 1.0,
                            mean = RH_mean, sd = RH_sd)
    e_RH_val <- rnorm(N, mean = 2.15, sd = 0.5)
    temporal_climate <- rnorm(N, mean = 0, sd = 0.5)
  }

  # === Forward simulate through rate equation ===
  T_K <- EHT + temporal_climate + 273.15

  log_rate <- a_eq +
    b_eq * (-Ea / (R_gas * T_K)) +
    c_eq * H2Ot +
    d_eq * (1.0 / T_K) +
    e_RH_val * (RH - 0.98)

  # Numerical stability (match Stan model)
  log_rate <- pmax(pmin(log_rate, 50), -50)
  rate <- exp(log_rate)

  # H2Ome prediction
  H2Ome_pred <- sqrt(rate * age_bp)
  H2Ome_rep  <- rnorm(N, mean = H2Ome_pred, sd = meas_sigma)

  # Convert age to AD
  age_ad <- 1950 - age_bp

  # Filter out physically impossible values (negative H2Ome, NaN)
  valid <- is.finite(H2Ome_pred) & H2Ome_pred > 0 &
           is.finite(H2Ome_rep) & is.finite(age_ad)

  tibble(
    tier     = sprintf("Tier %d", tier),
    age_bp   = age_bp[valid],
    age_ad   = age_ad[valid],
    H2Ome_pred = H2Ome_pred[valid],
    H2Ome_rep  = H2Ome_rep[valid]
  )
}

# Run all three tiers
results <- bind_rows(
  simulate_prior_predictive(1),
  simulate_prior_predictive(2),
  simulate_prior_predictive(3)
)

# --- Summary statistics ---
summary_stats <- results %>%
  group_by(tier) %>%
  summarise(
    n_valid = n(),
    age_ad_mean = mean(age_ad),
    age_ad_sd   = sd(age_ad),
    age_ad_lo   = quantile(age_ad, 0.025),
    age_ad_hi   = quantile(age_ad, 0.975),
    H2Ome_mean  = mean(H2Ome_rep),
    H2Ome_sd    = sd(H2Ome_rep),
    H2Ome_lo    = quantile(H2Ome_rep, 0.025),
    H2Ome_hi    = quantile(H2Ome_rep, 0.975),
    .groups = "drop"
  )

for (i in seq_len(nrow(summary_stats))) {
  row <- summary_stats[i, ]
  cat(sprintf("%s (n=%d valid of %d):\n", row$tier, row$n_valid, N_sim))
  cat(sprintf("  Age (AD):  mean=%.0f, SD=%.0f, 95%% CI=[%.0f, %.0f]\n",
              row$age_ad_mean, row$age_ad_sd, row$age_ad_lo, row$age_ad_hi))
  cat(sprintf("  H2Ome (wt%%): mean=%.4f, SD=%.4f, 95%% CI=[%.4f, %.4f]\n\n",
              row$H2Ome_mean, row$H2Ome_sd, row$H2Ome_lo, row$H2Ome_hi))
}

cat(sprintf("Observed H2Ome: mean=%.4f, SD=%.4f, range=[%.4f, %.4f]\n\n",
            mean(H2Ome_obs), sd(H2Ome_obs), min(H2Ome_obs), max(H2Ome_obs)))

# --- Output directory ---
output_dir <- here("output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Figure S3: Two-panel prior predictive check ---

# Color palette for tiers
tier_colors <- c("Tier 1" = "#1b9e77", "Tier 2" = "#d95f02", "Tier 3" = "#7570b3")

# Panel A: Prior predictive H2Ome density
obs_df <- tibble(H2Ome = H2Ome_obs)

panel_a <- ggplot() +
  geom_density(data = results, aes(x = H2Ome_rep, color = tier),
               linewidth = 0.7, key_glyph = "path") +
  geom_density(data = obs_df, aes(x = H2Ome, linetype = "Observed"),
               color = "black", linewidth = 1.0) +
  scale_color_manual(values = tier_colors, name = "Prior predictive") +
  scale_linetype_manual(values = c("Observed" = "solid"), name = "") +
  coord_cartesian(xlim = c(-0.05, 0.35)) +
  labs(
    x = expression(H[2]*O[me] ~ "(wt%)"),
    y = "Density",
    title = "A. Prior predictive H\u2082O\u2098\u2091 vs. observed data"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.75),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 12)
  )

# Panel B: Prior predictive age distribution (years AD)
panel_b <- ggplot() +
  annotate("rect",
           xmin = stevenson_age_min_ad, xmax = stevenson_age_max_ad,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.5) +
  annotate("text",
           x = mean(c(stevenson_age_min_ad, stevenson_age_max_ad)),
           y = Inf, vjust = 1.5,
           label = "Stevenson reported range",
           size = 3, color = "grey40") +
  geom_density(data = results, aes(x = age_ad, color = tier),
               linewidth = 0.7) +
  scale_color_manual(values = tier_colors, name = "Prior predictive") +
  labs(
    x = "Age (years AD)",
    y = "Density",
    title = "B. Prior predictive age distribution"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.75),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 12)
  )

# Combine panels
fig_s3 <- panel_a / panel_b +
  plot_annotation(
    caption = "Figure S3. Prior predictive checks for Tiers 1\u20133. (A) Prior predictive H\u2082O\u2098\u2091 density (colored)\noverlaid with observed data (black). Priors produce observables in the right range but with wider spread.\n(B) Prior predictive age distribution (years AD). Grey band shows Stevenson's reported date range.\nPriors encompass the full Rapa Nui settlement range without generating archaeologically implausible ages.",
    theme = theme(plot.caption = element_text(size = 9, hjust = 0))
  )

# Save output
ggsave(file.path(output_dir, "figure_s3_prior_predictive.png"),
       fig_s3, width = 8, height = 9, dpi = 300)

cat("Prior predictive check plot saved.\n")
cat(sprintf("  %s\n", file.path(output_dir, "figure_s3_prior_predictive.png")))
