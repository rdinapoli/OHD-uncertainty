# 08_measurement_error.R - Measurement error sensitivity and Stevenson assumptions test
# Purpose: Sweeps measurement precision levels and tests Stevenson's fixed-parameter assumptions
# Inputs: data/stevenson_2015.csv, stan/sensitivity/tier0_meas_error_sensitivity.stan,
#   stan/sensitivity/tier3_meas_error_sensitivity.stan,
#   stan/sensitivity/tier0_stevenson_assumptions.stan
# Outputs: output/tables/measurement_error_sweep_results.csv,
#   output/tables/measurement_error_floor_summary.csv
# Runtime: ~30 minutes

library(tidyverse)
library(rstan)
library(here)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(42)

# =============================================================================
# PART 1: Measurement Error Sensitivity Sweep
# =============================================================================
#
# What fraction of the irreducible floor comes from IR-PAS measurement
# precision vs. rate equation coefficient uncertainty?
#
# Sweeps sigma_meas from 0.03 wt% (calibrated value) down to 0.001 wt%
# (physically impossible precision) under Stevenson's exact assumptions.

# --- Load Data ---------------------------------------------------------------

data_path <- here("data/stevenson_2015.csv")
df <- read_csv(data_path, show_col_types = FALSE)

df <- df %>%
  mutate(study_area_idx = as.integer(factor(study_area, levels = c("SA1", "SA2", "SA3"))))

cat("Loaded", nrow(df), "samples\n\n")

# Rate equation coefficients
a_eq <- 23.8042281295
b_eq <- 0.3782846434
c_eq <- 0.0002071707
d_eq <- -6313.2075128583
R_gas <- 8.314

# Stevenson's exact EHT values (from his data)
EHT_fixed <- c(22.4, 21.8, 23.3)  # SA1, SA2, SA3

# --- Compile Model -----------------------------------------------------------

cat("Compiling tier0_meas_error_sensitivity.stan...\n")
model <- stan_model(here("stan/sensitivity/tier0_meas_error_sensitivity.stan"))
cat("  Model compiled successfully\n\n")

# --- Define Sweep Values -----------------------------------------------------

# sigma_meas values (wt%):
# - 0.03: Our calibrated value (from manuscript)
# - 0.01: Approximate IR-PAS instrumental precision
# - 0.005: Very optimistic precision
# - 0.001: Physically impossible (test asymptotic behavior)

sigma_values <- c(0.03, 0.01, 0.005, 0.001)

cat("sigma_meas sweep values (wt%):", paste(sigma_values, collapse = ", "), "\n\n")

# --- Analytical Predictions --------------------------------------------------

cat("Calculating analytical predictions...\n")

# For age = H2Ome^2 / rate:
#   d(age)/d(H2Ome) = 2*H2Ome / rate
#   sigma_age = |d(age)/d(H2Ome)| * sigma_H2Ome = (2*H2Ome/rate) * sigma_meas

# Use median values from data
median_H2Ome <- median(df$H2Ome_pct)
median_rate <- median(df$rate)

cat(sprintf("  Median H2Ome: %.4f wt%%\n", median_H2Ome))
cat(sprintf("  Median rate: %.2e (wt%%)^2/year\n", median_rate))

analytical_predictions <- tibble(
  sigma_meas = sigma_values,
  predicted_sigma_age = (2 * median_H2Ome / median_rate) * sigma_values
)

cat("\nAnalytical predictions:\n")
for (i in 1:nrow(analytical_predictions)) {
  cat(sprintf("  sigma_meas = %.3f wt%% -> sigma_age = %.0f years\n",
              analytical_predictions$sigma_meas[i],
              analytical_predictions$predicted_sigma_age[i]))
}
cat("\n")

# --- Run Sweep ---------------------------------------------------------------

results_list <- list()

for (i in seq_along(sigma_values)) {
  sigma <- sigma_values[i]
  cat(sprintf("[%d/%d] sigma_meas = %.3f wt%%...\n", i, length(sigma_values), sigma))

  # Stan data with Stevenson's FIXED values + current sigma_meas
  stan_data <- list(
    N = nrow(df),
    N_areas = 3,
    H2Ome = df$H2Ome_pct,
    study_area = df$study_area_idx,
    H2Ot_fixed = df$H2Ot_pct,       # Use Stevenson's measured values
    Ea_fixed = df$Ea_J_per_mol,     # Use Stevenson's measured values
    EHT_fixed = EHT_fixed,          # Stevenson's study-area constants
    RH_fixed = 0.98,                # Stevenson's universal constant
    sigma_meas = sigma,             # SWEPT VALUE
    R_gas = R_gas,
    a_eq = a_eq,
    b_eq = b_eq,
    c_eq = c_eq,
    d_eq = d_eq
  )

  fit <- sampling(
    model,
    data = stan_data,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    seed = 42,
    control = list(adapt_delta = 0.95),
    refresh = 0
  )

  # Extract results
  summary_fit <- summary(fit)$summary
  age_params <- grep("^age\\[", rownames(summary_fit), value = TRUE)
  age_summary <- summary_fit[age_params, ]

  # Compute uncertainties
  ci_widths <- age_summary[, "97.5%"] - age_summary[, "2.5%"]
  half_ci_widths <- ci_widths / 2
  posterior_sds <- age_summary[, "sd"]

  results_list[[i]] <- tibble(
    sigma_meas = sigma,
    median_uncertainty = median(half_ci_widths),
    mean_uncertainty = mean(half_ci_widths),
    median_posterior_sd = median(posterior_sds),
    min_uncertainty = min(half_ci_widths),
    max_uncertainty = max(half_ci_widths),
    ratio_to_claimed = median(half_ci_widths) / 30,
    rhat_max = max(age_summary[, "Rhat"], na.rm = TRUE),
    ess_min = min(age_summary[, "n_eff"], na.rm = TRUE)
  )

  cat(sprintf("  Empirical: +/-%.1f years (%.1fx claimed)\n",
              median(half_ci_widths), median(half_ci_widths) / 30))
  cat(sprintf("  Analytical prediction: +/-%.0f years\n",
              analytical_predictions$predicted_sigma_age[i]))
  cat(sprintf("  Convergence: Rhat max = %.3f, ESS min = %.0f\n",
              max(age_summary[, "Rhat"], na.rm = TRUE),
              min(age_summary[, "n_eff"], na.rm = TRUE)))
  cat("\n")

  rm(fit)
  gc()
}

# --- Combine Results ---------------------------------------------------------

results <- bind_rows(results_list)

# Add analytical predictions
results <- results %>%
  left_join(analytical_predictions, by = "sigma_meas") %>%
  mutate(
    empirical_vs_analytical = median_uncertainty / predicted_sigma_age,
    residual_uncertainty = median_uncertainty - predicted_sigma_age
  )

# Calculate what fraction of floor is measurement error vs residual
baseline_uncertainty <- results$median_uncertainty[results$sigma_meas == 0.03]
asymptotic_residual <- results$median_uncertainty[results$sigma_meas == 0.001]

meas_error_fraction <- (baseline_uncertainty - asymptotic_residual) / baseline_uncertainty * 100
residual_fraction <- asymptotic_residual / baseline_uncertainty * 100

cat(sprintf("Floor decomposition (at sigma_meas = 0.03 wt%%):\n"))
cat(sprintf("  Measurement error contribution: %.1f%%\n", meas_error_fraction))
cat(sprintf("  Residual (rate equation/other): %.1f%%\n", residual_fraction))
cat(sprintf("\n"))

cat(sprintf("Key insight:\n"))
cat(sprintf("  - At sigma_meas = 0.03 wt%%: +/-%.0f years\n", baseline_uncertainty))
cat(sprintf("  - At sigma_meas -> 0 (asymptotic): +/-%.0f years\n", asymptotic_residual))
cat(sprintf("  - The %.0f year residual represents rate equation coefficient uncertainty\n",
            asymptotic_residual))

# --- Save Results ------------------------------------------------------------

write_csv(results, here("output/tables/measurement_error_sweep_results.csv"))

# Summary for manuscript
summary_for_manuscript <- tibble(
  analysis = "Measurement Error Floor Decomposition",
  baseline_sigma_meas = 0.03,
  baseline_uncertainty = baseline_uncertainty,
  asymptotic_residual = asymptotic_residual,
  meas_error_fraction_pct = meas_error_fraction,
  residual_fraction_pct = residual_fraction,
  analytical_vs_empirical_ratio = results$empirical_vs_analytical[results$sigma_meas == 0.03]
)

write_csv(summary_for_manuscript, here("output/tables/measurement_error_floor_summary.csv"))

cat("\nPart 1 results saved.\n\n")

# =============================================================================
# PART 2: Stevenson's Exact Assumptions Test
# =============================================================================
#
# What is the MINIMUM uncertainty even if we accept ALL of Stevenson's assumptions?
# This tests whether the problem is:
#   (a) Our "diffuse" priors (Stevenson's likely counterargument), or
#   (b) Fundamental to the method (rate equation + measurement error)
#
# If uncertainty is still large with Stevenson's exact values, the problem is fundamental.

cat("Compiling model with Stevenson's exact assumptions...\n")
model_stevenson <- stan_model(here("stan/sensitivity/tier0_stevenson_assumptions.stan"))

# Stan data with Stevenson's FIXED values
stan_data_stevenson <- list(
  N = nrow(df),
  N_areas = 3,
  H2Ome = df$H2Ome_pct,
  study_area = df$study_area_idx,
  H2Ot_fixed = df$H2Ot_pct,       # Use Stevenson's measured values
  Ea_fixed = df$Ea_J_per_mol,     # Use Stevenson's measured values
  EHT_fixed = EHT_fixed,          # Stevenson's study-area constants
  RH_fixed = 0.98,                # Stevenson's universal constant
  R_gas = R_gas,
  a_eq = a_eq,
  b_eq = b_eq,
  c_eq = c_eq,
  d_eq = d_eq
)

cat("Running with ZERO parameter uncertainty (all fixed at Stevenson's values)...\n")
fit_stevenson <- sampling(
  model_stevenson,
  data = stan_data_stevenson,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = 42,
  control = list(adapt_delta = 0.95),
  refresh = 500
)

# Extract results
summary_fit <- summary(fit_stevenson)$summary
age_params <- grep("^age\\[", rownames(summary_fit), value = TRUE)
age_summary <- summary_fit[age_params, ]

# Compute uncertainties
ci_widths <- age_summary[, "97.5%"] - age_summary[, "2.5%"]
half_ci_widths <- ci_widths / 2

cat(sprintf("\nWith Stevenson's EXACT assumptions (no parameter uncertainty):\n"))
cat(sprintf("  Median uncertainty: +/-%.1f years\n", median(half_ci_widths)))
cat(sprintf("  Mean uncertainty: +/-%.1f years\n", mean(half_ci_widths)))
cat(sprintf("  Range: +/-%.1f to +/-%.1f years\n", min(half_ci_widths), max(half_ci_widths)))
cat(sprintf("  Ratio to claimed +/-30: %.1fx\n", median(half_ci_widths) / 30))

# Convergence
cat(sprintf("\nConvergence: Rhat max = %.3f, ESS min = %.0f\n",
            max(age_summary[, "Rhat"], na.rm = TRUE),
            min(age_summary[, "n_eff"], na.rm = TRUE)))

# Compare to our main analysis
cat(sprintf("\nCOMPARISON:\n"))
cat(sprintf("  Stevenson's claim:                    +/-30 years\n"))
cat(sprintf("  Irreducible floor (this test):       +/-%.0f years (%.1fx claimed)\n",
            median(half_ci_widths), median(half_ci_widths) / 30))
cat(sprintf("  Our full uncertainty (Tier 3):       +/-299 years (10x claimed)\n"))

improvement <- 299 / median(half_ci_widths)
cat(sprintf("\n  Even accepting ALL Stevenson's assumptions reduces uncertainty by only %.1fx\n", improvement))
cat(sprintf("  The claimed +/-30 years is %.0fx smaller than the irreducible floor!\n",
            median(half_ci_widths) / 30))

# Save results
results_stevenson <- tibble(
  test = "Stevenson_exact_assumptions",
  median_uncertainty = median(half_ci_widths),
  mean_uncertainty = mean(half_ci_widths),
  min_uncertainty = min(half_ci_widths),
  max_uncertainty = max(half_ci_widths),
  ratio_to_claimed = median(half_ci_widths) / 30
)
write_csv(results_stevenson, here("output/tables/stevenson_assumptions_results.csv"))

cat("\nPart 2 results saved.\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
