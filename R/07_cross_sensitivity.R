# 07_cross_sensitivity.R - Cross-parameter sensitivity and simultaneous constraint tests
# Purpose: Tests parameter interactions and simultaneous tightening of all priors
# Inputs: data/stevenson_2015.csv, stan/sensitivity/tier2_eht_sensitivity.stan,
#   stan/sensitivity/tier2_simultaneous_tight.stan
# Outputs: output/tables/cross_sensitivity_results.csv,
#   output/tables/cross_sensitivity_summary.csv,
#   output/tables/simultaneous_constraint_results.csv
# Runtime: ~30 minutes

library(tidyverse)
library(rstan)
library(here)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(42)

# =============================================================================
# PART 1: Cross-Sensitivity Sweep (Tier 0 vs Tier 3)
# =============================================================================
#
# What happens when measurement precision improves but environmental
# parameters remain uncertain? Does parameter uncertainty create an independent
# floor that persists regardless of sigma_meas?
#
# Design: 2x3 Factorial
#   - Tier 0 (params known): from measurement_error_sweep results
#   - Tier 3 (params uncertain): run here at 3 sigma_meas values

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
e_RH_mean <- 2.15
e_RH_sd <- 0.5

# --- Load Existing Tier 0 Results -------------------------------------------

cat("Loading existing Tier 0 results...\n")
tier0_results <- read_csv(here("output/tables/measurement_error_sweep_results.csv"),
                          show_col_types = FALSE)

# Filter to the 3 values we're using
tier0_subset <- tier0_results %>%
  filter(sigma_meas %in% c(0.03, 0.01, 0.001)) %>%
  select(sigma_meas, median_uncertainty) %>%
  rename(tier0_uncertainty = median_uncertainty)

cat("Tier 0 reference values:\n")
print(tier0_subset)
cat("\n")

# --- Compile Tier 3 Model ---------------------------------------------------

cat("Compiling tier3_meas_error_sensitivity.stan...\n")
model <- stan_model(here("stan/sensitivity/tier3_meas_error_sensitivity.stan"))
cat("  Model compiled successfully\n\n")

# --- Run Tier 3 Sweep -------------------------------------------------------

sigma_values <- c(0.03, 0.01, 0.001)
cat("sigma_meas sweep values (wt%):", paste(sigma_values, collapse = ", "), "\n\n")

results_list <- list()

for (i in seq_along(sigma_values)) {
  sigma <- sigma_values[i]
  cat(sprintf("[%d/%d] Tier 3 with sigma_meas = %.3f wt%%...\n", i, length(sigma_values), sigma))

  # Stan data for Tier 3 (full model with parameter uncertainty)
  stan_data <- list(
    N = nrow(df),
    N_areas = 3,
    H2Ome = df$H2Ome_pct,
    study_area = df$study_area_idx,
    sigma_meas = sigma,                  # SWEPT VALUE
    temporal_climate_sd = 0.5,           # +/-0.5 C temporal variation
    R_gas = R_gas,
    a_eq = a_eq,
    b_eq = b_eq,
    c_eq = c_eq,
    d_eq = d_eq,
    e_RH_mean = e_RH_mean,
    e_RH_sd = e_RH_sd
  )

  fit <- sampling(
    model,
    data = stan_data,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    seed = 42,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 500
  )

  # Extract results
  summary_fit <- summary(fit)$summary
  age_params <- grep("^age\\[", rownames(summary_fit), value = TRUE)
  age_summary <- summary_fit[age_params, ]

  # Compute uncertainties
  ci_widths <- age_summary[, "97.5%"] - age_summary[, "2.5%"]
  half_ci_widths <- ci_widths / 2

  results_list[[i]] <- tibble(
    sigma_meas = sigma,
    tier3_uncertainty = median(half_ci_widths),
    tier3_mean = mean(half_ci_widths),
    tier3_min = min(half_ci_widths),
    tier3_max = max(half_ci_widths),
    rhat_max = max(age_summary[, "Rhat"], na.rm = TRUE),
    ess_min = min(age_summary[, "n_eff"], na.rm = TRUE)
  )

  cat(sprintf("  Tier 3: +/-%.1f years (%.1fx claimed)\n",
              median(half_ci_widths), median(half_ci_widths) / 30))
  cat(sprintf("  Convergence: Rhat max = %.3f, ESS min = %.0f\n",
              max(age_summary[, "Rhat"], na.rm = TRUE),
              min(age_summary[, "n_eff"], na.rm = TRUE)))
  cat("\n")

  rm(fit)
  gc()
}

# --- Combine and Save -------------------------------------------------------

tier3_results <- bind_rows(results_list)

# Merge with Tier 0 results
combined <- tier3_results %>%
  left_join(tier0_subset, by = "sigma_meas") %>%
  mutate(
    gap = tier3_uncertainty - tier0_uncertainty,
    ratio = tier3_uncertainty / tier0_uncertainty
  )

# Key analysis: What happens as sigma_meas -> 0?
tier3_at_001 <- combined$tier3_uncertainty[combined$sigma_meas == 0.001]
tier0_at_001 <- combined$tier0_uncertainty[combined$sigma_meas == 0.001]
gap_at_001 <- combined$gap[combined$sigma_meas == 0.001]

cat("KEY FINDING: The divergence at sigma_meas = 0.001 wt%\n")
cat(sprintf("  Tier 0 (params known):     +/-%.0f years\n", tier0_at_001))
cat(sprintf("  Tier 3 (params uncertain): +/-%.0f years\n", tier3_at_001))
cat(sprintf("  GAP (parameter floor):     +%.0f years\n", gap_at_001))
cat("\n")

if (tier3_at_001 > 100) {
  cat("INTERPRETATION: Parameter uncertainty creates an independent floor.\n")
  cat("Even with perfect measurement, OHD cannot achieve sub-centennial precision\n")
  cat("because burial conditions (EHT, RH) remain unknowable.\n")
} else {
  cat("INTERPRETATION: Measurement precision is the dominant bottleneck.\n")
  cat("Improving IR-PAS precision could substantially improve OHD accuracy.\n")
}

write_csv(combined, here("output/tables/cross_sensitivity_results.csv"))

# Summary for manuscript
summary_for_manuscript <- tibble(
  analysis = "Cross-Sensitivity (Tier 0 vs Tier 3)",
  sigma_meas_values = "0.03, 0.01, 0.001",
  tier0_at_001 = tier0_at_001,
  tier3_at_001 = tier3_at_001,
  parameter_floor = gap_at_001,
  conclusion = ifelse(tier3_at_001 > 100,
                      "Dual bottleneck: measurement AND parameter uncertainty",
                      "Single bottleneck: measurement precision dominates")
)

write_csv(summary_for_manuscript, here("output/tables/cross_sensitivity_summary.csv"))

cat("Part 1 results saved.\n\n")

# =============================================================================
# PART 2: Simultaneous Constraint Test
# =============================================================================
#
# What happens when ALL parameters are constrained to their tightest values?
#   - EHT_sd = 0.1 C (30x tighter than baseline)
#   - H2Ot_sigma = 0.005 wt% (16x tighter than baseline)
#   - Ea_sigma = 100 J/mol (20x tighter than baseline)
#
# If this still produces +/-250+ years, even physically impossible joint
# precision cannot recover Stevenson's claims.

cat("Compiling simultaneous constraint model...\n")
model_simult <- stan_model(here("stan/sensitivity/tier2_simultaneous_tight.stan"))
cat("  Model compiled successfully\n\n")

# Define test configurations
tests <- list(
  list(
    name = "Baseline (our analysis)",
    EHT_sd = 3.0,
    H2Ot_sigma = 0.08,
    Ea_sigma = 2000
  ),
  list(
    name = "Individual tight (best single)",
    EHT_sd = 0.1,
    H2Ot_sigma = 0.08,
    Ea_sigma = 2000
  ),
  list(
    name = "Simultaneous tight (all three)",
    EHT_sd = 0.1,
    H2Ot_sigma = 0.005,
    Ea_sigma = 100
  ),
  list(
    name = "Extreme tight (physically impossible)",
    EHT_sd = 0.01,
    H2Ot_sigma = 0.001,
    Ea_sigma = 10
  )
)

results_list_simult <- list()

for (i in seq_along(tests)) {
  test <- tests[[i]]
  cat(sprintf("\n[%d/%d] %s\n", i, length(tests), test$name))
  cat(sprintf("  EHT_sd = %.2f C, H2Ot_sigma = %.4f wt%%, Ea_sigma = %.0f J/mol\n",
              test$EHT_sd, test$H2Ot_sigma, test$Ea_sigma))

  stan_data <- list(
    N = nrow(df),
    N_areas = 3,
    H2Ome = df$H2Ome_pct,
    study_area = df$study_area_idx,
    R_gas = R_gas,
    a_eq = a_eq,
    b_eq = b_eq,
    c_eq = c_eq,
    d_eq = d_eq,
    e_RH_mean = e_RH_mean,
    e_RH_sd = e_RH_sd,
    EHT_sd_prior = test$EHT_sd,
    H2Ot_sigma_prior = test$H2Ot_sigma,
    Ea_sigma_prior = test$Ea_sigma
  )

  fit <- sampling(
    model_simult,
    data = stan_data,
    chains = 4,
    iter = 4000,
    warmup = 2000,
    seed = 42,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 500
  )

  summary_fit <- summary(fit)$summary
  age_params <- grep("^age\\[", rownames(summary_fit), value = TRUE)
  age_summary <- summary_fit[age_params, ]

  ci_widths <- age_summary[, "97.5%"] - age_summary[, "2.5%"]
  half_ci_widths <- ci_widths / 2

  cat(sprintf("  Result: +/-%.1f years (%.1fx claimed)\n",
              median(half_ci_widths), median(half_ci_widths) / 30))
  cat(sprintf("  Convergence: Rhat max = %.3f, ESS min = %.0f\n",
              max(age_summary[, "Rhat"], na.rm = TRUE),
              min(age_summary[, "n_eff"], na.rm = TRUE)))

  results_list_simult[[i]] <- tibble(
    test_name = test$name,
    EHT_sd = test$EHT_sd,
    H2Ot_sigma = test$H2Ot_sigma,
    Ea_sigma = test$Ea_sigma,
    median_uncertainty = median(half_ci_widths),
    mean_uncertainty = mean(half_ci_widths),
    min_uncertainty = min(half_ci_widths),
    max_uncertainty = max(half_ci_widths),
    ratio_to_claimed = median(half_ci_widths) / 30,
    rhat_max = max(age_summary[, "Rhat"], na.rm = TRUE),
    ess_min = min(age_summary[, "n_eff"], na.rm = TRUE)
  )

  rm(fit)
  gc()
}

# --- Summary -----------------------------------------------------------------

results_simult <- bind_rows(results_list_simult)

baseline_uncertainty <- results_simult$median_uncertainty[1]
simultaneous_tight <- results_simult$median_uncertainty[results_simult$test_name == "Simultaneous tight (all three)"]
extreme_tight <- results_simult$median_uncertainty[results_simult$test_name == "Extreme tight (physically impossible)"]

cat(sprintf("\nBaseline (our analysis):              +/-%.0f years (%.1fx claimed)\n",
            baseline_uncertainty, baseline_uncertainty / 30))
cat(sprintf("Simultaneous tight constraints:       +/-%.0f years (%.1fx claimed)\n",
            simultaneous_tight, simultaneous_tight / 30))
cat(sprintf("Extreme tight (physically impossible): +/-%.0f years (%.1fx claimed)\n",
            extreme_tight, extreme_tight / 30))

cat(sprintf("\nReduction from baseline to simultaneous tight: %.1f%%\n",
            (1 - simultaneous_tight / baseline_uncertainty) * 100))
cat(sprintf("Reduction from baseline to extreme tight:      %.1f%%\n",
            (1 - extreme_tight / baseline_uncertainty) * 100))

if (simultaneous_tight > 250) {
  cat("\nEven with physically impossible joint precision on all parameters,\n")
  cat(sprintf("uncertainty remains +/-%.0f years, still %.1fx Stevenson's claim.\n",
              simultaneous_tight, simultaneous_tight / 30))
} else if (simultaneous_tight > 100) {
  cat(sprintf("\n'Perfect knowledge' achieves +/-%.0f years, still %.1fx claimed.\n",
              simultaneous_tight, simultaneous_tight / 30))
  cat("This is still insufficient for pre/post-contact discrimination (62 years).\n")
}

write_csv(results_simult, here("output/tables/simultaneous_constraint_results.csv"))

cat("\nPart 2 results saved.\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
