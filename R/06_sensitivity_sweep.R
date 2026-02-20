# 06_sensitivity_sweep.R - Individual parameter prior sensitivity analysis
# Purpose: Sweeps each parameter prior width independently to assess sensitivity
# Inputs: data/stevenson_2015.csv, stan/sensitivity/tier2_eht_sensitivity.stan,
#   stan/sensitivity/tier2_h2ot_sensitivity.stan, stan/sensitivity/tier2_ea_sensitivity.stan,
#   stan/sensitivity/tier2_rh_sensitivity.stan
# Outputs: output/tables/eht_sensitivity_results.csv,
#   output/tables/comprehensive_sensitivity_results.csv
# Runtime: ~2 hours

library(tidyverse)
library(rstan)
library(here)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(42)

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# 1. Load Data
# =============================================================================

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

# Base Stan data
stan_data_base <- list(
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
  e_RH_sd = e_RH_sd
)

# =============================================================================
# 2. Compile All Models
# =============================================================================

cat("Compiling models...\n")
model_eht <- stan_model(here("stan/sensitivity/tier2_eht_sensitivity.stan"))
model_h2ot <- stan_model(here("stan/sensitivity/tier2_h2ot_sensitivity.stan"))
model_ea <- stan_model(here("stan/sensitivity/tier2_ea_sensitivity.stan"))
model_rh <- stan_model(here("stan/sensitivity/tier2_rh_sensitivity.stan"))
cat("  All models compiled successfully\n\n")

# =============================================================================
# 3. Define Sensitivity Sweep Parameters
# =============================================================================

# EHT: 3.0C (our estimate) down to 0.1C
eht_values <- c(3.0, 1.0, 0.3, 0.1)

# H2Ot: ~0.08 wt% (from observed variation) down to 0.005 wt% (~5% CV at mean 0.1)
# Observed: range 0.05-0.43, SD ~0.08
h2ot_values <- c(0.08, 0.04, 0.02, 0.005)

# Ea: ~2000 J/mol (from observed variation) down to 100 J/mol (~0.1% CV)
# Observed: range 79642-85858, SD ~2000
ea_values <- c(2000, 1000, 500, 100)

# RH: 0.10 (realistic: reflects 2.1x rainfall gradient) down to 0.005 (implausible uniformity)
# Context: Rapa Nui rainfall 805-1,690 mm/year (2.1x range), but Stevenson used RH=0.98 for all
rh_values <- c(0.10, 0.05, 0.02, 0.005)

# =============================================================================
# 4. Helper Function to Run Sweep
# =============================================================================

run_sensitivity <- function(model, param_name, param_values, stan_data_base) {
  results_list <- list()

  for (i in seq_along(param_values)) {
    val <- param_values[i]
    cat(sprintf("  [%d/%d] %s = %.4f...\n", i, length(param_values), param_name, val))

    stan_data <- stan_data_base
    stan_data[[paste0(param_name, "_prior")]] <- val

    fit <- sampling(
      model,
      data = stan_data,
      chains = 4,
      iter = 2000,
      warmup = 1000,
      seed = 42,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      refresh = 0
    )

    summary_fit <- summary(fit)$summary
    age_params <- grep("^age\\[", rownames(summary_fit), value = TRUE)
    age_summary <- summary_fit[age_params, ]

    ci_widths <- age_summary[, "97.5%"] - age_summary[, "2.5%"]
    half_ci_widths <- ci_widths / 2

    results_list[[i]] <- tibble(
      parameter = param_name,
      prior_value = val,
      median_uncertainty = median(half_ci_widths),
      mean_uncertainty = mean(half_ci_widths),
      rhat_max = max(age_summary[, "Rhat"], na.rm = TRUE),
      ess_min = min(age_summary[, "n_eff"], na.rm = TRUE)
    )

    cat(sprintf("    Median uncertainty: +/-%.1f years\n", median(half_ci_widths)))

    rm(fit)
    gc()
  }

  bind_rows(results_list)
}

# =============================================================================
# 5. Run All Sweeps
# =============================================================================

cat("\nEHT sensitivity sweep\n")
results_eht <- run_sensitivity(model_eht, "EHT_sd", eht_values, stan_data_base)

cat("\nH2Ot sensitivity sweep\n")
results_h2ot <- run_sensitivity(model_h2ot, "H2Ot_sigma", h2ot_values, stan_data_base)

cat("\nEa sensitivity sweep\n")
results_ea <- run_sensitivity(model_ea, "Ea_sigma", ea_values, stan_data_base)

cat("\nRH sensitivity sweep\n")
results_rh <- run_sensitivity(model_rh, "RH_sd", rh_values, stan_data_base)

# =============================================================================
# 6. Combine Results
# =============================================================================

all_results <- bind_rows(results_eht, results_h2ot, results_ea, results_rh)

# Add normalized values and % change from baseline
all_results <- all_results %>%
  group_by(parameter) %>%
  mutate(
    baseline = first(median_uncertainty),
    pct_change = (median_uncertainty - baseline) / baseline * 100,
    fold_reduction = first(prior_value) / prior_value
  ) %>%
  ungroup()

# =============================================================================
# 7. Print Summary Table (for manuscript)
# =============================================================================

# Format for display
summary_table <- all_results %>%
  select(parameter, prior_value, fold_reduction, median_uncertainty, pct_change) %>%
  mutate(
    prior_display = case_when(
      parameter == "EHT_sd" ~ sprintf("%.1f C", prior_value),
      parameter == "H2Ot_sigma" ~ sprintf("%.3f wt%%", prior_value),
      parameter == "Ea_sigma" ~ sprintf("%.0f J/mol", prior_value),
      parameter == "RH_sd" ~ sprintf("%.3f", prior_value)
    ),
    fold_display = sprintf("%.0fx", fold_reduction),
    uncertainty_display = sprintf("+/-%.0f yr", median_uncertainty),
    change_display = sprintf("%.1f%%", pct_change)
  )

cat("\nParameter Sensitivity Summary:\n\n")

# EHT
cat("EHT (Temperature):\n")
cat("  Prior SD       | Tightening | Uncertainty | Change\n")
cat("  ---------------|------------|-------------|--------\n")
for (i in which(summary_table$parameter == "EHT_sd")) {
  cat(sprintf("  %-14s | %-10s | %-11s | %s\n",
              summary_table$prior_display[i],
              summary_table$fold_display[i],
              summary_table$uncertainty_display[i],
              summary_table$change_display[i]))
}

# H2Ot
cat("\nH2Ot (Structural Water):\n")
cat("  Prior SD       | Tightening | Uncertainty | Change\n")
cat("  ---------------|------------|-------------|--------\n")
for (i in which(summary_table$parameter == "H2Ot_sigma")) {
  cat(sprintf("  %-14s | %-10s | %-11s | %s\n",
              summary_table$prior_display[i],
              summary_table$fold_display[i],
              summary_table$uncertainty_display[i],
              summary_table$change_display[i]))
}

# Ea
cat("\nEa (Activation Energy):\n")
cat("  Prior SD       | Tightening | Uncertainty | Change\n")
cat("  ---------------|------------|-------------|--------\n")
for (i in which(summary_table$parameter == "Ea_sigma")) {
  cat(sprintf("  %-14s | %-10s | %-11s | %s\n",
              summary_table$prior_display[i],
              summary_table$fold_display[i],
              summary_table$uncertainty_display[i],
              summary_table$change_display[i]))
}

# RH
cat("\nRH (Relative Humidity):\n")
cat("  Prior SD       | Tightening | Uncertainty | Change\n")
cat("  ---------------|------------|-------------|--------\n")
for (i in which(summary_table$parameter == "RH_sd")) {
  cat(sprintf("  %-14s | %-10s | %-11s | %s\n",
              summary_table$prior_display[i],
              summary_table$fold_display[i],
              summary_table$uncertainty_display[i],
              summary_table$change_display[i]))
}

# =============================================================================
# 8. Key Finding
# =============================================================================

max_eht_reduction <- min(results_eht$pct_change)
max_h2ot_reduction <- min(results_h2ot$pct_change)
max_ea_reduction <- min(results_ea$pct_change)
max_rh_reduction <- min(results_rh$pct_change)

cat(sprintf("\nMaximum uncertainty reduction from tightening priors:\n"))
cat(sprintf("  EHT:  %.1f%% (from 30x tighter prior)\n", max_eht_reduction))
cat(sprintf("  H2Ot: %.1f%% (from 16x tighter prior)\n", max_h2ot_reduction))
cat(sprintf("  Ea:   %.1f%% (from 20x tighter prior)\n", max_ea_reduction))
cat(sprintf("  RH:   %.1f%% (from 20x tighter prior)\n", max_rh_reduction))

# =============================================================================
# 9. Save Results
# =============================================================================

write_csv(results_eht, here("output/tables/eht_sensitivity_results.csv"))
write_csv(all_results, here("output/tables/comprehensive_sensitivity_results.csv"))

# Create formatted table for manuscript
manuscript_table <- all_results %>%
  select(parameter, prior_value, fold_reduction, median_uncertainty, pct_change) %>%
  mutate(
    parameter = recode(parameter,
                       "EHT_sd" = "EHT (deg C)",
                       "H2Ot_sigma" = "H2Ot (wt%)",
                       "Ea_sigma" = "Ea (J/mol)",
                       "RH_sd" = "RH (fraction)")
  )
write_csv(manuscript_table, here("output/tables/manuscript_sensitivity_table.csv"))

cat("\nResults saved to output/tables/\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
