# 04_tier3_hybrid.R - Spatial vs. temporal variance decomposition
# Purpose: Compares within-site vs. between-site uncertainty components
# Inputs: output/fits/tier3_fit.rds
# Outputs: output/tables/variance_decomposition.csv, output/tables/variance_decomposition_summary.csv
# Runtime: ~5 minutes

library(tidyverse)
library(rstan)
library(here)

# ===== STEP 1: LOAD TIER 3 POSTERIOR =====
cat("Step 1: Loading Tier 3 posterior\n")

# Load Tier 3 Stan fit
tier3_fit <- readRDS(here("output/fits/tier3_fit.rds"))

# Extract posterior draws (1000 draws per parameter)
cat("Extracting posterior draws...\n")
n_draws <- 1000

# Extract all relevant parameters
# Note: With no-source model, H2Ot and Ea are per-sample (island-wide population)
posterior <- rstan::extract(tier3_fit,
                            pars = c("age", "EHT", "RH", "H2Ot", "Ea", "e_RH",
                                    "H2Ot_mu", "H2Ot_sigma", "Ea_mu", "Ea_sigma"),
                            permuted = TRUE)

# Sample 1000 draws (if we have more)
if (nrow(posterior$age) > n_draws) {
  draw_indices <- sample(1:nrow(posterior$age), n_draws, replace = FALSE)
  age_posterior <- posterior$age[draw_indices, ]
  EHT_posterior <- posterior$EHT[draw_indices, ]
  RH_posterior <- posterior$RH[draw_indices, ]
  H2Ot_posterior <- posterior$H2Ot[draw_indices, ]  # Per-sample (island-wide draws)
  Ea_posterior <- posterior$Ea[draw_indices, ]      # Per-sample (island-wide draws)
  e_RH_posterior <- posterior$e_RH[draw_indices]
} else {
  age_posterior <- posterior$age
  EHT_posterior <- posterior$EHT
  RH_posterior <- posterior$RH
  H2Ot_posterior <- posterior$H2Ot
  Ea_posterior <- posterior$Ea
  e_RH_posterior <- posterior$e_RH
  n_draws <- nrow(age_posterior)
}

N_samples <- ncol(age_posterior)

cat(sprintf("Extracted %d posterior draws for %d samples\n\n", n_draws, N_samples))

# ===== STEP 2: GENERATE TEMPERATURE SCENARIOS =====
cat("Step 2: Generating temperature scenarios\n")

# Rate equation constants (from RATE_EQUATION_FINAL.md)
a_eq <- 23.8042281295
b_eq <- 0.3782846434
c_eq <- 0.0002071707
d_eq <- -6313.2075128583
R_gas <- 8.314

# Burial depth parameters
burial_depth_range <- c(0.10, 0.70)  # 10-70 cm (screened matrix)
damping_depth <- 0.30  # 30 cm thermal diffusivity

# Temperature scenario generation function
generate_temperature_scenarios <- function(age_years, EHT_mean, burial_depth) {

  # Time vector (1-year steps)
  time_steps <- seq(0, age_years, by = 1)
  n_steps <- length(time_steps)

  # Burial depth thermal damping
  damping_factor <- exp(-burial_depth / damping_depth)

  # Scenario 1: BASELINE (stationary climate, no temporal variation)
  T_baseline <- rep(EHT_mean, n_steps)

  # Scenario 2: ENSO CYCLES ONLY (+/-1 C amplitude, 7-year periodicity, depth-damped)
  ENSO_amplitude <- 1.0 * damping_factor  # Depth-averaged
  ENSO_period <- 7  # years
  T_ENSO <- EHT_mean + ENSO_amplitude * sin(2 * pi * time_steps / ENSO_period)

  # Scenario 3: YEAR-TO-YEAR VARIATION ONLY (AR(1) process, rho=0.7, sd=0.5 C)
  # AR(1): T[t] = rho * T[t-1] + epsilon[t]
  rho <- 0.7
  sigma_epsilon <- 0.5 * damping_factor * sqrt(1 - rho^2)  # Stationary variance
  T_year_to_year <- numeric(n_steps)
  T_year_to_year[1] <- rnorm(1, 0, 0.5 * damping_factor)  # Initial value
  for (t in 2:n_steps) {
    T_year_to_year[t] <- rho * T_year_to_year[t-1] + rnorm(1, 0, sigma_epsilon)
  }
  T_year_to_year <- EHT_mean + T_year_to_year

  # Scenario 4: COMBINED REALISTIC (ENSO + year-to-year)
  T_combined <- EHT_mean +
    ENSO_amplitude * sin(2 * pi * time_steps / ENSO_period) +
    (T_year_to_year - EHT_mean)  # Add year-to-year deviation

  # Return all scenarios
  list(
    baseline = T_baseline,
    ENSO = T_ENSO,
    year_to_year = T_year_to_year,
    combined = T_combined,
    time_steps = time_steps
  )
}

cat("  Baseline: Stationary climate (EHT constant)\n")
cat("  ENSO: +/-1 C amplitude, 7-year period, depth-damped\n")
cat("  Year-to-year: AR(1) rho=0.7, sigma=0.5 C, depth-damped\n")
cat("  Combined: ENSO + year-to-year (most realistic)\n\n")

# ===== STEP 3: INTEGRATE HYDRATION RATE THROUGH TIME =====
cat("Step 3: Integrating hydration rates through time\n")

# Hydration rate calculation function WITH RH term
# Extended equation: log(rate) = a + b*(-Ea/(R*T)) + c*H2Ot + d*(1/T) + e_RH*(RH - 0.98)
calculate_rate <- function(T_celsius, H2Ot, Ea, RH, e_RH) {
  T_kelvin <- T_celsius + 273.15
  log_rate <- a_eq + b_eq * (-Ea / (R_gas * T_kelvin)) +
              c_eq * H2Ot + d_eq * (1.0 / T_kelvin) +
              e_RH * (RH - 0.98)  # RH term

  # Numerical stability
  log_rate <- pmin(pmax(log_rate, -50), 50)

  rate <- exp(log_rate)  # rate in (wt%)^2/year
  return(rate)
}

# Integrate hydration rate to calculate expected H2Ome
integrate_hydration <- function(T_scenario, H2Ot, Ea, RH, e_RH) {
  # T_scenario is a vector of temperatures over time
  # Returns sqrt(integral of rate(t) dt)

  rates <- calculate_rate(T_scenario, H2Ot, Ea, RH, e_RH)

  # Trapezoidal integration (1-year steps)
  integrated_rate <- sum(rates[-1] + rates[-length(rates)]) / 2

  # H2Ome = sqrt(integrated_rate)
  H2Ome_expected <- sqrt(integrated_rate)

  return(H2Ome_expected)
}

# ===== STEP 4: LOAD EMPIRICAL DATA =====
cat("Step 4: Loading empirical data\n")

data_2015 <- read_csv(here("data/stevenson_2015.csv"),
                      show_col_types = FALSE)

cat(sprintf("Loaded %d samples from Stevenson et al. (2015)\n\n", nrow(data_2015)))

# Extract study area assignments
study_area <- data_2015$study_area

# ===== STEP 5: RUN HYBRID BAYESIAN-MC ANALYSIS =====
cat("Step 5: Running hybrid Bayesian-MC analysis\n")
cat(sprintf("  %d samples x %d draws x 4 scenarios\n\n", N_samples, n_draws))

# Initialize result arrays
# Dimensions: [sample, draw, scenario]
scenarios <- c("baseline", "ENSO", "year_to_year", "combined")
H2Ome_temporal <- array(NA, dim = c(N_samples, n_draws, 4),
                        dimnames = list(NULL, NULL, scenarios))

# Progress tracking
start_time <- Sys.time()
pb <- txtProgressBar(min = 0, max = N_samples, style = 3)

for (i in 1:N_samples) {

  # Random burial depth for this sample (10-70 cm range)
  burial_depth <- runif(1, burial_depth_range[1], burial_depth_range[2])

  for (j in 1:n_draws) {

    # Extract posterior values for this draw
    age_j <- age_posterior[j, i]
    EHT_j <- EHT_posterior[j, i]
    RH_j <- RH_posterior[j, i]
    H2Ot_j <- H2Ot_posterior[j, i]  # Per-sample (island-wide, not source-specific)
    Ea_j <- Ea_posterior[j, i]      # Per-sample (island-wide, not source-specific)
    e_RH_j <- e_RH_posterior[j]

    # Generate temperature scenarios for this age/EHT/depth combination
    scenarios_j <- generate_temperature_scenarios(age_j, EHT_j, burial_depth)

    # Integrate hydration rate for each scenario (WITH RH)
    H2Ome_temporal[i, j, "baseline"] <- integrate_hydration(
      scenarios_j$baseline, H2Ot_j, Ea_j, RH_j, e_RH_j
    )

    H2Ome_temporal[i, j, "ENSO"] <- integrate_hydration(
      scenarios_j$ENSO, H2Ot_j, Ea_j, RH_j, e_RH_j
    )

    H2Ome_temporal[i, j, "year_to_year"] <- integrate_hydration(
      scenarios_j$year_to_year, H2Ot_j, Ea_j, RH_j, e_RH_j
    )

    H2Ome_temporal[i, j, "combined"] <- integrate_hydration(
      scenarios_j$combined, H2Ot_j, Ea_j, RH_j, e_RH_j
    )
  }

  setTxtProgressBar(pb, i)
}

close(pb)
end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat(sprintf("\nHybrid analysis complete in %.1f minutes\n\n", runtime))

# ===== STEP 6: CALCULATE VARIANCE DECOMPOSITION =====
cat("Step 6: Calculating variance decomposition\n")

# For each sample, calculate variance across draws and scenarios
variance_results <- data.frame(
  sample_id = 1:N_samples,
  study_area = study_area,
  var_baseline = numeric(N_samples),
  var_ENSO = numeric(N_samples),
  var_year_to_year = numeric(N_samples),
  var_combined = numeric(N_samples),
  var_spatial = numeric(N_samples),  # Tier 3 baseline
  var_temporal = numeric(N_samples)  # Combined - baseline
)

for (i in 1:N_samples) {
  # Variance for each scenario (across draws)
  variance_results$var_baseline[i] <- var(H2Ome_temporal[i, , "baseline"])
  variance_results$var_ENSO[i] <- var(H2Ome_temporal[i, , "ENSO"])
  variance_results$var_year_to_year[i] <- var(H2Ome_temporal[i, , "year_to_year"])
  variance_results$var_combined[i] <- var(H2Ome_temporal[i, , "combined"])

  # Spatial variance (from Tier 3 baseline scenario)
  variance_results$var_spatial[i] <- variance_results$var_baseline[i]

  # Temporal variance (difference between combined and baseline)
  variance_results$var_temporal[i] <- variance_results$var_combined[i] -
                                       variance_results$var_baseline[i]
}

# ===== STEP 7: SUMMARY STATISTICS =====
cat("Step 7: Summary statistics\n\n")

# Median variance contributions
median_var_spatial <- median(variance_results$var_spatial)
median_var_temporal <- median(variance_results$var_temporal)
median_var_total <- median(variance_results$var_combined)

# Handle negative temporal variance (can occur due to MC noise)
median_var_temporal <- max(median_var_temporal, 0)

pct_spatial <- 100 * median_var_spatial / median_var_total
pct_temporal <- 100 * median_var_temporal / median_var_total

cat("Variance Decomposition:\n")
cat(sprintf("  Median spatial variance:  %.6f (%.2f%%)\n",
            median_var_spatial, pct_spatial))
cat(sprintf("  Median temporal variance: %.6f (%.2f%%)\n",
            median_var_temporal, pct_temporal))
cat(sprintf("  Median total variance:    %.6f (100.0%%)\n\n", median_var_total))

if (pct_spatial > 95) {
  cat(sprintf("  Spatial uncertainty dominates at %.1f%%\n", pct_spatial))
  cat(sprintf("  Temporal climate contributes only %.2f%%\n\n", pct_temporal))
}

# ===== STEP 8: SAVE RESULTS =====
cat("Saving results\n")

dir.create(here("output/fits"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/figures"), recursive = TRUE, showWarnings = FALSE)

# Save H2Ome temporal array
saveRDS(H2Ome_temporal, here("output/fits/H2Ome_temporal.rds"))

# Save variance decomposition
write_csv(variance_results, here("output/tables/variance_decomposition.csv"))

# Save summary statistics
summary_stats <- data.frame(
  metric = c("median_var_spatial", "median_var_temporal", "median_var_total",
             "pct_spatial", "pct_temporal", "runtime_minutes"),
  value = c(median_var_spatial, median_var_temporal, median_var_total,
            pct_spatial, pct_temporal, runtime)
)
write_csv(summary_stats, here("output/tables/variance_decomposition_summary.csv"))

# ===== STEP 9: DIAGNOSTIC PLOTS =====
cat("Creating diagnostic plots...\n")

pdf(here("output/figures/tier3_hybrid_diagnostics.pdf"),
    width = 12, height = 8)

# 1. Variance decomposition scatter plot
p_variance <- ggplot(variance_results, aes(x = var_spatial, y = var_temporal, color = study_area)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Spatial vs Temporal Variance Contributions",
    subtitle = sprintf("Median: %.2f%% spatial, %.2f%% temporal", pct_spatial, pct_temporal),
    x = "Spatial Variance (from Tier 3 Bayesian)",
    y = "Temporal Variance (from MC climate scenarios)"
  ) +
  theme_minimal()
print(p_variance)

# 2. Percent contribution histogram
pct_temporal_samples <- 100 * pmax(variance_results$var_temporal, 0) /
                         (variance_results$var_spatial + pmax(variance_results$var_temporal, 0))

p_pct <- ggplot(data.frame(pct = pct_temporal_samples), aes(x = pct)) +
  geom_histogram(bins = 30, fill = "forestgreen", alpha = 0.7) +
  geom_vline(xintercept = median(pct_temporal_samples),
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Distribution of Temporal Variance Contribution",
    subtitle = sprintf("Median = %.2f%%", median(pct_temporal_samples)),
    x = "Temporal Contribution (%)",
    y = "Count"
  ) +
  theme_minimal()
print(p_pct)

# 3. Variance by study area
variance_long <- variance_results %>%
  select(sample_id, study_area, var_spatial, var_temporal) %>%
  pivot_longer(cols = c(var_spatial, var_temporal),
               names_to = "component", values_to = "variance") %>%
  mutate(component = factor(component,
                           levels = c("var_spatial", "var_temporal"),
                           labels = c("Spatial", "Temporal")))

p_by_area <- ggplot(variance_long, aes(x = study_area, y = variance, fill = component)) +
  geom_boxplot(alpha = 0.7) +
  scale_y_log10() +
  labs(
    title = "Variance Components by Study Area",
    subtitle = "Log scale",
    x = "Study Area",
    y = "Variance (log scale)",
    fill = "Component"
  ) +
  scale_fill_manual(values = c("Spatial" = "steelblue", "Temporal" = "forestgreen")) +
  theme_minimal()
print(p_by_area)

# 4. Scenario comparison
scenario_means <- data.frame(
  scenario = scenarios,
  mean_H2Ome = c(mean(H2Ome_temporal[, , "baseline"]),
                 mean(H2Ome_temporal[, , "ENSO"]),
                 mean(H2Ome_temporal[, , "year_to_year"]),
                 mean(H2Ome_temporal[, , "combined"])),
  sd_H2Ome = c(sd(H2Ome_temporal[, , "baseline"]),
               sd(H2Ome_temporal[, , "ENSO"]),
               sd(H2Ome_temporal[, , "year_to_year"]),
               sd(H2Ome_temporal[, , "combined"]))
)

p_scenarios <- ggplot(scenario_means, aes(x = scenario, y = sd_H2Ome)) +
  geom_bar(stat = "identity", fill = "forestgreen", alpha = 0.7) +
  labs(
    title = "H2Ome Variability by Temperature Scenario",
    subtitle = "Standard deviation of predicted H2Ome across all samples and draws",
    x = "Temperature Scenario",
    y = "SD of H2Ome (wt%)"
  ) +
  theme_minimal()
print(p_scenarios)

# 5. Pie chart of variance decomposition
var_data <- data.frame(
  component = c("Spatial\n(EHT + RH + composition)", "Temporal\n(climate variability)"),
  value = c(pct_spatial, pct_temporal)
)

p_pie <- ggplot(var_data, aes(x = "", y = value, fill = component)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(
    title = "Tier 3 Hybrid: Variance Decomposition",
    subtitle = "Spatial vs Temporal uncertainty attribution",
    fill = "Component"
  ) +
  scale_fill_manual(values = c("Spatial\n(EHT + RH + composition)" = "steelblue",
                               "Temporal\n(climate variability)" = "forestgreen")) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())
print(p_pie)

dev.off()

# Also save key plot as standalone PDF
ggsave(here("output/figures/variance_spatial_vs_temporal.pdf"),
       p_variance, width = 8, height = 6)

cat(sprintf("Tier 3 hybrid complete. Spatial: %.1f%%, Temporal: %.2f%%\n", pct_spatial, pct_temporal))
