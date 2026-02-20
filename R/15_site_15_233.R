# 15_site_15_233.R - Site 15-233 phase model analysis
# Purpose: Fits 6 Bayesian phase models (C14/OHD/combined x ordered/single) for Site 15-233
# Inputs: data/site_15_233_ohd.csv, data/site_15_233_radiocarbon.csv,
#   stan/radiocarbon_ordered.stan, stan/radiocarbon_single_phase.stan,
#   stan/ohd_ordered.stan, stan/ohd_single_phase.stan,
#   stan/combined_ordered.stan, stan/combined_single_phase.stan
# Outputs: output/fits/site15_*.rds (6 fit objects),
#   output/tables/site15_model_comparison.csv
# Runtime: ~20 minutes

library(rstan)
library(loo)
library(dplyr)
library(ggplot2)
library(Bchron)
library(here)

# Stan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Configuration
N_ITER <- 50000     # Production run
N_WARMUP <- 25000   # Production run
N_CHAINS <- 4
SEED <- 42

# Create output directories
dir.create(here("output", "fits"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)

# Source utilities
source(here("R", "utils.R"))

# Load and Prepare Data
cat("Loading Site 15-233 data...\n")

# Load OHD data
ohd_raw <- read.csv(here("data", "site_15_233_ohd.csv"))
cat(sprintf("Loaded %d OHD samples\n", nrow(ohd_raw)))

# Load C14 data
c14_raw <- read.csv(here("data", "site_15_233_radiocarbon.csv"))
cat(sprintf("Loaded %d C14 samples\n", nrow(c14_raw)))

# Fix Phase Assignments
# Data file has L1 = Phase 1 (wrong)
# We need: Phase 1 = OLDEST (deepest = L6)
# Formula: Phase = (max_level + 1) - level_number

# Extract level number from Level column (e.g., "L3" -> 3)
ohd_raw$Level_Num <- as.integer(gsub("L", "", ohd_raw$Level))
c14_raw$Level_Num <- as.integer(gsub("L", "", c14_raw$Level))

# Calculate correct phases: Phase 1 = deepest (L6), Phase 6 = shallowest (L1)
max_level <- 6
ohd_raw$Phase_Correct <- (max_level + 1) - ohd_raw$Level_Num
c14_raw$Phase_Correct <- (max_level + 1) - c14_raw$Level_Num

# Validate phase assignments
if (any(is.na(ohd_raw$Phase_Correct))) {
  stop("ERROR: Some OHD samples have undefined phase assignments")
}
if (any(is.na(c14_raw$Phase_Correct))) {
  stop("ERROR: Some C14 samples have undefined phase assignments")
}

# Convert AD dates to BP (BP = 1950 - AD)
ohd_raw$OHD_BP <- 1950 - ohd_raw$OHD_Date_AD

cat("Phase assignment: L1 (shallowest) -> Phase 6, L6 (deepest) -> Phase 1\n")

# Show OHD summary by level
ohd_summary <- ohd_raw %>%
  group_by(Level, Phase_Correct) %>%
  summarise(
    N = n(),
    Mean_AD = mean(OHD_Date_AD),
    Mean_BP = mean(OHD_BP),
    SD_BP = sd(OHD_BP),
    .groups = "drop"
  )
print(ohd_summary)

# Show C14 summary
c14_summary <- c14_raw %>%
  select(Lab_ID, Level, Phase_Correct, C14_Age, C14_SD)
print(c14_summary)

# Radiocarbon Calibration
cat("Calibrating radiocarbon dates...\n")

# Define calibration grid parameters
CAL_AGE_MIN <- 0      # AD 1950 in BP (can't go negative in Stan)
CAL_AGE_MAX <- 600    # AD 1350 in BP
N_GRID <- 501

# Create calibration grid
cal_grid <- seq(CAL_AGE_MIN, CAL_AGE_MAX, length.out = N_GRID)

# Calibrate each C14 sample using Bchron
cal_ages_list <- list()
cal_probs_list <- list()

for (i in 1:nrow(c14_raw)) {
  # Bchron calibration
  cal_result <- BchronCalibrate(
    ages = c14_raw$C14_Age[i],
    ageSds = c14_raw$C14_SD[i],
    calCurves = "shcal20",  # Southern Hemisphere for Rapa Nui (27 S)
    allowOutside = TRUE     # L6 date (100 BP) is slightly outside shcal20 range (118-50227 BP)
  )

  # Extract the calibrated distribution
  cal_df <- as.data.frame(cal_result[[1]])

  # Interpolate to our standard grid
  cal_interp <- approx(
    x = cal_df$ageGrid,
    y = cal_df$densities,
    xout = cal_grid,
    rule = 2
  )

  # Normalize to sum to 1
  cal_probs <- cal_interp$y
  cal_probs[cal_probs < 0] <- 0
  cal_probs <- cal_probs / sum(cal_probs)

  cal_ages_list[[i]] <- cal_grid
  cal_probs_list[[i]] <- cal_probs
}

cat(sprintf("Calibration complete. Grid: %d points, range: %d to %d BP\n",
            N_GRID, CAL_AGE_MIN, CAL_AGE_MAX))

# Define Age Bounds
ohd_age_range <- range(ohd_raw$OHD_BP)
c14_age_range <- c(CAL_AGE_MIN, CAL_AGE_MAX)

# Combined range with padding
AGE_MIN <- max(0, min(ohd_age_range[1], c14_age_range[1]) - 50)
AGE_MAX <- max(ohd_age_range[2], c14_age_range[2]) + 50

cat(sprintf("Age bounds: %.0f to %.0f BP\n", AGE_MIN, AGE_MAX))

# Compile Stan Models
cat("Compiling Stan models...\n")

# Single-phase models
stan_model_c14_single <- stan_model(here("stan", "radiocarbon_single_phase.stan"))
stan_model_ohd_single <- stan_model(here("stan", "ohd_single_phase.stan"))
stan_model_combined_single <- stan_model(here("stan", "combined_single_phase.stan"))

# Hard-ordered phase models
stan_model_c14_ordered <- stan_model(here("stan", "radiocarbon_ordered.stan"))
stan_model_ohd_ordered <- stan_model(here("stan", "ohd_ordered.stan"))
stan_model_combined_ordered <- stan_model(here("stan", "combined_ordered.stan"))

cat("All models compiled successfully.\n")

# Prepare Stan Data
N_PHASES <- 6

# OHD ordered data
stan_data_ohd_ordered <- list(
  N = nrow(ohd_raw),
  N_phases = N_PHASES,
  phase_id = ohd_raw$Phase_Correct,
  OHD_age_BP = ohd_raw$OHD_BP,
  OHD_SD = ohd_raw$OHD_SD,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# OHD single-phase data
stan_data_ohd_single <- list(
  N = nrow(ohd_raw),
  OHD_age_BP = ohd_raw$OHD_BP,
  OHD_SD = ohd_raw$OHD_SD,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# C14 ordered data
stan_data_c14_ordered <- list(
  N = nrow(c14_raw),
  N_phases = N_PHASES,
  phase_id = c14_raw$Phase_Correct,
  N_grid = N_GRID,
  cal_ages = cal_ages_list,
  cal_probs = cal_probs_list,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# C14 single-phase data
stan_data_c14_single <- list(
  N = nrow(c14_raw),
  N_grid = N_GRID,
  cal_ages = cal_ages_list,
  cal_probs = cal_probs_list,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# Combined ordered data
stan_data_combined_ordered <- list(
  N_c14 = nrow(c14_raw),
  N_ohd = nrow(ohd_raw),
  N_phases = N_PHASES,
  phase_id_c14 = c14_raw$Phase_Correct,
  phase_id_ohd = ohd_raw$Phase_Correct,
  N_grid = N_GRID,
  cal_ages = cal_ages_list,
  cal_probs = cal_probs_list,
  OHD_age_BP = ohd_raw$OHD_BP,
  OHD_SD = ohd_raw$OHD_SD,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# Combined single-phase data
stan_data_combined_single <- list(
  N_c14 = nrow(c14_raw),
  N_ohd = nrow(ohd_raw),
  N_grid = N_GRID,
  cal_ages = cal_ages_list,
  cal_probs = cal_probs_list,
  OHD_age_BP = ohd_raw$OHD_BP,
  OHD_SD = ohd_raw$OHD_SD,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

cat("Stan data prepared for all 6 models.\n")

# Fit Models
fit_model <- function(model, data, name) {
  cat(sprintf("Fitting Model %s...\n", name))
  start_time <- Sys.time()

  fit <- sampling(
    model,
    data = data,
    iter = N_ITER,
    warmup = N_WARMUP,
    chains = N_CHAINS,
    seed = SEED,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Completed in %.1f minutes\n", as.numeric(elapsed)))

  # Check convergence
  summary_df <- as.data.frame(summary(fit)$summary)
  max_rhat <- max(summary_df$Rhat, na.rm = TRUE)
  min_ess <- min(summary_df$n_eff, na.rm = TRUE)

  cat(sprintf("  Max Rhat: %.3f, Min ESS: %.0f\n", max_rhat, min_ess))

  # HMC diagnostics
  n_div <- sum(rstan::get_num_divergent(fit))
  n_tree <- sum(rstan::get_num_max_treedepth(fit))
  cat(sprintf("  Divergent transitions: %d\n", n_div))
  cat(sprintf("  Max treedepth exceedances: %d\n", n_tree))

  if (max_rhat > 1.05) {
    warning(sprintf("Model %s: High Rhat (%.3f) - may not have converged", name, max_rhat))
  }
  if (n_div > 0) warning(sprintf("Model %s: %d DIVERGENT TRANSITIONS", name, n_div))
  if (n_tree > 0) warning(sprintf("Model %s: %d MAX TREEDEPTH EXCEEDANCES", name, n_tree))

  fit
}

cat("\nFitting all 6 models...\n")

fit_1A <- fit_model(stan_model_c14_ordered, stan_data_c14_ordered, "1A (C14 Ordered)")
fit_1B <- fit_model(stan_model_c14_single, stan_data_c14_single, "1B (C14 Single)")
fit_1C <- fit_model(stan_model_ohd_ordered, stan_data_ohd_ordered, "1C (OHD Ordered)")
fit_1D <- fit_model(stan_model_ohd_single, stan_data_ohd_single, "1D (OHD Single)")
fit_1E <- fit_model(stan_model_combined_ordered, stan_data_combined_ordered, "1E (Combined Ordered)")
fit_1F <- fit_model(stan_model_combined_single, stan_data_combined_single, "1F (Combined Single)")

# Save fits
saveRDS(fit_1A, here("output", "fits", "site15_1A_c14_ordered.rds"))
saveRDS(fit_1B, here("output", "fits", "site15_1B_c14_single.rds"))
saveRDS(fit_1C, here("output", "fits", "site15_1C_ohd_ordered.rds"))
saveRDS(fit_1D, here("output", "fits", "site15_1D_ohd_single.rds"))
saveRDS(fit_1E, here("output", "fits", "site15_1E_combined_ordered.rds"))
saveRDS(fit_1F, here("output", "fits", "site15_1F_combined_single.rds"))

cat("All fits saved.\n")

# LOO-CV Model Comparison
cat("Computing LOO-CV...\n")

compute_loo <- function(fit, name) {
  log_lik <- extract_log_lik(fit, merge_chains = FALSE)
  r_eff <- relative_eff(exp(log_lik))
  loo_result <- loo(log_lik, r_eff = r_eff)

  cat(sprintf("\n%s:\n", name))
  print(loo_result)

  loo_result
}

loo_1A <- compute_loo(fit_1A, "Model 1A (C14 Ordered)")
loo_1B <- compute_loo(fit_1B, "Model 1B (C14 Single)")
loo_1C <- compute_loo(fit_1C, "Model 1C (OHD Ordered)")
loo_1D <- compute_loo(fit_1D, "Model 1D (OHD Single)")
loo_1E <- compute_loo(fit_1E, "Model 1E (Combined Ordered)")
loo_1F <- compute_loo(fit_1F, "Model 1F (Combined Single)")

# KEY COMPARISONS

# Comparison 1: C14 Ordered vs Single
cat("\nComparison 1A vs 1B: Does phase structure help C14?\n")
comp_c14 <- loo_compare(loo_1A, loo_1B)
print(comp_c14)

# Comparison 2: OHD Ordered vs Single (KEY TEST)
cat("\n*** KEY TEST: Comparison 1C vs 1D: Does phase structure help OHD? ***\n")
comp_ohd <- loo_compare(loo_1C, loo_1D)
print(comp_ohd)

# Interpret
elpd_diff <- comp_ohd[2, "elpd_diff"]
se_diff <- comp_ohd[2, "se_diff"]
z_score <- elpd_diff / se_diff

cat(sprintf("  ELPD difference: %.2f (SE: %.2f)\n", elpd_diff, se_diff))
cat(sprintf("  Z-score: %.2f\n", z_score))

if (abs(z_score) < 2) {
  cat("  Result: INDISTINGUISHABLE - OHD does not prefer phases over single population\n")
} else if (z_score > 2) {
  cat("  Result: Ordered model BETTER - OHD fits phase structure\n")
} else {
  cat("  Result: Single-phase model BETTER - OHD conflicts with stratigraphic ordering\n")
}

# Comparison 3: Combined Ordered vs Single
cat("\nComparison 1E vs 1F: Does phase structure help combined?\n")
comp_combined <- loo_compare(loo_1E, loo_1F)
print(comp_combined)

# Stacking Weights
cat("\nStacking weights:\n")

# C14 models
cat("C14 Models:\n")
stack_c14 <- loo_model_weights(list(loo_1A, loo_1B), method = "stacking")
names(stack_c14) <- c("1A (Ordered)", "1B (Single)")
print(round(stack_c14, 3))

# OHD models
cat("OHD Models (KEY):\n")
stack_ohd <- loo_model_weights(list(loo_1C, loo_1D), method = "stacking")
names(stack_ohd) <- c("1C (Ordered)", "1D (Single)")
print(round(stack_ohd, 3))

# Combined models
cat("Combined Models:\n")
stack_combined <- loo_model_weights(list(loo_1E, loo_1F), method = "stacking")
names(stack_combined) <- c("1E (Ordered)", "1F (Single)")
print(round(stack_combined, 3))

# Phase Boundary Estimates (Ordered Models)
cat("\nPhase boundary estimates:\n")

extract_boundaries <- function(fit, model_name) {
  posterior <- as.data.frame(fit)
  boundary_cols <- grep("^phase_boundary\\[", colnames(posterior), value = TRUE)

  if (length(boundary_cols) == 0) {
    cat(sprintf("%s: No phase boundaries found\n", model_name))
    return(NULL)
  }

  cat(sprintf("%s Phase Boundaries:\n", model_name))
  for (col in boundary_cols) {
    vals <- posterior[[col]]
    cat(sprintf("  %s: %.1f BP (95%% CI: %.1f - %.1f)\n",
                col,
                median(vals),
                quantile(vals, 0.025),
                quantile(vals, 0.975)))
  }
  cat("\n")
}

extract_boundaries(fit_1A, "Model 1A (C14 Ordered)")
extract_boundaries(fit_1C, "Model 1C (OHD Ordered)")
extract_boundaries(fit_1E, "Model 1E (Combined Ordered)")

# Summary Table
summary_table <- data.frame(
  Model = c("1A (C14 Ordered)", "1B (C14 Single)",
            "1C (OHD Ordered)", "1D (OHD Single)",
            "1E (Combined Ordered)", "1F (Combined Single)"),
  Data = c("C14", "C14", "OHD", "OHD", "Combined", "Combined"),
  Structure = c("Ordered", "Single", "Ordered", "Single", "Ordered", "Single"),
  N = c(nrow(c14_raw), nrow(c14_raw),
        nrow(ohd_raw), nrow(ohd_raw),
        nrow(c14_raw) + nrow(ohd_raw), nrow(c14_raw) + nrow(ohd_raw)),
  ELPD = c(loo_1A$estimates["elpd_loo", "Estimate"],
           loo_1B$estimates["elpd_loo", "Estimate"],
           loo_1C$estimates["elpd_loo", "Estimate"],
           loo_1D$estimates["elpd_loo", "Estimate"],
           loo_1E$estimates["elpd_loo", "Estimate"],
           loo_1F$estimates["elpd_loo", "Estimate"]),
  ELPD_SE = c(loo_1A$estimates["elpd_loo", "SE"],
              loo_1B$estimates["elpd_loo", "SE"],
              loo_1C$estimates["elpd_loo", "SE"],
              loo_1D$estimates["elpd_loo", "SE"],
              loo_1E$estimates["elpd_loo", "SE"],
              loo_1F$estimates["elpd_loo", "SE"])
)

print(summary_table)

# Save summary
write.csv(summary_table, here("output", "tables", "site15_model_comparison.csv"), row.names = FALSE)

cat("\nSite 15-233 analysis complete.\n")
