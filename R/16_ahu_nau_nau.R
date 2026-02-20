# 16_ahu_nau_nau.R - Ahu Nau Nau phase model analysis
# Purpose: Fits 6 Bayesian phase models (C14/OHD/combined x ordered/single) for Ahu Nau Nau
# Inputs: data/ahu_nau_nau_ohd.csv, data/ahu_nau_nau_radiocarbon.csv,
#   stan/radiocarbon_ordered.stan, stan/radiocarbon_single_phase.stan,
#   stan/ohd_ordered.stan, stan/ohd_single_phase.stan,
#   stan/combined_ordered.stan, stan/combined_single_phase.stan
# Outputs: output/fits/ann_*.rds (6 fit objects),
#   output/tables/ann_model_comparison.csv
# Runtime: ~15 minutes

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
cat("Loading Ahu Nau Nau data...\n")

# Load OHD data
ohd_raw <- read.csv(here("data", "ahu_nau_nau_ohd.csv"))
cat(sprintf("Loaded %d OHD samples (raw)\n", nrow(ohd_raw)))

# Load C14 data
c14_raw <- read.csv(here("data", "ahu_nau_nau_radiocarbon.csv"))
cat(sprintf("Loaded %d C14 samples (raw)\n", nrow(c14_raw)))

# Apply Exclusions (Following Stevenson's Own Assessment)
cat("Applying exclusions per Stevenson's assessment:\n")

# Exclude Level 1 OHD (UNRELIABLE - wind-blown sand, includes 2002 AD)
ohd_level1 <- ohd_raw[ohd_raw$Level == 1, ]
cat(sprintf("  Excluding %d Level 1 OHD samples (wind-blown sand, UNRELIABLE)\n", nrow(ohd_level1)))
cat(sprintf("    Level 1 range: %d - %d AD (includes impossible 2002 AD!)\n",
            min(ohd_level1$OHD_Date_AD), max(ohd_level1$OHD_Date_AD)))

# Exclude intrusive OHD samples
ohd_intrusive <- ohd_raw[ohd_raw$Intrusive == TRUE, ]
cat(sprintf("  Excluding %d 'intrusive' OHD samples per Stevenson\n", nrow(ohd_intrusive)))

# Apply OHD exclusions
ohd_data <- ohd_raw %>%
  filter(Level != 1, Intrusive == FALSE)
cat(sprintf("  OHD after exclusions: %d samples\n", nrow(ohd_data)))

# Exclude rat bone C14 samples (marine reservoir effect)
c14_rat <- c14_raw[c14_raw$Rat_Bone == TRUE, ]
cat(sprintf("  Excluding %d rat bone C14 samples (marine reservoir effect)\n", nrow(c14_rat)))

c14_data <- c14_raw %>%
  filter(Rat_Bone == FALSE)
cat(sprintf("  C14 after exclusions: %d samples\n", nrow(c14_data)))

# Calculate Phase Assignments
# Levels used: 2, 4, 5, 6 -> Map to phases 1-4
# Phase 1 = OLDEST (deepest = L6)
# Phase 4 = YOUNGEST (shallowest = L2)

level_order <- c(6, 5, 4, 2)  # deepest to shallowest
phase_map <- setNames(1:4, level_order)

ohd_data$Phase <- phase_map[as.character(ohd_data$Level)]
c14_data$Phase <- phase_map[as.character(c14_data$Level)]

# Validate phase assignments
if (any(is.na(ohd_data$Phase))) {
  stop("ERROR: Some OHD samples have undefined phase assignments")
}
if (any(is.na(c14_data$Phase))) {
  stop("ERROR: Some C14 samples have undefined phase assignments")
}

cat("Phase assignment: L6 -> Phase 1 (oldest), L5 -> Phase 2, L4 -> Phase 3, L2 -> Phase 4 (youngest)\n")

# Convert AD dates to BP (BP = 1950 - AD)
ohd_data$OHD_BP <- 1950 - ohd_data$OHD_Date_AD

# Show OHD summary by level
ohd_summary <- ohd_data %>%
  group_by(Level, Phase) %>%
  summarise(
    N = n(),
    Mean_AD = mean(OHD_Date_AD),
    Mean_BP = mean(OHD_BP),
    SD_BP = sd(OHD_BP),
    .groups = "drop"
  )
print(ohd_summary)

# Show C14 summary
c14_summary <- c14_data %>%
  select(Lab_ID, Level, Phase, C14_Age, C14_SD, Material)
print(c14_summary)

# Radiocarbon Calibration
cat("Calibrating radiocarbon dates...\n")

# Define calibration grid parameters (wider range for older site)
CAL_AGE_MIN <- 200    # AD 1750 in BP
CAL_AGE_MAX <- 900    # AD 1050 in BP
N_GRID <- 501

# Create calibration grid
cal_grid <- seq(CAL_AGE_MIN, CAL_AGE_MAX, length.out = N_GRID)

# Calibrate each C14 sample using Bchron
cal_ages_list <- list()
cal_probs_list <- list()

for (i in 1:nrow(c14_data)) {
  # Bchron calibration
  cal_result <- BchronCalibrate(
    ages = c14_data$C14_Age[i],
    ageSds = c14_data$C14_SD[i],
    calCurves = "shcal20"  # Southern Hemisphere for Rapa Nui (27 S)
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
ohd_age_range <- range(ohd_data$OHD_BP)
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
N_PHASES <- 4

# OHD ordered data
stan_data_ohd_ordered <- list(
  N = nrow(ohd_data),
  N_phases = N_PHASES,
  phase_id = ohd_data$Phase,
  OHD_age_BP = ohd_data$OHD_BP,
  OHD_SD = ohd_data$OHD_SD,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# OHD single-phase data
stan_data_ohd_single <- list(
  N = nrow(ohd_data),
  OHD_age_BP = ohd_data$OHD_BP,
  OHD_SD = ohd_data$OHD_SD,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# C14 ordered data
stan_data_c14_ordered <- list(
  N = nrow(c14_data),
  N_phases = N_PHASES,
  phase_id = c14_data$Phase,
  N_grid = N_GRID,
  cal_ages = cal_ages_list,
  cal_probs = cal_probs_list,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# C14 single-phase data
stan_data_c14_single <- list(
  N = nrow(c14_data),
  N_grid = N_GRID,
  cal_ages = cal_ages_list,
  cal_probs = cal_probs_list,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# Combined ordered data
stan_data_combined_ordered <- list(
  N_c14 = nrow(c14_data),
  N_ohd = nrow(ohd_data),
  N_phases = N_PHASES,
  phase_id_c14 = c14_data$Phase,
  phase_id_ohd = ohd_data$Phase,
  N_grid = N_GRID,
  cal_ages = cal_ages_list,
  cal_probs = cal_probs_list,
  OHD_age_BP = ohd_data$OHD_BP,
  OHD_SD = ohd_data$OHD_SD,
  age_min = AGE_MIN,
  age_max = AGE_MAX
)

# Combined single-phase data
stan_data_combined_single <- list(
  N_c14 = nrow(c14_data),
  N_ohd = nrow(ohd_data),
  N_grid = N_GRID,
  cal_ages = cal_ages_list,
  cal_probs = cal_probs_list,
  OHD_age_BP = ohd_data$OHD_BP,
  OHD_SD = ohd_data$OHD_SD,
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

fit_2A <- fit_model(stan_model_c14_ordered, stan_data_c14_ordered, "2A (C14 Ordered)")
fit_2B <- fit_model(stan_model_c14_single, stan_data_c14_single, "2B (C14 Single)")
fit_2C <- fit_model(stan_model_ohd_ordered, stan_data_ohd_ordered, "2C (OHD Ordered)")
fit_2D <- fit_model(stan_model_ohd_single, stan_data_ohd_single, "2D (OHD Single)")
fit_2E <- fit_model(stan_model_combined_ordered, stan_data_combined_ordered, "2E (Combined Ordered)")
fit_2F <- fit_model(stan_model_combined_single, stan_data_combined_single, "2F (Combined Single)")

# Save fits
saveRDS(fit_2A, here("output", "fits", "ann_2A_c14_ordered.rds"))
saveRDS(fit_2B, here("output", "fits", "ann_2B_c14_single.rds"))
saveRDS(fit_2C, here("output", "fits", "ann_2C_ohd_ordered.rds"))
saveRDS(fit_2D, here("output", "fits", "ann_2D_ohd_single.rds"))
saveRDS(fit_2E, here("output", "fits", "ann_2E_combined_ordered.rds"))
saveRDS(fit_2F, here("output", "fits", "ann_2F_combined_single.rds"))

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

loo_2A <- compute_loo(fit_2A, "Model 2A (C14 Ordered)")
loo_2B <- compute_loo(fit_2B, "Model 2B (C14 Single)")
loo_2C <- compute_loo(fit_2C, "Model 2C (OHD Ordered)")
loo_2D <- compute_loo(fit_2D, "Model 2D (OHD Single)")
loo_2E <- compute_loo(fit_2E, "Model 2E (Combined Ordered)")
loo_2F <- compute_loo(fit_2F, "Model 2F (Combined Single)")

# KEY COMPARISONS

# Comparison 1: C14 Ordered vs Single
cat("\nComparison 2A vs 2B: Does phase structure help C14?\n")
comp_c14 <- loo_compare(loo_2A, loo_2B)
print(comp_c14)

# Comparison 2: OHD Ordered vs Single (KEY TEST)
cat("\n*** KEY TEST: Comparison 2C vs 2D: Does phase structure help OHD? ***\n")
comp_ohd <- loo_compare(loo_2C, loo_2D)
print(comp_ohd)

# Interpret
elpd_diff <- comp_ohd[2, "elpd_diff"]
se_diff <- comp_ohd[2, "se_diff"]
z_score <- elpd_diff / se_diff

cat(sprintf("  ELPD difference: %.2f (SE: %.2f)\n", elpd_diff, se_diff))
cat(sprintf("  Z-score: %.2f\n", z_score))

if (abs(z_score) < 2) {
  cat("  Result: INDISTINGUISHABLE - OHD does not prefer phases over single population\n")
  cat("  At a STRATIFIED site, this proves OHD is uninformative about stratigraphy!\n")
} else if (z_score > 2) {
  cat("  Result: Ordered model BETTER - OHD fits phase structure\n")
} else {
  cat("  Result: Single-phase model BETTER - OHD conflicts with stratigraphic ordering\n")
}

# Comparison 3: Combined Ordered vs Single
cat("\nComparison 2E vs 2F: Does phase structure help combined?\n")
comp_combined <- loo_compare(loo_2E, loo_2F)
print(comp_combined)

# Cross-Method Comparison (C14 vs OHD Response to Phase Structure)
cat("\nCross-method comparison: C14 vs OHD response to phase structure\n")

c14_elpd_diff <- comp_c14[2, "elpd_diff"]
c14_se_diff <- comp_c14[2, "se_diff"]
c14_z <- c14_elpd_diff / c14_se_diff

ohd_elpd_diff <- comp_ohd[2, "elpd_diff"]
ohd_se_diff <- comp_ohd[2, "se_diff"]
ohd_z <- ohd_elpd_diff / ohd_se_diff

cat(sprintf("C14: ELPD diff = %+.2f (SE %.2f), Z = %+.2f\n", c14_elpd_diff, c14_se_diff, c14_z))
cat(sprintf("OHD: ELPD diff = %+.2f (SE %.2f), Z = %+.2f\n", ohd_elpd_diff, ohd_se_diff, ohd_z))

if (abs(c14_z) >= 2 && abs(ohd_z) < 2) {
  cat("C14 responds to stratigraphic phases but OHD does not.\n")
} else if (abs(c14_z) < 2 && abs(ohd_z) < 2) {
  cat("Neither method strongly responds to phases.\n")
} else {
  cat("Both methods respond similarly to phase structure.\n")
}

# Stacking Weights
cat("\nStacking weights:\n")

# C14 models
cat("C14 Models:\n")
stack_c14 <- loo_model_weights(list(loo_2A, loo_2B), method = "stacking")
names(stack_c14) <- c("2A (Ordered)", "2B (Single)")
print(round(stack_c14, 3))

# OHD models
cat("OHD Models (KEY):\n")
stack_ohd <- loo_model_weights(list(loo_2C, loo_2D), method = "stacking")
names(stack_ohd) <- c("2C (Ordered)", "2D (Single)")
print(round(stack_ohd, 3))

# Combined models
cat("Combined Models:\n")
stack_combined <- loo_model_weights(list(loo_2E, loo_2F), method = "stacking")
names(stack_combined) <- c("2E (Ordered)", "2F (Single)")
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

extract_boundaries(fit_2A, "Model 2A (C14 Ordered)")
extract_boundaries(fit_2C, "Model 2C (OHD Ordered)")
extract_boundaries(fit_2E, "Model 2E (Combined Ordered)")

# Summary Table
summary_table <- data.frame(
  Model = c("2A (C14 Ordered)", "2B (C14 Single)",
            "2C (OHD Ordered)", "2D (OHD Single)",
            "2E (Combined Ordered)", "2F (Combined Single)"),
  Data = c("C14", "C14", "OHD", "OHD", "Combined", "Combined"),
  Structure = c("Ordered", "Single", "Ordered", "Single", "Ordered", "Single"),
  N = c(nrow(c14_data), nrow(c14_data),
        nrow(ohd_data), nrow(ohd_data),
        nrow(c14_data) + nrow(ohd_data), nrow(c14_data) + nrow(ohd_data)),
  ELPD = c(loo_2A$estimates["elpd_loo", "Estimate"],
           loo_2B$estimates["elpd_loo", "Estimate"],
           loo_2C$estimates["elpd_loo", "Estimate"],
           loo_2D$estimates["elpd_loo", "Estimate"],
           loo_2E$estimates["elpd_loo", "Estimate"],
           loo_2F$estimates["elpd_loo", "Estimate"]),
  ELPD_SE = c(loo_2A$estimates["elpd_loo", "SE"],
              loo_2B$estimates["elpd_loo", "SE"],
              loo_2C$estimates["elpd_loo", "SE"],
              loo_2D$estimates["elpd_loo", "SE"],
              loo_2E$estimates["elpd_loo", "SE"],
              loo_2F$estimates["elpd_loo", "SE"])
)

print(summary_table)

# Save summary
write.csv(summary_table, here("output", "tables", "ann_model_comparison.csv"), row.names = FALSE)

cat("\nAhu Nau Nau analysis complete.\n")
