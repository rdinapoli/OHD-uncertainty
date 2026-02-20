# 10_variance_attribution.R - Sensitivity-weighted variance decomposition
# Purpose: Delta-method variance attribution by parameter source (Table 2, Table S2)
# Inputs: output/fits/tier3_fit.rds
# Outputs: output/tables/sensitivity_decomposition_results.csv,
#   output/tables/sensitivity_decomposition_summary.csv
# Runtime: ~5 minutes

library(tidyverse)
library(rstan)
library(here)

options(mc.cores = parallel::detectCores())

# ===== STEP 1: LOAD TIER 3 POSTERIOR =====
cat("Step 1: Loading Tier 3 posterior\n")

tier3_fit <- readRDS(here("output/fits/tier3_fit.rds"))
cat("  Loaded Tier 3 fit object\n\n")

# Extract all parameter posteriors
posterior <- rstan::extract(tier3_fit)

# Get dimensions
n_draws <- dim(posterior$age)[1]  # Number of posterior draws
n_samples <- dim(posterior$age)[2]  # Number of OHD samples (242)

cat(sprintf("Posterior dimensions: %d draws x %d samples\n\n", n_draws, n_samples))

# ===== STEP 2: DEFINE RATE EQUATION =====
cat("Step 2: Define OHD rate equation\n")

# Final rate equation coefficients (from rate equation audit)
a_eq <- 23.8042281295
b_eq <- 0.3782846434
c_eq <- 0.0002071707
d_eq <- -6313.2075128583
R_gas <- 8.314

# RH coefficient (will use posterior if available, else fixed)
e_RH_default <- 2.15

# Function to compute age given parameters
compute_age <- function(H2Ome, H2Ot, Ea, EHT, RH, e_RH = e_RH_default) {
  T_K <- EHT + 273.15
  log_rate <- a_eq + b_eq * (-Ea / (R_gas * T_K)) +
              c_eq * H2Ot + d_eq * (1 / T_K) +
              e_RH * (RH - 0.98)

  # Clamp for numerical stability
  log_rate <- pmax(pmin(log_rate, 50), -50)
  rate <- exp(log_rate)

  age <- H2Ome^2 / rate
  return(age)
}

cat("  log(rate) = a + b*(-Ea/(R*T)) + c*H2Ot + d*(1/T) + e_RH*(RH - 0.98)\n")
cat("  age = H2Ome^2 / rate\n\n")

# ===== STEP 3: COMPUTE NUMERICAL PARTIAL DERIVATIVES =====
cat("Step 3: Computing numerical partial derivatives\n")

# Load empirical data for H2Ome values
empirical_data <- read_csv(here("data/stevenson_2015.csv"),
                           show_col_types = FALSE)
H2Ome <- empirical_data$H2Ome_pct

# Step sizes for numerical differentiation (central differences)
eps_H2Ot <- 0.001   # wt%
eps_Ea <- 100       # J/mol
eps_EHT <- 0.1      # C
eps_RH <- 0.01      # fraction

# Initialize arrays for sensitivities
# We'll compute at the posterior mean for each sample
sens_H2Ot <- matrix(NA, n_draws, n_samples)
sens_Ea <- matrix(NA, n_draws, n_samples)
sens_EHT <- matrix(NA, n_draws, n_samples)
sens_RH <- matrix(NA, n_draws, n_samples)

cat("Computing partial derivatives for all posterior draws...\n")
cat("This may take a few minutes...\n\n")

# Get e_RH posterior if available
if ("e_RH" %in% names(posterior)) {
  e_RH_draws <- posterior$e_RH
} else {
  e_RH_draws <- rep(e_RH_default, n_draws)
}

# Subsample draws for efficiency (use every 10th draw)
draw_indices <- seq(1, n_draws, by = 10)
n_subdraws <- length(draw_indices)

cat(sprintf("Using %d subsampled draws for efficiency\n\n", n_subdraws))

pb <- txtProgressBar(min = 0, max = n_subdraws, style = 3)

for (d_idx in seq_along(draw_indices)) {
  d <- draw_indices[d_idx]

  for (i in 1:n_samples) {
    # Get current parameter values
    H2Ot_i <- posterior$H2Ot[d, i]
    Ea_i <- posterior$Ea[d, i]
    EHT_i <- posterior$EHT[d, i]
    RH_i <- posterior$RH[d, i]
    e_RH_i <- e_RH_draws[d]
    H2Ome_i <- H2Ome[i]

    # Central difference for dage/dH2Ot
    age_plus <- compute_age(H2Ome_i, H2Ot_i + eps_H2Ot, Ea_i, EHT_i, RH_i, e_RH_i)
    age_minus <- compute_age(H2Ome_i, H2Ot_i - eps_H2Ot, Ea_i, EHT_i, RH_i, e_RH_i)
    sens_H2Ot[d, i] <- (age_plus - age_minus) / (2 * eps_H2Ot)

    # Central difference for dage/dEa
    age_plus <- compute_age(H2Ome_i, H2Ot_i, Ea_i + eps_Ea, EHT_i, RH_i, e_RH_i)
    age_minus <- compute_age(H2Ome_i, H2Ot_i, Ea_i - eps_Ea, EHT_i, RH_i, e_RH_i)
    sens_Ea[d, i] <- (age_plus - age_minus) / (2 * eps_Ea)

    # Central difference for dage/dEHT
    age_plus <- compute_age(H2Ome_i, H2Ot_i, Ea_i, EHT_i + eps_EHT, RH_i, e_RH_i)
    age_minus <- compute_age(H2Ome_i, H2Ot_i, Ea_i, EHT_i - eps_EHT, RH_i, e_RH_i)
    sens_EHT[d, i] <- (age_plus - age_minus) / (2 * eps_EHT)

    # Central difference for dage/dRH
    age_plus <- compute_age(H2Ome_i, H2Ot_i, Ea_i, EHT_i, RH_i + eps_RH, e_RH_i)
    age_minus <- compute_age(H2Ome_i, H2Ot_i, Ea_i, EHT_i, RH_i - eps_RH, e_RH_i)
    sens_RH[d, i] <- (age_plus - age_minus) / (2 * eps_RH)
  }

  setTxtProgressBar(pb, d_idx)
}
close(pb)

cat("\nComputed partial derivatives\n\n")

# ===== STEP 4: COMPUTE VARIANCE CONTRIBUTIONS =====
cat("Step 4: Computing variance contributions\n")

# For each sample, compute the variance contribution from each parameter
# Using: Var_from_theta = (dage/dtheta)^2 * Var(theta)

# Get parameter variances across posterior (for each sample)
var_H2Ot <- apply(posterior$H2Ot, 2, var)
var_Ea <- apply(posterior$Ea, 2, var)
var_EHT <- apply(posterior$EHT, 2, var)
var_RH <- apply(posterior$RH, 2, var)

# Get mean sensitivities (averaged across subsampled draws)
mean_sens_H2Ot <- apply(sens_H2Ot[draw_indices, ], 2, mean, na.rm = TRUE)
mean_sens_Ea <- apply(sens_Ea[draw_indices, ], 2, mean, na.rm = TRUE)
mean_sens_EHT <- apply(sens_EHT[draw_indices, ], 2, mean, na.rm = TRUE)
mean_sens_RH <- apply(sens_RH[draw_indices, ], 2, mean, na.rm = TRUE)

# Compute variance contributions
var_contrib_H2Ot <- mean_sens_H2Ot^2 * var_H2Ot
var_contrib_Ea <- mean_sens_Ea^2 * var_Ea
var_contrib_EHT <- mean_sens_EHT^2 * var_EHT
var_contrib_RH <- mean_sens_RH^2 * var_RH

# Total from diagonal terms
var_total_diag <- var_contrib_H2Ot + var_contrib_Ea + var_contrib_EHT + var_contrib_RH

# Actual total variance from age posterior
var_age_actual <- apply(posterior$age, 2, var)

cat("Summary of mean sensitivities (dage/dtheta):\n")
cat(sprintf("  dage/dH2Ot: %.1f years per 1 wt%% (median)\n", median(mean_sens_H2Ot)))
cat(sprintf("  dage/dEa:   %.4f years per 1 J/mol (median)\n", median(mean_sens_Ea)))
cat(sprintf("  dage/dEHT:  %.1f years per 1C (median)\n", median(mean_sens_EHT)))
cat(sprintf("  dage/dRH:   %.1f years per 0.01 RH (median)\n", median(mean_sens_RH)))
cat("\n")

# ===== STEP 5: COMPUTE COVARIANCE CONTRIBUTIONS =====
cat("Step 5: Computing covariance contributions\n")

# Covariance terms: 2 * (dage/dtheta_i)(dage/dtheta_j) * Cov(theta_i, theta_j)
# These capture the confounding between parameters

# For each sample, compute covariances
cov_H2Ot_Ea <- sapply(1:n_samples, function(i) cov(posterior$H2Ot[, i], posterior$Ea[, i]))
cov_H2Ot_EHT <- sapply(1:n_samples, function(i) cov(posterior$H2Ot[, i], posterior$EHT[, i]))
cov_H2Ot_RH <- sapply(1:n_samples, function(i) cov(posterior$H2Ot[, i], posterior$RH[, i]))
cov_Ea_EHT <- sapply(1:n_samples, function(i) cov(posterior$Ea[, i], posterior$EHT[, i]))
cov_Ea_RH <- sapply(1:n_samples, function(i) cov(posterior$Ea[, i], posterior$RH[, i]))
cov_EHT_RH <- sapply(1:n_samples, function(i) cov(posterior$EHT[, i], posterior$RH[, i]))

# Covariance contributions to variance
cov_contrib_H2Ot_Ea <- 2 * mean_sens_H2Ot * mean_sens_Ea * cov_H2Ot_Ea
cov_contrib_H2Ot_EHT <- 2 * mean_sens_H2Ot * mean_sens_EHT * cov_H2Ot_EHT
cov_contrib_H2Ot_RH <- 2 * mean_sens_H2Ot * mean_sens_RH * cov_H2Ot_RH
cov_contrib_Ea_EHT <- 2 * mean_sens_Ea * mean_sens_EHT * cov_Ea_EHT
cov_contrib_Ea_RH <- 2 * mean_sens_Ea * mean_sens_RH * cov_Ea_RH
cov_contrib_EHT_RH <- 2 * mean_sens_EHT * mean_sens_RH * cov_EHT_RH

# Total covariance contribution
cov_total <- cov_contrib_H2Ot_Ea + cov_contrib_H2Ot_EHT + cov_contrib_H2Ot_RH +
             cov_contrib_Ea_EHT + cov_contrib_Ea_RH + cov_contrib_EHT_RH

# Full variance from Taylor expansion
var_total_taylor <- var_total_diag + cov_total

cat("Median covariances (showing confounding):\n")
cat(sprintf("  Cov(H2Ot, Ea):  %.6f\n", median(cov_H2Ot_Ea)))
cat(sprintf("  Cov(H2Ot, EHT): %.6f\n", median(cov_H2Ot_EHT)))
cat(sprintf("  Cov(H2Ot, RH):  %.6f\n", median(cov_H2Ot_RH)))
cat(sprintf("  Cov(Ea, EHT):   %.6f\n", median(cov_Ea_EHT)))
cat(sprintf("  Cov(Ea, RH):    %.6f\n", median(cov_Ea_RH)))
cat(sprintf("  Cov(EHT, RH):   %.6f\n", median(cov_EHT_RH)))
cat("\n")

# ===== STEP 6: VARIANCE DECOMPOSITION RESULTS =====
cat("Step 6: Variance decomposition results\n\n")

# Compute percentages (using absolute values for proportion calculation)
total_abs_contrib <- abs(var_contrib_H2Ot) + abs(var_contrib_Ea) +
                     abs(var_contrib_EHT) + abs(var_contrib_RH) + abs(cov_total)

pct_H2Ot <- 100 * median(abs(var_contrib_H2Ot) / total_abs_contrib)
pct_Ea <- 100 * median(abs(var_contrib_Ea) / total_abs_contrib)
pct_EHT <- 100 * median(abs(var_contrib_EHT) / total_abs_contrib)
pct_RH <- 100 * median(abs(var_contrib_RH) / total_abs_contrib)
pct_cov <- 100 * median(abs(cov_total) / total_abs_contrib)

# Group into composition vs environment
pct_composition <- pct_H2Ot + pct_Ea
pct_environment <- pct_EHT + pct_RH

cat("Individual Parameter Contributions:\n")
cat(sprintf("  H2Ot (structural water):    %5.1f%%\n", pct_H2Ot))
cat(sprintf("  Ea (activation energy):     %5.1f%%\n", pct_Ea))
cat(sprintf("  EHT (temperature):          %5.1f%%\n", pct_EHT))
cat(sprintf("  RH (relative humidity):     %5.1f%%\n", pct_RH))
cat(sprintf("  Covariance terms:           %5.1f%%\n", pct_cov))
cat("\n")

cat("Grouped Contributions:\n")
cat(sprintf("  COMPOSITION (H2Ot + Ea):    %5.1f%%\n", pct_composition))
cat(sprintf("  ENVIRONMENT (EHT + RH):     %5.1f%%\n", pct_environment))
cat(sprintf("  COVARIANCE (confounding):   %5.1f%%\n", pct_cov))
cat("\n")

# Verify Taylor approximation quality
cat("Taylor Approximation Quality:\n")
cat(sprintf("  Median actual age variance:  %.1f years^2\n", median(var_age_actual)))
cat(sprintf("  Median Taylor variance:      %.1f years^2\n", median(var_total_taylor)))
cat(sprintf("  Ratio (Taylor/Actual):       %.2f\n", median(var_total_taylor) / median(var_age_actual)))
cat("\n")

# ===== STEP 7: INTERPRETATION =====
cat("Step 7: Key findings\n\n")

cat(sprintf("1. A 1C change in EHT changes age by ~%.0f years\n", abs(median(mean_sens_EHT))))
cat(sprintf("2. A 0.01 wt%% change in H2Ot changes age by ~%.0f years\n", abs(median(mean_sens_H2Ot)) * 0.01))
cat(sprintf("3. Covariance terms contribute %.1f%% of total variance\n\n", pct_cov))

# ===== STEP 8: SAVE RESULTS =====
cat("Step 8: Saving results\n")

# Create results dataframe
sensitivity_results <- data.frame(
  sample = 1:n_samples,
  lab_id = empirical_data$lab_id,
  study_area = empirical_data$study_area,
  sens_H2Ot = mean_sens_H2Ot,
  sens_Ea = mean_sens_Ea,
  sens_EHT = mean_sens_EHT,
  sens_RH = mean_sens_RH,
  var_contrib_H2Ot = var_contrib_H2Ot,
  var_contrib_Ea = var_contrib_Ea,
  var_contrib_EHT = var_contrib_EHT,
  var_contrib_RH = var_contrib_RH,
  cov_total = cov_total,
  var_taylor = var_total_taylor,
  var_actual = var_age_actual
)

write_csv(sensitivity_results,
          here("output/tables/sensitivity_decomposition_results.csv"))
cat("  Saved per-sample sensitivity results\n")

# Save summary
summary_stats <- data.frame(
  parameter = c("H2Ot", "Ea", "EHT", "RH", "Covariance", "Composition", "Environment"),
  pct_contribution = c(pct_H2Ot, pct_Ea, pct_EHT, pct_RH, pct_cov, pct_composition, pct_environment),
  median_sensitivity = c(median(mean_sens_H2Ot), median(mean_sens_Ea),
                         median(mean_sens_EHT), median(mean_sens_RH), NA, NA, NA),
  units = c("years/wt%", "years/(J/mol)", "years/C", "years/RH", NA, NA, NA)
)

write_csv(summary_stats,
          here("output/tables/sensitivity_decomposition_summary.csv"))
cat("  Saved summary statistics\n")

# ===== STEP 9: CREATE DIAGNOSTIC PLOTS =====
cat("\nCreating diagnostic plots...\n")

pdf(here("output/figures/sensitivity_diagnostics.pdf"),
    width = 12, height = 10)

# 1. Pie chart of variance contributions
par(mfrow = c(2, 2))

# Panel 1: Individual contributions
contributions <- c(pct_H2Ot, pct_Ea, pct_EHT, pct_RH, pct_cov)
labels <- c(sprintf("H2Ot\n%.1f%%", pct_H2Ot),
            sprintf("Ea\n%.1f%%", pct_Ea),
            sprintf("EHT\n%.1f%%", pct_EHT),
            sprintf("RH\n%.1f%%", pct_RH),
            sprintf("Cov\n%.1f%%", pct_cov))
pie(contributions, labels = labels,
    main = "Sensitivity-Weighted Variance Decomposition\n(Individual Parameters)",
    col = c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#7f7f7f"))

# Panel 2: Grouped contributions
grouped <- c(pct_composition, pct_environment, pct_cov)
labels_grouped <- c(sprintf("Composition\n%.1f%%", pct_composition),
                    sprintf("Environment\n%.1f%%", pct_environment),
                    sprintf("Covariance\n%.1f%%", pct_cov))
pie(grouped, labels = labels_grouped,
    main = "Grouped Variance Decomposition",
    col = c("#1f77b4", "#ff7f0e", "#7f7f7f"))

# Panel 3: Sensitivity distributions
boxplot(list(H2Ot = mean_sens_H2Ot,
             EHT = mean_sens_EHT),
        main = "Sensitivities by Parameter\n(years per unit)",
        ylab = "Sensitivity (dage/dtheta)",
        col = c("#1f77b4", "#ff7f0e"))

# Panel 4: Taylor vs Actual variance
plot(var_age_actual, var_total_taylor,
     xlab = "Actual Age Variance (years^2)",
     ylab = "Taylor Approximation Variance (years^2)",
     main = "Taylor Approximation Quality",
     pch = 16, col = rgb(0, 0, 0, 0.3))
abline(0, 1, col = "red", lwd = 2)

dev.off()
cat("  Saved diagnostic plots\n")

# ===== FINAL SUMMARY =====
cat("\nFinal variance attribution:\n")
cat(sprintf("  COMPOSITION (H2Ot + Ea):    %5.1f%%\n", pct_composition))
cat(sprintf("  ENVIRONMENT (EHT + RH):     %5.1f%%\n", pct_environment))
cat(sprintf("  COVARIANCE (confounding):   %5.1f%%\n", pct_cov))
cat(sprintf("\n  dage/dEHT = %.1f years/C (median)\n", median(mean_sens_EHT)))
cat(sprintf("  dage/dH2Ot = %.0f years/wt%% (median)\n", median(mean_sens_H2Ot)))
cat("\n")
