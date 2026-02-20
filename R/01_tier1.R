# 01_tier1.R - Fit Tier 1 (composition-only) OHD model
# Purpose: Fits Stan model with H2Ome measurement uncertainty only
# Inputs: data/stevenson_2015.csv, stan/tier1_composition.stan
# Outputs: output/fits/tier1_fit.rds, output/tables/tier1_results.csv, output/tables/tier1_summary.csv
# Runtime: ~5 minutes

library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(here)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ===== STEP 1: LOAD EMPIRICAL DATA =====
cat("Step 1: Loading Stevenson 2015 empirical data\n")

empirical_data <- read_csv(
  here("data/stevenson_2015.csv"),
  show_col_types = FALSE
)

cat(sprintf("Loaded %d samples from Stevenson et al. (2015) PNAS study\n\n", nrow(empirical_data)))

# ===== STEP 2: PREPARE DATA FOR STAN =====
cat("Step 2: Preparing data for Tier 1 analysis\n")

# Map study areas to integers (for reference, not used in Tier 1)
study_area_map <- c("SA1" = 1, "SA2" = 2, "SA3" = 3)
empirical_data <- empirical_data %>%
  mutate(study_area_id = study_area_map[study_area])

cat("Study area distribution:\n")
print(table(empirical_data$study_area))
cat("\n")

# Rate equation coefficients from Stevenson (2015) supplementary material
a_eq <- 23.8042281295      # Intercept
b_eq <- 0.3782846434       # Coefficient of (-Ea/(R×T))
c_eq <- 0.0002071707       # Coefficient of H₂Ot
d_eq <- -6313.2075128583   # Coefficient of (1/T)
R_gas <- 8.314             # J/(mol·K)

# ===== TIER 1: COMPOSITION UNCERTAINTY ONLY =====
cat("Fitting Tier 1: Composition uncertainty (H2Ot, Ea)\n\n")

# Prepare Stan data - NO SOURCE ASSIGNMENTS
tier1_data <- list(
  N = nrow(empirical_data),
  H2Ome = empirical_data$H2Ome_pct,
  EHT = empirical_data$eht_c,          # Stevenson's study-area constants
  R_gas = R_gas,
  a_eq = a_eq,
  b_eq = b_eq,
  c_eq = c_eq,
  d_eq = d_eq
)

cat("Compiling Tier 1 Stan model...\n")
tier1_model <- stan_model(here("stan/tier1_composition.stan"))

cat("Sampling from Tier 1 model...\n")
# adapt_delta = 0.99: increased from default 0.80 to reduce divergent transitions in hierarchical model
# max_treedepth = 15: increased from default 10 for complex posterior geometries
tier1_fit <- sampling(
  tier1_model,
  data = tier1_data,
  chains = 4,
  iter = 10000,
  warmup = 5000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 42
)

cat("Tier 1 sampling complete\n\n")

# Extract posteriors
tier1_ages <- rstan::extract(tier1_fit, "age")[[1]]
tier1_age_summary <- apply(tier1_ages, 2, function(x) {
  c(mean = mean(x), sd = sd(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)))
})

# Create results dataframe
tier1_results <- data.frame(
  lab_id = empirical_data$lab_id,
  study_area = empirical_data$study_area,
  stevenson_age = empirical_data$date_bp,
  stevenson_sd = empirical_data$sd,
  post_mean = tier1_age_summary["mean", ],
  post_sd = tier1_age_summary["sd", ],
  q025 = tier1_age_summary["q025", ],
  q975 = tier1_age_summary["q975", ]
) %>%
  mutate(
    ci_width = q975 - q025,
    age_diff = post_mean - stevenson_age,
    exceeds_claimed = ci_width > 60  # Stevenson claimed ±30 years (60 year CI width)
  )

cat("Tier 1 Results:\n")
cat(sprintf("  Median posterior SD: +/-%0.1f years\n", median(tier1_results$post_sd)))
cat(sprintf("  Median 95%% CI width: +/-%0.1f years\n", median(tier1_results$ci_width)/2))
cat(sprintf("  Ratio: %.1fx claimed precision\n", (median(tier1_results$ci_width)/2)/30))
cat(sprintf("  Samples exceeding +/-30 year claim: %.1f%%\n\n",
            100 * mean(tier1_results$exceeds_claimed)))

# Convergence diagnostics
tier1_summary <- summary(tier1_fit)$summary
tier1_rhat_max <- max(tier1_summary[, "Rhat"], na.rm = TRUE)
tier1_ess_min <- min(tier1_summary[, "n_eff"], na.rm = TRUE)

cat(sprintf("  Convergence: max Rhat = %.3f, min ESS = %.0f\n", tier1_rhat_max, tier1_ess_min))

# HMC diagnostics
n_divergent <- sum(rstan::get_num_divergent(tier1_fit))
n_max_treedepth <- sum(rstan::get_num_max_treedepth(tier1_fit))
cat(sprintf("  Divergent transitions: %d\n", n_divergent))
cat(sprintf("  Max treedepth exceedances: %d\n\n", n_max_treedepth))

if (tier1_rhat_max > 1.01) {
  cat("WARNING: R-hat > 1.01 suggests convergence issues\n")
  cat("  Consider increasing iterations or adjusting priors\n\n")
}
if (n_divergent > 0) warning("DIVERGENT TRANSITIONS DETECTED")
if (n_max_treedepth > 0) warning("MAX TREEDEPTH EXCEEDANCES DETECTED")

# ===== SAVE RESULTS =====
cat("Saving results\n")

# Create output directories
dir.create(here("output/fits"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/figures"), recursive = TRUE, showWarnings = FALSE)

# Save fit object
saveRDS(tier1_fit, here("output/fits/tier1_fit.rds"))

# Save results table
write_csv(tier1_results, here("output/tables/tier1_results.csv"))

# Save summary statistics
tier1_summary_stats <- data.frame(
  tier = "Tier 1",
  median_sd = median(tier1_results$post_sd),
  median_ci_width = median(tier1_results$ci_width),
  pct_exceeding_claimed = 100 * mean(tier1_results$exceeds_claimed),
  rhat_max = tier1_rhat_max,
  ess_min = tier1_ess_min,
  n_divergent = n_divergent,
  n_max_treedepth = n_max_treedepth,
  n_samples = nrow(empirical_data)
)
write_csv(tier1_summary_stats, here("output/tables/tier1_summary.csv"))

# ===== CREATE DIAGNOSTIC PLOTS =====
cat("Creating diagnostic plots...\n")

pdf(here("output/figures/tier1_diagnostics.pdf"), width = 12, height = 8)

# 1. Comparison with Stevenson's reported ages
p1 <- ggplot(tier1_results, aes(x = stevenson_age, y = post_mean)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Tier 1: Posterior Ages vs Stevenson's Reported Ages",
       x = "Stevenson's Reported Age (BP)",
       y = "Tier 1 Posterior Mean Age (BP)") +
  theme_minimal()
print(p1)

# 2. CI width distribution
p2 <- ggplot(tier1_results, aes(x = ci_width/2)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 30, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 30, y = Inf, label = "Stevenson's +/-30 years",
           vjust = 2, hjust = -0.1, color = "red") +
  labs(title = "Tier 1: Distribution of 95% CI Half-Widths",
       subtitle = sprintf("Median: +/-%.1f years (%.1fx claimed)",
                          median(tier1_results$ci_width)/2,
                          (median(tier1_results$ci_width)/2)/30),
       x = "95% CI Half-Width (years)", y = "Count") +
  theme_minimal()
print(p2)

# 3. CI width by study area
p3 <- ggplot(tier1_results, aes(x = study_area, y = ci_width/2)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 30, color = "red", linetype = "dashed") +
  labs(title = "Tier 1: CI Half-Width by Study Area",
       subtitle = "Composition uncertainty only (EHT/RH known)",
       x = "Study Area", y = "95% CI Half-Width (years)") +
  theme_minimal()
print(p3)

dev.off()

cat(sprintf("Tier 1 complete. Median uncertainty: +/-%0.1f years (%.1fx Stevenson's +/-30)\n",
            median(tier1_results$ci_width)/2,
            (median(tier1_results$ci_width)/2)/30))
