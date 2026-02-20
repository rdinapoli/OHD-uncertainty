# 02_tier2.R - Fit Tier 2 (composition + environment) OHD model
# Purpose: Adds hierarchical environmental parameters (EHT, RH) to Tier 1
# Inputs: data/stevenson_2015.csv, stan/tier2_environment.stan
# Outputs: output/fits/tier2_fit.rds, output/tables/tier2_results.csv, output/tables/tier2_summary.csv
# Runtime: ~10 minutes

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
cat("Step 2: Preparing data for Tier 2 analysis\n")

# Map study areas to integers
study_area_map <- c("SA1" = 1, "SA2" = 2, "SA3" = 3)
empirical_data <- empirical_data %>%
  mutate(study_area_id = study_area_map[study_area])

cat("Study area distribution:\n")
print(table(empirical_data$study_area))
cat("\n")

# Final rate equation coefficients
a_eq <- 23.8042281295
b_eq <- 0.3782846434
c_eq <- 0.0002071707
d_eq <- -6313.2075128583
R_gas <- 8.314

# ===== TIER 2: + ENVIRONMENTAL UNCERTAINTY =====
cat("Fitting Tier 2: Composition + environmental uncertainty\n\n")

# RH coefficient parameters (from Mazer et al. 1991 - Easter Island specific)
e_RH_mean <- 2.15
e_RH_sd <- 0.5

# Prepare Stan data - NO SOURCE ASSIGNMENTS
tier2_data <- list(
  N = nrow(empirical_data),
  N_areas = 3,
  H2Ome = empirical_data$H2Ome_pct,
  study_area = empirical_data$study_area_id,
  R_gas = R_gas,
  a_eq = a_eq,
  b_eq = b_eq,
  c_eq = c_eq,
  d_eq = d_eq,
  e_RH_mean = e_RH_mean,
  e_RH_sd = e_RH_sd
)

cat("Compiling Tier 2 Stan model...\n")
tier2_model <- stan_model(here("stan/tier2_environment.stan"))

cat("Sampling from Tier 2 model...\n")
# adapt_delta = 0.99: increased from default 0.80 to reduce divergent transitions in hierarchical model
# max_treedepth = 15: increased from default 10 for complex posterior geometries
tier2_fit <- sampling(
  tier2_model,
  data = tier2_data,
  chains = 4,
  iter = 10000,
  warmup = 5000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 42
)

cat("Tier 2 sampling complete\n\n")

# Extract posteriors
tier2_ages <- rstan::extract(tier2_fit, "age")[[1]]
tier2_age_summary <- apply(tier2_ages, 2, function(x) {
  c(mean = mean(x), sd = sd(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)))
})

# Create results dataframe
tier2_results <- data.frame(
  lab_id = empirical_data$lab_id,
  study_area = empirical_data$study_area,
  stevenson_age = empirical_data$date_bp,
  stevenson_sd = empirical_data$sd,
  post_mean = tier2_age_summary["mean", ],
  post_sd = tier2_age_summary["sd", ],
  q025 = tier2_age_summary["q025", ],
  q975 = tier2_age_summary["q975", ]
) %>%
  mutate(
    ci_width = q975 - q025,
    age_diff = post_mean - stevenson_age,
    exceeds_claimed = ci_width > 60
  )

cat("Tier 2 Results:\n")
cat(sprintf("  Median posterior SD: +/-%0.1f years\n", median(tier2_results$post_sd)))
cat(sprintf("  Median 95%% CI width: +/-%0.1f years\n", median(tier2_results$ci_width)/2))
cat(sprintf("  Ratio: %.1fx claimed precision\n", (median(tier2_results$ci_width)/2)/30))
cat(sprintf("  Samples exceeding +/-30 year claim: %.1f%%\n\n",
            100 * mean(tier2_results$exceeds_claimed)))

# Load Tier 1 results for comparison
tier1_results <- read_csv(here("output/tables/tier1_results.csv"),
                          show_col_types = FALSE)

cat("Comparison with Tier 1:\n")
cat(sprintf("  Tier 1 median CI: +/-%0.1f years\n", median(tier1_results$ci_width)/2))
cat(sprintf("  Tier 2 median CI: +/-%0.1f years\n", median(tier2_results$ci_width)/2))
cat(sprintf("  Increase: +%.1f years (%.0f%% increase)\n\n",
            median(tier2_results$ci_width)/2 - median(tier1_results$ci_width)/2,
            100 * ((median(tier2_results$ci_width)/2) / (median(tier1_results$ci_width)/2) - 1)))

# Convergence diagnostics
tier2_summary <- summary(tier2_fit)$summary
tier2_rhat_max <- max(tier2_summary[, "Rhat"], na.rm = TRUE)
tier2_ess_min <- min(tier2_summary[, "n_eff"], na.rm = TRUE)

cat(sprintf("  Convergence: max Rhat = %.3f, min ESS = %.0f\n", tier2_rhat_max, tier2_ess_min))

# HMC diagnostics
n_divergent <- sum(rstan::get_num_divergent(tier2_fit))
n_max_treedepth <- sum(rstan::get_num_max_treedepth(tier2_fit))
cat(sprintf("  Divergent transitions: %d\n", n_divergent))
cat(sprintf("  Max treedepth exceedances: %d\n\n", n_max_treedepth))

if (tier2_rhat_max > 1.01) {
  cat("WARNING: R-hat > 1.01 suggests convergence issues\n\n")
}
if (n_divergent > 0) warning("DIVERGENT TRANSITIONS DETECTED")
if (n_max_treedepth > 0) warning("MAX TREEDEPTH EXCEEDANCES DETECTED")

# ===== SAVE RESULTS =====
cat("Saving results\n")

dir.create(here("output/fits"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/figures"), recursive = TRUE, showWarnings = FALSE)

saveRDS(tier2_fit, here("output/fits/tier2_fit.rds"))

write_csv(tier2_results, here("output/tables/tier2_results.csv"))

tier2_summary_stats <- data.frame(
  tier = "Tier 2",
  median_sd = median(tier2_results$post_sd),
  median_ci_width = median(tier2_results$ci_width),
  pct_exceeding_claimed = 100 * mean(tier2_results$exceeds_claimed),
  rhat_max = tier2_rhat_max,
  ess_min = tier2_ess_min,
  n_divergent = n_divergent,
  n_max_treedepth = n_max_treedepth,
  n_samples = nrow(empirical_data)
)
write_csv(tier2_summary_stats, here("output/tables/tier2_summary.csv"))

# ===== CREATE DIAGNOSTIC PLOTS =====
cat("Creating diagnostic plots...\n")

pdf(here("output/figures/tier2_diagnostics.pdf"), width = 12, height = 8)

# 1. Tier 1 vs Tier 2 comparison
comparison_data <- data.frame(
  lab_id = tier1_results$lab_id,
  tier1_ci = tier1_results$ci_width,
  tier2_ci = tier2_results$ci_width,
  study_area = tier1_results$study_area
)

p1 <- ggplot(comparison_data, aes(x = tier1_ci/2, y = tier2_ci/2)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Tier 2 vs Tier 1: Impact of Environmental Uncertainty",
       subtitle = sprintf("Median increase: +%.1f years",
                          median(tier2_results$ci_width)/2 - median(tier1_results$ci_width)/2),
       x = "Tier 1 CI Half-Width (years)", y = "Tier 2 CI Half-Width (years)") +
  theme_minimal()
print(p1)

# 2. CI width distribution comparison
comparison_long <- bind_rows(
  tier1_results %>% mutate(tier = "Tier 1", ci_half = ci_width/2),
  tier2_results %>% mutate(tier = "Tier 2", ci_half = ci_width/2)
)

p2 <- ggplot(comparison_long, aes(x = ci_half, fill = tier)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 30, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 30, y = Inf, label = "Stevenson's +/-30 years",
           vjust = 2, hjust = -0.1, color = "red") +
  labs(title = "Tier 1 vs Tier 2: Distribution of Uncertainties",
       x = "95% CI Half-Width (years)", y = "Count") +
  scale_fill_manual(values = c("Tier 1" = "steelblue", "Tier 2" = "darkorange")) +
  theme_minimal()
print(p2)

# 3. CI width by study area (Tier 2)
p3 <- ggplot(tier2_results, aes(x = study_area, y = ci_width/2)) +
  geom_boxplot(fill = "darkorange", alpha = 0.7) +
  geom_hline(yintercept = 30, color = "red", linetype = "dashed") +
  labs(title = "Tier 2: CI Half-Width by Study Area",
       subtitle = "Composition + environment uncertainty",
       x = "Study Area", y = "95% CI Half-Width (years)") +
  theme_minimal()
print(p3)

# 4. Inferred EHT distributions
eht_samples <- rstan::extract(tier2_fit, "EHT_mean")[[1]]
eht_df <- data.frame(
  SA1 = eht_samples[, 1],
  SA2 = eht_samples[, 2],
  SA3 = eht_samples[, 3]
) %>% pivot_longer(everything(), names_to = "study_area", values_to = "EHT")

p4 <- ggplot(eht_df, aes(x = EHT, fill = study_area)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = c(21.8, 22.4, 23.3), linetype = "dashed") +
  labs(title = "Tier 2: Posterior EHT Distributions",
       subtitle = "Dashed lines: Stevenson's study-area constants",
       x = "Effective Hydration Temperature (C)", y = "Density") +
  theme_minimal()
print(p4)

# 5. Inferred RH distributions
rh_samples <- rstan::extract(tier2_fit, "RH_mean")[[1]]
rh_df <- data.frame(
  SA1 = rh_samples[, 1],
  SA2 = rh_samples[, 2],
  SA3 = rh_samples[, 3]
) %>% pivot_longer(everything(), names_to = "study_area", values_to = "RH")

p5 <- ggplot(rh_df, aes(x = RH, fill = study_area)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0.98, linetype = "dashed") +
  labs(title = "Tier 2: Posterior RH Distributions",
       subtitle = "Dashed line: Stevenson's universal constant (0.98)",
       x = "Relative Humidity (fraction)", y = "Density") +
  theme_minimal()
print(p5)

dev.off()

cat(sprintf("Tier 2 complete. Median uncertainty: +/-%0.1f years (%.1fx Stevenson's +/-30)\n",
            median(tier2_results$ci_width)/2,
            (median(tier2_results$ci_width)/2)/30))
