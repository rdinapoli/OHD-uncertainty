# 03_tier3.R - Fit Tier 3 (complete) OHD model
# Purpose: Full model with composition, environment, and temporal climate uncertainty
# Inputs: data/stevenson_2015.csv, stan/tier3_complete.stan
# Outputs: output/fits/tier3_fit.rds, output/tables/tier3_results.csv, output/tables/tier3_summary.csv
# Runtime: ~15 minutes

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
cat("Step 2: Preparing data for Tier 3 analysis\n")

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

# ===== TIER 3: COMPLETE MODEL =====
cat("Fitting Tier 3: Complete model (composition + environment + temporal climate)\n\n")

# RH coefficient parameters (from Mazer et al. 1991 - Easter Island specific)
e_RH_mean <- 2.15
e_RH_sd <- 0.5

# Temporal climate variation (±0.5 C over artifact lifetime)
# Conservative estimate based on Holocene climate variability
temporal_climate_sd <- 0.5

# Prepare Stan data - NO SOURCE ASSIGNMENTS
tier3_data <- list(
  N = nrow(empirical_data),
  N_areas = 3,
  H2Ome = empirical_data$H2Ome_pct,
  study_area = empirical_data$study_area_id,
  temporal_climate_sd = temporal_climate_sd,
  R_gas = R_gas,
  a_eq = a_eq,
  b_eq = b_eq,
  c_eq = c_eq,
  d_eq = d_eq,
  e_RH_mean = e_RH_mean,
  e_RH_sd = e_RH_sd
)

cat("Compiling Tier 3 Stan model...\n")
tier3_model <- stan_model(here("stan/tier3_complete.stan"))

cat("Sampling from Tier 3 model...\n")
# adapt_delta = 0.99: increased from default 0.80 to reduce divergent transitions in hierarchical model
# max_treedepth = 15: increased from default 10 for complex posterior geometries
tier3_fit <- sampling(
  tier3_model,
  data = tier3_data,
  chains = 4,
  iter = 10000,
  warmup = 5000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 42
)

cat("Tier 3 sampling complete\n\n")

# Extract posteriors
tier3_ages <- rstan::extract(tier3_fit, "age")[[1]]
tier3_age_summary <- apply(tier3_ages, 2, function(x) {
  c(mean = mean(x), sd = sd(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)))
})

# Create results dataframe
tier3_results <- data.frame(
  lab_id = empirical_data$lab_id,
  study_area = empirical_data$study_area,
  stevenson_age = empirical_data$date_bp,
  stevenson_sd = empirical_data$sd,
  post_mean = tier3_age_summary["mean", ],
  post_sd = tier3_age_summary["sd", ],
  q025 = tier3_age_summary["q025", ],
  q975 = tier3_age_summary["q975", ]
) %>%
  mutate(
    ci_width = q975 - q025,
    age_diff = post_mean - stevenson_age,
    exceeds_claimed = ci_width > 60
  )

cat("Tier 3 Results:\n")
cat(sprintf("  Median posterior SD: +/-%0.1f years\n", median(tier3_results$post_sd)))
cat(sprintf("  Median 95%% CI width: +/-%0.1f years\n", median(tier3_results$ci_width)/2))
cat(sprintf("  Ratio: %.1fx claimed precision\n", (median(tier3_results$ci_width)/2)/30))
cat(sprintf("  Samples exceeding +/-30 year claim: %.1f%%\n\n",
            100 * mean(tier3_results$exceeds_claimed)))

# Load Tier 1 and Tier 2 results for comparison
tier1_results <- read_csv(here("output/tables/tier1_results.csv"),
                          show_col_types = FALSE)
tier2_results <- read_csv(here("output/tables/tier2_results.csv"),
                          show_col_types = FALSE)

cat("Comparison with Tier 1 and Tier 2:\n")
cat(sprintf("  Tier 1 median CI: +/-%0.1f years (composition only)\n", median(tier1_results$ci_width)/2))
cat(sprintf("  Tier 2 median CI: +/-%0.1f years (+ environment)\n", median(tier2_results$ci_width)/2))
cat(sprintf("  Tier 3 median CI: +/-%0.1f years (+ temporal climate)\n", median(tier3_results$ci_width)/2))
cat(sprintf("  Increase from Tier 2: +%.1f years (%.0f%% increase)\n\n",
            median(tier3_results$ci_width)/2 - median(tier2_results$ci_width)/2,
            100 * ((median(tier3_results$ci_width)/2) / (median(tier2_results$ci_width)/2) - 1)))

# Convergence diagnostics
tier3_summary <- summary(tier3_fit)$summary
tier3_rhat_max <- max(tier3_summary[, "Rhat"], na.rm = TRUE)
tier3_ess_min <- min(tier3_summary[, "n_eff"], na.rm = TRUE)

cat(sprintf("  Convergence: max Rhat = %.3f, min ESS = %.0f\n", tier3_rhat_max, tier3_ess_min))

# HMC diagnostics
n_divergent <- sum(rstan::get_num_divergent(tier3_fit))
n_max_treedepth <- sum(rstan::get_num_max_treedepth(tier3_fit))
cat(sprintf("  Divergent transitions: %d\n", n_divergent))
cat(sprintf("  Max treedepth exceedances: %d\n\n", n_max_treedepth))

if (tier3_rhat_max > 1.01) {
  cat("WARNING: R-hat > 1.01 suggests convergence issues\n\n")
}
if (n_divergent > 0) warning("DIVERGENT TRANSITIONS DETECTED")
if (n_max_treedepth > 0) warning("MAX TREEDEPTH EXCEEDANCES DETECTED")

# ===== SAVE RESULTS =====
cat("Saving results\n")

dir.create(here("output/fits"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/figures"), recursive = TRUE, showWarnings = FALSE)

saveRDS(tier3_fit, here("output/fits/tier3_fit.rds"))

write_csv(tier3_results, here("output/tables/tier3_results.csv"))

tier3_summary_stats <- data.frame(
  tier = "Tier 3",
  median_sd = median(tier3_results$post_sd),
  median_ci_width = median(tier3_results$ci_width),
  pct_exceeding_claimed = 100 * mean(tier3_results$exceeds_claimed),
  rhat_max = tier3_rhat_max,
  ess_min = tier3_ess_min,
  n_divergent = n_divergent,
  n_max_treedepth = n_max_treedepth,
  n_samples = nrow(empirical_data)
)
write_csv(tier3_summary_stats, here("output/tables/tier3_summary.csv"))

# ===== CREATE DIAGNOSTIC PLOTS =====
cat("Creating diagnostic plots...\n")

pdf(here("output/figures/tier3_diagnostics.pdf"), width = 12, height = 8)

# 1. Tier comparison (all three tiers)
comparison_long <- bind_rows(
  tier1_results %>% mutate(tier = "Tier 1", ci_half = ci_width/2),
  tier2_results %>% mutate(tier = "Tier 2", ci_half = ci_width/2),
  tier3_results %>% mutate(tier = "Tier 3", ci_half = ci_width/2)
)

p1 <- ggplot(comparison_long, aes(x = ci_half, fill = tier)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 30, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 30, y = Inf, label = "Stevenson's +/-30 years",
           vjust = 2, hjust = -0.1, color = "red") +
  labs(title = "Tier Comparison: Distribution of Uncertainties",
       subtitle = "Tier 1: Composition | Tier 2: + Environment | Tier 3: + Temporal Climate",
       x = "95% CI Half-Width (years)", y = "Count") +
  scale_fill_manual(values = c("Tier 1" = "steelblue", "Tier 2" = "darkorange", "Tier 3" = "forestgreen")) +
  theme_minimal()
print(p1)

# 2. Tier 2 vs Tier 3 scatter comparison
comparison_data <- data.frame(
  lab_id = tier2_results$lab_id,
  tier2_ci = tier2_results$ci_width,
  tier3_ci = tier3_results$ci_width,
  study_area = tier2_results$study_area
)

p2 <- ggplot(comparison_data, aes(x = tier2_ci/2, y = tier3_ci/2, color = study_area)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dashed") +
  labs(title = "Tier 3 vs Tier 2: Impact of Temporal Climate Variation",
       subtitle = sprintf("Median increase: +%.1f years",
                          median(tier3_results$ci_width)/2 - median(tier2_results$ci_width)/2),
       x = "Tier 2 CI Half-Width (years)", y = "Tier 3 CI Half-Width (years)") +
  theme_minimal()
print(p2)

# 3. CI width by study area (Tier 3)
p3 <- ggplot(tier3_results, aes(x = study_area, y = ci_width/2)) +
  geom_boxplot(fill = "forestgreen", alpha = 0.7) +
  geom_hline(yintercept = 30, color = "red", linetype = "dashed") +
  labs(title = "Tier 3: CI Half-Width by Study Area",
       subtitle = "Complete model (composition + environment + temporal climate)",
       x = "Study Area", y = "95% CI Half-Width (years)") +
  theme_minimal()
print(p3)

# 4. Posterior ages vs Stevenson's reported ages
p4 <- ggplot(tier3_results, aes(x = stevenson_age, y = post_mean, color = study_area)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_errorbar(aes(ymin = q025, ymax = q975), alpha = 0.2, width = 0) +
  labs(title = "Tier 3: Posterior Ages vs Stevenson's Reported Ages",
       subtitle = "Error bars show 95% credible intervals",
       x = "Stevenson's Reported Age (BP)",
       y = "Tier 3 Posterior Mean Age (BP)") +
  theme_minimal()
print(p4)

# 5. Temporal climate effect distribution
temp_samples <- rstan::extract(tier3_fit, "temporal_climate_effect")[[1]]
temp_df <- data.frame(
  effect = as.vector(temp_samples[sample(nrow(temp_samples), 1000), ])
)

p5 <- ggplot(temp_df, aes(x = effect)) +
  geom_histogram(bins = 50, fill = "forestgreen", alpha = 0.7) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
  labs(title = "Tier 3: Distribution of Temporal Climate Effects",
       subtitle = "Temperature variation over artifact burial lifetime",
       x = "Temperature Deviation (C)", y = "Count") +
  theme_minimal()
print(p5)

dev.off()

# ===== MODEL COMPARISON (LOO) =====
cat("Computing LOO for model comparison...\n")

tier3_loo <- loo(tier3_fit)
saveRDS(tier3_loo, here("output/fits/tier3_loo.rds"))

# Load Tier 1 and Tier 2 LOO if available
tier1_loo_file <- here("output/fits/tier1_loo.rds")
tier2_loo_file <- here("output/fits/tier2_loo.rds")

if (file.exists(tier1_loo_file) && file.exists(tier2_loo_file)) {
  tier1_loo <- readRDS(tier1_loo_file)
  tier2_loo <- readRDS(tier2_loo_file)

  cat("Model Comparison (LOO-CV):\n")
  loo_comparison <- loo_compare(tier1_loo, tier2_loo, tier3_loo)
  print(loo_comparison)

  write.csv(as.data.frame(loo_comparison),
            here("output/tables/loo_comparison.csv"))
}

cat(sprintf("Tier 3 complete. Median uncertainty: +/-%0.1f years (%.1fx Stevenson's +/-30)\n",
            median(tier3_results$ci_width)/2,
            (median(tier3_results$ci_width)/2)/30))
