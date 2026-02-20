# 09_counterfactual.R - Environment-only counterfactual analysis
# Purpose: Fits model varying only environmental parameters to isolate their contribution
# Inputs: data/stevenson_2015.csv, stan/counterfactual_environment_only.stan
# Outputs: output/fits/counterfactual_env_only_fit.rds,
#   output/tables/counterfactual_env_only_results.csv,
#   output/tables/counterfactual_env_only_summary.csv
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
cat("Step 2: Preparing data\n")

# Map study areas to integers
study_area_map <- c("SA1" = 1, "SA2" = 2, "SA3" = 3)
empirical_data <- empirical_data %>%
  mutate(study_area_id = study_area_map[study_area])

cat("Study area distribution:\n")
print(table(empirical_data$study_area))
cat("\n")

# Fixed composition values (island-wide means from Stevenson 2015)
H2Ot_fixed <- mean(empirical_data$H2Ot_pct)  # ~0.18 wt%
Ea_fixed <- mean(empirical_data$Ea_J_per_mol)  # ~84,000 J/mol

cat(sprintf("Fixed composition values (population means):\n"))
cat(sprintf("  H2Ot = %.4f wt%% (island-wide mean)\n", H2Ot_fixed))
cat(sprintf("  Ea = %.0f J/mol (island-wide mean)\n\n", Ea_fixed))

# Final rate equation coefficients
a_eq <- 23.8042281295
b_eq <- 0.3782846434
c_eq <- 0.0002071707
d_eq <- -6313.2075128583
R_gas <- 8.314

# RH coefficient parameters
e_RH_mean <- 2.15
e_RH_sd <- 0.5

# ===== RUN COUNTERFACTUAL MODEL =====
cat("Running counterfactual: environment only\n\n")

cat("Model assumptions:\n")
cat("  - H2Ome measured (IR-PAS spectroscopy)\n")
cat(sprintf("  - H2Ot FIXED at %.4f wt%% (population mean)\n", H2Ot_fixed))
cat(sprintf("  - Ea FIXED at %.0f J/mol (population mean)\n", Ea_fixed))
cat("  - EHT UNCERTAIN (hierarchical by study area)\n")
cat("  - RH UNCERTAIN (varying by study area)\n\n")

# Prepare Stan data
env_only_data <- list(
  N = nrow(empirical_data),
  N_areas = 3,
  H2Ome = empirical_data$H2Ome_pct,
  study_area = empirical_data$study_area_id,
  H2Ot_fixed = H2Ot_fixed,
  Ea_fixed = Ea_fixed,
  R_gas = R_gas,
  a_eq = a_eq,
  b_eq = b_eq,
  c_eq = c_eq,
  d_eq = d_eq,
  e_RH_mean = e_RH_mean,
  e_RH_sd = e_RH_sd
)

cat("Compiling counterfactual Stan model...\n")
env_only_model <- stan_model(here("stan/counterfactual_environment_only.stan"))

cat("Sampling from environment-only model (this may take 25-35 minutes with 10000 iterations)...\n")
env_only_fit <- sampling(
  env_only_model,
  data = env_only_data,
  chains = 4,
  iter = 10000,
  warmup = 5000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 42
)

cat("Environment-only sampling complete\n\n")

# Extract posteriors
env_only_ages <- rstan::extract(env_only_fit, "age")[[1]]
env_only_age_summary <- apply(env_only_ages, 2, function(x) {
  c(mean = mean(x), sd = sd(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)))
})

# Create results dataframe
env_only_results <- data.frame(
  lab_id = empirical_data$lab_id,
  study_area = empirical_data$study_area,
  stevenson_age = empirical_data$date_bp,
  stevenson_sd = empirical_data$sd,
  post_mean = env_only_age_summary["mean", ],
  post_sd = env_only_age_summary["sd", ],
  q025 = env_only_age_summary["q025", ],
  q975 = env_only_age_summary["q975", ]
) %>%
  mutate(
    ci_width = q975 - q025,
    age_diff = post_mean - stevenson_age,
    exceeds_claimed = ci_width > 60
  )

cat("Environment-Only Counterfactual Results:\n")
cat(sprintf("  Median posterior SD: +/-%.1f years\n", median(env_only_results$post_sd)))
cat(sprintf("  Median 95%% CI width: +/-%.1f years\n", median(env_only_results$ci_width)/2))
cat(sprintf("  Compare to Stevenson's claim: +/-30 years\n"))
cat(sprintf("  Ratio: %.1fx claimed precision\n", (median(env_only_results$ci_width)/2)/30))
cat(sprintf("  Samples exceeding +/-30 year claim: %.1f%%\n\n",
            100 * mean(env_only_results$exceeds_claimed)))

# Load Tier 1 results for comparison
tier1_results <- read_csv(here("output/tables/tier1_results.csv"),
                          show_col_types = FALSE)

cat("Comparison with Tier 1 (composition only, environment fixed):\n")
cat(sprintf("  Tier 1 (comp only): +/-%.1f years\n", median(tier1_results$ci_width)/2))
cat(sprintf("  This (env only):    +/-%.1f years\n", median(env_only_results$ci_width)/2))
cat(sprintf("  Difference: %.1f years\n\n",
            median(env_only_results$ci_width)/2 - median(tier1_results$ci_width)/2))

# Convergence diagnostics
env_only_summary <- summary(env_only_fit)$summary
env_only_rhat_max <- max(env_only_summary[, "Rhat"], na.rm = TRUE)
env_only_ess_min <- min(env_only_summary[, "n_eff"], na.rm = TRUE)

cat(sprintf("  Convergence: max Rhat = %.3f, min ESS = %.0f\n\n", env_only_rhat_max, env_only_ess_min))

# ===== SAVE RESULTS =====
cat("Saving results\n")

saveRDS(env_only_fit, here("output/fits/counterfactual_env_only_fit.rds"))
cat("  Saved Stan fit object\n")

write_csv(env_only_results, here("output/tables/counterfactual_env_only_results.csv"))
cat("  Saved results table\n")

env_only_summary_stats <- data.frame(
  model = "Environment Only (Counterfactual)",
  H2Ot_fixed = H2Ot_fixed,
  Ea_fixed = Ea_fixed,
  median_sd = median(env_only_results$post_sd),
  median_ci_width = median(env_only_results$ci_width),
  pct_exceeding_claimed = 100 * mean(env_only_results$exceeds_claimed),
  rhat_max = env_only_rhat_max,
  ess_min = env_only_ess_min,
  n_samples = nrow(empirical_data)
)
write_csv(env_only_summary_stats, here("output/tables/counterfactual_env_only_summary.csv"))
cat("  Saved summary statistics\n")

# ===== DIAGNOSTIC PLOTS =====
cat("\nCreating diagnostic plots...\n")

pdf(here("output/figures/counterfactual_env_only_diagnostics.pdf"), width = 12, height = 8)

# 1. Comparison with Tier 1
comparison_long <- bind_rows(
  tier1_results %>% mutate(model = "Tier 1\n(Comp. Only)", ci_half = ci_width/2),
  env_only_results %>% mutate(model = "Counterfactual\n(Env. Only)", ci_half = ci_width/2)
)

p1 <- ggplot(comparison_long, aes(x = ci_half, fill = model)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 30, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 30, y = Inf, label = "Stevenson's +/-30 years",
           vjust = 2, hjust = -0.1, color = "red") +
  labs(title = "Counterfactual Comparison: Composition vs Environment",
       subtitle = "Which uncertainty source contributes more?",
       x = "95% CI Half-Width (years)", y = "Count") +
  scale_fill_manual(values = c("Tier 1\n(Comp. Only)" = "steelblue",
                               "Counterfactual\n(Env. Only)" = "darkorange")) +
  theme_minimal()
print(p1)

# 2. CI width by study area
p2 <- ggplot(env_only_results, aes(x = study_area, y = ci_width/2)) +
  geom_boxplot(fill = "darkorange", alpha = 0.7) +
  geom_hline(yintercept = 30, color = "red", linetype = "dashed") +
  labs(title = "Environment-Only: CI Half-Width by Study Area",
       subtitle = "Composition fixed at population means",
       x = "Study Area", y = "95% CI Half-Width (years)") +
  theme_minimal()
print(p2)

# 3. Posterior ages vs Stevenson
p3 <- ggplot(env_only_results, aes(x = stevenson_age, y = post_mean)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Environment-Only: Posterior Ages vs Stevenson's Reported Ages",
       x = "Stevenson's Reported Age (BP)",
       y = "Environment-Only Posterior Mean Age (BP)") +
  theme_minimal()
print(p3)

# 4. EHT posteriors by study area
eht_samples <- rstan::extract(env_only_fit, "EHT_mean")[[1]]
eht_df <- data.frame(
  SA1 = eht_samples[, 1],
  SA2 = eht_samples[, 2],
  SA3 = eht_samples[, 3]
) %>% pivot_longer(everything(), names_to = "study_area", values_to = "EHT")

p4 <- ggplot(eht_df, aes(x = EHT, fill = study_area)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = c(21.8, 22.4, 23.3), linetype = "dashed") +
  labs(title = "Environment-Only: Posterior EHT Distributions",
       subtitle = "Dashed lines: Stevenson's study-area constants",
       x = "Effective Hydration Temperature (C)", y = "Density") +
  theme_minimal()
print(p4)

dev.off()
cat("  Saved diagnostic plots\n")

# ===== SUMMARY =====
cat("\nKey Findings:\n")
cat(sprintf("  Environment-only uncertainty: +/-%.1f years\n",
            median(env_only_results$ci_width)/2))
cat(sprintf("  Composition-only uncertainty: +/-%.1f years (Tier 1)\n",
            median(tier1_results$ci_width)/2))

if (median(env_only_results$ci_width) < median(tier1_results$ci_width)) {
  cat("\n  COMPOSITION DOMINATES: Even with perfect environmental knowledge,\n")
  cat("    compositional uncertainty produces larger age uncertainty.\n\n")
} else {
  cat("\n  ENVIRONMENT DOMINATES: Even with perfect compositional knowledge,\n")
  cat("    environmental uncertainty produces larger age uncertainty.\n\n")
}
