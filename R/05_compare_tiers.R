# 05_compare_tiers.R - LOO-CV comparison between tier models
# Purpose: Compares predictive performance of Tiers 1-3 using LOO-CV
# Inputs: output/fits/tier1_fit.rds, output/fits/tier2_fit.rds, output/fits/tier3_fit.rds
# Outputs: output/tables/tier_loo_comparison.csv, output/tables/pareto_k_summary.csv
# Runtime: ~3 minutes

library(tidyverse)
library(rstan)
library(loo)
library(here)

options(mc.cores = parallel::detectCores())

cat("LOO-CV comparison between Tiers 1-3\n\n")

dir.create(here("output/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/fits"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output/figures"), recursive = TRUE, showWarnings = FALSE)

# ===== LOAD FITS AND COMPUTE LOO =====

compute_loo_from_fit <- function(fit, tier_name) {
  cat(sprintf("Computing LOO for %s...\n", tier_name))

  log_lik <- loo::extract_log_lik(fit, merge_chains = FALSE)
  r_eff <- loo::relative_eff(exp(log_lik))
  loo_result <- loo::loo(log_lik, r_eff = r_eff)

  cat(sprintf("  %s: ELPD = %.1f (SE = %.1f)\n", tier_name,
              loo_result$estimates["elpd_loo", "Estimate"],
              loo_result$estimates["elpd_loo", "SE"]))

  # Pareto k diagnostics
  k_vals <- loo_result$diagnostics$pareto_k
  cat(sprintf("  Pareto k: max = %.3f, mean = %.3f\n", max(k_vals), mean(k_vals)))
  cat(sprintf("  k > 0.5: %d (%.1f%%); k > 0.7: %d (%.1f%%); k > 1.0: %d (%.1f%%)\n",
              sum(k_vals > 0.5), 100 * mean(k_vals > 0.5),
              sum(k_vals > 0.7), 100 * mean(k_vals > 0.7),
              sum(k_vals > 1.0), 100 * mean(k_vals > 1.0)))
  cat("\n")

  return(loo_result)
}

# Load fits
cat("Loading tier fit objects...\n\n")

tier1_fit <- readRDS(here("output/fits/tier1_fit.rds"))
loo1 <- compute_loo_from_fit(tier1_fit, "Tier 1")
rm(tier1_fit); gc(verbose = FALSE)

tier2_fit <- readRDS(here("output/fits/tier2_fit.rds"))
loo2 <- compute_loo_from_fit(tier2_fit, "Tier 2")
rm(tier2_fit); gc(verbose = FALSE)

# Check if Tier 3 LOO already exists
tier3_loo_path <- here("output/fits/tier3_loo.rds")
if (file.exists(tier3_loo_path)) {
  cat("Loading existing Tier 3 LOO from file...\n")
  loo3 <- readRDS(tier3_loo_path)
  k_vals <- loo3$diagnostics$pareto_k
  cat(sprintf("  Tier 3: ELPD = %.1f (SE = %.1f)\n",
              loo3$estimates["elpd_loo", "Estimate"],
              loo3$estimates["elpd_loo", "SE"]))
  cat(sprintf("  Pareto k: max = %.3f, mean = %.3f\n", max(k_vals), mean(k_vals)))
  cat(sprintf("  k > 0.5: %d (%.1f%%); k > 0.7: %d (%.1f%%); k > 1.0: %d (%.1f%%)\n\n",
              sum(k_vals > 0.5), 100 * mean(k_vals > 0.5),
              sum(k_vals > 0.7), 100 * mean(k_vals > 0.7),
              sum(k_vals > 1.0), 100 * mean(k_vals > 1.0)))
} else {
  tier3_fit <- readRDS(here("output/fits/tier3_fit.rds"))
  loo3 <- compute_loo_from_fit(tier3_fit, "Tier 3")
  rm(tier3_fit); gc(verbose = FALSE)
}

# Save LOO objects
saveRDS(loo1, here("output/fits/tier1_loo.rds"))
saveRDS(loo2, here("output/fits/tier2_loo.rds"))
saveRDS(loo3, here("output/fits/tier3_loo.rds"))

# ===== MODEL COMPARISON =====

cat("Model comparison (loo_compare)\n")

comp <- loo::loo_compare(list("Tier 1" = loo1, "Tier 2" = loo2, "Tier 3" = loo3))
print(comp)

# Save comparison
comp_df <- as.data.frame(comp)
comp_df$tier <- rownames(comp_df)
comp_df <- comp_df %>% select(tier, everything())
write_csv(comp_df, here("output/tables/tier_loo_comparison.csv"))
cat("\nSaved: tier_loo_comparison.csv\n")

# ===== STACKING WEIGHTS =====

cat("\nStacking weights\n")

# Extract pointwise ELPD for stacking
lpd_point <- cbind(
  loo1$pointwise[, "elpd_loo"],
  loo2$pointwise[, "elpd_loo"],
  loo3$pointwise[, "elpd_loo"]
)

weights <- loo::stacking_weights(lpd_point)
cat(sprintf("  Tier 1: %.3f\n", weights[1]))
cat(sprintf("  Tier 2: %.3f\n", weights[2]))
cat(sprintf("  Tier 3: %.3f\n", weights[3]))

# ===== PARETO K SUMMARY TABLE =====

cat("\nPareto k diagnostic summary\n")

pareto_summary <- data.frame(
  tier = c("Tier 1", "Tier 2", "Tier 3"),
  elpd_loo = c(loo1$estimates["elpd_loo", "Estimate"],
               loo2$estimates["elpd_loo", "Estimate"],
               loo3$estimates["elpd_loo", "Estimate"]),
  elpd_se = c(loo1$estimates["elpd_loo", "SE"],
              loo2$estimates["elpd_loo", "SE"],
              loo3$estimates["elpd_loo", "SE"]),
  k_max = c(max(loo1$diagnostics$pareto_k),
            max(loo2$diagnostics$pareto_k),
            max(loo3$diagnostics$pareto_k)),
  k_mean = c(mean(loo1$diagnostics$pareto_k),
             mean(loo2$diagnostics$pareto_k),
             mean(loo3$diagnostics$pareto_k)),
  n_k_above_05 = c(sum(loo1$diagnostics$pareto_k > 0.5),
                    sum(loo2$diagnostics$pareto_k > 0.5),
                    sum(loo3$diagnostics$pareto_k > 0.5)),
  n_k_above_07 = c(sum(loo1$diagnostics$pareto_k > 0.7),
                    sum(loo2$diagnostics$pareto_k > 0.7),
                    sum(loo3$diagnostics$pareto_k > 0.7)),
  n_k_above_10 = c(sum(loo1$diagnostics$pareto_k > 1.0),
                    sum(loo2$diagnostics$pareto_k > 1.0),
                    sum(loo3$diagnostics$pareto_k > 1.0)),
  stacking_weight = as.numeric(weights)
)

# Add delta-ELPD relative to best
best_elpd <- max(pareto_summary$elpd_loo)
pareto_summary$delta_elpd <- pareto_summary$elpd_loo - best_elpd

write_csv(pareto_summary, here("output/tables/pareto_k_summary.csv"))
cat("Saved: pareto_k_summary.csv\n\n")

# Print markdown table
cat("| Tier | ELPD | SE | delta-ELPD | Weight | k_max | k > 0.7 |\n")
cat("|------|------|----|------------|--------|-------|---------|\n")
for (i in seq_len(nrow(pareto_summary))) {
  r <- pareto_summary[i, ]
  cat(sprintf("| %s | %.1f | %.1f | %.1f | %.3f | %.3f | %d |\n",
              r$tier, r$elpd_loo, r$elpd_se, r$delta_elpd,
              r$stacking_weight, r$k_max, r$n_k_above_07))
}

# ===== PARETO K PLOT =====

cat("\nGenerating Pareto k diagnostic plot...\n")

k_data <- data.frame(
  sample = rep(1:242, 3),
  k = c(loo1$diagnostics$pareto_k,
        loo2$diagnostics$pareto_k,
        loo3$diagnostics$pareto_k),
  tier = rep(c("Tier 1", "Tier 2", "Tier 3"), each = 242)
)

p <- ggplot(k_data, aes(x = sample, y = k, color = tier)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  annotate("text", x = 250, y = 0.52, label = "k = 0.5", hjust = 1, color = "orange", size = 3) +
  annotate("text", x = 250, y = 0.72, label = "k = 0.7", hjust = 1, color = "red", size = 3) +
  facet_wrap(~ tier, ncol = 1) +
  scale_color_manual(values = c("Tier 1" = "#1b9e77", "Tier 2" = "#d95f02", "Tier 3" = "#7570b3")) +
  labs(title = "Pareto k Diagnostics by Tier",
       subtitle = "k < 0.7 indicates reliable LOO-CV estimates",
       x = "Sample index", y = "Pareto k") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(here("output/figures/pareto_k_diagnostics.pdf"), p, width = 10, height = 8)
ggsave(here("output/figures/pareto_k_diagnostics.png"), p, width = 10, height = 8, dpi = 150)

cat("LOO-CV tier comparison complete\n")
