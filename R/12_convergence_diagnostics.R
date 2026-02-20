# 12_convergence_diagnostics.R - MCMC convergence diagnostics
# Purpose: Generates trace plots and pairs plots for Tier 3 model
# Inputs: output/fits/tier3_fit.rds
# Outputs: output/figures/figure_s1_trace.png, output/figures/figure_s2_pairs.png
# Runtime: ~1 minute

library(tidyverse)
library(rstan)
library(bayesplot)
library(here)

options(mc.cores = parallel::detectCores())

# Output directory
out_dir <- here("output/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load saved Tier 3 fit
cat("Loading Tier 3 fit...\n")
tier3_fit <- readRDS(here("output/fits/tier3_fit.rds"))
cat("Fit loaded successfully.\n")

# Extract posterior draws as array (iterations x chains x parameters)
posterior_array <- as.array(tier3_fit)

# --- Identify a representative sample (median-age sample) ---
age_summary <- summary(tier3_fit, pars = "age")$summary
median_age_idx <- which.min(abs(age_summary[, "50%"] - median(age_summary[, "50%"])))
cat(sprintf("Representative sample: index %d (median posterior age = %.0f BP)\n",
            median_age_idx, age_summary[median_age_idx, "50%"]))

# --- Figure S1: Trace plots ---
cat("Generating trace plots...\n")

# Parameters: EHT_mean[1], EHT_sd[1], RH_mean[1], e_RH, and one representative age
trace_params <- c(
  "EHT_mean[1]",
  "EHT_sd[1]",
  "RH_mean[1]",
  "e_RH",
  sprintf("age[%d]", median_age_idx)
)

color_scheme_set("mix-brightblue-gray")

p_trace <- mcmc_trace(posterior_array, pars = trace_params,
                       facet_args = list(ncol = 1, strip.position = "left")) +
  theme(strip.text = element_text(size = 9),
        axis.text = element_text(size = 8))

ggsave(file.path(out_dir, "figure_s1_trace.png"), p_trace,
       width = 10, height = 10, dpi = 300)
cat("Trace plots saved.\n")

# --- Figure S2: Pairs plot ---
cat("Generating pairs plot...\n")

# For the representative sample, show joint posteriors of key parameters
pairs_params <- c(
  sprintf("age[%d]", median_age_idx),
  sprintf("H2Ot[%d]", median_age_idx),
  sprintf("Ea[%d]", median_age_idx),
  sprintf("EHT[%d]", median_age_idx),
  sprintf("RH[%d]", median_age_idx)
)

# Extract number of divergent transitions per chain
np <- nuts_params(tier3_fit)

p_pairs <- mcmc_pairs(posterior_array, pars = pairs_params, np = np,
                       off_diag_args = list(size = 0.3, alpha = 0.15),
                       diag_fun = "dens")

ggsave(file.path(out_dir, "figure_s2_pairs.png"), p_pairs,
       width = 10, height = 10, dpi = 300)
cat("Pairs plot saved.\n")

# Report diagnostics
n_divergent <- sum(rstan::get_num_divergent(tier3_fit))
n_max_treedepth <- sum(rstan::get_num_max_treedepth(tier3_fit))
cat(sprintf("Divergent transitions: %d\n", n_divergent))
cat(sprintf("Max treedepth exceedances: %d\n", n_max_treedepth))

cat("Done.\n")
cat(sprintf("  %s\n", file.path(out_dir, "figure_s1_trace.png")))
cat(sprintf("  %s\n", file.path(out_dir, "figure_s2_pairs.png")))
