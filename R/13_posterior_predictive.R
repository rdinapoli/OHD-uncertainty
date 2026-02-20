# 13_posterior_predictive.R - Posterior predictive checks
# Purpose: Generates posterior predictive check plots for model validation
# Inputs: output/fits/tier3_fit.rds, data/stevenson_2015.csv
# Outputs: output/figures/figure_s10_ppc.png
# Runtime: ~1 minute

library(tidyverse)
library(rstan)
library(bayesplot)
library(patchwork)
library(here)

# Load empirical data
empirical_data <- read_csv(
  here("data/stevenson_2015.csv"),
  show_col_types = FALSE
)
y_obs <- empirical_data$H2Ome_pct

# Output directory
output_dir <- here("output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Function to extract H2Ome_rep and generate PPC plots
generate_ppc <- function(fit_path, tier_name) {
  cat(sprintf("Loading %s fit...\n", tier_name))

  if (!file.exists(fit_path)) {
    cat(sprintf("  Fit file not found: %s. Skipping.\n\n", fit_path))
    return(NULL)
  }

  fit <- readRDS(fit_path)

  # Extract replicated observations (H2Ome_rep if available, otherwise generate from H2Ome_pred)
  fit_pars <- fit@sim$pars_oi
  if ("H2Ome_rep" %in% fit_pars) {
    y_rep <- rstan::extract(fit, "H2Ome_rep")[[1]]
  } else {
    cat("  H2Ome_rep not in saved fit; generating from H2Ome_pred + N(0, 0.03)\n")
    y_pred <- rstan::extract(fit, "H2Ome_pred")[[1]]
    y_rep <- y_pred + matrix(rnorm(length(y_pred), 0, 0.03), nrow = nrow(y_pred))
  }
  cat(sprintf("  Extracted %d posterior draws x %d samples\n", nrow(y_rep), ncol(y_rep)))

  # Subsample for plotting (bayesplot recommends ~200 draws)
  n_draws <- min(200, nrow(y_rep))
  draw_idx <- sample(nrow(y_rep), n_draws)
  y_rep_sub <- y_rep[draw_idx, ]

  # PPC density overlay
  p1 <- ppc_dens_overlay(y_obs, y_rep_sub) +
    labs(x = "H2Ome (wt%)")

  # PPC stat (mean)
  stat_mean <- apply(y_rep, 1, mean)
  p2 <- ppc_stat(y_obs, y_rep, stat = "mean") +
    labs(x = "Mean H2Ome (wt%)")

  # PPC stat (SD)
  p3 <- ppc_stat(y_obs, y_rep, stat = "sd") +
    labs(x = "SD of H2Ome (wt%)")

  # PPC intervals
  p4 <- ppc_intervals(y_obs, y_rep_sub, x = seq_along(y_obs), prob = 0.5, prob_outer = 0.95) +
    labs(x = "Sample index", y = "H2Ome (wt%)") +
    theme(axis.text.x = element_blank())

  list(density = p1, stat_mean = p2, stat_sd = p3, intervals = p4)
}

# Generate PPC for all three tiers
tier1_path <- here("output/fits/tier1_fit.rds")
tier2_path <- here("output/fits/tier2_fit.rds")
tier3_path <- here("output/fits/tier3_fit.rds")

ppc_tier1 <- generate_ppc(tier1_path, "Tier 1 (Composition)")
ppc_tier2 <- generate_ppc(tier2_path, "Tier 2 (+ Environment)")
ppc_tier3 <- generate_ppc(tier3_path, "Tier 3 (Complete)")

# Save composite PNG: 3x4 grid (rows = tiers, columns = plot types)
all_tiers <- list(
  "Tier 1: Composition" = ppc_tier1,
  "Tier 2: + Environment" = ppc_tier2,
  "Tier 3: Complete" = ppc_tier3
)
all_tiers <- all_tiers[!sapply(all_tiers, is.null)]

if (length(all_tiers) > 0) {
  col_names <- c("Density overlay", "Mean statistic", "SD statistic", "Per-sample intervals")
  plot_keys <- c("density", "stat_mean", "stat_sd", "intervals")

  # Build grid row by row
  rows <- list()
  for (i in seq_along(all_tiers)) {
    tier_name <- names(all_tiers)[i]
    tier_plots <- all_tiers[[i]]
    row_plots <- lapply(plot_keys, function(k) tier_plots[[k]] + theme(legend.position = "none"))
    # Add tier label to leftmost plot
    row_plots[[1]] <- row_plots[[1]] +
      labs(tag = tier_name) +
      theme(plot.tag = element_text(size = 10, face = "bold"),
            plot.tag.position = "left")
    rows[[i]] <- wrap_plots(row_plots, nrow = 1)
  }

  # Add column headers as a top annotation row
  header_plots <- lapply(col_names, function(lbl) {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = lbl, fontface = "bold", size = 4) +
      theme_void()
  })
  header_row <- wrap_plots(header_plots, nrow = 1)

  composite <- wrap_plots(c(list(header_row), rows), ncol = 1,
                          heights = c(0.06, rep(1, length(rows))))
  ggsave(file.path(output_dir, "figure_s10_ppc.png"), composite,
         width = 20, height = 15, dpi = 300)
}

cat("Posterior predictive check plots saved.\n")
cat(sprintf("  %s\n", file.path(output_dir, "figure_s10_ppc.png")))
