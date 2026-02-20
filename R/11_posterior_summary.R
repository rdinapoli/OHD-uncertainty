# 11_posterior_summary.R - Posterior parameter summaries and correlations
# Purpose: Extracts posterior parameter estimates, contraction ratios, and correlation structure
# Inputs: output/fits/tier3_fit.rds
# Outputs: output/tables/posterior_parameter_summary.csv,
#   output/tables/prior_posterior_contraction.csv,
#   output/tables/posterior_correlation_matrix.csv,
#   output/figures/figure_s5_prior_posterior.png,
#   output/figures/figure_s6_correlation.png
# Runtime: ~2 minutes

library(tidyverse)
library(rstan)
library(here)

options(mc.cores = parallel::detectCores())

# ===== HELPER FUNCTIONS =====

summarize_param <- function(draws, param_name, prior_str) {
  data.frame(
    parameter = param_name,
    prior = prior_str,
    post_mean = mean(draws),
    post_sd = sd(draws),
    q025 = unname(quantile(draws, 0.025)),
    q975 = unname(quantile(draws, 0.975)),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# PART 1: Posterior Parameter Summary
# =============================================================================

# --- Load Fits ---------------------------------------------------------------

cat("Loading tier fit objects (this may take a minute)...\n")

tier1_fit <- readRDS(here("output/fits/tier1_fit.rds"))
cat("  Tier 1 loaded\n")

tier2_fit <- readRDS(here("output/fits/tier2_fit.rds"))
cat("  Tier 2 loaded\n")

tier3_fit <- readRDS(here("output/fits/tier3_fit.rds"))
cat("  Tier 3 loaded\n\n")

# --- Tier 1 Parameters ------------------------------------------------------

cat("Extracting Tier 1 posteriors...\n")

t1_draws <- rstan::extract(tier1_fit)

tier1_table <- bind_rows(
  summarize_param(t1_draws$H2Ot_mu, "mu_H2Ot (wt%)", "N(0.18, 0.08)"),
  summarize_param(t1_draws$H2Ot_sigma, "sigma_H2Ot (wt%)", "Exp(10)"),
  summarize_param(t1_draws$Ea_mu, "mu_Ea (J/mol)", "N(84000, 2000)"),
  summarize_param(t1_draws$Ea_sigma, "sigma_Ea (J/mol)", "Exp(0.0005)")
)
tier1_table$tier <- "Tier 1"

# Age summary
t1_ages <- t1_draws$age
tier1_age <- data.frame(
  parameter = "age (years BP)",
  prior = "N(500, 200)",
  post_mean = mean(apply(t1_ages, 2, mean)),
  post_sd = mean(apply(t1_ages, 2, sd)),
  q025 = mean(apply(t1_ages, 2, quantile, 0.025)),
  q975 = mean(apply(t1_ages, 2, quantile, 0.975)),
  tier = "Tier 1",
  stringsAsFactors = FALSE
)
tier1_table <- bind_rows(tier1_table, tier1_age)

cat("  Tier 1 done\n")

# --- Tier 2 Parameters ------------------------------------------------------

cat("Extracting Tier 2 posteriors...\n")

t2_draws <- rstan::extract(tier2_fit)

# Composition (same structure as Tier 1)
tier2_table <- bind_rows(
  summarize_param(t2_draws$H2Ot_mu, "mu_H2Ot (wt%)", "N(0.18, 0.08)"),
  summarize_param(t2_draws$H2Ot_sigma, "sigma_H2Ot (wt%)", "Exp(10)"),
  summarize_param(t2_draws$Ea_mu, "mu_Ea (J/mol)", "N(84000, 2000)"),
  summarize_param(t2_draws$Ea_sigma, "sigma_Ea (J/mol)", "Exp(0.0005)")
)

# Environmental parameters (hierarchical by study area)
sa_labels <- c("SA1", "SA2", "SA3")
for (k in 1:3) {
  tier2_table <- bind_rows(
    tier2_table,
    summarize_param(t2_draws$EHT_mean[, k],
                    sprintf("mu_EHT[%s] (C)", sa_labels[k]),
                    "N(22.5, 2.5)"),
    summarize_param(t2_draws$EHT_sd[, k],
                    sprintf("sigma_EHT[%s] (C)", sa_labels[k]),
                    "Exp(0.33)")
  )
}
for (k in 1:3) {
  tier2_table <- bind_rows(
    tier2_table,
    summarize_param(t2_draws$RH_mean[, k],
                    sprintf("mu_RH[%s]", sa_labels[k]),
                    "N(0.85, 0.1)"),
    summarize_param(t2_draws$RH_sd[, k],
                    sprintf("sigma_RH[%s]", sa_labels[k]),
                    "Exp(5)")
  )
}

# RH sensitivity coefficient
tier2_table <- bind_rows(
  tier2_table,
  summarize_param(t2_draws$e_RH, "e_RH", "N(2.15, 0.5)")
)

# Age summary
t2_ages <- t2_draws$age
tier2_age <- data.frame(
  parameter = "age (years BP)",
  prior = "N(500, 200)",
  post_mean = mean(apply(t2_ages, 2, mean)),
  post_sd = mean(apply(t2_ages, 2, sd)),
  q025 = mean(apply(t2_ages, 2, quantile, 0.025)),
  q975 = mean(apply(t2_ages, 2, quantile, 0.975)),
  stringsAsFactors = FALSE
)
tier2_table <- bind_rows(tier2_table, tier2_age)
tier2_table$tier <- "Tier 2"

cat("  Tier 2 done\n")

# --- Tier 3 Parameters ------------------------------------------------------

cat("Extracting Tier 3 posteriors...\n")

t3_draws <- rstan::extract(tier3_fit)

# Same structure as Tier 2 plus temporal climate
tier3_table <- bind_rows(
  summarize_param(t3_draws$H2Ot_mu, "mu_H2Ot (wt%)", "N(0.18, 0.08)"),
  summarize_param(t3_draws$H2Ot_sigma, "sigma_H2Ot (wt%)", "Exp(10)"),
  summarize_param(t3_draws$Ea_mu, "mu_Ea (J/mol)", "N(84000, 2000)"),
  summarize_param(t3_draws$Ea_sigma, "sigma_Ea (J/mol)", "Exp(0.0005)")
)

for (k in 1:3) {
  tier3_table <- bind_rows(
    tier3_table,
    summarize_param(t3_draws$EHT_mean[, k],
                    sprintf("mu_EHT[%s] (C)", sa_labels[k]),
                    "N(22.5, 2.5)"),
    summarize_param(t3_draws$EHT_sd[, k],
                    sprintf("sigma_EHT[%s] (C)", sa_labels[k]),
                    "Exp(0.33)")
  )
}
for (k in 1:3) {
  tier3_table <- bind_rows(
    tier3_table,
    summarize_param(t3_draws$RH_mean[, k],
                    sprintf("mu_RH[%s]", sa_labels[k]),
                    "N(0.85, 0.1)"),
    summarize_param(t3_draws$RH_sd[, k],
                    sprintf("sigma_RH[%s]", sa_labels[k]),
                    "Exp(5)")
  )
}

tier3_table <- bind_rows(
  tier3_table,
  summarize_param(t3_draws$e_RH, "e_RH", "N(2.15, 0.5)")
)

# Age summary
t3_ages <- t3_draws$age
tier3_age <- data.frame(
  parameter = "age (years BP)",
  prior = "N(500, 200)",
  post_mean = mean(apply(t3_ages, 2, mean)),
  post_sd = mean(apply(t3_ages, 2, sd)),
  q025 = mean(apply(t3_ages, 2, quantile, 0.025)),
  q975 = mean(apply(t3_ages, 2, quantile, 0.975)),
  stringsAsFactors = FALSE
)
tier3_table <- bind_rows(tier3_table, tier3_age)
tier3_table$tier <- "Tier 3"

cat("  Tier 3 done\n\n")

# --- Combine and Save -------------------------------------------------------

all_params <- bind_rows(tier1_table, tier2_table, tier3_table)

# Round for presentation
all_params <- all_params %>%
  mutate(
    post_mean = case_when(
      grepl("Ea|J/mol", parameter) ~ round(post_mean, 0),
      grepl("age", parameter) ~ round(post_mean, 1),
      TRUE ~ round(post_mean, 4)
    ),
    post_sd = case_when(
      grepl("Ea|J/mol", parameter) ~ round(post_sd, 0),
      grepl("age", parameter) ~ round(post_sd, 1),
      TRUE ~ round(post_sd, 4)
    ),
    q025 = case_when(
      grepl("Ea|J/mol", parameter) ~ round(q025, 0),
      grepl("age", parameter) ~ round(q025, 1),
      TRUE ~ round(q025, 4)
    ),
    q975 = case_when(
      grepl("Ea|J/mol", parameter) ~ round(q975, 0),
      grepl("age", parameter) ~ round(q975, 1),
      TRUE ~ round(q975, 4)
    )
  )

write_csv(all_params, here("output/tables/posterior_parameter_summary.csv"))
cat("Saved: posterior_parameter_summary.csv\n")

# Print markdown table for SI
cat("\nTier 3 posterior summary (complete model):\n\n")

t3_out <- all_params %>% filter(tier == "Tier 3", parameter != "age (years BP)")

cat("| Parameter | Prior | Posterior Mean | Posterior SD | 95% CI |\n")
cat("|-----------|-------|---------------|-------------|--------|\n")
for (i in seq_len(nrow(t3_out))) {
  row <- t3_out[i, ]
  ci_str <- sprintf("[%s, %s]", format(row$q025, nsmall = ifelse(grepl("Ea|J/mol", row$parameter), 0, 4)),
                    format(row$q975, nsmall = ifelse(grepl("Ea|J/mol", row$parameter), 0, 4)))
  cat(sprintf("| %s | %s | %s | %s | %s |\n",
              row$parameter, row$prior,
              format(row$post_mean, nsmall = ifelse(grepl("Ea|J/mol", row$parameter), 0, 4)),
              format(row$post_sd, nsmall = ifelse(grepl("Ea|J/mol", row$parameter), 0, 4)),
              ci_str))
}

# --- Prior-Posterior Comparison Figure (Figure S5) ---------------------------

cat("\nGenerating prior-posterior comparison (Figure S5)...\n")

# Helper: plot prior vs posterior density
plot_prior_post <- function(draws, prior_fn, param_label, xlim_range = NULL) {
  if (is.null(xlim_range)) {
    xlim_range <- range(draws)
    # Extend to show prior
    xlim_range <- c(xlim_range[1] - 0.3 * diff(xlim_range),
                    xlim_range[2] + 0.3 * diff(xlim_range))
  }
  d <- density(draws)
  x_seq <- seq(xlim_range[1], xlim_range[2], length.out = 500)
  prior_vals <- prior_fn(x_seq)

  y_max <- max(c(d$y, prior_vals), na.rm = TRUE) * 1.1

  plot(d$x, d$y, type = "l", col = "steelblue", lwd = 2,
       xlim = xlim_range, ylim = c(0, y_max),
       xlab = param_label, ylab = "Density", main = param_label)
  lines(x_seq, prior_vals, col = "firebrick", lwd = 2, lty = 2)
  legend("topright", legend = c("Posterior", "Prior"),
         col = c("steelblue", "firebrick"), lwd = 2, lty = c(1, 2), cex = 0.8)
}

# Save as PNG
png(here("output/figures/figure_s5_prior_posterior.png"), width = 1400, height = 1000, res = 100)

par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

plot_prior_post(t3_draws$H2Ot_mu,
                function(x) dnorm(x, 0.18, 0.08),
                expression(mu[H[2]*O[t]] ~ "(wt%)"),
                xlim_range = c(-0.05, 0.45))

plot_prior_post(t3_draws$H2Ot_sigma,
                function(x) dexp(x, 10),
                expression(sigma[H[2]*O[t]] ~ "(wt%)"),
                xlim_range = c(0, 0.5))

plot_prior_post(t3_draws$Ea_mu,
                function(x) dnorm(x, 84000, 2000),
                expression(mu[E[a]] ~ "(J/mol)"),
                xlim_range = c(78000, 90000))

plot_prior_post(t3_draws$Ea_sigma,
                function(x) dexp(x, 0.0005),
                expression(sigma[E[a]] ~ "(J/mol)"),
                xlim_range = c(0, 8000))

plot_prior_post(t3_draws$EHT_mean[, 1],
                function(x) dnorm(x, 22.5, 2.5),
                expression(mu[EHT] * "[SA1] (C)"),
                xlim_range = c(15, 30))

plot_prior_post(t3_draws$EHT_mean[, 2],
                function(x) dnorm(x, 22.5, 2.5),
                expression(mu[EHT] * "[SA2] (C)"),
                xlim_range = c(15, 30))

plot_prior_post(t3_draws$RH_mean[, 1],
                function(x) dnorm(x, 0.85, 0.1),
                expression(mu[RH] * "[SA1]"),
                xlim_range = c(0.5, 1.1))

plot_prior_post(t3_draws$e_RH,
                function(x) dnorm(x, 2.15, 0.5),
                expression(e[RH]),
                xlim_range = c(0.5, 4.0))

avg_age_per_draw <- rowMeans(t3_draws$age)
plot_prior_post(avg_age_per_draw,
                function(x) dnorm(x, 500, 200),
                "Mean age (years BP)",
                xlim_range = c(0, 1000))

dev.off()

cat("Saved: figure_s5_prior_posterior.png\n")

# --- Prior-Posterior Contraction Ratios --------------------------------------

cat("\nComputing prior-posterior contraction ratios...\n")

# Contraction ratio = 1 - (posterior SD / prior SD)
# 0 = no contraction (no learning), 1 = complete contraction (fully determined)
contraction <- data.frame(
  parameter = c("mu_H2Ot", "sigma_H2Ot", "mu_Ea", "sigma_Ea",
                "mu_EHT[SA1]", "mu_EHT[SA2]", "mu_EHT[SA3]",
                "mu_RH[SA1]", "mu_RH[SA2]", "mu_RH[SA3]",
                "e_RH"),
  prior_sd = c(0.08, 1/10, 2000, 1/0.0005,
               2.5, 2.5, 2.5,
               0.1, 0.1, 0.1,
               0.5),
  post_sd = c(sd(t3_draws$H2Ot_mu), sd(t3_draws$H2Ot_sigma),
              sd(t3_draws$Ea_mu), sd(t3_draws$Ea_sigma),
              sd(t3_draws$EHT_mean[, 1]), sd(t3_draws$EHT_mean[, 2]),
              sd(t3_draws$EHT_mean[, 3]),
              sd(t3_draws$RH_mean[, 1]), sd(t3_draws$RH_mean[, 2]),
              sd(t3_draws$RH_mean[, 3]),
              sd(t3_draws$e_RH)),
  stringsAsFactors = FALSE
)
contraction$ratio <- 1 - contraction$post_sd / contraction$prior_sd
contraction$ratio <- round(contraction$ratio, 3)

cat("| Parameter | Prior SD | Posterior SD | Contraction |\n")
cat("|-----------|----------|-------------|-------------|\n")
for (i in seq_len(nrow(contraction))) {
  r <- contraction[i, ]
  cat(sprintf("| %s | %.4g | %.4g | %.1f%% |\n",
              r$parameter, r$prior_sd, r$post_sd, r$ratio * 100))
}

write_csv(contraction, here("output/tables/prior_posterior_contraction.csv"))
cat("\nSaved: prior_posterior_contraction.csv\n")

# =============================================================================
# PART 2: Posterior Correlation Matrix
# =============================================================================

cat("\nComputing posterior parameter correlations...\n\n")

# For each posterior draw, we have 242 sample-level parameters.
# To show the population-level correlation structure, compute the
# cross-sample correlation for each draw, then report median correlation.

N <- ncol(t3_draws$age)  # 242 samples
n_draws <- nrow(t3_draws$age)

# Subsample draws for efficiency
set.seed(42)
draw_idx <- sort(sample(n_draws, min(1000, n_draws)))

cor_matrices <- array(NA, dim = c(5, 5, length(draw_idx)))

for (d_idx in seq_along(draw_idx)) {
  d <- draw_idx[d_idx]
  mat <- cbind(
    age = t3_draws$age[d, ],
    H2Ot = t3_draws$H2Ot[d, ],
    Ea = t3_draws$Ea[d, ],
    EHT = t3_draws$EHT[d, ],
    RH = t3_draws$RH[d, ]
  )
  cor_matrices[, , d_idx] <- cor(mat)
}

# Median correlation across draws
median_cor <- apply(cor_matrices, c(1, 2), median)
rownames(median_cor) <- colnames(median_cor) <- c("age", "H2Ot", "Ea", "EHT", "RH")

# IQR for uncertainty
q25_cor <- apply(cor_matrices, c(1, 2), quantile, 0.25)
q75_cor <- apply(cor_matrices, c(1, 2), quantile, 0.75)

cat("Median posterior correlation matrix:\n\n")
print(round(median_cor, 3))

# Save correlation matrix as CSV
cor_df <- as.data.frame(round(median_cor, 4))
cor_df$parameter <- rownames(cor_df)
cor_df <- cor_df %>% select(parameter, everything())
write_csv(cor_df, here("output/tables/posterior_correlation_matrix.csv"))
cat("\nSaved: posterior_correlation_matrix.csv\n")

# --- Correlation Heatmap (Figure S6) ----------------------------------------

cat("\nGenerating correlation heatmap (Figure S6)...\n")

cor_long <- expand.grid(
  row = c("age", "H2Ot", "Ea", "EHT", "RH"),
  col = c("age", "H2Ot", "Ea", "EHT", "RH"),
  stringsAsFactors = FALSE
)
cor_long$value <- as.vector(median_cor)
cor_long$iqr_low <- as.vector(q25_cor)
cor_long$iqr_high <- as.vector(q75_cor)

# Format labels for heatmap
cor_long$label <- sprintf("%.2f", cor_long$value)
# Don't show diagonal
cor_long$label[cor_long$row == cor_long$col] <- ""

# Order factors
param_order <- c("age", "H2Ot", "Ea", "EHT", "RH")
cor_long$row <- factor(cor_long$row, levels = rev(param_order))
cor_long$col <- factor(cor_long$col, levels = param_order)

p <- ggplot(cor_long, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 4.5, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Correlation") +
  labs(title = "Posterior Correlation Matrix (Tier 3)",
       subtitle = "Median across posterior draws; cross-sample correlations reveal parameter confounding",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13),
    axis.text.y = element_text(size = 13),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave(here("output/figures/figure_s6_correlation.png"), p, width = 8, height = 7, dpi = 150)

cat("Saved: figure_s6_correlation.png\n")

# --- Pairs Plot (subset) ----------------------------------------------------

cat("\nGenerating pairs plot for key parameter pairs...\n")

# For a representative draw, show scatter of age vs each parameter
set.seed(42)
rep_draw <- sample(n_draws, 1)

pairs_data <- data.frame(
  age = t3_draws$age[rep_draw, ],
  H2Ot = t3_draws$H2Ot[rep_draw, ],
  Ea = t3_draws$Ea[rep_draw, ],
  EHT = t3_draws$EHT[rep_draw, ],
  RH = t3_draws$RH[rep_draw, ]
)

# More informative: show the across-draw correlation for a single sample
# Pick a representative sample (median-age sample)
median_ages <- apply(t3_draws$age, 2, median)
rep_sample <- which.min(abs(median_ages - median(median_ages)))

across_draw_data <- data.frame(
  age = t3_draws$age[draw_idx, rep_sample],
  H2Ot = t3_draws$H2Ot[draw_idx, rep_sample],
  Ea = t3_draws$Ea[draw_idx, rep_sample],
  EHT = t3_draws$EHT[draw_idx, rep_sample],
  RH = t3_draws$RH[draw_idx, rep_sample]
)

pdf(here("output/figures/posterior_pairs_plot.pdf"), width = 12, height = 12)
pairs(across_draw_data,
      labels = c("age (BP)", expression(H[2]*O[t]), expression(E[a]),
                 "EHT (C)", "RH"),
      main = sprintf("Posterior Pairs Plot (sample %d, representative)", rep_sample),
      pch = ".", col = adjustcolor("steelblue", alpha.f = 0.3),
      gap = 0.3)
dev.off()

cat("Saved: posterior_pairs_plot.pdf\n")

# --- Notable Correlations ---------------------------------------------------

cat("\nNotable correlations (median across draws):\n")
for (i in 1:4) {
  for (j in (i + 1):5) {
    r <- median_cor[i, j]
    if (abs(r) > 0.1) {
      cat(sprintf("  %s ~ %s: r = %.3f [IQR: %.3f, %.3f]\n",
                  rownames(median_cor)[i], colnames(median_cor)[j],
                  r, q25_cor[i, j], q75_cor[i, j]))
    }
  }
}

cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
