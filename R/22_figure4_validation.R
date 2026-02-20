# 22_figure4_validation.R - Phase model comparison and dating uncertainty (Figure 4)
# Purpose: Shows ΔELPD for phase model comparisons and tier uncertainty bars
# Inputs: output/tables/site15_model_comparison.csv, output/tables/ann_model_comparison.csv,
#   output/tables/tier1_summary.csv, output/tables/tier2_summary.csv, output/tables/tier3_summary.csv
# Outputs: output/figures/figure4_validation.png
# Runtime: ~1 minute

library(tidyverse)
library(patchwork)
library(here)

base_dir <- here::here()
output_dir <- file.path(base_dir, "output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Shared style

theme_fig <- theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

col_claimed <- "grey25"
col_actual  <- "firebrick"

# Panel A: ELPD Difference Point-Range Plot

cat("Loading phase model comparison data...\n")

site233 <- read_csv(
  file.path(base_dir, "output/tables/site15_model_comparison.csv"),
  show_col_types = FALSE
)
ann <- read_csv(
  file.path(base_dir, "output/tables/ann_model_comparison.csv"),
  show_col_types = FALSE
)

# Compute ELPD differences for key comparisons
# Convention: ΔELPD = ELPD(comparison) - ELPD(reference)
# Negative = comparison fits worse than reference

compute_delta <- function(df, model_comp, model_ref) {
  elpd_comp <- df$ELPD[grepl(model_comp, df$Model)]
  elpd_ref  <- df$ELPD[grepl(model_ref, df$Model)]
  se_comp   <- df$ELPD_SE[grepl(model_comp, df$Model)]
  se_ref    <- df$ELPD_SE[grepl(model_ref, df$Model)]
  # Approximate SE of difference (conservative)
  se_diff <- sqrt(se_comp^2 + se_ref^2)
  list(delta = elpd_comp - elpd_ref, se = se_diff)
}

# Key comparisons from the manuscript:
comparisons <- tribble(
  ~site, ~label, ~model_comp, ~model_ref,
  # Site 15-233
  "Site 15-233", "OHD: Ordered vs Single", "1C", "1D",
  "Site 15-233", "C14: Ordered vs Single", "1A", "1B",
  "Site 15-233", "Combined: Ordered vs Single", "1E", "1F",
  # Ahu Nau Nau
  "Ahu Nau Nau", "OHD: Ordered vs Single", "2C", "2D",
  "Ahu Nau Nau", "C14: Ordered vs Single", "2A", "2B",
  "Ahu Nau Nau", "Combined vs C14 only\n(ordered models)", "2E", "2A",
)

delta_results <- comparisons %>%
  rowwise() %>%
  mutate(
    df = list(if (grepl("15-233", site)) site233 else ann),
    result = list(compute_delta(df, model_comp, model_ref)),
    delta_elpd = result$delta,
    se = result$se
  ) %>%
  ungroup() %>%
  select(site, label, delta_elpd, se)

# Order for display
delta_results$label <- factor(delta_results$label,
  levels = rev(c(
    "C14: Ordered vs Single",
    "OHD: Ordered vs Single",
    "Combined: Ordered vs Single",
    "Combined vs C14 only\n(ordered models)"
  ))
)

delta_results$site <- factor(delta_results$site,
  levels = c("Site 15-233", "Ahu Nau Nau")
)

# Color by data type
delta_results <- delta_results %>%
  mutate(data_type = case_when(
    grepl("C14:", label) ~ "Radiocarbon",
    grepl("OHD:", label) ~ "OHD",
    grepl("Combined", label) ~ "Combined"
  ))

panel_a <- ggplot(delta_results, aes(x = delta_elpd, y = label, color = data_type)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey70", linewidth = 0.3) +
  geom_pointrange(aes(xmin = delta_elpd - 2 * se, xmax = delta_elpd + 2 * se),
                  size = 0.5, linewidth = 0.7) +
  facet_wrap(~site, ncol = 2) +
  scale_color_manual(
    values = c("Radiocarbon" = "steelblue", "OHD" = col_actual, "Combined" = "grey35"),
    name = NULL
  ) +
  labs(
    x = expression(Delta*ELPD ~ "(\u00b12 SE)"),
    y = NULL,
    tag = "A"
  ) +
  theme_fig +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.15),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 10),
    plot.tag = element_text(face = "bold", size = 12)
  )

# Panel B: Tier Uncertainty Bar Chart

cat("Loading tier summary data...\n")

tier3 <- read_csv(file.path(base_dir, "output/tables/tier3_summary.csv"),
                  show_col_types = FALSE)

tier_data <- bind_rows(
  tibble(tier = "Claimed\n(\u00b130 yr)", uncertainty = 30, type = "claimed"),
  tibble(tier = "Tier 3\n(Complete model)", uncertainty = tier3$median_ci_width / 2, type = "actual")
)

tier_data$tier <- factor(tier_data$tier,
  levels = c("Claimed\n(\u00b130 yr)", "Tier 3\n(Complete model)")
)

panel_b <- ggplot(tier_data, aes(x = tier, y = uncertainty, fill = type)) +
  # 520-year sequence reference (drawn first so bars layer on top only if taller)
  geom_hline(yintercept = 520 / 2, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  # Bars
  geom_col(width = 0.55) +
  # Reference line label -- left-aligned to avoid occlusion by Tier 3 bar
  annotate("text", x = 0.55, y = 520 / 2 + 8,
           label = "Rapa Nui sequence / 2 = \u00b1260 yr",
           size = 2.7, color = "grey30", hjust = 0, vjust = 0) +
  # Value labels
  geom_text(aes(label = paste0("\u00b1", round(uncertainty), " yr")),
            vjust = -0.5, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c("claimed" = "steelblue", "actual" = col_actual),
                    guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.4)) +
  scale_y_continuous(limits = c(0, 340), breaks = seq(0, 300, 50),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Dating Uncertainty\n(\u00b1 years, 95% CI half-width)",
    tag = "B"
  ) +
  theme_fig +
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    axis.title.y = element_text(margin = margin(r = 2))
  )

# Combine and Save

cat("Compositing Figure 4...\n")

fig4 <- panel_a + panel_b +
  plot_layout(widths = c(1.4, 1))

ggsave(file.path(output_dir, "figure4_validation.pdf"), fig4,
       width = 10, height = 4.5, device = cairo_pdf)
ggsave(file.path(output_dir, "figure4_validation.png"), fig4,
       width = 10, height = 4.5, dpi = 300)

cat("Figure 4 saved to:", output_dir, "\n")
