# 24_figure6_classification.R - Classification accuracy degradation (Figure 6)
# Purpose: Shows how dating uncertainty degrades pre/post-contact classification
# Inputs: output/tables/tier3_summary.csv
# Outputs: output/figures/figure6_classification.png
# Runtime: ~1 minute

library(tidyverse)
library(patchwork)
library(here)

output_dir <- here("output/figures")
drafts_dir <- file.path(output_dir, "drafts")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(drafts_dir, recursive = TRUE, showWarnings = FALSE)

CONTACT_YEAR <- 1722

# =============================================================================
# Analytical classification accuracy
# =============================================================================
# For a true year and given SE, the probability of correct classification:
#   pre-contact (true < 1722):  P(obs < 1722) = pnorm(1722, true_year, se)
#   post-contact (true >= 1722): P(obs >= 1722) = 1 - pnorm(1722, true_year, se)

classify_accuracy_analytical <- function(true_year, se) {
  ifelse(true_year < CONTACT_YEAR,
         pnorm(CONTACT_YEAR, true_year, se) * 100,
         (1 - pnorm(CONTACT_YEAR, true_year, se)) * 100)
}

# =============================================================================
# Shared theme
# =============================================================================

theme_fig1 <- theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

# Colors for the two uncertainty levels
col_claimed <- "grey25"
col_actual  <- "firebrick"

# =============================================================================
# Option A: Two-Line Sigmoid
# =============================================================================

cat("Generating Option A: Two-Line Sigmoid...\n")

years_a <- seq(1000, 2000, by = 1)
df_a <- bind_rows(
  tibble(year = years_a, se = 30, accuracy = classify_accuracy_analytical(years_a, 30),
         label = "\u00b130 yr (claimed)"),
  tibble(year = years_a, se = 299, accuracy = classify_accuracy_analytical(years_a, 299),
         label = "\u00b1299 yr (actual)")
)

# For the ribbon between curves
df_ribbon <- tibble(
  year = years_a,
  acc_30  = classify_accuracy_analytical(years_a, 30),
  acc_299 = classify_accuracy_analytical(years_a, 299)
)

# Accuracy at AD 1660 for annotations
acc_1660_30  <- classify_accuracy_analytical(1660, 30)
acc_1660_299 <- classify_accuracy_analytical(1660, 299)

p_a <- ggplot() +
  # Ribbon showing lost resolving power
  geom_ribbon(data = df_ribbon, aes(x = year, ymin = acc_299, ymax = acc_30),
              fill = "firebrick", alpha = 0.12) +
  # Curves
  geom_line(data = df_a, aes(x = year, y = accuracy, linetype = label, color = label),
            linewidth = 0.9) +
  scale_color_manual(values = c("\u00b130 yr (claimed)" = col_claimed,
                                "\u00b1299 yr (actual)" = col_actual),
                     name = "Uncertainty") +
  scale_linetype_manual(values = c("\u00b130 yr (claimed)" = "solid",
                                   "\u00b1299 yr (actual)" = "solid"),
                        name = "Uncertainty") +
  # Contact year
  geom_vline(xintercept = CONTACT_YEAR, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  annotate("text", x = CONTACT_YEAR + 8, y = 97,
           label = "European Contact (AD 1722)",
           color = "grey40", size = 2.8, hjust = 0) +
  # Random guessing
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  annotate("text", x = 1005, y = 52, label = "random guessing",
           color = "grey50", size = 2.5, hjust = 0) +
  # AD 1660 annotation
  geom_vline(xintercept = 1660, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  annotate("segment", x = 1660, xend = 1660,
           y = acc_1660_299, yend = acc_1660_30,
           arrow = arrow(ends = "both", length = unit(0.06, "inches")),
           color = col_actual, linewidth = 0.5) +
  annotate("text", x = 1660 - 10, y = (acc_1660_30 + acc_1660_299) / 2,
           label = paste0(round(acc_1660_30), "% \u2192 ", round(acc_1660_299), "%"),
           color = col_actual, size = 3, hjust = 1, fontface = "bold") +
  annotate("text", x = 1660, y = acc_1660_299 - 3,
           label = "AD 1660\n(SA1 claim)",
           color = "grey40", size = 2.5, hjust = 0.5, vjust = 1) +
  # Labels
  labs(x = "True Year (AD)", y = "Classification Accuracy (%)") +
  scale_y_continuous(limits = c(45, 100), breaks = seq(50, 100, by = 10)) +
  scale_x_continuous(limits = c(1000, 2000), breaks = seq(1000, 2000, by = 200)) +
  theme_fig1 +
  theme(legend.position = c(0.15, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.key.width = unit(1.2, "cm"))

ggsave(file.path(drafts_dir, "figure6_alt_A_sigmoid.pdf"), p_a,
       width = 6, height = 4, device = cairo_pdf)
cat("  Saved drafts/figure6_alt_A_sigmoid.pdf\n")


# =============================================================================
# Option B: Degradation Curve
# =============================================================================

cat("Generating Option B: Degradation Curve...\n")

ses_b <- seq(1, 350, by = 1)
true_years_b <- c(1660, 1600, 1500)
year_labels <- c("1660" = "AD 1660 (SA1 claim)",
                 "1600" = "AD 1600",
                 "1500" = "AD 1500")

df_b <- expand.grid(se = ses_b, true_year = true_years_b) %>%
  as_tibble() %>%
  mutate(
    accuracy = classify_accuracy_analytical(true_year, se),
    year_label = year_labels[as.character(true_year)]
  )

p_b <- ggplot(df_b, aes(x = se, y = accuracy, color = year_label, linetype = year_label)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("AD 1660 (SA1 claim)" = col_actual,
                                "AD 1600" = "grey35",
                                "AD 1500" = "grey60"),
                     name = "True Year") +
  scale_linetype_manual(values = c("AD 1660 (SA1 claim)" = "solid",
                                   "AD 1600" = "longdash",
                                   "AD 1500" = "longdash"),
                        name = "True Year") +
  # Reference lines at +/-30 and +/-299
  geom_vline(xintercept = 30, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  annotate("text", x = 33, y = 52, label = "\u00b130 yr\n(claimed)", color = "grey40",
           size = 2.5, hjust = 0) +
  geom_vline(xintercept = 299, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  annotate("text", x = 302, y = 52, label = "\u00b1299 yr\n(actual)", color = "grey40",
           size = 2.5, hjust = 0) +
  # Random guessing
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  # Labels
  labs(x = "Dating Uncertainty (\u00b1 years)", y = "Classification Accuracy (%)") +
  scale_y_continuous(limits = c(45, 100), breaks = seq(50, 100, by = 10)) +
  scale_x_continuous(limits = c(0, 350), breaks = seq(0, 350, by = 50)) +
  theme_fig1 +
  theme(legend.position = c(0.78, 0.82),
        legend.background = element_rect(fill = "white", color = "grey80"))

ggsave(file.path(drafts_dir, "figure6_alt_B_degradation.pdf"), p_b,
       width = 6, height = 4, device = cairo_pdf)
cat("  Saved drafts/figure6_alt_B_degradation.pdf\n")


# =============================================================================
# Option C: Density Pair
# =============================================================================

cat("Generating Option C: Density Pair...\n")

true_year_c <- 1660
x_range <- seq(800, 2500, by = 1)

make_density_panel <- function(se, panel_label) {
  dens <- dnorm(x_range, mean = true_year_c, sd = se)
  df <- tibble(x = x_range, density = dens)

  # Classification regions
  df <- df %>%
    mutate(region = ifelse(x < CONTACT_YEAR, "Correct\n(pre-contact)", "Misclassified\n(post-contact)"))

  acc <- pnorm(CONTACT_YEAR, true_year_c, se) * 100

  ggplot(df, aes(x = x, y = density, fill = region)) +
    geom_area(alpha = 0.6, position = "identity") +
    geom_line(aes(x = x, y = density), inherit.aes = FALSE, color = "black", linewidth = 0.4) +
    scale_fill_manual(values = c("Correct\n(pre-contact)" = "grey70",
                                 "Misclassified\n(post-contact)" = col_actual),
                      name = NULL) +
    geom_vline(xintercept = CONTACT_YEAR, linetype = "dashed", color = "grey30", linewidth = 0.5) +
    geom_vline(xintercept = true_year_c, linetype = "dotted", color = "grey50", linewidth = 0.4) +
    # Accuracy annotation
    annotate("text", x = true_year_c - se * 0.3, y = max(dens) * 0.85,
             label = paste0(round(acc), "% correct"),
             color = "grey20", size = 3.5, fontface = "bold") +
    # Misclassified annotation
    annotate("text", x = CONTACT_YEAR + se * 0.5,
             y = max(dens) * 0.5,
             label = paste0(round(100 - acc), "%\nmisclassified"),
             color = col_actual, size = 3, fontface = "bold") +
    labs(title = panel_label,
         x = "Year (AD)", y = "Probability Density") +
    scale_x_continuous(limits = c(1000, 2200), breaks = seq(1000, 2200, by = 200)) +
    theme_fig1 +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          plot.title = element_text(size = 10))
}

p_c1 <- make_density_panel(30, paste0("OHD date: AD 1660 \u00b1 30 yr"))
p_c2 <- make_density_panel(299, paste0("OHD date: AD 1660 \u00b1 299 yr"))

p_c <- p_c1 + p_c2 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(drafts_dir, "figure6_alt_C_density_pair.pdf"), p_c,
       width = 7, height = 3.5, device = cairo_pdf)
cat("  Saved drafts/figure6_alt_C_density_pair.pdf\n")


# =============================================================================
# Option D: Composite (C on top, A below)
# =============================================================================

cat("Generating composite figure...\n")

# Rebuild A and C without individual legends for the composite
p_a_comp <- p_a +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold", size = 12))

p_c1_comp <- p_c1 + labs(tag = "A")
p_c2_comp <- p_c2 + labs(tag = "")

p_c_comp <- (p_c1_comp + p_c2_comp) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 12))

p_d <- p_c_comp / p_a_comp +
  plot_layout(heights = c(1, 1.1))

ggsave(file.path(output_dir, "figure6_classification.pdf"), p_d,
       width = 7, height = 6, device = cairo_pdf)
ggsave(file.path(output_dir, "figure6_classification.png"), p_d,
       width = 7, height = 6, dpi = 300)
cat("Saved figure6_classification.pdf and .png\n")
