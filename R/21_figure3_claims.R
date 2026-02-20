# 21_figure3_claims.R - Stevenson's published OHD chronology (Figure 3)
# Purpose: Plots 242 published dates with ±30-year error bars and SPDs by study area
# Inputs: data/stevenson_2015.csv
# Outputs: output/figures/figure3_claims.png
# Runtime: ~1 minute

library(tidyverse)
library(patchwork)
library(here)

base_dir <- here::here()
output_dir <- file.path(base_dir, "output/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data

dat <- read_csv(file.path(base_dir, "data/stevenson_2015.csv"),
                show_col_types = FALSE)

sa_labels <- c(
  "SA1" = "SA1 (n = 127)",
  "SA2" = "SA2 (n = 50)",
  "SA3" = "SA3 (n = 65)"
)

dat <- dat %>%
  mutate(sa_label = sa_labels[study_area],
         sa_label = factor(sa_label, levels = sa_labels))

cat(sprintf("Loaded %d samples\n", nrow(dat)))
cat(sprintf("Date range: AD %.0f to AD %.0f\n", min(dat$date_ad), max(dat$date_ad)))

# Timeline
colonization <- 1200
contact <- 1722

# Shared theme

theme_pub <- theme_bw(base_size = 9) +
  theme(
    text = element_text(family = "sans"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

# Panel A: Dot plot with +/-30 yr error bars

cat("Generating Panel A: Date estimates with claimed uncertainty\n")

# Sort within each SA by date for cleaner visual
dat <- dat %>%
  group_by(study_area) %>%
  mutate(rank = rank(date_ad, ties.method = "first")) %>%
  ungroup()

p_a <- ggplot(dat, aes(x = date_ad, y = rank)) +
  # Colonization and contact shading
  annotate("rect", xmin = -Inf, xmax = colonization,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.5) +
  annotate("rect", xmin = contact, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.5) +
  # Error bars (+/-30 yr)
  geom_errorbarh(aes(xmin = date_ad - 30, xmax = date_ad + 30),
                 height = 0, linewidth = 0.2, color = "grey50", alpha = 0.6) +
  # Points
  geom_point(size = 0.4, color = "grey20", alpha = 0.8) +
  # Contact line
  geom_vline(xintercept = contact, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  # Colonization line
  geom_vline(xintercept = colonization, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  facet_grid(sa_label ~ ., scales = "free_y",
             space = "free_y") +
  scale_x_continuous(breaks = seq(1200, 1800, by = 100),
                     limits = c(1100, 1800)) +
  labs(
    x = "Calculated date (AD)",
    y = "Sample (ranked by date)"
  ) +
  theme_pub +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

# Panel B: Summed Probability Distributions (SPDs)

cat("Generating Panel B: Summed probability distributions\n")

# Create SPD: sum of normal densities, each with sd = 30
time_grid <- seq(1100, 1800, by = 1)

spd_by_sa <- dat %>%
  group_by(study_area, sa_label) %>%
  summarize(
    spd = list({
      densities <- sapply(date_ad, function(mu) dnorm(time_grid, mu, 30))
      rowSums(densities)
    }),
    .groups = "drop"
  ) %>%
  mutate(data = map2(sa_label, spd, ~tibble(
    year = time_grid,
    density = .y,
    sa_label = .x
  ))) %>%
  select(data) %>%
  unnest(data)

sa_line_labels <- c(
  "SA1 (n = 127)" = "SA1",
  "SA2 (n = 50)" = "SA2",
  "SA3 (n = 65)" = "SA3"
)

p_b <- ggplot() +
  # Colonization and contact shading
  annotate("rect", xmin = -Inf, xmax = colonization,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.5) +
  annotate("rect", xmin = contact, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.5) +
  # Per-SA SPDs as lines
  geom_line(data = spd_by_sa,
            aes(x = year, y = density, linetype = sa_label),
            linewidth = 0.4, color = "grey40") +
  # Contact line
  geom_vline(xintercept = contact, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = colonization, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  # Labels
  annotate("text", x = colonization - 5, y = Inf,
           label = "~AD 1200", size = 1.8,
           hjust = 1, vjust = 1.5, color = "grey40") +
  annotate("text", x = contact + 5, y = Inf,
           label = "AD 1722", size = 1.8,
           hjust = 0, vjust = 1.5, color = "grey40") +
  scale_x_continuous(breaks = seq(1200, 1800, by = 100),
                     limits = c(1100, 1800)) +
  scale_linetype_manual(
    values = c("solid", "longdash", "dotted"),
    labels = sa_line_labels,
    name = NULL
  ) +
  labs(
    x = "Calculated date (AD)",
    y = "Summed probability"
  ) +
  theme_pub +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.12, 0.55),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.8, "lines"),
    legend.key.width = unit(1.2, "lines"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA)
  )

# Compose

cat("Composing figure\n")

fig <- (p_a + ggtitle("A")) /
  (p_b + ggtitle("B")) +
  plot_layout(heights = c(3, 2))

ggsave(file.path(output_dir, "figure3_claims.png"),
       fig, width = 4, height = 7, dpi = 300)
ggsave(file.path(output_dir, "figure3_claims.pdf"),
       fig, width = 4, height = 7, device = cairo_pdf)

cat("Figure 3 saved to:", output_dir, "\n")
