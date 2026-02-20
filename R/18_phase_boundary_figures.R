# 18_phase_boundary_figures.R - OxCal-style phase boundary figures
# Purpose: Generates posterior density plots for all 12 phase models (Figures S8a-f, S9a-f)
# Inputs: output/fits/site15_*.rds, output/fits/ann_*.rds, R/utils.R
# Outputs: output/figures/figure_s8a_site15_c14_ordered.png through figure_s8f_site15_combined_single.png,
#   output/figures/figure_s9a_ann_c14_ordered.png through figure_s9f_ann_combined_single.png
# Runtime: ~2 minutes

library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rcarbon)
library(here)

# Source utilities
source(here("R", "utils.R"))

# Create output directory
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)

# Style constants
POSTERIOR_FILL <- "#2171b5"
POSTERIOR_ALPHA <- 0.7
POSTERIOR_OUTLINE <- "black"
UNMODELED_FILL <- "grey80"
UNMODELED_ALPHA <- 0.3
UNMODELED_OUTLINE <- "grey60"
POLYGON_LINEWIDTH <- 0.3
POLYGON_HEIGHT <- 0.85

bp_to_ad <- function(bp) 1950 - bp

# Publication theme
theme_oxcal <- theme_minimal() +
  theme(
    text = element_text(size = 10),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0),
    legend.position = "none"
  )

# Polygon helpers

create_density_polygon <- function(mean_val, sd_val, y_baseline, height = POLYGON_HEIGHT,
                                   n_points = 100, range_sd = 3.5) {
  x_seq <- seq(mean_val - range_sd * sd_val, mean_val + range_sd * sd_val, length.out = n_points)
  density_vals <- dnorm(x_seq, mean_val, sd_val)
  density_norm <- density_vals / max(density_vals) * height
  n <- length(x_seq)
  data.frame(
    x = c(x_seq[1], x_seq, x_seq[n]),
    y = c(y_baseline, y_baseline + density_norm, y_baseline)
  )
}

create_posterior_density <- function(samples, y_baseline, height = POLYGON_HEIGHT, n_points = 100) {
  dens <- density(samples, n = n_points, adjust = 1)
  density_norm <- dens$y / max(dens$y) * height
  n <- length(dens$x)
  data.frame(
    x = c(dens$x[1], dens$x, dens$x[n]),
    y = c(y_baseline, y_baseline + density_norm, y_baseline)
  )
}

calibrate_c14_polygon <- function(c14_age, c14_sd, y_baseline, height = POLYGON_HEIGHT) {
  cal <- rcarbon::calibrate(
    x = c14_age, errors = c14_sd,
    calCurves = "shcal20", verbose = FALSE
  )
  cal_grid <- cal$grids[[1]]
  cal_ad <- bp_to_ad(cal_grid$calBP)
  prob_norm <- cal_grid$PrDens / max(cal_grid$PrDens) * height
  n <- length(cal_ad)
  data.frame(
    x = c(cal_ad[1], cal_ad, cal_ad[n]),
    y = c(y_baseline, y_baseline + prob_norm, y_baseline)
  )
}

# Core plotting: build an OxCal-style figure from polygon data

build_oxcal_plot <- function(unmodeled_polys, posterior_polys, plot_data, title,
                             x_range = c(1000, 2000)) {
  # Determine axis limits from data
  all_x <- c(unmodeled_polys$x, posterior_polys$x)
  x_min <- floor(min(all_x, na.rm = TRUE) / 50) * 50 - 50
  x_max <- ceiling(max(all_x, na.rm = TRUE) / 50) * 50 + 50
  x_min <- max(x_range[1], x_min)
  x_max <- min(x_range[2], x_max)

  p <- ggplot() +
    geom_polygon(data = unmodeled_polys,
                 aes(x = x, y = y, group = sample_idx),
                 fill = UNMODELED_FILL, color = UNMODELED_OUTLINE,
                 alpha = UNMODELED_ALPHA, linewidth = POLYGON_LINEWIDTH) +
    geom_polygon(data = posterior_polys,
                 aes(x = x, y = y, group = sample_idx),
                 fill = POSTERIOR_FILL, color = POSTERIOR_OUTLINE,
                 alpha = POSTERIOR_ALPHA, linewidth = POLYGON_LINEWIDTH) +
    scale_y_continuous(
      breaks = plot_data$y_pos + 0.4,
      labels = plot_data$sample_label,
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    scale_x_continuous(
      breaks = seq(x_range[1], x_range[2], 100),
      limits = c(x_min, x_max)
    ) +
    labs(title = title, x = "Calendar Date (AD)") +
    theme_oxcal

  # Level separator lines
  level_changes <- which(diff(plot_data$Level_num) != 0)
  for (lc in level_changes) {
    p <- p + geom_hline(yintercept = lc + 0.5, color = "grey30",
                        linewidth = 0.5, linetype = "solid")
  }

  # Right-side level labels
  level_groups <- plot_data %>%
    group_by(Level_num) %>%
    summarise(y_mid = mean(y_pos) + 0.4, .groups = "drop")

  p <- p + annotate("text", x = x_max - 10, y = level_groups$y_mid,
                     label = paste0("L", level_groups$Level_num),
                     size = 3.5, fontface = "bold", hjust = 1)

  return(p)
}

# Figure function: C14-only models (A/B)

make_oxcal_c14 <- function(fit, c14_data, title) {
  posterior <- rstan::extract(fit)

  # Get true_age parameter (true_age for standalone, sample_age as fallback)
  if ("true_age" %in% names(posterior)) {
    true_age_post <- posterior$true_age
  } else if ("sample_age" %in% names(posterior)) {
    true_age_post <- posterior$sample_age
  } else {
    stop("No true_age or sample_age parameter found")
  }

  # Sort: oldest (highest Level_num) at bottom, youngest at top
  plot_data <- c14_data %>%
    mutate(orig_idx = row_number()) %>%
    arrange(desc(Level_num), C14_Age) %>%
    mutate(
      y_pos = row_number(),
      sample_label = paste0(Lab_ID, " (L", Level_num, ")")
    )

  n_samples <- nrow(plot_data)

  # Unmodeled: SHCal20 calibration curves
  unmodeled_polys <- lapply(1:n_samples, function(i) {
    poly <- calibrate_c14_polygon(
      c14_age = plot_data$C14_Age[i],
      c14_sd = plot_data$C14_SD[i],
      y_baseline = plot_data$y_pos[i]
    )
    poly$sample_idx <- i
    poly
  }) %>% bind_rows()

  # Posterior
  posterior_polys <- lapply(1:n_samples, function(i) {
    post_AD <- bp_to_ad(true_age_post[, plot_data$orig_idx[i]])
    poly <- create_posterior_density(samples = post_AD, y_baseline = plot_data$y_pos[i])
    poly$sample_idx <- i
    poly
  }) %>% bind_rows()

  build_oxcal_plot(unmodeled_polys, posterior_polys, plot_data, title)
}

# Figure function: OHD-only models (C/D)

make_oxcal_ohd <- function(fit, ohd_data, title) {
  posterior <- rstan::extract(fit)
  true_age_post <- posterior$true_age

  # Sort: oldest (highest Level_num) at bottom, youngest at top
  plot_data <- ohd_data %>%
    mutate(orig_idx = row_number()) %>%
    arrange(desc(Level_num), OHD_Date_AD) %>%
    mutate(
      y_pos = row_number(),
      sample_label = paste0(Sample_ID, " (L", Level_num, ")")
    )

  n_samples <- nrow(plot_data)

  # Unmodeled: normal distributions from claimed OHD dates
  unmodeled_polys <- lapply(1:n_samples, function(i) {
    poly <- create_density_polygon(
      mean_val = plot_data$OHD_Date_AD[i],
      sd_val = plot_data$OHD_SD[i],
      y_baseline = plot_data$y_pos[i]
    )
    poly$sample_idx <- i
    poly
  }) %>% bind_rows()

  # Posterior
  posterior_polys <- lapply(1:n_samples, function(i) {
    post_AD <- bp_to_ad(true_age_post[, plot_data$orig_idx[i]])
    poly <- create_posterior_density(samples = post_AD, y_baseline = plot_data$y_pos[i])
    poly$sample_idx <- i
    poly
  }) %>% bind_rows()

  build_oxcal_plot(unmodeled_polys, posterior_polys, plot_data, title)
}

# Figure function: Combined models (E/F)

make_oxcal_combined <- function(fit, c14_data, ohd_data, title) {
  posterior <- rstan::extract(fit)
  true_age_c14 <- posterior$true_age_c14  # draws x N_c14
  true_age_ohd <- posterior$true_age_ohd  # draws x N_ohd

  # Build a merged sample list with type information
  c14_info <- c14_data %>%
    mutate(
      orig_idx = row_number(),
      type = "C14",
      sample_label_raw = Lab_ID,
      date_AD = bp_to_ad(C14_Age),
      date_sd = C14_SD,
      c14_age = C14_Age
    ) %>%
    select(orig_idx, type, Level_num, sample_label_raw, date_AD, date_sd, c14_age)

  ohd_info <- ohd_data %>%
    mutate(
      orig_idx = row_number(),
      type = "OHD",
      sample_label_raw = Sample_ID,
      date_AD = OHD_Date_AD,
      date_sd = OHD_SD,
      c14_age = NA_real_
    ) %>%
    select(orig_idx, type, Level_num, sample_label_raw, date_AD, date_sd, c14_age)

  # Merge and sort: oldest (highest Level_num) at bottom
  plot_data <- bind_rows(c14_info, ohd_info) %>%
    arrange(desc(Level_num), date_AD) %>%
    mutate(
      y_pos = row_number(),
      type_tag = ifelse(type == "C14", "\u00b9\u2074C", "OHD"),
      sample_label = paste0(sample_label_raw, " [", type_tag, "] (L", Level_num, ")")
    )

  n_samples <- nrow(plot_data)

  # Unmodeled polygons: calibrated C14 or normal OHD depending on type
  unmodeled_polys <- lapply(1:n_samples, function(i) {
    if (plot_data$type[i] == "C14") {
      poly <- calibrate_c14_polygon(
        c14_age = plot_data$c14_age[i],
        c14_sd = plot_data$date_sd[i],
        y_baseline = plot_data$y_pos[i]
      )
    } else {
      poly <- create_density_polygon(
        mean_val = plot_data$date_AD[i],
        sd_val = plot_data$date_sd[i],
        y_baseline = plot_data$y_pos[i]
      )
    }
    poly$sample_idx <- i
    poly
  }) %>% bind_rows()

  # Posterior polygons: use true_age_c14 or true_age_ohd by type
  posterior_polys <- lapply(1:n_samples, function(i) {
    if (plot_data$type[i] == "C14") {
      post_AD <- bp_to_ad(true_age_c14[, plot_data$orig_idx[i]])
    } else {
      post_AD <- bp_to_ad(true_age_ohd[, plot_data$orig_idx[i]])
    }
    poly <- create_posterior_density(samples = post_AD, y_baseline = plot_data$y_pos[i])
    poly$sample_idx <- i
    poly
  }) %>% bind_rows()

  build_oxcal_plot(unmodeled_polys, posterior_polys, plot_data, title)
}

# Load data
cat("Loading data...\n")

# Site 15-233
site15_ohd <- read.csv(here("data", "site_15_233_ohd.csv"))
site15_ohd$OHD_BP <- 1950 - site15_ohd$OHD_Date_AD
site15_ohd$Level_num <- as.numeric(gsub("L", "", site15_ohd$Level))

site15_c14 <- read.csv(here("data", "site_15_233_radiocarbon.csv"))
site15_c14$Level_num <- as.numeric(gsub("L", "", site15_c14$Level))

# Ahu Nau Nau (with exclusions matching the analysis scripts)
ann_ohd_raw <- read.csv(here("data", "ahu_nau_nau_ohd.csv"))
ann_ohd <- ann_ohd_raw %>%
  filter(Level != 1, Intrusive == FALSE)
ann_ohd$OHD_BP <- 1950 - ann_ohd$OHD_Date_AD
ann_ohd$Level_num <- ann_ohd$Level

ann_c14_raw <- read.csv(here("data", "ahu_nau_nau_radiocarbon.csv"))
ann_c14 <- ann_c14_raw %>%
  filter(Rat_Bone == FALSE)
ann_c14$Level_num <- ann_c14$Level

cat(sprintf("  Site 15-233: %d C14, %d OHD\n", nrow(site15_c14), nrow(site15_ohd)))
cat(sprintf("  Ahu Nau Nau: %d C14, %d OHD\n", nrow(ann_c14), nrow(ann_ohd)))

# Generate all 12 figures

# --- Site 15-233 (Figures S8a-f) ---
cat("\nGenerating Site 15-233 figures...\n")

# S8a: C14 Ordered
cat("  Figure S8a (C14 Ordered)...\n")
fit <- readRDS(here("output", "fits", "site15_1A_c14_ordered.rds"))
p <- make_oxcal_c14(fit, site15_c14, "Model 1A: Radiocarbon, Phase-Ordered (6 phases)")
ggsave(here("output", "figures", "figure_s8a_site15_c14_ordered.png"), p, width = 10, height = 6, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S8b: C14 Single
cat("  Figure S8b (C14 Single)...\n")
fit <- readRDS(here("output", "fits", "site15_1B_c14_single.rds"))
p <- make_oxcal_c14(fit, site15_c14, "Model 1B: Radiocarbon, Single Phase")
ggsave(here("output", "figures", "figure_s8b_site15_c14_single.png"), p, width = 10, height = 6, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S8c: OHD Ordered
cat("  Figure S8c (OHD Ordered)...\n")
fit <- readRDS(here("output", "fits", "site15_1C_ohd_ordered.rds"))
p <- make_oxcal_ohd(fit, site15_ohd, "Model 1C: OHD, Phase-Ordered (6 phases)")
ggsave(here("output", "figures", "figure_s8c_site15_ohd_ordered.png"), p, width = 10, height = 10, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S8d: OHD Single
cat("  Figure S8d (OHD Single)...\n")
fit <- readRDS(here("output", "fits", "site15_1D_ohd_single.rds"))
p <- make_oxcal_ohd(fit, site15_ohd, "Model 1D: OHD, Single Phase")
ggsave(here("output", "figures", "figure_s8d_site15_ohd_single.png"), p, width = 10, height = 10, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S8e: Combined Ordered
cat("  Figure S8e (Combined Ordered)...\n")
fit <- readRDS(here("output", "fits", "site15_1E_combined_ordered.rds"))
p <- make_oxcal_combined(fit, site15_c14, site15_ohd,
                         "Model 1E: Combined (C14 + OHD), Phase-Ordered (6 phases)")
ggsave(here("output", "figures", "figure_s8e_site15_combined_ordered.png"), p, width = 10, height = 12, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S8f: Combined Single
cat("  Figure S8f (Combined Single)...\n")
fit <- readRDS(here("output", "fits", "site15_1F_combined_single.rds"))
p <- make_oxcal_combined(fit, site15_c14, site15_ohd,
                         "Model 1F: Combined (C14 + OHD), Single Phase")
ggsave(here("output", "figures", "figure_s8f_site15_combined_single.png"), p, width = 10, height = 12, dpi = 300)
rm(fit); gc(verbose = FALSE)

# --- Ahu Nau Nau (Figures S9a-f) ---
cat("\nGenerating Ahu Nau Nau figures...\n")

# S9a: C14 Ordered
cat("  Figure S9a (C14 Ordered)...\n")
fit <- readRDS(here("output", "fits", "ann_2A_c14_ordered.rds"))
p <- make_oxcal_c14(fit, ann_c14, "Model 2A: Radiocarbon, Phase-Ordered (4 phases)")
ggsave(here("output", "figures", "figure_s9a_ann_c14_ordered.png"), p, width = 10, height = 8, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S9b: C14 Single
cat("  Figure S9b (C14 Single)...\n")
fit <- readRDS(here("output", "fits", "ann_2B_c14_single.rds"))
p <- make_oxcal_c14(fit, ann_c14, "Model 2B: Radiocarbon, Single Phase")
ggsave(here("output", "figures", "figure_s9b_ann_c14_single.png"), p, width = 10, height = 8, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S9c: OHD Ordered
cat("  Figure S9c (OHD Ordered)...\n")
fit <- readRDS(here("output", "fits", "ann_2C_ohd_ordered.rds"))
p <- make_oxcal_ohd(fit, ann_ohd, "Model 2C: OHD, Phase-Ordered (4 phases)")
ggsave(here("output", "figures", "figure_s9c_ann_ohd_ordered.png"), p, width = 10, height = 10, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S9d: OHD Single
cat("  Figure S9d (OHD Single)...\n")
fit <- readRDS(here("output", "fits", "ann_2D_ohd_single.rds"))
p <- make_oxcal_ohd(fit, ann_ohd, "Model 2D: OHD, Single Phase")
ggsave(here("output", "figures", "figure_s9d_ann_ohd_single.png"), p, width = 10, height = 10, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S9e: Combined Ordered
cat("  Figure S9e (Combined Ordered)...\n")
fit <- readRDS(here("output", "fits", "ann_2E_combined_ordered.rds"))
p <- make_oxcal_combined(fit, ann_c14, ann_ohd,
                         "Model 2E: Combined (C14 + OHD), Phase-Ordered (4 phases)")
ggsave(here("output", "figures", "figure_s9e_ann_combined_ordered.png"), p, width = 10, height = 12, dpi = 300)
rm(fit); gc(verbose = FALSE)

# S9f: Combined Single
cat("  Figure S9f (Combined Single)...\n")
fit <- readRDS(here("output", "fits", "ann_2F_combined_single.rds"))
p <- make_oxcal_combined(fit, ann_c14, ann_ohd,
                         "Model 2F: Combined (C14 + OHD), Single Phase")
ggsave(here("output", "figures", "figure_s9f_ann_combined_single.png"), p, width = 10, height = 12, dpi = 300)
rm(fit); gc(verbose = FALSE)

cat("\nAll 12 phase model figures generated.\n")
