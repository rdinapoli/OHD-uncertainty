# 25_spd_reconstruction.R - Summed probability distribution reconstruction
# Purpose: Computes SPD envelopes under claimed and corrected uncertainty
# Inputs: data/stevenson_2015.csv, output/tables/tier3_summary.csv
# Outputs: output/tables/spd_metrics.csv, output/tables/spd_envelopes.csv
# Runtime: ~2 minutes

library(tidyverse)
library(here)

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

data_path <- here("data/stevenson_2015.csv")
output_dir <- here("output/tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
cat("Loading Stevenson 2015 data...\n")
ohd_data <- read_csv(data_path, show_col_types = FALSE)

cat(sprintf("Loaded %d samples\n", nrow(ohd_data)))

# Our uncertainty estimate from Tier 3
OUR_SE <- 299  # +/-299 years (10x Stevenson's claimed +/-30)
STEVENSON_SE <- 30  # Stevenson's claimed SE

# -----------------------------------------------------------------------------
# Function: Generate SPD from dates and uncertainties
# -----------------------------------------------------------------------------

generate_spd <- function(dates_ad, se, year_range = c(1000, 2000), resolution = 1) {
  years <- seq(year_range[1], year_range[2], by = resolution)
  spd <- rep(0, length(years))

  if (length(se) == 1) {
    se <- rep(se, length(dates_ad))
  }

  for (i in seq_along(dates_ad)) {
    if (!is.na(dates_ad[i]) && !is.na(se[i])) {
      prob <- dnorm(years, mean = dates_ad[i], sd = se[i])
      spd <- spd + prob
    }
  }

  spd <- spd / sum(spd * resolution)

  return(data.frame(year = years, density = spd))
}

# -----------------------------------------------------------------------------
# 4.1a: Stevenson's Analysis (+/-30 years)
# -----------------------------------------------------------------------------

cat("Generating SPDs with Stevenson's +/-30 year SE...\n")

spd_30_sa1 <- generate_spd(
  ohd_data$date_ad[ohd_data$study_area == "SA1"],
  se = STEVENSON_SE
) %>% mutate(scenario = "Stevenson (\u00b130 yr)", study_area = "SA1")

spd_30_sa2 <- generate_spd(
  ohd_data$date_ad[ohd_data$study_area == "SA2"],
  se = STEVENSON_SE
) %>% mutate(scenario = "Stevenson (\u00b130 yr)", study_area = "SA2")

spd_30_sa3 <- generate_spd(
  ohd_data$date_ad[ohd_data$study_area == "SA3"],
  se = STEVENSON_SE
) %>% mutate(scenario = "Stevenson (\u00b130 yr)", study_area = "SA3")

cat("  SA1 peak year:", spd_30_sa1$year[which.max(spd_30_sa1$density)], "\n")
cat("  SA2 peak year:", spd_30_sa2$year[which.max(spd_30_sa2$density)], "\n")
cat("  SA3 peak year:", spd_30_sa3$year[which.max(spd_30_sa3$density)], "\n")

# -----------------------------------------------------------------------------
# 4.1b: Stevenson's dates with our uncertainty (+/-299 years)
# -----------------------------------------------------------------------------

cat("Generating SPDs with corrected +/-299 year SE...\n")

spd_299_sa1 <- generate_spd(
  ohd_data$date_ad[ohd_data$study_area == "SA1"],
  se = OUR_SE
) %>% mutate(scenario = "Corrected (\u00b1299 yr)", study_area = "SA1")

spd_299_sa2 <- generate_spd(
  ohd_data$date_ad[ohd_data$study_area == "SA2"],
  se = OUR_SE
) %>% mutate(scenario = "Corrected (\u00b1299 yr)", study_area = "SA2")

spd_299_sa3 <- generate_spd(
  ohd_data$date_ad[ohd_data$study_area == "SA3"],
  se = OUR_SE
) %>% mutate(scenario = "Corrected (\u00b1299 yr)", study_area = "SA3")

cat("  SA1 peak year:", spd_299_sa1$year[which.max(spd_299_sa1$density)], "\n")
cat("  SA2 peak year:", spd_299_sa2$year[which.max(spd_299_sa2$density)], "\n")
cat("  SA3 peak year:", spd_299_sa3$year[which.max(spd_299_sa3$density)], "\n")

# -----------------------------------------------------------------------------
# Combine all SPDs
# -----------------------------------------------------------------------------

spd_all <- bind_rows(
  spd_30_sa1, spd_30_sa2, spd_30_sa3,
  spd_299_sa1, spd_299_sa2, spd_299_sa3
)

# Set factor levels for plotting
spd_all$scenario <- factor(spd_all$scenario,
                           levels = c("Stevenson (\u00b130 yr)",
                                      "Corrected (\u00b1299 yr)"))

spd_all$study_area <- factor(spd_all$study_area,
                             levels = c("SA1", "SA2", "SA3"),
                             labels = c("SA1 (Te Niu)",
                                        "SA2 (Maunga O'Koro)",
                                        "SA3 (Anakena West)"))

# Save envelopes
write_csv(spd_all, file.path(output_dir, "spd_envelopes.csv"))
cat("Saved: spd_envelopes.csv\n")

# -----------------------------------------------------------------------------
# Quantitative Analysis: Curve characteristics
# -----------------------------------------------------------------------------

cat("Computing SPD metrics...\n")

calc_spd_metrics <- function(spd_df) {
  spd_df %>%
    group_by(scenario, study_area) %>%
    summarize(
      peak_year = year[which.max(density)],
      peak_density = max(density),
      halfmax = max(density) / 2,
      fwhm_low = year[which(density >= halfmax)[1]],
      fwhm_high = year[tail(which(density >= halfmax), 1)],
      fwhm = fwhm_high - fwhm_low,
      cumprob = list(cumsum(density) / sum(density)),
      ci95_low = year[which(cumsum(density)/sum(density) >= 0.025)[1]],
      ci95_high = year[which(cumsum(density)/sum(density) >= 0.975)[1]],
      ci95_width = ci95_high - ci95_low,
      .groups = "drop"
    ) %>%
    select(-cumprob)
}

spd_metrics <- calc_spd_metrics(spd_all)

cat("\nSPD Metrics Summary:\n")
print(spd_metrics, n = Inf)

write_csv(spd_metrics, file.path(output_dir, "spd_metrics.csv"))
cat("Saved: spd_metrics.csv\n")

# -----------------------------------------------------------------------------
# Key Finding: Can we distinguish decline timing?
# -----------------------------------------------------------------------------

stevenson_claims <- data.frame(
  study_area = c("SA1 (Te Niu)", "SA2 (Maunga O'Koro)", "SA3 (Anakena West)"),
  claimed_decline = c(1660, 1705, 1850),
  claimed_confidence = c("99.36%", "95.00%", "32.8%")
)

cat("\nStevenson's claims vs. corrected 95% CI widths:\n")
for (sa in unique(spd_metrics$study_area)) {
  cat(sprintf("\n%s:\n", sa))

  for (scen in unique(spd_metrics$scenario)) {
    m <- spd_metrics %>% filter(study_area == sa, scenario == scen)
    if (nrow(m) > 0) {
      cat(sprintf("  %s: Peak at %d, 95%% CI width = %d years\n",
                  scen, m$peak_year, m$ci95_width))
    }
  }

  claim <- stevenson_claims$claimed_decline[stevenson_claims$study_area == sa]
  cat(sprintf("  Stevenson claimed decline: AD %d\n", claim))

  diff_from_contact <- abs(claim - 1722)
  corrected_ci <- spd_metrics %>%
    filter(study_area == sa, scenario == "Corrected (\u00b1299 yr)") %>%
    pull(ci95_width)

  if (length(corrected_ci) > 0) {
    if (diff_from_contact < corrected_ci / 2) {
      cat(sprintf("  Cannot distinguish AD %d from AD 1722 (diff=%d < half-CI=%d)\n",
                  claim, diff_from_contact, round(corrected_ci/2)))
    } else {
      cat(sprintf("  Potentially distinguishable (diff=%d >= half-CI=%d)\n",
                  diff_from_contact, round(corrected_ci/2)))
    }
  }
}

cat("\nSPD reconstruction complete.\n")
