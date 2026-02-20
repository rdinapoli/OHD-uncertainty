# 17_ordering_test.R - Stratigraphic ordering test
# Purpose: Extracts ordering violation percentages from phase model posteriors (Table S9)
# Inputs: output/fits/site15_*.rds, output/fits/ann_*.rds
# Outputs: output/tables/ordering_test_results.csv (printed to console)
# Runtime: ~1 minute

# DEPRECATION NOTE (2025-12-11):
# This script was designed for an "unconstrained" model that allows phases
# but doesn't enforce ordering. The current model suite uses:
#   - fit_2D_ohd_single.rds (single-phase, no boundaries)
#   - fit_2C_ohd_ordered.rds (ordered phases, forced correct ordering)
#
# The "0.3% correct ordering" analysis is replaced by ELPD-based evidence:
#   - DELTA_ELPD = -98 (Site 15-233): OHD ordered 98 worse than single
#   - DELTA_ELPD = -16 (Ahu Nau Nau): OHD ordered 16 worse than single
#
# To revive this analysis, an unconstrained model would need to be created.

stop("DEPRECATED: This script requires fit_model_2D_unconstrained.rds which no longer exists.
See deprecation note above for alternative ELPD-based evidence.")

library(rstan)
library(tidyverse)
library(here)

# Load the unconstrained OHD model fit
fit <- readRDS(here("output", "fits", "ann_2D_ohd_single.rds"))

# Extract posteriors
posterior <- rstan::extract(fit)

# The model has phase_boundary (3 boundaries between 4 phases)
# Boundaries are: boundary[1] between Phase 1 and 2, boundary[2] between 2 and 3, etc.
# Correct stratigraphic order: boundary[1] < boundary[2] < boundary[3]

cat("Available parameters:\n")
for (p in names(posterior)) {
  cat(sprintf("  %s: dim = %s\n", p, paste(dim(posterior[[p]]), collapse=" x ")))
}

phase_boundaries <- posterior$phase_boundary  # 8000 x 3
n_samples <- nrow(phase_boundaries)

cat(sprintf("\nPhase boundaries dimensions: %d samples x %d boundaries\n",
            nrow(phase_boundaries), ncol(phase_boundaries)))

# Check the boundaries_ordered indicator from the model
cat("\nModel's boundaries_ordered indicator:\n")
ordered_pct_model <- mean(posterior$boundaries_ordered) * 100
cat(sprintf("Percentage with correct boundary ordering: %.1f%%\n", ordered_pct_model))

# Analyze what orderings DO occur
# Convert boundaries to ordering string
get_boundary_order <- function(bounds) {
  paste(order(bounds), collapse="-")
}

orderings <- apply(phase_boundaries, 1, get_boundary_order)
ordering_counts <- table(orderings)
ordering_df <- data.frame(
  ordering = names(ordering_counts),
  count = as.numeric(ordering_counts),
  percentage = as.numeric(ordering_counts) / n_samples * 100
) %>%
  arrange(desc(count))

cat("\nBoundary ordering '1-2-3' means: boundary1 < boundary2 < boundary3 (correct stratigraphic order)\n")
cat("Phases: Phase 1 (Level 2), Phase 2 (Level 4), Phase 3 (Level 5), Phase 4 (Level 6)\n")

# The correct ordering is "1-2-3"
correct_ordering <- "1-2-3"
correct_count <- sum(orderings == correct_ordering)
correct_pct <- correct_count / n_samples * 100

cat(sprintf("\nCorrect ordering (1-2-3): %.1f%% (%d of %d samples)\n",
            correct_pct, correct_count, n_samples))

cat("\nAll orderings:\n")
print(ordering_df)

# Interpret the most common ordering
most_common <- ordering_df[1, ]
cat(sprintf("\nMost common ordering: %s (%.1f%%)\n",
            most_common$ordering, most_common$percentage))

# Parse the ordering
ord <- as.numeric(strsplit(most_common$ordering, "-")[[1]])

cat("\nInterpretation of most common ordering:\n")

boundary_ages <- c("youngest", "middle", "oldest")
cat("Boundaries sorted by age:\n")
for (i in 1:3) {
  cat(sprintf("  %s: Boundary %d (between Phase %d and Phase %d)\n",
              boundary_ages[i], ord[i], ord[i], ord[i]+1))
}

# Key result
cat(sprintf("\nCorrect stratigraphic ordering: %.1f%% of posterior samples\n", correct_pct))
cat(sprintf("Most common ordering: '%s' (%.1f%% of samples)\n",
            most_common$ordering, most_common$percentage))

if (most_common$ordering != correct_ordering) {
  cat("The most common ordering is INCORRECT.\n")
}

# Save results
results <- list(
  n_samples = n_samples,
  correct_percentage = correct_pct,
  ordering_df = ordering_df,
  most_common = most_common
)

saveRDS(results, here("output", "fits", "ordering_percentages.rds"))
write.csv(ordering_df, here("output", "tables", "ordering_test_results.csv"), row.names = FALSE)

cat("\nResults saved.\n")
