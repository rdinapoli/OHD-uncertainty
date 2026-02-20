#!/usr/bin/env bash
# run_all.sh - Run all analyses in dependency order
# Total runtime: ~4-5 hours on a modern desktop (4+ cores, 16 GB RAM)
#
# Usage: bash run_all.sh
#
# Scripts are organized into groups by dependency. Within each group,
# scripts are run sequentially. Groups 2 and 3 could run in parallel
# if desired (they have no mutual dependencies).

set -e  # Exit on first error

echo "=== Group 1: Tier model fits (sequential) ==="
Rscript R/01_tier1.R
Rscript R/02_tier2.R
Rscript R/03_tier3.R

echo "=== Group 2: Post-fit analyses (require tier fits) ==="
Rscript R/04_tier3_hybrid.R
Rscript R/05_compare_tiers.R
Rscript R/09_counterfactual.R
Rscript R/10_variance_attribution.R
Rscript R/11_posterior_summary.R
Rscript R/12_convergence_diagnostics.R
Rscript R/13_posterior_predictive.R
Rscript R/14_prior_predictive.R

echo "=== Group 3: Sensitivity analyses (independent) ==="
Rscript R/06_sensitivity_sweep.R
Rscript R/07_cross_sensitivity.R
Rscript R/08_measurement_error.R

echo "=== Group 4: Phase model validation ==="
Rscript R/15_site_15_233.R
Rscript R/16_ahu_nau_nau.R
# R/17_ordering_test.R is deprecated (requires unconstrained model that no longer exists;
# see ELPD-based alternative documented in the script header)
Rscript R/18_phase_boundary_figures.R

echo "=== Group 5: Main text figures ==="
# R/19_figure1_map.R requires TIF raster files not included in the repo
# (pre-generated figure1_map.png is provided in output/figures/)
# Rscript R/19_figure1_map.R
Rscript R/20_figure2_sensitivity.R
Rscript R/21_figure3_claims.R
Rscript R/22_figure4_validation.R
Rscript R/23_figure5_bottleneck.R
Rscript R/24_figure6_classification.R
Rscript R/25_spd_reconstruction.R
Rscript R/26_figure7_spd.R

echo "=== Group 6: SI analyses and figures ==="
Rscript R/27_discrimination.R
Rscript R/28_correlation.R
Rscript R/29_figure_s4_variance.R
Rscript R/30_figure_s7_sensitivity.R

echo "=== Group 7: Session info ==="
Rscript R/31_session_info.R

echo "=== All analyses complete ==="
