# Fundamental Uncertainties Invalidate the Precision of Obsidian Hydration Dating: A Rapa Nui Case Study

Bayesian reanalysis of obsidian hydration dating (OHD) uncertainty using the Stevenson et al. (2015) Rapa Nui dataset. Propagating realistic parameter uncertainty through the hydration rate equation yields dating precision of approximately ±300 years (95% CI ~600 years), compared to the claimed ±30 years. This uncertainty exceeds the island's entire prehistoric sequence (~520 years), rendering OHD dates uninformative for the chronological questions they were used to address.

**Authors**: Beau DiNapoli and Carl P. Lipo

## Requirements

- **R** >= 4.4.0
- **C++ toolchain** for Stan model compilation (g++ or clang++)
- ~16 GB RAM recommended (Stan models use 4 chains in parallel)

### R package installation

```r
install.packages(c(
  "rstan", "loo", "bayesplot", "tidyverse", "patchwork",
  "here", "gridExtra", "Bchron", "rcarbon", "truncnorm",
  "ggridges", "scales", "terra", "stringr"
))
```

## Reproduction

All scripts use `here::here()` for portable paths. Run from the project root directory.

### Quick start (pre-generated outputs)

All figures and summary tables are included in `output/figures/` and `output/tables/`. The manuscript can be reviewed without running any code.

### Full reproduction

Run all analyses in dependency order:

```bash
bash run_all.sh
```

Or run individual script groups:

```bash
# Group 1: Tier model fits (sequential, each depends on previous)
Rscript R/01_tier1.R
Rscript R/02_tier2.R
Rscript R/03_tier3.R

# Group 2: Post-fit analyses (require tier fits)
Rscript R/04_tier3_hybrid.R
Rscript R/05_compare_tiers.R
Rscript R/09_counterfactual.R
Rscript R/10_variance_attribution.R
Rscript R/11_posterior_summary.R
Rscript R/12_convergence_diagnostics.R
Rscript R/13_posterior_predictive.R
Rscript R/14_prior_predictive.R

# Group 3: Sensitivity analyses (independent of tier fits)
Rscript R/06_sensitivity_sweep.R
Rscript R/07_cross_sensitivity.R
Rscript R/08_measurement_error.R

# Group 4: Phase model validation
Rscript R/15_site_15_233.R
Rscript R/16_ahu_nau_nau.R
# R/17_ordering_test.R is deprecated (see script header for details)
Rscript R/18_phase_boundary_figures.R

# Group 5: Main text figures
# R/19_figure1_map.R requires TIF raster files not included in repo
# (pre-generated figure1_map.png provided in output/figures/)
# Rscript R/19_figure1_map.R
Rscript R/20_figure2_sensitivity.R
Rscript R/21_figure3_claims.R
Rscript R/22_figure4_validation.R
Rscript R/23_figure5_bottleneck.R
Rscript R/24_figure6_classification.R
Rscript R/25_spd_reconstruction.R
Rscript R/26_figure7_spd.R

# Group 6: SI analyses and figures
Rscript R/27_discrimination.R
Rscript R/28_correlation.R
Rscript R/29_figure_s4_variance.R
Rscript R/30_figure_s7_sensitivity.R

# Group 7: Session info
Rscript R/31_session_info.R
```

### Estimated runtimes

| Script group | Runtime |
|-------------|---------|
| Tier 1-3 model fits | ~30 min |
| Post-fit analyses | ~15 min |
| Sensitivity sweeps (06-08) | ~3 hours |
| Phase model validation (15-18) | ~40 min |
| Figure generation (19-30) | ~10 min |
| **Total** | **~4-5 hours** |

Runtimes measured on a modern desktop with 4+ cores and 16 GB RAM. Stan models run 4 chains in parallel by default.

## Directory layout

```
OHD-uncertainty/
├── R/                    # Analysis scripts (31 numbered + utils.R)
├── stan/                 # Stan model files
│   └── sensitivity/      # Sensitivity analysis variants
├── data/                 # Input datasets
├── output/
│   ├── fits/             # Stan fit objects (.rds, gitignored)
│   ├── figures/          # All figures (7 main text + 20 SI)
│   └── tables/           # All CSV summary tables
├── run_all.sh            # Orchestration script
├── references.bib        # BibTeX references
└── LICENSE               # MIT License
```

## Contact

Beau DiNapoli, Department of Anthropology, Binghamton University
dinapoli@binghamton.edu
