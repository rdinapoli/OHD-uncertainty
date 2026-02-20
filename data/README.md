# Data Dictionary

## stevenson_2015.csv

Primary dataset: 242 obsidian hydration samples from three study areas on Rapa Nui (Easter Island), extracted from the supplementary tables of Stevenson et al. (2015).

**Source**: Stevenson, C. M., Puleston, C. O., Vitousek, P. M., Chadwick, O. A., Haoa, S., & Ladefoged, T. N. (2015). Variation in Rapa Nui (Easter Island) land use indicates production and population peaks prior to European contact. *Proceedings of the National Academy of Sciences*, 112(4), 1025-1030.

| Column | Description | Units | Notes |
|--------|-------------|-------|-------|
| `study_area` | Study area identifier | SA1, SA2, SA3 | SA1: Te Niu (n=127), SA2: Maunga O'Koro (n=50), SA3: Anakena West (n=65) |
| `lab_id` | Laboratory sample identifier | — | DHR-series numbers |
| `provenience` | Archaeological provenience | — | Site-unit designation |
| `1630_cm1` | IR-PAS absorbance at 1630 cm⁻¹ | absorbance | Raw spectroscopic measurement |
| `H2Ome_pct` | Molecular water content | wt% | Measured by IR-PAS; primary observable for dating |
| `3570_cm1` | IR-PAS absorbance at 3570 cm⁻¹ | absorbance | Raw spectroscopic measurement |
| `H2Ot_pct` | Structural (total) water content | wt% | Measured by IR-PAS; rate equation parameter |
| `eht_c` | Effective hydration temperature | °C | Study-area constants (SA1=22.4, SA2=21.8, SA3=23.3) |
| `RH_fraction` | Relative humidity | fraction | Universal constant (0.98 for all samples) |
| `rate` | Hydration rate | μm²/day | Calculated from EHT, RH, H2Ot, Ea |
| `H2Ome2_per_day` | Squared molecular water per day | (wt%)²/day | Intermediate calculation |
| `Ea_J_per_mol` | Activation energy | J/mol | From Arrhenius fit to composition data |
| `date_bp` | Calculated OHD date | years BP | Derived from assumed EHT/RH (circular) |
| `date_ad` | Calculated OHD date | years AD | Derived from assumed EHT/RH (circular) |
| `sd` | Claimed uncertainty | years | ±30 for all samples |

**Important**: The `date_bp`, `date_ad`, `sd`, `eht_c`, and `RH_fraction` columns are not independent measurements. The dates are calculated from the other columns using assumed (not measured) EHT and RH values. Our Bayesian analysis treats EHT and RH as uncertain parameters rather than known constants.

## Validation datasets (Stevenson 2013)

Four files from Stevenson et al. (2013) used for independent phase model validation:

**Source**: Stevenson, C. M., Ladefoged, T. N., & Novak, S. W. (2013). Prehistoric settlement chronology on Rapa Nui, Chile: Obsidian hydration dating using infrared photoacoustic spectroscopy. *Journal of Archaeological Science*, 40(7), 3021-3030.

| File | Description | n |
|------|-------------|---|
| `site_15_233_ohd.csv` | OHD dates for Site 15-233 | 23 |
| `site_15_233_radiocarbon.csv` | Radiocarbon dates for Site 15-233 | 3 |
| `ahu_nau_nau_ohd.csv` | OHD dates for Ahu Nau Nau | 29 |
| `ahu_nau_nau_radiocarbon.csv` | Radiocarbon dates for Ahu Nau Nau | 9 |
