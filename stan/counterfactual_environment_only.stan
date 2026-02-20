// Counterfactual: Environment Only
//
// This model asks: "What would uncertainty be if we knew composition perfectly?"
//
// FIXED at population means:
//   - H₂Ot = 0.18 wt% (island-wide mean)
//   - Ea = 84,000 J/mol (island-wide mean)
//
// UNCERTAIN (hierarchical):
//   - EHT (hierarchical by study area)
//   - RH (varying by study area)
//
// Purpose: Quantify the IRREDUCIBLE environmental contribution to OHD uncertainty
//          independent of compositional confounding
//
// Expected: Lower uncertainty than Tier 1/2 because composition is a major driver
//
// UPDATED 2025-12-11: Fixed truncation using inverse CDF method
//   - Replaced manual clamping with proper truncated normal via unit interval transform
//   - Eliminates boundary pile-up and gradient discontinuities
//
// Date: 2025-12-03

functions {
  // Transform unit interval to truncated normal using inverse CDF
  real unit_to_truncnorm(real u, real mu, real sigma, real lower, real upper) {
    real Phi_lower = Phi((lower - mu) / sigma);
    real Phi_upper = Phi((upper - mu) / sigma);
    real p = Phi_lower + u * (Phi_upper - Phi_lower);

    // Clamp p to avoid numerical issues at boundaries
    if (p < 1e-10) p = 1e-10;
    if (p > 1 - 1e-10) p = 1 - 1e-10;

    return mu + sigma * inv_Phi(p);
  }
}

data {
  int<lower=1> N;              // Number of samples
  int<lower=1> N_areas;        // Number of study areas (3)

  // Observed data
  vector<lower=0>[N] H2Ome;    // Molecular water content (wt%)

  // Study area assignments
  array[N] int<lower=1,upper=N_areas> study_area;   // Study area ID per sample

  // FIXED composition parameters (population means)
  real H2Ot_fixed;             // Fixed at 0.18 wt%
  real Ea_fixed;               // Fixed at 84000 J/mol

  // Constants
  real R_gas;                  // Gas constant (8.314 J/(mol·K))

  // Final rate equation coefficients (R²=1.000 fit to Stevenson 2015)
  real a_eq;                   // 23.8042281295 (Intercept)
  real b_eq;                   // 0.3782846434 (Coefficient of -Ea/(R×T))
  real c_eq;                   // 0.0002071707 (Coefficient of H₂Ot)
  real d_eq;                   // -6313.2075128583 (Coefficient of 1/T)

  // RH coefficient (from Mazer et al. 1991)
  real e_RH_mean;              // 2.15 (Easter Island specific)
  real<lower=0> e_RH_sd;       // 0.5 (uncertainty)
}

transformed data {
  // Truncation bounds
  real EHT_lower = 18;
  real EHT_upper = 28;
  real RH_lower = 0.7;
  real RH_upper = 1.0;
}

parameters {
  // Age parameters
  vector<lower=100,upper=900>[N] age;  // Age in years BP

  // Environmental parameters ONLY (hierarchical by study area)
  vector<lower=18,upper=28>[N_areas] EHT_mean;  // Mean EHT by study area
  vector<lower=0.1>[N_areas] EHT_sd;            // SD EHT by study area

  vector<lower=0.7,upper=1.0>[N_areas] RH_mean;  // Mean RH by study area
  vector<lower=0.01>[N_areas] RH_sd;             // SD RH by study area

  // NON-CENTERED: Environment on unit interval
  vector<lower=0.001,upper=0.999>[N] EHT_unit;
  vector<lower=0.001,upper=0.999>[N] RH_unit;

  // RH coefficient with uncertainty
  real<lower=0> e_RH;  // Sampled RH sensitivity coefficient
}

transformed parameters {
  // Environmental parameters (uncertain)
  vector[N] EHT;
  vector[N] RH;
  vector[N] T_K;

  for (i in 1:N) {
    // Transform unit interval to truncated normal
    EHT[i] = unit_to_truncnorm(EHT_unit[i], EHT_mean[study_area[i]], EHT_sd[study_area[i]], EHT_lower, EHT_upper);
    RH[i] = unit_to_truncnorm(RH_unit[i], RH_mean[study_area[i]], RH_sd[study_area[i]], RH_lower, RH_upper);

    T_K[i] = EHT[i] + 273.15;  // Convert to Kelvin
  }
}

model {
  // Priors on study-area level environmental parameters
  EHT_mean ~ normal(22.5, 2.5);   // Wide prior for study-area means
  EHT_sd ~ exponential(0.33);      // Allows ~3°C within-area SD

  // RH priors per study area (rainfall gradient)
  for (k in 1:N_areas) {
    RH_mean[k] ~ normal(0.85, 0.1);  // Independent priors per area
  }
  RH_sd ~ exponential(5);  // Allows ~0.2 within-area variation

  // Implicit uniform prior on unit interval parameters
  // The inverse CDF transformation gives proper truncated normal

  // Age priors (weakly informative for Rapa Nui range)
  age ~ normal(500, 200);

  // RH sensitivity coefficient prior
  e_RH ~ normal(e_RH_mean, e_RH_sd);

  // Likelihood: H₂Ome measurements given ages and parameters
  // NOTE: H₂Ot and Ea are FIXED at population means
  for (i in 1:N) {
    // Rate equation with FIXED composition
    real log_rate = a_eq + b_eq * (-Ea_fixed / (R_gas * T_K[i])) +
                    c_eq * H2Ot_fixed + d_eq * (1.0 / T_K[i]) +
                    e_RH * (RH[i] - 0.98);

    // Numerical stability
    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;

    real rate = exp(log_rate);

    // Hydration model: H2Ome = sqrt(rate * t)
    real H2Ome_pred = sqrt(rate * age[i]);

    // Measurement model
    H2Ome[i] ~ normal(H2Ome_pred, 0.03);
  }
}

generated quantities {
  // Posterior predictive checks
  vector[N] H2Ome_pred;
  vector[N] log_lik;

  // Variance components (environment only in this model)
  real var_environment;

  for (i in 1:N) {
    real log_rate = a_eq + b_eq * (-Ea_fixed / (R_gas * T_K[i])) +
                    c_eq * H2Ot_fixed + d_eq * (1.0 / T_K[i]) +
                    e_RH * (RH[i] - 0.98);

    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;

    real rate = exp(log_rate);
    H2Ome_pred[i] = sqrt(rate * age[i]);
    log_lik[i] = normal_lpdf(H2Ome[i] | H2Ome_pred[i], 0.03);
  }

  // Posterior variance in environmental parameters across all samples (sum of
  // EHT and RH variances), representing the environmental contribution to age
  // uncertainty when composition is held fixed at population means
  var_environment = variance(EHT) + variance(RH);
}
