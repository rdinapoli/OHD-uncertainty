// Tier 3 with Measurement Error Sensitivity: Cross-Sensitivity Analysis
//
// Purpose: What happens when measurement precision improves but environmental
// parameters remain uncertain? This model tests whether parameter uncertainty
// creates an independent floor that persists regardless of sigma_meas.
//
// This is the COMPLETE model (Tier 3) with sigma_meas passed as data for sweeping.
// Compare results to tier0_meas_error_sensitivity.stan (all params fixed).
//
// Expected behavior:
//   - If sigma_meas is the only bottleneck: uncertainty should scale with sigma_meas
//   - If parameter uncertainty also matters: uncertainty plateaus as sigma_meas -> 0
//   - The GAP between Tier 0 and Tier 3 at small sigma_meas = irreducible parameter floor
//
// From variance decomposition:
//   - EHT: -59.8 years/C x sigma_EHT ~2.5C -> ~+/-150 years
//   - RH: -853 years/RH x sigma_RH ~0.05 -> ~+/-43 years
//   - Expected parameter floor: ~+/-150-200 years
//
// Date: 2025-12-04
// Purpose: Cross-sensitivity analysis to determine which bottleneck dominates
// UPDATED 2026-02-19: Standardized truncation to use unit_to_truncnorm(),
//   matching tier1/tier2/tier3 and counterfactual models

functions {
  // Transform unit interval to truncated normal using inverse CDF
  // u: value in (0,1)
  // mu: mean of underlying normal
  // sigma: sd of underlying normal
  // lower, upper: truncation bounds
  real unit_to_truncnorm(real u, real mu, real sigma, real lower, real upper) {
    real Phi_lower = Phi((lower - mu) / sigma);
    real Phi_upper = Phi((upper - mu) / sigma);
    real p = Phi_lower + u * (Phi_upper - Phi_lower);

    // Numerical stability: prevents log(0) at distribution boundaries
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

  // MEASUREMENT ERROR - NOW PASSED AS DATA FOR SENSITIVITY SWEEP
  real<lower=0> sigma_meas;    // IR-PAS measurement error (wt%)

  // Temporal climate parameters
  real<lower=0> temporal_climate_sd;  // +/-0.5C temperature variation over artifact lifetime

  // Constants
  real R_gas;                  // Gas constant (8.314 J/(mol*K))

  // Final rate equation coefficients (R^2=1.000 fit to Stevenson 2015)
  real a_eq;                   // 23.8042281295 (Intercept)
  real b_eq;                   // 0.3782846434 (Coefficient of -Ea/(R*T))
  real c_eq;                   // 0.0002071707 (Coefficient of H2Ot)
  real d_eq;                   // -6313.2075128583 (Coefficient of 1/T)

  // RH coefficient (from Mazer et al. 1991)
  real e_RH_mean;              // 2.15 (Easter Island specific)
  real<lower=0> e_RH_sd;       // 0.5 (uncertainty)
}

transformed data {
  // Truncation bounds (matching tier3_complete.stan)
  real H2Ot_lower = 0.05;
  real H2Ot_upper = 0.43;
  real Ea_lower = 75000;
  real Ea_upper = 92000;
  real EHT_lower = 18;
  real EHT_upper = 28;
  real RH_lower = 0.7;
  real RH_upper = 1.0;
}

parameters {
  // Age parameters
  vector<lower=100,upper=900>[N] age;  // Age in years BP (narrower Rapa Nui range)

  // Composition parameters (ISLAND-WIDE population, not source-specific)
  real<lower=0.05,upper=0.5> H2Ot_mu;    // Population mean H2Ot
  real<lower=0.01> H2Ot_sigma;            // Population SD H2Ot

  real<lower=79000,upper=86000> Ea_mu;   // Population mean Ea
  real<lower=100> Ea_sigma;               // Population SD Ea

  // NON-CENTERED: Composition on unit interval
  vector<lower=0.001,upper=0.999>[N] H2Ot_unit;
  vector<lower=0.001,upper=0.999>[N] Ea_unit;

  // Environmental parameters (hierarchical by study area)
  vector<lower=18,upper=28>[N_areas] EHT_mean;  // Mean EHT by study area
  vector<lower=0.1>[N_areas] EHT_sd;            // SD EHT by study area (within-area variation)

  vector<lower=0.7,upper=1.0>[N_areas] RH_mean;  // Mean RH by study area
  vector<lower=0.01>[N_areas] RH_sd;             // SD RH by study area

  // NON-CENTERED: Environment on unit interval
  vector<lower=0.001,upper=0.999>[N] EHT_unit;
  vector<lower=0.001,upper=0.999>[N] RH_unit;

  // Temporal climate effects (unbounded - no truncation needed)
  vector[N] temporal_climate_effect;  // Temperature deviation over artifact lifetime

  // RH coefficient with uncertainty
  real<lower=0> e_RH;  // Sampled RH sensitivity coefficient
}

transformed parameters {
  // Transform unit interval to truncated normal
  vector[N] H2Ot;
  vector[N] Ea;
  vector[N] EHT;
  vector[N] RH;

  // Effective temperature includes temporal variation
  vector[N] T_eff;
  vector[N] T_K;

  for (i in 1:N) {
    // Composition (island-wide)
    H2Ot[i] = unit_to_truncnorm(H2Ot_unit[i], H2Ot_mu, H2Ot_sigma, H2Ot_lower, H2Ot_upper);
    Ea[i] = unit_to_truncnorm(Ea_unit[i], Ea_mu, Ea_sigma, Ea_lower, Ea_upper);

    // Environment (hierarchical by study area)
    EHT[i] = unit_to_truncnorm(EHT_unit[i], EHT_mean[study_area[i]], EHT_sd[study_area[i]], EHT_lower, EHT_upper);
    RH[i] = unit_to_truncnorm(RH_unit[i], RH_mean[study_area[i]], RH_sd[study_area[i]], RH_lower, RH_upper);

    // Add temporal climate effect to burial temperature
    T_eff[i] = EHT[i] + temporal_climate_effect[i];
    T_K[i] = T_eff[i] + 273.15;  // Convert to Kelvin
  }
}

model {
  // Priors on island-wide composition parameters
  H2Ot_mu ~ normal(0.18, 0.08);
  H2Ot_sigma ~ exponential(10);
  Ea_mu ~ normal(84000, 2000);
  Ea_sigma ~ exponential(0.0005);

  // Priors on study-area level environmental parameters
  EHT_mean ~ normal(22.5, 2.5);
  EHT_sd ~ exponential(0.33);

  // Allow RH to vary substantially between areas (rainfall gradient!)
  for (k in 1:N_areas) {
    RH_mean[k] ~ normal(0.85, 0.1);
  }
  RH_sd ~ exponential(5);

  // Implicit uniform prior on unit interval parameters
  // The inverse CDF transformation gives proper truncated normal

  // Temporal climate variation prior
  temporal_climate_effect ~ normal(0, temporal_climate_sd);

  // Age priors (weakly informative for Rapa Nui range)
  age ~ normal(500, 200);

  // RH sensitivity coefficient prior
  e_RH ~ normal(e_RH_mean, e_RH_sd);

  // Likelihood: H2Ome measurements given ages and parameters
  for (i in 1:N) {
    // EXTENDED rate equation WITH RH term
    real log_rate = a_eq + b_eq * (-Ea[i] / (R_gas * T_K[i])) +
                    c_eq * H2Ot[i] + d_eq * (1.0 / T_K[i]) + e_RH * (RH[i] - 0.98);

    // Numerical stability
    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;

    real rate = exp(log_rate);

    // Hydration model: H2Ome = sqrt(rate * t)
    real H2Ome_pred = sqrt(rate * age[i]);

    // Measurement model - sigma_meas now from data, not hardcoded
    H2Ome[i] ~ normal(H2Ome_pred, sigma_meas);
  }
}

generated quantities {
  // Posterior predictive checks
  vector[N] H2Ome_pred;
  vector[N] log_lik;

  // Store sigma_meas used for this run
  real sigma_meas_used = sigma_meas;

  for (i in 1:N) {
    // EXTENDED rate equation WITH RH term
    real log_rate = a_eq + b_eq * (-Ea[i] / (R_gas * T_K[i])) +
                    c_eq * H2Ot[i] + d_eq * (1.0 / T_K[i]) + e_RH * (RH[i] - 0.98);

    // Numerical stability
    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;

    real rate = exp(log_rate);
    H2Ome_pred[i] = sqrt(rate * age[i]);
    log_lik[i] = normal_lpdf(H2Ome[i] | H2Ome_pred[i], sigma_meas);
  }
}
