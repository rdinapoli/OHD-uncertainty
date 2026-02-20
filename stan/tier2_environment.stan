// Tier 2: Composition + Environmental Uncertainty
//
// Adds EHT and RH uncertainty to Tier 1 composition model
// Assumptions:
//   - H2Ome measured perfectly (given)
//   - H2Ot, Ea uncertain (island-wide population distributions)
//   - EHT uncertain (hierarchical by study area)
//   - RH uncertain (varying by study area)
//
// NO SOURCE STRUCTURE - Rationale:
//   1. Stevenson 2015 has NO per-sample XRF (source assignments would be assumptions)
//   2. Orito and Rano Kau I are geochemically nearly identical
//   3. Source mixing contributed 0% variance in testing
//   4. Island-wide distributions are more honest about what we know
//
// UPDATED 2025-12-11: Fixed truncation using inverse CDF method (Option B)
//   - Replaced manual clamping with proper truncated normal via unit interval transform
//   - Eliminates boundary pile-up and gradient discontinuities
//
// Purpose: Quantify environmental contribution to OHD uncertainty
// Expected: ~+/-298 years (9.9x Stevenson's claimed +/-30 years)

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

  // Constants
  real R_gas;                  // Gas constant (8.314 J/(mol*K))

  // Final rate equation coefficients (R^2=1.000 fit to Stevenson 2015)
  // Rate equation coefficients from Stevenson (2015) supplementary material
  real a_eq;                   // 23.8042281295 (Intercept)
  real b_eq;                   // 0.3782846434 (Coefficient of -Ea/(R*T))
  real c_eq;                   // 0.0002071707 (Coefficient of H2Ot)
  real d_eq;                   // -6313.2075128583 (Coefficient of 1/T)

  // RH coefficient (from Mazer et al. 1991)
  // RH coefficient derived from Mazer et al. (1991) Easter Island basalt data
  real e_RH_mean;              // 2.15 (Easter Island specific)
  real<lower=0> e_RH_sd;       // 0.5 (uncertainty)
}

transformed data {
  // Truncation bounds
  // H2Ot: [0.05, 0.43] wt% — matches observed data range (0.052-0.432)
  real H2Ot_lower = 0.05;
  real H2Ot_upper = 0.43;
  // Ea: [75000, 92000] J/mol — extends beyond observed range (79642-85858) to avoid boundary effects
  real Ea_lower = 75000;
  real Ea_upper = 92000;
  // EHT: [18, 28] C — encompasses study-area means (21.8-23.3C) +/- ~5C for burial depth variation
  real EHT_lower = 18;
  real EHT_upper = 28;
  // RH: [0.7, 1.0] — subsurface range; lower bound from driest plausible conditions
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

  // RH coefficient with uncertainty
  real<lower=0> e_RH;  // Sampled RH sensitivity coefficient
}

transformed parameters {
  // Transform unit interval to truncated normal
  vector[N] H2Ot;
  vector[N] Ea;
  vector[N] EHT;
  vector[N] RH;
  vector[N] T_K;

  for (i in 1:N) {
    // Composition (island-wide)
    H2Ot[i] = unit_to_truncnorm(H2Ot_unit[i], H2Ot_mu, H2Ot_sigma, H2Ot_lower, H2Ot_upper);
    Ea[i] = unit_to_truncnorm(Ea_unit[i], Ea_mu, Ea_sigma, Ea_lower, Ea_upper);

    // Environment (hierarchical by study area)
    EHT[i] = unit_to_truncnorm(EHT_unit[i], EHT_mean[study_area[i]], EHT_sd[study_area[i]], EHT_lower, EHT_upper);
    RH[i] = unit_to_truncnorm(RH_unit[i], RH_mean[study_area[i]], RH_sd[study_area[i]], RH_lower, RH_upper);

    T_K[i] = EHT[i] + 273.15;  // Convert to Kelvin
  }
}

model {
  // Priors on island-wide composition parameters
  // Based on observed range in Stevenson 2015: H2Ot 0.05-0.43 wt%, Ea 79.6-85.9 kJ/mol
  H2Ot_mu ~ normal(0.18, 0.08);      // Centered on observed mean (0.052-0.432 range)
  H2Ot_sigma ~ exponential(10);       // Mean 0.1 wt%; Rogers & Stevenson (2022): within-source CV 20-40%
  Ea_mu ~ normal(84000, 2000);        // Centered on observed mean (79642-85858 range)
  Ea_sigma ~ exponential(0.0005);     // Mean 2000 J/mol; weakly informative (observed SD = 532)

  // Priors on study-area level environmental parameters
  EHT_mean ~ normal(22.5, 2.5);   // Wider prior (allow 18-27C range across areas)
  EHT_sd ~ exponential(0.33);      // Mean 3C; combined burial depth, spatial, and paleoclimate uncertainty

  // Allow RH to vary substantially between areas (rainfall gradient!)
  for (k in 1:N_areas) {
    RH_mean[k] ~ normal(0.85, 0.1);  // Independent priors per area (wide)
  }
  RH_sd ~ exponential(5);  // Mean 0.2; 2.1x rainfall gradient across study areas

  // Implicit uniform prior on unit interval parameters
  // The inverse CDF transformation gives proper truncated normal

  // Age priors (weakly informative for Rapa Nui range)
  age ~ normal(500, 200);

  // RH sensitivity coefficient prior
  e_RH ~ normal(e_RH_mean, e_RH_sd);  // 2.15 +/- 0.5 (Easter Island specific)

  // Likelihood: H2Ome measurements given ages and parameters
  for (i in 1:N) {
    // EXTENDED rate equation WITH RH term
    // log(rate) = a + b*(-Ea/(R*T)) + c*H2Ot + d*(1/T) + e_RH*(RH - 0.98)
    // When RH = 0.98, RH term = 0, preserving R^2=1.000 fit to Stevenson's data
    real log_rate = a_eq + b_eq * (-Ea[i] / (R_gas * T_K[i])) + c_eq * H2Ot[i] + d_eq * (1.0 / T_K[i]) + e_RH * (RH[i] - 0.98);

    // Numerical stability
    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;

    real rate = exp(log_rate);

    // Hydration model: H2Ome = sqrt(rate * t)
    real H2Ome_pred = sqrt(rate * age[i]);

    // Measurement model
    H2Ome[i] ~ normal(H2Ome_pred, 0.03);  // +/-0.03 wt%
  }
}

generated quantities {
  // Posterior predictive checks
  vector[N] H2Ome_pred;  // Deterministic mean prediction (for log_lik)
  vector[N] H2Ome_rep;   // Replicated observations with measurement noise (for PPC)
  vector[N] log_lik;

  for (i in 1:N) {
    // EXTENDED rate equation WITH RH term
    real log_rate = a_eq + b_eq * (-Ea[i] / (R_gas * T_K[i])) + c_eq * H2Ot[i] + d_eq * (1.0 / T_K[i]) + e_RH * (RH[i] - 0.98);

    // Numerical stability
    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;

    real rate = exp(log_rate);
    H2Ome_pred[i] = sqrt(rate * age[i]);
    H2Ome_rep[i] = normal_rng(H2Ome_pred[i], 0.03);
    log_lik[i] = normal_lpdf(H2Ome[i] | H2Ome_pred[i], 0.03);
  }
  // Variance decomposition performed in sensitivity_variance_decomposition.R (delta method)
}
