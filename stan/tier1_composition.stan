// Tier 1: Composition Uncertainty Only
//
// Infers ages from H2Ome measurements given uncertain H2Ot and Ea
// Assumptions:
//   - H2Ome measured perfectly (given)
//   - H2Ot uncertain (island-wide population distribution)
//   - Ea uncertain (island-wide population distribution)
//   - EHT, RH known exactly (UNREALISTIC, but isolates composition effect)
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
// Purpose: Quantify composition contribution to OHD uncertainty
// Expected: ~+/-205 years (6.8x Stevenson's claimed +/-30 years)

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

  // Observed data
  vector<lower=0>[N] H2Ome;    // Molecular water content (wt%)

  // Known parameters (for Tier 1 only - treats Stevenson's values as truth)
  vector[N] EHT;               // Effective hydration temperature (C)

  // Constants
  real R_gas;                  // Gas constant (8.314 J/(mol*K))

  // Rate equation coefficients from Stevenson (2015) supplementary material
  real a_eq;                   // 23.8042281295 (Intercept)
  real b_eq;                   // 0.3782846434 (Coefficient of -Ea/(R*T))
  real c_eq;                   // 0.0002071707 (Coefficient of H2Ot)
  real d_eq;                   // -6313.2075128583 (Coefficient of 1/T)
}

transformed data {
  vector[N] T_K = EHT + 273.15;  // Convert to Kelvin

  // Truncation bounds
  // H2Ot: [0.05, 0.43] wt% — matches observed data range (0.052-0.432)
  real H2Ot_lower = 0.05;
  real H2Ot_upper = 0.43;
  // Ea: [75000, 92000] J/mol — extends beyond observed range (79642-85858) to avoid boundary effects
  real Ea_lower = 75000;
  real Ea_upper = 92000;
}

parameters {
  // Age parameters
  vector<lower=100,upper=900>[N] age;  // Age in years BP (narrower Rapa Nui range)

  // Composition parameters (ISLAND-WIDE population, not source-specific)
  real<lower=0.05,upper=0.5> H2Ot_mu;    // Population mean H2Ot
  real<lower=0.01> H2Ot_sigma;            // Population SD H2Ot

  real<lower=79000,upper=86000> Ea_mu;   // Population mean Ea
  real<lower=100> Ea_sigma;               // Population SD Ea

  // NON-CENTERED: Raw parameters on unit interval (proper truncation via inverse CDF)
  vector<lower=0.001,upper=0.999>[N] H2Ot_unit;
  vector<lower=0.001,upper=0.999>[N] Ea_unit;
}

transformed parameters {
  // Transform unit interval to truncated normal
  vector[N] H2Ot;
  vector[N] Ea;

  for (i in 1:N) {
    H2Ot[i] = unit_to_truncnorm(H2Ot_unit[i], H2Ot_mu, H2Ot_sigma, H2Ot_lower, H2Ot_upper);
    Ea[i] = unit_to_truncnorm(Ea_unit[i], Ea_mu, Ea_sigma, Ea_lower, Ea_upper);
  }
}

model {
  // Priors on island-wide composition parameters
  // Based on observed range in Stevenson 2015: H2Ot 0.05-0.43 wt%, Ea 79.6-85.9 kJ/mol
  H2Ot_mu ~ normal(0.18, 0.08);      // Centered on observed mean (0.052-0.432 range)
  H2Ot_sigma ~ exponential(10);       // Mean 0.1 wt%; Rogers & Stevenson (2022): within-source CV 20-40%
  Ea_mu ~ normal(84000, 2000);        // Centered on observed mean (79642-85858 range)
  Ea_sigma ~ exponential(0.0005);     // Mean 2000 J/mol; weakly informative (observed SD = 532)

  // Implicit uniform prior on unit interval parameters
  // The inverse CDF transformation gives proper truncated normal

  // Age priors (weakly informative for Rapa Nui range)
  age ~ normal(500, 200);

  // Likelihood: H2Ome measurements given ages and parameters
  for (i in 1:N) {
    // Final rate equation (R^2=1.000 fit to Stevenson 2015)
    // log(rate) = a + b*(-Ea/(R*T)) + c*H2Ot + d*(1/T)
    // NOTE: No RH term in Tier 1 (EHT/RH treated as known, using Stevenson's RH=0.98)
    real log_rate = a_eq + b_eq * (-Ea[i] / (R_gas * T_K[i])) + c_eq * H2Ot[i] + d_eq * (1.0 / T_K[i]);

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
    real log_rate = a_eq + b_eq * (-Ea[i] / (R_gas * T_K[i])) + c_eq * H2Ot[i] + d_eq * (1.0 / T_K[i]);

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
