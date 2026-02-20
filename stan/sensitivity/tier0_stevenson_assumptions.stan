// Tier 0: Stevenson's Exact Assumptions (Irreducible Floor Test)
//
// This model uses Stevenson's EXACT parameter values with NO uncertainty:
//   - EHT: Fixed at study-area constants (22.4, 21.8, 23.3°C)
//   - RH: Fixed at 0.98 for all samples
//   - H₂Ot: Fixed at measured values per sample
//   - Ea: Fixed at measured values per sample
//
// Purpose: Test the irreducible uncertainty floor even under Stevenson's assumptions.
// If uncertainty is still large, the problem is fundamental (rate equation + measurement).
//
// CREATED 2025-12-03: Prior sensitivity analysis

data {
  int<lower=1> N;              // Number of samples
  int<lower=1> N_areas;        // Number of study areas (3)

  // Observed data
  vector<lower=0>[N] H2Ome;    // Molecular water content (wt%)

  // Study area assignments
  array[N] int<lower=1,upper=N_areas> study_area;   // Study area ID per sample

  // FIXED VALUES FROM STEVENSON (passed as data, not parameters)
  vector[N] H2Ot_fixed;        // Stevenson's measured H₂Ot per sample
  vector[N] Ea_fixed;          // Stevenson's measured Ea per sample
  vector[N_areas] EHT_fixed;   // Stevenson's EHT per study area (22.4, 21.8, 23.3)
  real RH_fixed;               // Stevenson's universal RH (0.98)

  // Constants
  real R_gas;                  // Gas constant (8.314 J/(mol·K))

  // Final rate equation coefficients
  real a_eq;
  real b_eq;
  real c_eq;
  real d_eq;
}

parameters {
  // ONLY age is unknown - everything else is fixed
  vector<lower=100,upper=900>[N] age;  // Age in years BP
}

model {
  // Age prior (weakly informative)
  age ~ normal(500, 200);

  // Likelihood
  for (i in 1:N) {
    // Use FIXED values from Stevenson
    real H2Ot = H2Ot_fixed[i];
    real Ea = Ea_fixed[i];
    real EHT = EHT_fixed[study_area[i]];
    real RH = RH_fixed;
    real T_K = EHT + 273.15;

    // Rate equation (no RH term since RH = 0.98 exactly)
    real log_rate = a_eq + b_eq * (-Ea / (R_gas * T_K)) + c_eq * H2Ot + d_eq * (1.0 / T_K);

    // Numerical stability
    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;

    real rate = exp(log_rate);
    real H2Ome_pred = sqrt(rate * age[i]);

    // Measurement model - this is the ONLY source of uncertainty
    H2Ome[i] ~ normal(H2Ome_pred, 0.03);
  }
}

generated quantities {
  vector[N] H2Ome_pred;
  vector[N] log_lik;
  vector[N] age_residual;  // Difference from Stevenson's reported ages

  for (i in 1:N) {
    real H2Ot = H2Ot_fixed[i];
    real Ea = Ea_fixed[i];
    real EHT = EHT_fixed[study_area[i]];
    real T_K = EHT + 273.15;

    real log_rate = a_eq + b_eq * (-Ea / (R_gas * T_K)) + c_eq * H2Ot + d_eq * (1.0 / T_K);

    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;

    real rate = exp(log_rate);
    H2Ome_pred[i] = sqrt(rate * age[i]);
    log_lik[i] = normal_lpdf(H2Ome[i] | H2Ome_pred[i], 0.03);
  }
}
