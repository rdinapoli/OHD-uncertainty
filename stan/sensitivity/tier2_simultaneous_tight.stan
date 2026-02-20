// Tier 2 with ALL Parameters Simultaneously Constrained
//
// The coup de grâce test: What happens when we constrain ALL parameters
// to their tightest values simultaneously?
//
// If this still produces ±250+ years, even physically impossible joint
// precision cannot recover Stevenson's claims.
//
// Parameters fixed at:
//   - EHT_sd = 0.1°C (30× tighter than baseline)
//   - H2Ot_sigma = 0.005 wt% (16× tighter than baseline)
//   - Ea_sigma = 100 J/mol (20× tighter than baseline)
//
// This tests the "fix the calibration" response: if better measurements
// of controlling variables could fix OHD, this model would show it.

data {
  int<lower=1> N;
  int<lower=1> N_areas;
  vector<lower=0>[N] H2Ome;
  array[N] int<lower=1,upper=N_areas> study_area;
  real R_gas;
  real a_eq;
  real b_eq;
  real c_eq;
  real d_eq;
  real e_RH_mean;
  real<lower=0> e_RH_sd;

  // ALL THREE SENSITIVITY PARAMETERS PASSED AS DATA
  real<lower=0> EHT_sd_prior;
  real<lower=0> H2Ot_sigma_prior;
  real<lower=0> Ea_sigma_prior;
}

parameters {
  vector<lower=100,upper=900>[N] age;

  // Population means (SDs are fixed)
  real<lower=0.05,upper=0.5> H2Ot_mu;
  real<lower=79000,upper=86000> Ea_mu;

  // Raw parameters for non-centered parameterization
  vector[N] H2Ot_raw;
  vector[N] Ea_raw;

  // Environmental parameters
  vector<lower=18,upper=28>[N_areas] EHT_mean;
  vector<lower=0.7,upper=1.0>[N_areas] RH_mean;
  vector<lower=0>[N_areas] RH_sd;

  vector[N] EHT_raw;
  vector[N] RH_raw;

  real<lower=0> e_RH;
}

transformed parameters {
  vector<lower=0.05,upper=0.43>[N] H2Ot;
  vector<lower=75000,upper=92000>[N] Ea;
  vector<lower=18,upper=28>[N] EHT;
  vector<lower=0.7,upper=1.0>[N] RH;
  vector[N] T_K;

  for (i in 1:N) {
    // ALL use fixed sigmas from data
    H2Ot[i] = H2Ot_mu + H2Ot_sigma_prior * H2Ot_raw[i];
    Ea[i] = Ea_mu + Ea_sigma_prior * Ea_raw[i];
    EHT[i] = EHT_mean[study_area[i]] + EHT_sd_prior * EHT_raw[i];
    RH[i] = RH_mean[study_area[i]] + RH_sd[study_area[i]] * RH_raw[i];

    // Bounds
    if (H2Ot[i] < 0.05) H2Ot[i] = 0.05;
    if (H2Ot[i] > 0.43) H2Ot[i] = 0.43;
    if (Ea[i] < 75000) Ea[i] = 75000;
    if (Ea[i] > 92000) Ea[i] = 92000;
    if (EHT[i] < 18) EHT[i] = 18;
    if (EHT[i] > 28) EHT[i] = 28;
    if (RH[i] < 0.7) RH[i] = 0.7;
    if (RH[i] > 1.0) RH[i] = 1.0;

    T_K[i] = EHT[i] + 273.15;
  }
}

model {
  H2Ot_mu ~ normal(0.18, 0.08);
  Ea_mu ~ normal(84000, 2000);

  H2Ot_raw ~ std_normal();
  Ea_raw ~ std_normal();

  EHT_mean ~ normal(22.5, 2.5);

  for (k in 1:N_areas) {
    RH_mean[k] ~ normal(0.85, 0.1);
  }
  RH_sd ~ exponential(5);

  EHT_raw ~ std_normal();
  RH_raw ~ std_normal();

  age ~ normal(500, 200);
  e_RH ~ normal(e_RH_mean, e_RH_sd);

  for (i in 1:N) {
    real log_rate = a_eq + b_eq * (-Ea[i] / (R_gas * T_K[i])) + c_eq * H2Ot[i] + d_eq * (1.0 / T_K[i]) + e_RH * (RH[i] - 0.98);
    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;
    real rate = exp(log_rate);
    real H2Ome_pred = sqrt(rate * age[i]);
    H2Ome[i] ~ normal(H2Ome_pred, 0.03);
  }
}

generated quantities {
  vector[N] H2Ome_pred;
  vector[N] log_lik;

  // Store all fixed values for output
  real EHT_sd_used = EHT_sd_prior;
  real H2Ot_sigma_used = H2Ot_sigma_prior;
  real Ea_sigma_used = Ea_sigma_prior;

  for (i in 1:N) {
    real log_rate = a_eq + b_eq * (-Ea[i] / (R_gas * T_K[i])) + c_eq * H2Ot[i] + d_eq * (1.0 / T_K[i]) + e_RH * (RH[i] - 0.98);
    if (log_rate > 50) log_rate = 50;
    if (log_rate < -50) log_rate = -50;
    real rate = exp(log_rate);
    H2Ome_pred[i] = sqrt(rate * age[i]);
    log_lik[i] = normal_lpdf(H2Ome[i] | H2Ome_pred[i], 0.03);
  }
}
