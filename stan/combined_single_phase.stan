// combined_single_phase.stan
//
// Purpose: Single population model combining C14 and OHD dates
// All samples assumed to come from one temporal population.
// Uses correct calibration density interpolation for C14 likelihood.
//
// Phase Convention:
//   N/A - single population model, no phases
//
// Date: 2025-12-10
// Reference: stan_model_rewrite_instructions.md

functions {
  // Interpolate calibration density at a given calendar age
  real interp_calibration_density(vector cal_ages, vector cal_probs, real query_age) {
    int N = num_elements(cal_ages);

    if (query_age <= cal_ages[1]) {
      return fmax(cal_probs[1], 1e-10);
    }
    if (query_age >= cal_ages[N]) {
      return fmax(cal_probs[N], 1e-10);
    }

    int lo = 1;
    int hi = N;
    while (hi - lo > 1) {
      int mid = (lo + hi) / 2;
      if (cal_ages[mid] <= query_age) {
        lo = mid;
      } else {
        hi = mid;
      }
    }

    real t = (query_age - cal_ages[lo]) / (cal_ages[hi] - cal_ages[lo]);
    real density = cal_probs[lo] + t * (cal_probs[hi] - cal_probs[lo]);

    return fmax(density, 1e-10);
  }
}

data {
  // C14 data
  int<lower=0> N_c14;                          // Number of C14 samples
  int<lower=1> N_grid;                         // Calibration grid length
  array[N_c14] vector[N_grid] cal_ages;        // Calendar ages for calibration (BP)
  array[N_c14] vector[N_grid] cal_probs;       // Calibration probability densities

  // OHD data
  int<lower=0> N_ohd;                          // Number of OHD samples
  vector[N_ohd] OHD_age_BP;                    // OHD point estimates (BP)
  vector<lower=0>[N_ohd] OHD_SD;               // OHD uncertainties (1 SD)

  // Age bounds
  real<lower=0> age_min;
  real<lower=age_min> age_max;
}

parameters {
  real<lower=age_min, upper=age_max> population_mean;
  real<lower=0> population_sd;
  vector<lower=age_min, upper=age_max>[N_c14] true_age_c14;
  vector<lower=age_min, upper=age_max>[N_ohd] true_age_ohd;
}

model {
  // Priors
  population_mean ~ uniform(age_min, age_max);
  population_sd ~ normal(0, 200);

  // Hierarchical prior: all true ages from same population
  true_age_c14 ~ normal(population_mean, population_sd);
  true_age_ohd ~ normal(population_mean, population_sd);

  // C14 likelihood: calibration density interpolation
  for (i in 1:N_c14) {
    real density = interp_calibration_density(cal_ages[i], cal_probs[i], true_age_c14[i]);
    target += log(density);
  }

  // OHD likelihood: normal
  OHD_age_BP ~ normal(true_age_ohd, OHD_SD);
}

generated quantities {
  vector[N_c14 + N_ohd] log_lik;

  // C14 log-likelihoods
  for (i in 1:N_c14) {
    real density = interp_calibration_density(cal_ages[i], cal_probs[i], true_age_c14[i]);
    log_lik[i] = log(density);
  }

  // OHD log-likelihoods
  for (i in 1:N_ohd) {
    log_lik[N_c14 + i] = normal_lpdf(OHD_age_BP[i] | true_age_ohd[i], OHD_SD[i]);
  }
}
