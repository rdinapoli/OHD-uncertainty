// radiocarbon_single_phase.stan
//
// Purpose: Single population model for radiocarbon dates only
// All samples assumed to come from one temporal population.
// Uses correct calibration density interpolation for likelihood.
//
// Phase Convention:
//   N/A - single population model, no phases
//
// Date: 2025-12-10
// Reference: stan_model_rewrite_instructions.md

functions {
  // Interpolate calibration density at a given calendar age
  // Returns the probability density (not log) at the query age
  // cal_ages should be sorted in ascending order (youngest BP to oldest BP)
  real interp_calibration_density(vector cal_ages, vector cal_probs, real query_age) {
    int N = num_elements(cal_ages);

    // Handle out-of-bounds: return small but non-zero density
    if (query_age <= cal_ages[1]) {
      return fmax(cal_probs[1], 1e-10);
    }
    if (query_age >= cal_ages[N]) {
      return fmax(cal_probs[N], 1e-10);
    }

    // Binary search for bracketing interval
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

    // Linear interpolation
    real t = (query_age - cal_ages[lo]) / (cal_ages[hi] - cal_ages[lo]);
    real density = cal_probs[lo] + t * (cal_probs[hi] - cal_probs[lo]);

    // Ensure positive density for numerical stability
    return fmax(density, 1e-10);
  }
}

data {
  int<lower=1> N;                           // Number of C14 samples
  int<lower=1> N_grid;                      // Calibration grid length
  array[N] vector[N_grid] cal_ages;         // Calendar ages for calibration (BP, ascending)
  array[N] vector[N_grid] cal_probs;        // Calibration probability densities (normalized)
  real<lower=0> age_min;                    // Youngest possible age (smallest BP)
  real<lower=age_min> age_max;              // Oldest possible age (largest BP)
}

parameters {
  real<lower=age_min, upper=age_max> population_mean;  // Population mean age
  real<lower=0> population_sd;                          // Population SD
  vector<lower=age_min, upper=age_max>[N] true_age;    // Latent true calendar ages
}

model {
  // Priors
  population_mean ~ uniform(age_min, age_max);
  population_sd ~ normal(0, 200);  // Weakly informative

  // Hierarchical prior: true ages come from population
  true_age ~ normal(population_mean, population_sd);

  // Likelihood: calibration density evaluated at true_age
  for (i in 1:N) {
    real density = interp_calibration_density(cal_ages[i], cal_probs[i], true_age[i]);
    target += log(density);
  }
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    real density = interp_calibration_density(cal_ages[i], cal_probs[i], true_age[i]);
    log_lik[i] = log(density);
  }
}
