// combined_ordered.stan
//
// Purpose: Hard-ordered phase model combining C14 and OHD dates
// Samples constrained to stratigraphic phases with ordered boundaries.
// Uses correct calibration density interpolation for C14 likelihood.
//
// Phase Convention:
//   Phase 1 = OLDEST (deepest stratigraphy, highest BP values)
//   Phase N_phases = YOUNGEST (shallowest stratigraphy, lowest BP values)
//   phase_boundary[1] separates Phase 1 from Phase 2
//   phase_boundary[p] separates Phase p from Phase p+1
//
// Date: 2025-12-10
// Reference: stan_model_rewrite_instructions.md

functions {
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
  int<lower=0> N_c14;
  int<lower=1> N_grid;
  array[N_c14] vector[N_grid] cal_ages;
  array[N_c14] vector[N_grid] cal_probs;
  array[N_c14] int<lower=1> phase_id_c14;

  // OHD data
  int<lower=0> N_ohd;
  vector[N_ohd] OHD_age_BP;
  vector<lower=0>[N_ohd] OHD_SD;
  array[N_ohd] int<lower=1> phase_id_ohd;

  // Phase structure
  int<lower=2> N_phases;

  // Age bounds
  real<lower=0> age_min;
  real<lower=age_min> age_max;
}

parameters {
  positive_ordered[N_phases-1] boundary_raw;
  vector<lower=0, upper=1>[N_c14] event_position_c14;
  vector<lower=0, upper=1>[N_ohd] event_position_ohd;
}

transformed parameters {
  vector[N_phases-1] phase_boundary;
  vector[N_c14] true_age_c14;
  vector[N_ohd] true_age_ohd;

  // Scale raw boundaries to [age_min, age_max] range
  {
    real max_raw = boundary_raw[N_phases-1] + 1.0;
    for (p in 1:(N_phases-1)) {
      real normalized = boundary_raw[N_phases - p] / max_raw;
      phase_boundary[p] = age_min + normalized * (age_max - age_min);
    }
  }

  // Compute true_age for C14 samples
  for (i in 1:N_c14) {
    int p = phase_id_c14[i];
    real phase_start;
    real phase_end;

    if (p == 1) {
      phase_start = age_max;
      phase_end = phase_boundary[1];
    } else if (p == N_phases) {
      phase_start = phase_boundary[N_phases-1];
      phase_end = age_min;
    } else {
      phase_start = phase_boundary[p-1];
      phase_end = phase_boundary[p];
    }

    true_age_c14[i] = phase_start - event_position_c14[i] * (phase_start - phase_end);
  }

  // Compute true_age for OHD samples
  for (i in 1:N_ohd) {
    int p = phase_id_ohd[i];
    real phase_start;
    real phase_end;

    if (p == 1) {
      phase_start = age_max;
      phase_end = phase_boundary[1];
    } else if (p == N_phases) {
      phase_start = phase_boundary[N_phases-1];
      phase_end = age_min;
    } else {
      phase_start = phase_boundary[p-1];
      phase_end = phase_boundary[p];
    }

    true_age_ohd[i] = phase_start - event_position_ohd[i] * (phase_start - phase_end);
  }
}

model {
  // Priors
  boundary_raw ~ exponential(0.1);
  event_position_c14 ~ uniform(0, 1);
  event_position_ohd ~ uniform(0, 1);

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
  vector[N_phases] phase_duration;

  // C14 log-likelihoods
  for (i in 1:N_c14) {
    real density = interp_calibration_density(cal_ages[i], cal_probs[i], true_age_c14[i]);
    log_lik[i] = log(density);
  }

  // OHD log-likelihoods
  for (i in 1:N_ohd) {
    log_lik[N_c14 + i] = normal_lpdf(OHD_age_BP[i] | true_age_ohd[i], OHD_SD[i]);
  }

  // Phase durations
  phase_duration[1] = age_max - phase_boundary[1];
  for (p in 2:(N_phases-1)) {
    phase_duration[p] = phase_boundary[p-1] - phase_boundary[p];
  }
  phase_duration[N_phases] = phase_boundary[N_phases-1] - age_min;
}
