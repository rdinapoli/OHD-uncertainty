// radiocarbon_ordered.stan
//
// Purpose: Hard-ordered phase model for radiocarbon dates only
// Samples constrained to stratigraphic phases with ordered boundaries.
// Uses correct calibration density interpolation for likelihood.
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
  // Interpolate calibration density at a given calendar age
  // Returns the probability density (not log) at the query age
  // cal_ages should be sorted in ascending order
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
  int<lower=1> N;                                      // Number of C14 samples
  int<lower=2> N_phases;                               // Number of phases
  array[N] int<lower=1, upper=N_phases> phase_id;      // Phase assignment for each sample
  int<lower=1> N_grid;                                 // Calibration grid length
  array[N] vector[N_grid] cal_ages;                    // Calendar ages for calibration (BP)
  array[N] vector[N_grid] cal_probs;                   // Calibration probability densities
  real<lower=0> age_min;                               // Youngest possible age (smallest BP)
  real<lower=age_min> age_max;                         // Oldest possible age (largest BP)
}

parameters {
  // Boundaries in DECREASING order (older to younger in BP)
  positive_ordered[N_phases-1] boundary_raw;

  // Position of each event within its phase (0 = oldest end, 1 = youngest end)
  vector<lower=0, upper=1>[N] event_position;
}

transformed parameters {
  vector[N_phases-1] phase_boundary;
  vector[N] true_age;

  // Scale raw boundaries to [age_min, age_max] range
  {
    real max_raw = boundary_raw[N_phases-1] + 1.0;
    for (p in 1:(N_phases-1)) {
      real normalized = boundary_raw[N_phases - p] / max_raw;
      phase_boundary[p] = age_min + normalized * (age_max - age_min);
    }
  }

  // Compute true_age for each sample based on phase and event_position
  for (i in 1:N) {
    int p = phase_id[i];
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

    true_age[i] = phase_start - event_position[i] * (phase_start - phase_end);
  }
}

model {
  // Priors on raw boundaries
  boundary_raw ~ exponential(0.1);

  // Uniform prior on event positions
  event_position ~ uniform(0, 1);

  // Likelihood: calibration density evaluated at true_age
  for (i in 1:N) {
    real density = interp_calibration_density(cal_ages[i], cal_probs[i], true_age[i]);
    target += log(density);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N_phases] phase_duration;

  for (i in 1:N) {
    real density = interp_calibration_density(cal_ages[i], cal_probs[i], true_age[i]);
    log_lik[i] = log(density);
  }

  phase_duration[1] = age_max - phase_boundary[1];
  for (p in 2:(N_phases-1)) {
    phase_duration[p] = phase_boundary[p-1] - phase_boundary[p];
  }
  phase_duration[N_phases] = phase_boundary[N_phases-1] - age_min;
}
