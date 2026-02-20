// ohd_ordered.stan
//
// Purpose: Hard-ordered phase model for OHD dates only
// Samples constrained to stratigraphic phases with ordered boundaries.
//
// Phase Convention:
//   Phase 1 = OLDEST (deepest stratigraphy, highest BP values)
//   Phase N_phases = YOUNGEST (shallowest stratigraphy, lowest BP values)
//   phase_boundary[1] separates Phase 1 from Phase 2
//   phase_boundary[p] separates Phase p from Phase p+1
//
// Date: 2025-12-10
// Reference: stan_model_rewrite_instructions.md

data {
  int<lower=1> N;                                      // Number of OHD samples
  int<lower=2> N_phases;                               // Number of phases
  array[N] int<lower=1, upper=N_phases> phase_id;      // Phase assignment for each sample
  vector[N] OHD_age_BP;                                // OHD point estimates (BP)
  vector<lower=0>[N] OHD_SD;                           // OHD uncertainties (1 SD)
  real<lower=0> age_min;                               // Youngest possible age (smallest BP)
  real<lower=age_min> age_max;                         // Oldest possible age (largest BP)
}

parameters {
  // Boundaries in DECREASING order (older to younger in BP)
  // Using positive_ordered which gives increasing values, then we reverse
  positive_ordered[N_phases-1] boundary_raw;

  // Position of each event within its phase (0 = oldest end, 1 = youngest end)
  vector<lower=0, upper=1>[N] event_position;
}

transformed parameters {
  vector[N_phases-1] phase_boundary;
  vector[N] true_age;

  // Scale raw boundaries to [age_min, age_max] range
  // positive_ordered gives values in increasing order
  // We want phase_boundary[1] > phase_boundary[2] > ... (decreasing in BP)
  {
    real max_raw = boundary_raw[N_phases-1] + 1.0;  // Ensure normalization doesn't collapse
    for (p in 1:(N_phases-1)) {
      // Map boundary_raw[N_phases-p] (largest to smallest) to phase_boundary[p]
      real normalized = boundary_raw[N_phases - p] / max_raw;
      phase_boundary[p] = age_min + normalized * (age_max - age_min);
    }
  }

  // Compute true_age for each sample based on phase and event_position
  for (i in 1:N) {
    int p = phase_id[i];
    real phase_start;  // Older (larger BP) end of phase
    real phase_end;    // Younger (smaller BP) end of phase

    if (p == 1) {
      // Phase 1 (oldest): from age_max to phase_boundary[1]
      phase_start = age_max;
      phase_end = phase_boundary[1];
    } else if (p == N_phases) {
      // Phase N (youngest): from phase_boundary[N_phases-1] to age_min
      phase_start = phase_boundary[N_phases-1];
      phase_end = age_min;
    } else {
      // Middle phases: from phase_boundary[p-1] to phase_boundary[p]
      phase_start = phase_boundary[p-1];
      phase_end = phase_boundary[p];
    }

    // Interpolate within phase: event_position=0 -> phase_start, event_position=1 -> phase_end
    true_age[i] = phase_start - event_position[i] * (phase_start - phase_end);
  }
}

model {
  // Priors on raw boundaries (weakly informative)
  boundary_raw ~ exponential(0.1);

  // Uniform prior on event positions within phases
  event_position ~ uniform(0, 1);

  // Likelihood: OHD observations given true ages
  OHD_age_BP ~ normal(true_age, OHD_SD);
}

generated quantities {
  vector[N] log_lik;
  vector[N_phases] phase_duration;

  // Log-likelihood for LOO-CV
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(OHD_age_BP[i] | true_age[i], OHD_SD[i]);
  }

  // Phase durations
  phase_duration[1] = age_max - phase_boundary[1];
  for (p in 2:(N_phases-1)) {
    phase_duration[p] = phase_boundary[p-1] - phase_boundary[p];
  }
  phase_duration[N_phases] = phase_boundary[N_phases-1] - age_min;
}
