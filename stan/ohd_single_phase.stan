// ohd_single_phase.stan
//
// Purpose: Single population model for OHD dates only
// All samples assumed to come from one temporal population.
//
// Phase Convention:
//   N/A - single population model, no phases
//
// Date: 2025-12-10
// Reference: stan_model_rewrite_instructions.md

data {
  int<lower=1> N;                        // Number of OHD samples
  vector[N] OHD_age_BP;                  // OHD point estimates (BP)
  vector<lower=0>[N] OHD_SD;             // OHD uncertainties (1 SD)
  real<lower=0> age_min;                 // Youngest possible age (smallest BP)
  real<lower=age_min> age_max;           // Oldest possible age (largest BP)
}

parameters {
  real<lower=age_min, upper=age_max> population_mean;  // Population mean age
  real<lower=0> population_sd;                          // Population SD
  vector<lower=age_min, upper=age_max>[N] true_age;    // Latent true ages
}

model {
  // Priors
  population_mean ~ uniform(age_min, age_max);
  population_sd ~ normal(0, 200);  // Weakly informative, allows wide populations

  // Hierarchical prior: true ages come from population
  true_age ~ normal(population_mean, population_sd);

  // Likelihood: OHD observations given true ages
  OHD_age_BP ~ normal(true_age, OHD_SD);
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = normal_lpdf(OHD_age_BP[i] | true_age[i], OHD_SD[i]);
  }
}
