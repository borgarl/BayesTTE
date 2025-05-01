data {
  int<lower=0> N; // Number of patients
  int<lower=0,upper=1> t[N]; // Treatment arm (0 = control, 1 = treatment)
  vector[N] y; // Survival times
  int<lower=0,upper=1> v[N]; // Event status (1 = event, 0 = censored)
  real mu_slope0; // Was mu_beta0
  real sigma_slope0; // Was sigma_beta0
  real mu_intercept0; // Was mu_baseline_hazard0
  real sigma_intercept0; // Was sigma_baseline_hazard0
}

parameters {
  real beta; // Corresponds to slope
  real baseline_hazard; // Corresponds to intercept
}

model {
  // Priors
  beta ~ normal(mu_slope0, sigma_slope0); // Updated variable names
  baseline_hazard ~ normal(mu_intercept0, sigma_intercept0); // Updated variable names

  // Likelihood
  for (i in 1:N) {
    real hazard = exp(baseline_hazard + t[i] * beta);
    if (v[i] == 1) {
      target += log(hazard) - hazard * y[i];
    } else {
      target += - hazard * y[i];
    }
  }
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    real hazard = exp(baseline_hazard + t[i] * beta);
    if (v[i] == 1) {
      log_lik[i] = log(hazard) - hazard * y[i];
    } else {
      log_lik[i] = -hazard * y[i];
    }
  }
}

