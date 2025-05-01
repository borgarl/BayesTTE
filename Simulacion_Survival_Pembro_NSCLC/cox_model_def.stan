data {
  int<lower=0> N; // Nº patients
  int<lower=0,upper=1> t[N]; // Binary treatment variable
  vector[N] y;    // Survival times
  vector[N] v;    // Event status: 1 for event, 0 for censored
  real mu_intercept0;
  real sigma_intercept0;
  real mu_slope0;
  real sigma_slope0;
}

parameters {
  real intercept;
  real slope;
}

transformed parameters {
  vector[N] linear_predictor;
  for (i in 1:N) {
    linear_predictor[i] = intercept + t[i] * slope;
  }
}

model {
  // Priors
  intercept ~ normal(mu_intercept0, sigma_intercept0);
  slope ~ normal(mu_slope0, sigma_slope0);

  // Likelihood
  for (i in 1:N) {
    if (v[i] == 1) {
      target += linear_predictor[i] - log_sum_exp(linear_predictor);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  
  for (i in 1:N) {
    if (v[i] == 1) {
      log_lik[i] = linear_predictor[i] - log_sum_exp(linear_predictor);
    } else {
      log_lik[i] = 0;
    }
  }
}
