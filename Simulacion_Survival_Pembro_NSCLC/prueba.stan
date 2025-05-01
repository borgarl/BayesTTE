data {
  int<lower=0> N; // Number of patients
  int<lower=0,upper=1> t[N]; // Binary treatment variable
  vector[N] y;    // vector of times
  vector[N] v;    // vector of event status
  real<lower=0> alpha_shape[2]; // Shape for Gamma prior on shape parameters
  real<lower=0> beta_shape[2];  // Rate for Gamma prior on shape parameters
  real<lower=0> alpha_scale[2]; // Shape for Gamma prior on scale parameters
  real<lower=0> beta_scale[2];  // Rate for Gamma prior on scale parameters
}

parameters {
  real<lower=0> shape[2]; // Shape parameters for control and treatment
  real<lower=0> scale[2]; // Scale parameters for control and treatment
}

model {
  // Priors
  shape[1] ~ gamma(alpha_shape[1], beta_shape[1]);
  scale[1] ~ gamma(alpha_scale[1], beta_scale[1]);
  shape[2] ~ gamma(alpha_shape[2], beta_shape[2]);
  scale[2] ~ gamma(alpha_scale[2], beta_scale[2]);

  // Likelihood
  for (i in 1:N) {
    if (v[i] == 1) {
      target += weibull_lpdf(y[i] | shape[t[i]+1], scale[t[i]+1]);
    } else {
      target += weibull_lccdf(y[i] | shape[t[i]+1], scale[t[i]+1]);
    }
  }
}

generated quantities {
  real hr;
  hr = (shape[2] / shape[1]) * (scale[1] / scale[2]);
}


