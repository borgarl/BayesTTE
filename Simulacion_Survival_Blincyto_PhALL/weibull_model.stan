data {
  int<lower=0> N; // Nº pacientes
  int<lower=0,upper=1> t[N]; // Binary treatment variable
  vector[N] y;    // vector de tiempos
  vector[N] v;    // vector status: evento o censura
  real<lower=0> mu_shape0;
  real<lower=0> sigma_shape0;
  real mu_scale0;
  real sigma_scale0;
  real<lower=0> mu_intercept0;
  real<lower=0> sigma_intercept0;
  real mu_slope0;
  real sigma_slope0;
}

parameters {
  real shape;
  real intercept;
  real slope;
  real scale;
}

model {
  // Priors
  shape ~ normal(mu_shape0, sigma_shape0);
  intercept ~ normal(mu_intercept0, sigma_intercept0);
  slope ~ normal(mu_slope0, sigma_slope0);
  scale ~ normal(mu_scale0, sigma_scale0);

  for (i in 1:N) {
    real linear_predictor = intercept + t[i] * slope;

    if (v[i] == 1) {
      target += weibull_lpdf(y[i] | shape, exp(linear_predictor) * scale);
    } else {
      target += weibull_lccdf(y[i] | shape, exp(linear_predictor) * scale);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    real linear_predictor = intercept + t[i] * slope;

    if (v[i] == 1) {
      log_lik[i] = weibull_lpdf(y[i] | shape, exp(linear_predictor) * scale);
    } else {
      log_lik[i] = weibull_lccdf(y[i] | shape, exp(linear_predictor) * scale);
    }
  }
}

