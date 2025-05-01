data {
  int<lower=0> N; // Nº pacientes
  int<lower=0,upper=1> t[N]; // Binary treatment variable
  vector[N] y;    // vector de tiempos
  vector[N] v;    // vector status: evento o censura
  real<lower=0> mu_shape0;
  real<lower=0> sigma_shape0;
  real mu_scale0[2];
  real sigma_scale0[2];
}
transformed data {
  real<lower=0> alpha_shape0;  // Parameters for the Gamma prior on shape
  real<lower=0> beta_shape0;
  real alpha_scale0[2];  // Parameters for the Gamma prior on scale
  real beta_scale0[2];
  // Transform means and SDs to shape and rate parameters
  alpha_shape0 = (mu_shape0 / sigma_shape0)^2;
  beta_shape0 = mu_shape0 / (sigma_shape0^2);
  for (i in 1:2) {
    alpha_scale0[i] = (mu_scale0[i] / sigma_scale0[i])^2;
    beta_scale0[i] = mu_scale0[i] / (sigma_scale0[i]^2);
  }
}
parameters {
  real<lower=0> shape;
  real<lower=0> scale[2];
}
transformed parameters {
  vector[N] linear_predictor;
  for (i in 1:N) {
    linear_predictor[i] = shape / scale[t[i]+1];  // lambda = k / lambda
  }
}
model {
  // Priors
  shape ~ gamma(alpha_shape0, beta_shape0);
  scale ~ gamma(alpha_scale0, beta_scale0);
  for (i in 1:N) {
    if (v[i] == 1) {
      target += weibull_lpdf(y[i] | shape, linear_predictor[i]);
    } else {
      target += weibull_lccdf(y[i] | shape, linear_predictor[i]);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  real hr;
  for (i in 1:N) {
    if (v[i] == 1) {
      log_lik[i] = weibull_lpdf(y[i] | shape, linear_predictor[i]);
    } else {
      log_lik[i] = weibull_lccdf(y[i] | shape, linear_predictor[i]);
    }
  }
  hr = pow((scale[2] / scale[1]), shape);  // Hazard ratio
}

