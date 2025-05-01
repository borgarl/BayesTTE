data {
  int<lower=0> N; // Nº pacientes
  int<lower=0,upper=1> t[N]; // Binary treatment variable
  vector[N] y;    // vector de tiempos
  vector[N] v;    // vector status: evento o censura
  real<lower=0> mu_shape0;
  real<lower=0> sigma_shape0;
  real<lower=0> prior_param1_arm1; // Gamma prior parameter for control arm
  real<lower=0> prior_param2_arm1; // Gamma prior parameter for control arm
  real<lower=0> prior_param1_arm2; // Uniform prior lower bound for experimental arm
  real<lower=0> prior_param2_arm2; // Uniform prior upper bound for experimental arm
}

transformed data {
  real<lower=0> alpha_shape0;  // Parameters for the Gamma prior on shape
  real<lower=0> beta_shape0;

  // Transform means and SDs to shape and rate parameters for the shape prior
  alpha_shape0 = (mu_shape0 / sigma_shape0)^2;
  beta_shape0 = mu_shape0 / (sigma_shape0^2);
}

parameters {
  real<lower=0> shape;
  real<lower=0> scale[2];
}

model {
  // Priors
  shape ~ gamma(alpha_shape0, beta_shape0);
  scale[1] ~ gamma(prior_param1_arm1, prior_param2_arm1); // Gamma prior for control arm
  scale[2] ~ uniform(prior_param1_arm2, prior_param2_arm2); // Uniform prior for treatment arm
  for (i in 1:N) {
    if (v[i] == 1) {
      target += weibull_lpdf(y[i] | shape, y[i] / scale[t[i]+1]);
    } else {
      target += weibull_lccdf(y[i] | shape, y[i] / scale[t[i]+1]);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  real hr;
  real scale_log;
  real intercept;
  for (i in 1:N) {
    if (v[i] == 1) {
      log_lik[i] = weibull_lpdf(y[i] | shape, y[i] / scale[t[i]+1]);
    } else {
      log_lik[i] = weibull_lccdf(y[i] | shape, y[i] / scale[t[i]+1]);
    }
  }
  hr = exp((log(scale[1]) - log(scale[2]))*shape);  // Hazard ratio
  scale_log = -log(shape);
  intercept = log(scale[1]);
}

