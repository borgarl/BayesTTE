data {
  int<lower=0> N; // Nº pacientes
  int<lower=0,upper=1> t[N]; // Binary treatment variable
  vector[N] y;    // vector de tiempos
  vector[N] v;    // vector status: evento o censura
  real<lower=0> mu_shape0;
  real<lower=0> sigma_shape0;
  real prior_param1_arm1; 
  real prior_param2_arm1;
  real prior_param1_arm2;
  real prior_param2_arm2;
}

parameters {
  real<lower=0> shape;
  real log_scale1;
  real log_scale2;
}

transformed parameters {
  vector[N] linear_predictor;
  for (i in 1:N) {
    if(t[i] == 0)
      linear_predictor[i] = shape / exp(log_scale1);
    else
      linear_predictor[i] = shape / exp(log_scale2);
  }
}

model {
  // Priors
  shape ~ normal(mu_shape0, sigma_shape0);
  log_scale1 ~ normal(prior_param1_arm1, prior_param2_arm1);
  log_scale2 ~ normal(prior_param1_arm2, prior_param2_arm2);
  
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
  real scale_log;
  real intercept;
  for (i in 1:N) {
    if (v[i] == 1) {
      log_lik[i] = weibull_lpdf(y[i] | shape, linear_predictor[i]);
    } else {
      log_lik[i] = weibull_lccdf(y[i] | shape, linear_predictor[i]);
    }
  }
  hr = pow((exp(log_scale2) / exp(log_scale1)), shape);  // Hazard ratio
  scale_log = log(1 / shape);
  intercept = log_scale1;
}
