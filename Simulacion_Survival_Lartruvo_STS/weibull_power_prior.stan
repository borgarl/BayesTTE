data {
  int<lower=0> N; 
  int<lower=0> N_external; 
  int<lower=0,upper=1> t[N];
  int<lower=0,upper=1> t_external[N_external]; 
  vector[N] y; 
  vector[N_external] y_external; 
  vector[N] v; 
  vector[N_external] v_external; 
  real<lower=0> alpha_shape0; 
  real<lower=0> beta_shape0;
  real<lower=0> alpha_scale0_control;
  real<lower=0> beta_scale0_control;
  real<lower=0> alpha_scale0_treatment;
  real<lower=0> beta_scale0_treatment;
  real<lower=0, upper=1> a0;  // a0 weight power prior
}

parameters {
  real<lower=0> shape;
  real<lower=0> scale_control;
  real<lower=0> scale_treatment;
}

model {
  // Priors
  shape ~ uniform(0.1, 4);
  scale_control ~ gamma(alpha_scale0_control, beta_scale0_control);
  scale_treatment ~ gamma(alpha_scale0_treatment, beta_scale0_treatment);

  // Likelihood for the external control data, weighted by a0
  for (j in 1:N_external) {
    if (v_external[j] == 1) {
      target += a0 * weibull_lpdf(y_external[j] | shape, scale_control);
    } else {
      target += a0 * weibull_lccdf(y_external[j] | shape, scale_control);
    }
  }

  // Likelihood for the current trial data
  for (i in 1:N) {
    if (t[i] == 0) {  // Control group
      if (v[i] == 1) {
        target += weibull_lpdf(y[i] | shape, scale_control);
      } else {
        target += weibull_lccdf(y[i] | shape, scale_control);
      }
    } else {  // Treatment group
      if (v[i] == 1) {
        target += weibull_lpdf(y[i] | shape, scale_treatment);
      } else {
        target += weibull_lccdf(y[i] | shape, scale_treatment);
      }
    }
  }
}

generated quantities {
  real HR;
  vector[N] log_lik;
  real used_a0 = a0; 

  HR = (scale_control / scale_treatment) ^ shape;

  for (i in 1:N) {
    if (t[i] == 0) {  // Control group
      if (v[i] == 1) {
        log_lik[i] = weibull_lpdf(y[i] | shape, scale_control);
      } else {
        log_lik[i] = weibull_lccdf(y[i] | shape, scale_control);
      }
    } else {  // Treatment group
      if (v[i] == 1) {
        log_lik[i] = weibull_lpdf(y[i] | shape, scale_treatment);
      } else {
        log_lik[i] = weibull_lccdf(y[i] | shape, scale_treatment);
      }
    }
  }
}
