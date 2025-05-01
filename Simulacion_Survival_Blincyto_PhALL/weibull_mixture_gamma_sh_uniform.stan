data {
  int<lower=0> N; 
  int<lower=0,upper=1> t[N]; 
  vector[N] y; 
  vector[N] v; 
  real<lower=0> alpha_shape0; 
  real<lower=0> beta_shape0;
  real<lower=0> alpha_scale0_control;
  real<lower=0> beta_scale0_control;
  real<lower=0> alpha_scale0_treatment;
  real<lower=0> beta_scale0_treatment;
  real<lower=0, upper=1> w;  
}

parameters {
  real<lower=0> shape;
  real<lower=0> scale_control;
  real<lower=0> scale_treatment;
}

model {
  // Priors
  shape ~ uniform(0.1, 4);
  
  // Mixture prior for scale_control
  target += log_mix(w,
                    gamma_lpdf(scale_control | alpha_scale0_control, beta_scale0_control),
                    gamma_lpdf(scale_control | 1, 0.0625)); 
  
  scale_treatment ~ gamma(alpha_scale0_treatment, beta_scale0_treatment);
  
  // Likelihood
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
