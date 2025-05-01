data {
  int<lower=0> N; // number of patients
  int<lower=0,upper=1> t[N]; // treatment assignment
  vector<lower=0>[N] y; // survival times
  int<lower=0,upper=1> v[N]; // censoring status
  real<lower=0> prior_param1_arm1; // control arm scale prior mean
  real<lower=0> prior_param2_arm1; // control arm scale prior sd
  real<lower=0> prior_param1_arm2; // experimental arm scale prior mean
  real<lower=0> prior_param2_arm2; // experimental arm scale prior sd
  real<lower=0> mu_shape0; // shape prior mean
  real<lower=0> sigma_shape0; // shape prior sd
}

parameters {
  real<lower=0> shape;
  real<lower=0> scale[2];
}

model {
  shape ~ normal(mu_shape0, sigma_shape0) T[0,];
  scale[1] ~ normal(prior_param1_arm1, prior_param2_arm1) T[0,]; // control arm
  scale[2] ~ normal(prior_param1_arm2, prior_param2_arm2) T[0,]; // experimental arm
  for (i in 1:N) {
    if (v[i] == 1) {
      target += weibull_lpdf(y[i] | shape, exp(t[i]) * scale[t[i]+1]);
    } else {
      target += weibull_lccdf(y[i] | shape, exp(t[i]) * scale[t[i]+1]);
    }
  }
}

generated quantities {
  real HR;
  HR = scale[1] / scale[2]; // hazard ratio, control / experimental
}


