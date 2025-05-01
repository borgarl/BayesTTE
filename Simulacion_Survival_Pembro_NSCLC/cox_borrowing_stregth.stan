data {
  int<lower=0> N; // Total number of observations
  int<lower=0, upper=1> status[N]; // Event indicator (1=event, 0=censored)
  vector[N] time; // Survival times
  int<lower=0, upper=1> dataset[N]; // Indicator for dataset (0=current, 1=external)
}

parameters {
  real beta; // Coefficient for the dataset
  real<lower=0> baseline_hazard; // Baseline hazard rate
}

model {
  // Priors
  beta ~ normal(0, 10); // Weakly informative prior for the dataset coefficient
  baseline_hazard ~ normal(0, 10); // Weakly informative prior for the baseline hazard

  // Likelihood
  for (i in 1:N) {
    // Cox proportional hazards model
    real lambda = baseline_hazard * exp(beta * dataset[i]); // Hazard function
    real S = exp(-lambda * time[i]); // Survival function

    // If the event is observed
    if (status[i] == 1) {
      target += log(lambda) - lambda * time[i]; // Log-likelihood for an observed event
    } else {
      target += log(S); // Log-likelihood for a censored observation
    }
  }
}

generated quantities {
  real HR = exp(-beta); // Hazard ratio for the external dataset vs. current
}
