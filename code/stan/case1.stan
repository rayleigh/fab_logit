
data { 
  int<lower=1> N;
  array[N] int y;
  array[N] int n_sample;
  int<lower=0> K;
  matrix[N, K] X;
  int<lower=1> N_pop;
  int<lower=0> K_pop;
  matrix[N_pop, K_pop] X_pop;
  real<lower=0> sens;
  real<lower=0> spec;
}

parameters { 
  real Intercept;
  vector[K] beta;
}

transformed parameters { 
  vector<lower=0, upper=1>[N] p = inv_logit(Intercept + X * beta);
  vector<lower=0, upper=1>[N] p_sample = p * sens + (1 - p) * (1 - spec);
}

model { 
  y ~ binomial(n_sample, p_sample);
  Intercept ~ normal(0, 5);
  beta[1] ~ normal(0, 3);
  beta[2] ~ normal(0, 3);
  beta[3] ~ normal(0, 3);
}

generated quantities { 
  vector<lower=0, upper=1>[N_pop] p_pop = inv_logit(Intercept + X_pop * beta);
  array[N] int<lower = 0> y_rep = binomial_rng(n_sample, p_sample);
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_lpmf(y[n] | n_sample[n], p_sample[n]);
  }
}
  
