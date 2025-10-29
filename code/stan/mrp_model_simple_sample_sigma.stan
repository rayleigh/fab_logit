
data { 
  int<lower=1> N;
  array[N] int y;
  array[N] int n_sample;
  int<lower=0> K;
  matrix[N, K] X;
  int<lower=1> N_race;
  array[N] int<lower=1, upper=N_race> J_race;
  int<lower=1> N_zip;
  array[N] int<lower=1, upper=N_zip> J_zip;
  int<lower=1> N_age;
  array[N] int<lower=1, upper=N_age> J_age;
  real<lower=0> sens;
  real<lower=0> spec;
  matrix[N_zip, N_zip] zip_code_distance;
  real<lower=0> l;
  //real<lower=0> sigma;
}

transformed data {
  matrix[N_zip, N_zip] K_m =
    exp(-square(zip_code_distance) / (2 * l^2));
  matrix[N_zip, N_zip] chol_K;

  for (n in 1:N_zip) {
    K_m[n, n] = K_m[n, n] + 1e-3;
  }
  chol_K = cholesky_decompose(K_m);
}

parameters { 
  real Intercept;
  vector[K] beta;
  real<lower=0> lambda_race;
  vector[N_race] z_race;
  real<lower=0> lambda_zip;
  vector[N_zip] z_zip;
  real<lower=0> lambda_age;
  vector[N_age] z_age;
  real<lower=0> sigma;
}

transformed parameters { 
  real<lower=0> scaled_lambda_race = lambda_race;
  vector[N_race] a_race = z_race * scaled_lambda_race;
  real<lower=0> scaled_lambda_zip = lambda_zip;
  // vector[N_zip] a_zip = z_zip * scaled_lambda_zip;
  vector[N_zip] a_zip = sigma * chol_K * z_zip;
  real<lower=0> scaled_lambda_age = lambda_age;
  vector[N_age] a_age = z_age * scaled_lambda_age;
  vector<lower=0, upper=1>[N] p = inv_logit(Intercept + X * beta + a_race[J_race] + a_zip[J_zip] + a_age[J_age]);
  vector<lower=0, upper=1>[N] p_sample = p * sens + (1 - p) * (1 - spec);
}

model { 
  y ~ binomial(n_sample, p_sample);
  Intercept ~ normal(0, 5);
  beta[1] ~ normal(0, 3);
  beta[2] ~ normal(0, 3);
  beta[3] ~ normal(0, 3);
  beta[4] ~ normal(0, 3);
  beta[5] ~ normal(0, 3);
  beta[6] ~ normal(0, 3);
  beta[7] ~ normal(0, 3);
  z_race ~ std_normal();
  z_zip ~ std_normal();
  z_age ~ std_normal();
  lambda_race ~ normal(0, 3);
  lambda_zip ~ normal(0, 3);
  lambda_age ~ normal(0, 3);
  sigma ~ cauchy(0, 1);
}

generated quantities { 
  array[N] int<lower = 0> y_rep = binomial_rng(n_sample, p_sample);
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_lpmf(y[n] | n_sample[n], p_sample[n]);
  }
}
  
