//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> num_obs;
  int<lower=0> num_vars;
  int<lower=0> num_loc;
  array[num_obs, num_vars] int X;
  // array[num_vars] num_levels;
  int max_num_level;
  array[num_obs] int loc_info_v;
  array[num_obs] int infected_v;
  array[num_obs] int total_test_v;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real intercept;
  matrix[max_num_level, num_vars] params;
  vector[num_loc] loc_param;
  real<lower=0> tau;
  vector<lower=0>[num_loc] sigma_v;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[num_obs] mu = rep_vector(0, num_obs);
  
  intercept ~ normal(0, 5);
  to_vector(params) ~ normal(0, 5);
  
  tau ~ cauchy(0, 1);
  sigma_v ~ cauchy(0, 1);
  loc_param ~ normal(0, sigma_v * tau);

  for (i in 1:num_vars) {
    mu += params[X[:, i],i];
  }
  mu += loc_param[loc_info_v];
  
  infected_v ~ binomial_logit(total_test_v, mu);
}

