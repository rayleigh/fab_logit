library(tidyverse)
library(brms)
library(boot)
library(parallel)

prior_binomial_y <- function(vals, y, n, prior_mean, prior_sd) {
  return(dbinom(y, n, inv.logit(vals)) * dnorm(vals, prior_mean, prior_sd))
}

calc_coverage_risk <- function(opt_v, theta, alpha, n, prior_mean, prior_sd) {

  s_j <- opt_v[1]
  interval_info <- derive_intervals_theta(theta, alpha, s_j, n)
  prop_y <- 0:n / n
  interested_y <- which(prop_y >= interval_info[1] & prop_y <= interval_info[2]) - 1
  sum(sapply(interested_y, function(y)
    integrate(prior_binomial_y, -Inf, Inf, y = y, n = n,
              prior_mean = prior_mean, prior_sd = prior_sd)$val))
}

optim_deriv_intervals <- function(s_j, mean_prop, alpha, n) {
  ci <- derive_intervals_theta(mean_prop, alpha, s_j, n)
  return(ci[2] - ci[1])
}

get_coverage_risk_boundary <- function(
  start_s_j, end_s_j, theta, alpha, n, prior_mean, prior_sd, min_coverage_risk) {

  while(abs(start_s_j - end_s_j) > 1e-6) {
    try_s_j <- (start_s_j + end_s_j) / 2
    tmp_risk <-
      calc_coverage_risk(try_s_j, theta = theta, alpha = alpha,
                         n = n, prior_mean = prior_mean,
                         prior_sd = prior_sd)
    if (tmp_risk > min_coverage_risk) {
      end_s_j = try_s_j
    } else {
      start_s_j = try_s_j
    }
  }
  return(start_s_j)
}

get_coverage_risk_region_s_j <- function(theta, alpha, n, prior_mean, prior_sd) {
  learn_region <- seq(0, 1, by = 0.01)
  coverage_risk_region <-
    sapply(learn_region, function(s_j) {
      calc_coverage_risk(s_j, theta = theta, alpha = alpha,
                         n = n, prior_mean = prior_mean,
                         prior_sd = prior_sd)
    })
  min_risk = min(coverage_risk_region)
  coverage_risk_region_ind <- which(coverage_risk_region == min_risk)
  min_ind = min(coverage_risk_region_ind)
  max_ind = max(coverage_risk_region_ind)
  lower_bnd = 1e-9
  upper_bnd = 1 - 1e-9
  if (min_ind != 1) {
    lower_bnd = get_coverage_risk_boundary(
      learn_region[min_ind], learn_region[min_ind - 1],
      theta, alpha, n, prior_mean, prior_sd, min_risk)
  }
  if (max_ind != length(learn_region)) {
    upper_bnd  = get_coverage_risk_boundary(
      learn_region[max_ind], learn_region[max_ind + 1],
      theta, alpha, n, prior_mean, prior_sd, min_risk)
  }
  return(c(lower_bnd, upper_bnd))
}

get_opt_s_j <- function(theta, alpha, n, prior_mean, prior_sd) {

  opt_risk_region <-
    get_coverage_risk_region_s_j(theta, alpha, n, prior_mean, prior_sd)
  optim(mean(opt_risk_region),
        optim_deriv_intervals, mean_prop = theta, alpha = alpha, n = n,
        method = "L-BFGS-B", lower = opt_risk_region[1],
        upper = opt_risk_region[2])$par
}

get_s_j_function <- function(alpha, n, prior_mean, prior_sd, spline = T) {

  theta_list <- seq(0, 1, by = 0.01)
  s_j_info <- sapply(theta_list, function(theta) {
    get_opt_s_j(theta, alpha, n, prior_mean, prior_sd)
  })
  if (spline) {
    return(smooth.spline(theta_list, s_j_info, cv = T))
  }
  s_j_info
}

refine_end_point <- function(
  init_end_point, y_prop, alpha, n, s_j_spline, interested_ind) {

  prev_end_point <- init_end_point
  tmp <- predict(s_j_spline, prev_end_point)$y
  tmp = min(tmp, 1 - 1e-9)
  tmp = max(tmp, 1e-9)
  next_end_point <-
    derive_intervals(y_prop, alpha, tmp, n)[interested_ind]
  iter = 0
  while (abs(next_end_point - prev_end_point) < 1e-3 && (iter < 1000)) {
    prev_end_point <- next_end_point
    tmp <- predict(s_j_spline, prev_end_point)$y
    tmp = min(tmp, 1 - 1e-9)
    tmp = max(tmp, 1e-9)
    next_end_point <-
      derive_intervals(y_prop, alpha, tmp, n)[interested_ind]
    iter = iter + 1
  }
  return(next_end_point)
}

refine_end_point_theta <- function(
  int_point, ext_point, y_prop, alpha, n, s_j_spline, interested_ind) {

  while(abs(int_point - ext_point) > 1e-6) {
    try_point = (int_point + ext_point) / 2
    tmp <- predict(s_j_spline, try_point)$y
    tmp <- min(tmp, 1- 1e-9)
    tmp <- max(tmp, 1e-9)
    ci <- derive_intervals_theta(try_point, alpha, tmp, n)
    if (ci[1] <= y_prop & y_prop <= ci[2]) {
      int_point = try_point
    } else {
      ext_point = try_point
    }
  }
  return(int_point)
}

refine_end_point_no_spline_theta <- function(
  int_point, ext_point, y_prop, alpha, n, interested_ind, prior_mean, prior_sd) {

  while(abs(int_point - ext_point) > 1e-6) {
    try_point = (int_point + ext_point) / 2
    opt_s_j <- get_opt_s_j(try_point, alpha, n, prior_mean, prior_sd)
    ci <- derive_intervals_theta(try_point, alpha, opt_s_j, n)
    if (ci[1] <= y_prop & y_prop <= ci[2]) {
      int_point = try_point
    } else {
      ext_point = try_point
    }
  }
  return(int_point)
}

derive_lower_intervals <- function(mean_prop, alpha, s_j, n) {
  s_j = min(s_j, 1 - 1e-9)
  s_j = max(s_j, 1e-9)
  prev_theta = 0.5
  next_theta =
    mean_prop - qnorm(alpha * ( 1 - s_j)) * sqrt(1 / n * prev_theta * (1 - prev_theta))
  iter = 0
  if (next_theta > 1) {
    next_theta = 1
  }
  while (abs(next_theta - prev_theta) < 1e-3 && (iter < 1000)) {
    prev_theta = next_theta
    next_theta =
      mean_prop - qnorm(alpha * ( 1 - s_j)) * sqrt(1 / n * prev_theta * (1 - prev_theta))
    if (next_theta > 1) {
      next_theta = 1
    }
    iter = iter + 1
  }
  return(next_theta)
}

derive_upper_intervals <- function(mean_prop, alpha, s_j, n) {
  s_j = min(s_j, 1 - 1e-9)
  s_j = max(s_j, 1e-9)
  prev_theta = 0.5
  next_theta =
    mean_prop - qnorm(1 - alpha * s_j) * sqrt(1 / n * prev_theta * (1 - prev_theta))
  iter = 0
  if (next_theta < 0) {
    next_theta = 0
  }
  while (abs(next_theta - prev_theta) < 1e-3 && (iter < 1000)) {
    prev_theta = next_theta
    next_theta =
      mean_prop - qnorm(1 - alpha * s_j) * sqrt(1 / n * prev_theta * (1 - prev_theta))
    if (next_theta < 0) {
      next_theta = 0
    }
    iter = iter + 1
  }
  return(next_theta)
}

derive_intervals <- function(mean_prop, alpha, s_j, n) {
  return(c(min(derive_lower_intervals(mean_prop, alpha, s_j, n), 1),
           max(derive_upper_intervals(mean_prop, alpha, s_j, n), 0)))
}

derive_intervals_theta <- function(theta, alpha, s_j, n) {

  s_j = min(s_j, 1 - 1e-9)
  s_j = max(s_j, 1e-9)

  sd <- sqrt(1 / n * theta * (1 - theta))
  ci <- c(theta + qnorm(alpha * ( 1 - s_j)) * sqrt(1 / n * theta * (1 - theta)),
            theta + qnorm(1 - alpha * s_j) * sqrt(1 / n * theta * (1 - theta)))
  ci <- pmax(ci, 0)
  ci <- pmin(ci, 1)

  return(ci)

}

#From https://stackoverflow.com/questions/71594430/how-to-find-where-the-interval-of-continuous-numbers-starts-and-ends
find_longest_cont_interval <- function(interested_inds) {
  
  group_id <- cumsum(c(1, diff(interested_inds)) != 1)
  split_inds <- split(interested_inds, group_id)
  return(split_inds[[which.max(sapply(split_inds, length))]])
}

derive_fab_intervals <- function(y_prop, alpha, n, prior_mean, prior_sd) {

  theta_list <- seq(0, 1, by = 0.01)
  s_j_spline <- get_s_j_function(alpha, n, prior_mean, prior_sd)
  s_j_vals <- predict(s_j_spline, theta_list)$y
  s_j_vals[s_j_vals > 1 - 1e-9] = 1 - 1e-9
  s_j_vals[s_j_vals < 1e-9] = 1e-9
  interval_vals <- sapply(1:length(s_j_vals), function(i) {
    derive_intervals_theta(theta_list[i], alpha, s_j_vals[i], n)
  })
  interested_inds <-
    which(interval_vals[1,] <= y_prop &
            y_prop <= interval_vals[2,])
  interested_inds <-
    find_longest_cont_interval(interested_inds)
  ci <- c(0, 1)
  ci <- c(0, 1)
  ind_info <- c(min(interested_inds), max(interested_inds))
  if (ind_info[1] > 1) {
    ci[1] <-
      refine_end_point_theta(
        theta_list[ind_info[1]], theta_list[ind_info[1] - 1],
        y_prop, alpha, n, s_j_spline, 1)
  }
  if (ind_info[2] < length(seq(0, 1, by = 0.01))) {
    ci[2] <-
      refine_end_point_theta(
        theta_list[ind_info[2]], theta_list[ind_info[2] + 1],
        y_prop, alpha, n, s_j_spline, 2)
  }
  return(ci)
}

derive_fab_intervals_no_spline <- function(y_prop, alpha, n, prior_mean, prior_sd) {

  theta_list <- seq(0, 1, by = 0.01)
  s_j_vals <- get_s_j_function(alpha, n, prior_mean, prior_sd, F)
  interval_vals <- sapply(1:length(s_j_vals), function(i) {
    derive_intervals_theta(theta_list[i], alpha, s_j_vals[i], n)
  })
  interested_inds <-
    which(interval_vals[1,] <= y_prop &
            y_prop <= interval_vals[2,])
  interested_inds <-
    find_longest_cont_interval(interested_inds)
  ci <- c(0, 1)
  ind_info <- c(min(interested_inds), max(interested_inds))
  if (ind_info[1] > 1) {
    ci[1] <-
      refine_end_point_no_spline_theta(
        theta_list[ind_info[1]], theta_list[ind_info[1] - 1],
        y_prop, alpha, n, 2, prior_mean, prior_sd)
  }
  if (ind_info[2] < length(theta_list)) {
    ci[2] <-
      refine_end_point_no_spline_theta(
        theta_list[ind_info[2]], theta_list[ind_info[2] + 1],
        y_prop, alpha, n, 2, prior_mean, prior_sd)
  }
  return(ci)
}

get_mrp_N <- function(gender, race, age, zip_code, acs_data, crosswalk_data) {

  interested_inds <- which(crosswalk_data$zip == zip_code)
  geoid_list <- as.numeric(crosswalk_data$geoid[interested_inds])
  interested_col <- paste(gender, race, gsub("-", ".", age), sep = "_")
  interested_data <-
    acs_data[acs_data$GEOID %in% geoid_list,] %>%
    arrange(GEOID) %>% select(interested_col)

  ratio_list <- crosswalk_data$tot_ratio[interested_inds]
  ratio_list <- ratio_list[order(geoid_list)]
  ratio_list <- ratio_list[which(sort(geoid_list) %in% acs_data$GEOID)]

  return(sum(ratio_list * acs_data[acs_data$GEOID %in% geoid_list, interested_col]) / sum(ratio_list))

}

brms_get_zip_code_draws <- function(new_stan_results, sim_reduced_covid) {

  zip_code_v <- sort(unique(sim_reduced_covid$zip))
  post_prob <- inv.logit(posterior_linpred(new_stan_results))
  zip_draws <- matrix(0, nrow = nrow(post_prob), ncol = length(zip_code_v))
  for (i in 1:length(zip_code_v)) {

    interested_inds <- which(zip_code_v == zip_code_v[i])
    zip_draws[, i] <-
      rowSums(post_prob[, interested_inds, drop = F] *
        sim_reduced_covid_data$total[interested_inds]) /
          sum(sim_reduced_covid_data$total[interested_inds])
  }
  zip_draws[zip_draws < 1e-9] = 1e-9
  zip_draws[zip_draws > 1 - 1e-9] = 1 - 1e-9
  return(zip_draws)
}

post_draws_get_zip_code_draws <- function(draws_df, sim_reduced_covid) {

  zip_code_v <- sort(unique(sim_reduced_covid$zip))
  zip_draws <- matrix(0, nrow = nrow(draws_df), ncol = length(zip_code_v))
  for (i in 1:length(zip_code_v)) {

    interested_inds <- which(zip_code_v == zip_code_v[i])
    zip_draws[, i] <-
      rowSums(inv.logit(draws_df[, interested_inds, drop = F]) *
                sim_reduced_covid_data$total[interested_inds]) /
      sum(sim_reduced_covid_data$total[interested_inds])
  }
  return(logit(zip_draws))
}

aggregate_sim_data <- function(sim_reduced_covid) {

   tmp <- sim_reduced_covid %>% group_by(zip) %>%
     summarize(pos = sum(positive), tot = sum(total))
   tmp$positive = tmp$pos
   tmp$total = tmp$tot
   return(tmp)
}

rstanarm_results_get_coverage <- function(new_stan_results, sim_reduced_covid_data) {

  post_count_by_zip <- apply(posterior_predict(new_stan_results), 1, function(row) {
    tmp <- as.data.frame(cbind(sim_reduced_covid_data$zip, row))
    colnames(tmp) <- c("zip", "count")
    tmp %>% group_by(zip) %>% summarize(case_count = sum(count))
  })
  post_count_by_zip <- sapply(post_count_by_zip, function(count_info) {
    count_info$case_count
  })
  post_count_credible_interval <-
    apply(post_count_by_zip, 1, function(row) quantile(row, probs = c(0.025, 0.975)))
  sim_zip_count <- sim_reduced_covid_data %>%
    group_by(zip) %>% summarize(case_count = sum(positive))
  return(
    (sim_zip_count$case_count > post_count_credible_interval[1,] &
      sim_zip_count$case_count < post_count_credible_interval[2,]) * 1)
}

post_draws_get_coverage <- function(post_draws, sim_reduced_covid_data) {

  post_count_by_zip <- apply(post_draws, 1, function(row) {
    tmp <- as.data.frame(cbind(sim_reduced_covid_data$zip, row))
    colnames(tmp) <- c("zip", "count")
    tmp %>% group_by(zip) %>% summarize(case_count = sum(count))
  })
  post_count_by_zip <- sapply(post_count_by_zip, function(count_info) {
    count_info$case_count
  })
  post_count_credible_interval <-
    apply(post_count_by_zip, 1, function(row) quantile(row, probs = c(0.025, 0.975)))
  sim_zip_count <- sim_reduced_covid_data %>%
    group_by(zip) %>% summarize(case_count = sum(positive))
  return(
    (sim_zip_count$case_count > post_count_credible_interval[1,] &
       sim_zip_count$case_count < post_count_credible_interval[2,]) * 1)
}

rstanarm_results_get_all_coverage <- function(new_stan_results, sim_reduced_covid_data) {

  post_count_credible_interval <-
    apply(posterior_predict(new_stan_results), 2, function(row) quantile(row, probs = c(0.025, 0.975)))
  sim_zip_count <- sim_reduced_covid_data %>%
    group_by(zip) %>% summarize(case_count = sum(positive))
  return(
    (sim_reduced_covid_data$positive > post_count_credible_interval[1,] &
       sim_reduced_covid_data$positive < post_count_credible_interval[2,]) * 1)
}

post_draws_get_all_coverage <- function(post_draws, sim_reduced_covid_data) {

  post_count_credible_interval <-
    apply(post_draws, 2, function(row) quantile(row, probs = c(0.025, 0.975)))
  return(
    (sim_reduced_covid_data$positive > post_count_credible_interval[1,] &
       sim_reduced_covid_data$positive < post_count_credible_interval[2,]) * 1)
}

brms_get_prop_coverage <- function(new_stan_results, gen_prop_val) {

  ci <- apply(inv.logit(posterior_linpred(new_stan_results)), 2, function(col) {
    quantile(col, probs = c(0.025, 0.975))
  })
  return(
    list((gen_prop_val > ci[1,] & gen_prop_val < ci[2,]) * 1,
	 ci[2,] - ci[1,]))
}

brms_get_prop_coverage_adjusted <- function(new_stan_results, gen_prop_val, sim_reduced_covid) {

  post_prob <- inv.logit(posterior_linpred(new_stan_results))
  post_prob[post_prob < 1e-9] = 1e-9
  post_prob[post_prob > 1 - 1e-9] = 1 - 1e-9
  post_mean = colMeans(logit(post_prob))
  post_sd <- apply(logit(post_prob), 2, sd)
  obs_prop <- sim_reduced_covid$positive / sim_reduced_covid$total
  ci <- sapply(1:length(post_mean), function(i) {
      print(i)
      if (sim_reduced_covid$total[i] < 3) {
         obs_prop[i] = inv.logit(post_mean[i])
      }
      derive_fab_intervals_no_spline(obs_prop[i], 0.05, sim_reduced_covid$total[i], post_mean[i], post_sd[i])
  })
  
  return(list((gen_prop_val > ci[1,] & gen_prop_val < ci[2,]) * 1,
	      ci[2,] - ci[1,]))

}

post_draws_get_prop_coverage <- function(post_draws, gen_prop_val) {

  ci <- apply(post_draws, 2, function(col) {
    quantile(col, probs = c(0.025, 0.975))
  })
  return(
    list((gen_prop_val > ci[1,] & gen_prop_val < ci[2,]) * 1,
	 ci[2,] - ci[1,]))
}

post_draws_get_prop_coverage_adjusted <- function(post_draws, gen_prop_val, sim_reduced_covid) {

  post_mean = colMeans(post_draws)
  post_sd <- apply(post_draws, 2, sd)
  obs_prop <- sim_reduced_covid$positive / sim_reduced_covid$total
  ci <- sapply(1:length(post_mean), function(i) {
    print(i)
    if (sim_reduced_covid$total[i] < 3) {
       obs_prop[i] = inv.logit(post_mean[i]) 
    }
    derive_fab_intervals_no_spline(obs_prop[i], 0.05, sim_reduced_covid$total[i], post_mean[i], post_sd[i])
  })
  return(
    list((gen_prop_val > ci[1,] & gen_prop_val < ci[2,]) * 1,
	 ci[2,] - ci[1,]))
}

get_direct_coverage_for_sim <- function(mean_prop, n, obs_prop) {

  ci <- sapply(1:length(mean_prop), function(i) {
    derive_intervals(mean_prop[i], 0.05, 0.5, n[i])
  })
  return(list((ci[2,] < obs_prop & obs_prop < ci[1,]) * 1,
         ci[1,] - ci[2,]))

}

get_direct_coverage_m <- function(reduced_covid_data, sim_prop_m, gen_prop_v) {

  coverage_ci_m <- matrix(0, nrow = ncol(sim_prop_m), ncol = nrow(sim_prop_m))
  coverage_ci_length_m <- matrix(0, nrow = ncol(sim_prop_m), ncol = nrow(sim_prop_m))

  for (i in 1:ncol(sim_prop_m)) {
    ci_info <-
      get_direct_coverage_for_sim(sim_prop_m[,i] / reduced_covid_data$total,
                                  reduced_covid_data$total, gen_prop_v)
    coverage_ci_m[i,] <- ci_info[[1]]
    coverage_ci_length_m[i,] <- ci_info[[2]]
  }

  return(list(coverage_ci_m, coverage_ci_length_m))
}

get_zip_info <- function(reduced_covid_data, zip_code_db, sim_data_100, all_draws_post_prop) {
  county_zip_code_info <-
    sapply(reduced_covid_data$zip, function(zip_code) {
      zip_code_db[zip_code_db$zipcode == zip_code,]$county
    })
  positive_prop_info <- data.frame("prop" = all_draws_post_prop[1,],
                                   "total" = reduced_covid_data$total,
                                   "total_mrp" = reduced_covid_data$mrpN,
                                   "zip" = reduced_covid_data$zip)
  mrpN_county_positive_prop_info <-
    positive_prop_info %>% group_by(zip) %>%
    summarize(pos_prop = sum(prop * total_mrp) / sum(total_mrp))
  
  county_reduced_covid_data <- reduced_covid_data
  tmp <- aggregate_sim_data(county_reduced_covid_data)
  county_prop_m <- matrix(0, nrow = nrow(tmp), ncol = ncol(sim_data_100))
  for (i in 1:100) {
    county_reduced_covid_data$positive <- sim_data_100[,i]
    tmp <- aggregate_sim_data(county_reduced_covid_data)
    county_prop_m[,i] <- tmp$positive
    
  }
  county_reduced_covid_data <- aggregate_sim_data(county_reduced_covid_data)
  
  return(list(mrpN_county_positive_prop_info, county_prop_m, county_reduced_covid_data))
}

