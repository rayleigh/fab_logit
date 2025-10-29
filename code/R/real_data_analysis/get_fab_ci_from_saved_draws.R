source("code/R/gen_functions/small_area_estimate_helper_functions.R")

library(tidycensus)
acs_data <- get_acs(
      geography = "tract", 
      state = c("IN", "IL"),
      variable = "B02001_001",
      year = 2020)

covid_by_county <- merge(spike_in_covid_data, in_county_info, by.x = "LOCATION_ID", by.y = "fips_state_county_cd")

zip_code_draws <- get_poststratified_draws_table("result_files/peak_indiana_covid_results_simple_sample_sigma_7.5.RDS", preprocess_data, poststrat_info)
zip_code_draws_pop <- get_pop_poststratified_draws("result_files/peak_indiana_covid_results_sample_sigma.RDS", new_data, new_data)

zip_observed_input_data <- get_observed_county_poststratified_draws_table(preprocess_data, F)
zip_code_input_data <- aggregate_sim_data(preprocess_data)

merge_info <- merge(crosswalk_info, acs_data, by.x = "geoid", by.y = "GEOID")

# input_data <- preprocess_data
zip_codes <- unique(preprocess_data$zip)
#Ignore 60421
zip_codes <- zip_codes[zip_codes != 60421]
all_counties <- unique(as.numeric(merge_info$geoid[merge_info$zip %in% zip_codes]) %/% 1e6)
mrp_county_draws <- matrix(0, nrow = nrow(zip_code_draws), ncol = length(all_counties))
# mrp_county_draws_pop <- matrix(0, nrow = nrow(zip_code_draws), ncol = length(all_counties))
mrp_input_data <- matrix(0, nrow = length(all_counties), ncol = 3)

zip_adj_ci_intervals <- matrix(0, nrow = 2, ncol = length(all_counties))
zip_adj_norm_prop_ci_intervals <- matrix(0, nrow = 2, ncol = length(all_counties))
zip_adj_ac_ci_intervals <- matrix(0, nrow = 2, ncol = length(all_counties))

zip_obs_adj_ci_intervals <- matrix(0, nrow = 2, ncol = length(all_counties))
zip_obs_adj_norm_prop_ci_intervals <- matrix(0, nrow = 2, ncol = length(all_counties))
zip_obs_adj_ac_ci_intervals <- matrix(0, nrow = 2, ncol = length(all_counties))

for (i in 1:length(all_counties)) {
  county = all_counties[i]
  # interested_zip_ind <- which(merge_info$fips_state_county_cd == county)
  interested_zip_ind <- which(as.numeric(merge_info$geoid) %/% 1e6 == county)
  # weights <- merge_info$tot_ratio[interested_zip_ind] * merge_info$pop_size[interested_zip_ind]
  weights <- merge_info$tot_ratio[interested_zip_ind] * merge_info$estimate[interested_zip_ind]
  ratio_weights <- merge_info$tot_ratio[interested_zip_ind]
  ratio_info <- data.frame("zip" = merge_info$zip[interested_zip_ind], 
                           "weights" = merge_info$tot_ratio[interested_zip_ind] * merge_info$estimate[interested_zip_ind],
                           "ratio_weights" = merge_info$tot_ratio[interested_zip_ind])
  ratio_info <- ratio_info %>% group_by(zip) %>% 
    summarize("weights" = sum(weights),
              "ratio_weights" = sum(ratio_weights))
  match_zip <- rep(NA, nrow(ratio_info))
  for (zip in ratio_info$zip) {
    print(zip)
    # if(any(input_data$zip_raw == zip)) {
    if(any(preprocess_data$zip == zip)) {
      match_zip[which(ratio_info$zip == zip)] <-
        # unique(preprocess_data$zip[input_data$zip_raw == zip])
        unique(preprocess_data$zip[preprocess_data$zip == zip])
    }
  }
  match_zip[ratio_info$zip == 60421] = NA
  weights <- ratio_info$weights
  ratio_weights <- ratio_info$ratio_weights
  match_zip <- sapply(match_zip, function(zip) {
    if (is.na(zip)) {
      return(NA)
    }
    return(unique(preprocess_data$zip_level[preprocess_data$zip == zip]))
  })
  print(mean(is.na(match_zip)))
  mrp_county_draws[,i] <- 
    rowSums(sweep(zip_code_draws[, match_zip[!is.na(match_zip)], drop = F], 2, 
         weights[!is.na(match_zip)], "*")) / sum(weights[!is.na(match_zip)])
  # mrp_county_draws_pop[,i]  <- 
  #   rowSums(sweep(inv.logit(zip_code_draws_pop[, match_zip[!is.na(match_zip)], drop = F]), 2, 
  #                 weights[!is.na(match_zip)], "*")) / sum(weights[!is.na(match_zip)])
  mrp_input_data[i,] <-
    c(i, sum(zip_observed_input_data$pos_prop[match_zip[!is.na(match_zip)]] * ratio_weights[!is.na(match_zip)]) / 
        sum(ratio_weights[!is.na(match_zip)]),
      sum(zip_observed_input_data$total[match_zip[!is.na(match_zip)]] * ratio_weights[!is.na(match_zip)]))
  
  zip_adj_ci_intervals[,i] <-
    rowSums(sweep(zip_code_post_wilson_intervals[, match_zip[!is.na(match_zip)], drop = F], 2,
                  ratio_weights[!is.na(match_zip)], "*")) /
    sum(ratio_weights[!is.na(match_zip)])
  zip_adj_norm_prop_ci_intervals[,i] <-
    rowSums(sweep(zip_code_post_wald_intervals[, match_zip[!is.na(match_zip)], drop = F], 2,
                  ratio_weights[!is.na(match_zip)], "*")) /
    sum(ratio_weights[!is.na(match_zip)])
  zip_adj_ac_ci_intervals[,i] <-
    rowSums(sweep(zip_code_post_ac_intervals[, match_zip[!is.na(match_zip)], drop = F], 2,
                  ratio_weights[!is.na(match_zip)], "*")) /
    sum(ratio_weights[!is.na(match_zip)])
  
  zip_obs_adj_ci_intervals[,i] <-
    rowSums(sweep(zip_code_obs_wilson_intervals[, match_zip[!is.na(match_zip)], drop = F], 2,
                  ratio_weights[!is.na(match_zip)], "*")) /
    sum(ratio_weights[!is.na(match_zip)])
  zip_obs_adj_norm_prop_ci_intervals[,i] <-
    rowSums(sweep(zip_code_obs_wald_intervals[, match_zip[!is.na(match_zip)], drop = F], 2,
                  ratio_weights[!is.na(match_zip)], "*")) /
    sum(ratio_weights[!is.na(match_zip)])
  zip_obs_adj_ac_ci_intervals[,i] <-
    rowSums(sweep(zip_code_obs_ac_intervals[, match_zip[!is.na(match_zip)], drop = F], 2,
                  ratio_weights[!is.na(match_zip)], "*")) /
    sum(ratio_weights[!is.na(match_zip)])
  
  # mrp_input_data[i,2:3] <- mrp_input_data[i,2:3] / sum(ratio_weights[!is.na(match_zip)])
  # print(sum(ratio_weights[!is.na(match_zip)]))
}
mrp_input_data <- as.data.frame(mrp_input_data)
mrp_input_data$county <- all_counties[mrp_input_data$V1]
mrp_county_draws <- logit(mrp_county_draws)
mrp_county_draws_pop <- logit(mrp_county_draws_pop)

# obs_point <- mrp_input_data$V2 / mrp_input_data$V3
obs_point <- mrp_input_data$V2
post_estimates <- colMeans(inv.logit(mrp_county_draws))

interested_inds <- which(round(mrp_input_data$V3) >= 20 & round(mrp_input_data$V3) <= 1500)
interested_inds <- interested_inds[interested_inds != 67]

post_means <- colMeans(mrp_county_draws)
post_sd <- apply(mrp_county_draws, 2, sd)
zip_adj_ci_intervals <- 
  do.call(cbind, lapply(interested_inds, function(i) {
    # sapply(1:ncol(zip_code_draws), function(i) {
    print(i)
    # derive_fab_intervals_multiple_y(c(inv.logit(post_means[i]),obs_point[i]), 0.05, round(mrp_input_data$V3[i]), post_means[i], post_sd[i])
    derive_fab_intervals_multiple_y(c(post_estimates[i],obs_point[i]), 0.05, round(mrp_input_data$V3[i]), post_means[i], post_sd[i], spline = F, all_in = T)
  }))

zip_adj_norm_prop_ci_intervals <- 
  do.call(cbind, lapply(interested_inds, function(i) {
    # sapply(1:ncol(zip_code_draws), function(i) {
    print(i)
    derive_fab_intervals_normal_multiple_y(c(post_estimates[i],obs_point[i]), 0.05, round(mrp_input_data$V3[i]), post_means[i], post_sd[i], spline = F)
    # derive_fab_intervals_normal_multiple_y(c(inv.logit(post_means[i]),obs_point[i]), 0.05, round(mrp_input_data$V3[i]), post_means[i], post_sd[i])
  }))

zip_adj_ac_ci_intervals <- 
  do.call(cbind, lapply(interested_inds, function(i) {
    # sapply(1:ncol(zip_code_draws), function(i) {
    print(i)
    mod_obs_val = (obs_point[i] * mrp_input_data$V3[i] + 2) / 
      (round(mrp_input_data$V3[i]) + 4)
    mod_post_val = (post_estimates[i] * mrp_input_data$V3[i] + 2) / 
      (round(mrp_input_data$V3[i]) + 4)
    derive_fab_intervals_ac_multiple_y(c(post_estimates[i],obs_point[i]), 0.05, round(mrp_input_data$V3[i]), post_means[i], post_sd[i], spline = F, all_in = T)
  }))
