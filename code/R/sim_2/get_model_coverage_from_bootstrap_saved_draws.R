library(tidyverse)
library(brms)
library(cmdstanr)
library(boot)

args = commandArgs(trailingOnly=TRUE)
sim_data_file = args[1]
draw_file = args[2]
result_file = args[3]
iter_num = as.numeric(args[4])

load("data_files/stan_indiana_covid_bootstrap_sim_data.Rdata")
load("data_files/zip_code_db.Rdata")
source("code/R/gen_functions/small_area_estimate_helper_functions.R")

load(sim_data_file)
#sim_data_100 <- spatial_sim_prop_m
print(sim_data_file)
print(result_file)


county_zip_code_info <-
  sapply(reduced_covid_data$zip, function(zip_code) {
     zip_code_db[zip_code_db$zipcode == zip_code,]$county
  })
  positive_prop_info <- data.frame("prop" = all_draws_post_prop[1,],
				   "total" = reduced_covid_data$total,
				   "total_mrp" = reduced_covid_data$mrpN,
                                   "zip" = reduced_covid_data$zip)
  county_positive_prop_info <-
    positive_prop_info %>% group_by(zip) %>%
    summarize(pos_prop = sum(prop * total) / sum(total),
	      zip_total = sum(total))

mrpN_county_positive_prop_info <-
  positive_prop_info %>% group_by(zip) %>%
  summarize(pos_prop = sum(prop * total_mrp) / sum(total_mrp))

ref_data <- reduced_covid_data
reduced_covid_data <- reduced_covid_data_list[[iter_num]]

observed_cells <- sort(which(ref_data$combined %in% reduced_covid_data$combined))
all_draws_post_prop <- all_draws_post_prop[,observed_cells]

observed_zip <- sort(unique(as.numeric(reduced_covid_data$zip)))
county_positive_prop_info <- county_positive_prop_info[observed_zip,]
mrpN_county_positive_prop_info[observed_zip,] <-
  mrpN_county_positive_prop_info[observed_zip,]

obs_positive_prop_info <- data.frame(
  "prop" = reduced_covid_data$positive / reduced_covid_data$total,
  "total" = reduced_covid_data$total,
  "total_mrp" = reduced_covid_data$mrpN,
  "zip" = reduced_covid_data$zip)
mrpN_obs_county_positive_prop_info <-
  obs_positive_prop_info %>% group_by(zip) %>%
  #summarize(pos_prop = sum(prop * total_mrp) / sum(total_mrp),
  summarize(pos_prop = sum(prop * total) / sum(total),
	    total = sum(total)) %>%
  arrange(zip)  

#rm(new_stan_results)

#mod <- cmdstan_model("default_brms_stan_code_gp.stan")
coverage_m <- matrix(0, nrow = 100, ncol = 30)
prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
prop_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
prop_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
direct_post_prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
direct_post_prop_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_norm_prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_norm_prop_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_norm_prop_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_ac_prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_ac_prop_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_ac_prop_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_obs_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_obs_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_obs_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_post_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_post_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_post_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_norm_prop_obs_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_norm_prop_obs_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_norm_prop_obs_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_norm_prop_post_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_norm_prop_post_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_norm_prop_post_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_ac_prop_obs_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_ac_prop_obs_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_ac_prop_obs_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_ac_prop_post_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_ac_prop_post_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_ac_prop_post_bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
zip_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
zip_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
zip_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_zip_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_zip_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_zip_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
direct_post_zip_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
direct_post_zip_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
mrp_zip_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
mrp_zip_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
mrp_zip_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_obs_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_obs_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_obs_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_post_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_post_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_post_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_norm_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_norm_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_norm_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_norm_obs_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_norm_obs_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_norm_obs_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_norm_post_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_norm_post_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_norm_post_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_ac_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_ac_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_ac_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_ac_obs_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_ac_obs_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_ac_obs_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_ac_post_coverage_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_ac_post_coverage_length_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
adj_mrp_zip_ac_post_bias_m <- matrix(0, nrow = 100, ncol = length(observed_zip))
#prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
#adj_prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
#for (i in 1:100) {
  i = iter_num
  print(paste0("Iter: ", i))
  sim_reduced_covid_data <- reduced_covid_data
  #sim_reduced_covid_data$positive <- sim_data_100[,i]
  sim_reduced_covid_data$negative <-
    sim_reduced_covid_data$total - sim_reduced_covid_data$positive
  sim_reduced_covid_data$pos_prop <-
    sim_reduced_covid_data$positive / sim_reduced_covid_data$total

  print(sim_reduced_covid_data)
  print(dim(all_draws_post_prop))

  new_stan_results <- readRDS(draw_file)

  #new_stan_results <- mod$sample(
  #  data = data_list, adapt_delta = 0.95, max_treedepth = 12,
  #  parallel_chains = 4)
  #print(colMeans(inv.logit(draws_df)))

  # fit_summary_info <- new_stan_results$summary(
  #   variables = "post_draws",
  #   posterior::default_summary_measures(),
  #   extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
  # )

  draws_df <- new_stan_results$draws(format = "df")
  draws_df <- data.matrix(draws_df[, grep("post_draws", colnames(draws_df))])
  #coverage_m[i,] <-
  #  post_draws_get_coverage(draws_df, sim_reduced_covid_data)
  bias_m[i,] <- all_draws_post_prop[1,] - colMeans(draws_df)  



  draws_df <- new_stan_results$draws(format = "df")
  draws_df <- data.matrix(draws_df[, grep("mu", colnames(draws_df))])
  prop_coverage_info <-
    post_draws_get_prop_coverage(inv.logit(draws_df), all_draws_post_prop[1,])
  prop_coverage_m[i,] <- prop_coverage_info[[1]]
  prop_coverage_length_m[i,] <- prop_coverage_info[[2]]
  prop_bias_m[i,] <- colMeans(inv.logit(draws_df)) - all_draws_post_prop[1,]
  prop_coverage_info <-
    post_draws_get_direct_post_prop_coverage(draws_df, sim_reduced_covid_data$total, all_draws_post_prop[1,])
  direct_post_prop_coverage_m[i,] <- prop_coverage_info[[1]]
  direct_post_prop_coverage_length_m[i,] <- prop_coverage_info[[2]]
  #prop_coverage_m[i,] <-
  #  post_draws_get_prop_coverage(inv.logit(draws_df), all_draws_post_prop[1,]) 
  #adj_prop_coverage_m[i,] <-
  #  post_draws_get_prop_coverage_adjusted(draws_df, all_draws_post_prop[1,], sim_reduced_covid_data) 
  adj_prop_coverage_info <-
    post_draws_get_prop_coverage_adjusted_multiple_non_parallel(draws_df, all_draws_post_prop[1,], sim_reduced_covid_data)
  adj_prop_coverage_m[i,] <-
    adj_prop_coverage_info[[1]][1:ncol(all_draws_post_prop) * 3]
  adj_prop_coverage_length_m[i,] <-
    adj_prop_coverage_info[[2]][1:ncol(all_draws_post_prop) * 3]
  adj_prop_bias_m[i,] <-
    adj_prop_coverage_info[[3]][1:ncol(all_draws_post_prop) * 3]
  adj_prop_obs_coverage_m[i,] <-
    adj_prop_coverage_info[[1]][1:ncol(all_draws_post_prop) * 3 - 2]
  adj_prop_obs_coverage_length_m[i,] <-
    adj_prop_coverage_info[[2]][1:ncol(all_draws_post_prop) * 3 - 2]
  adj_prop_obs_bias_m[i,] <-
    adj_prop_coverage_info[[3]][1:ncol(all_draws_post_prop) * 3 - 2]
  adj_prop_post_coverage_m[i,] <-
     adj_prop_coverage_info[[1]][1:ncol(all_draws_post_prop) * 3 - 1]
  adj_prop_post_coverage_length_m[i,] <-
     adj_prop_coverage_info[[2]][1:ncol(all_draws_post_prop) * 3 - 1]
  adj_prop_post_bias_m[i,] <-
     adj_prop_coverage_info[[3]][1:ncol(all_draws_post_prop) * 3 - 1]

  adj_prop_coverage_info <-
    post_draws_get_norm_prop_coverage_adjusted_multiple_non_parallel(draws_df, all_draws_post_prop[1,], sim_reduced_covid_data)
  adj_norm_prop_coverage_m[i,] <-
    adj_prop_coverage_info[[1]][1:ncol(all_draws_post_prop) * 3]
  adj_norm_prop_coverage_length_m[i,] <-
    adj_prop_coverage_info[[2]][1:ncol(all_draws_post_prop) * 3]
  adj_norm_prop_bias_m[i,] <-
    adj_prop_coverage_info[[3]][1:ncol(all_draws_post_prop) * 3]
  adj_norm_prop_obs_coverage_m[i,] <-
    adj_prop_coverage_info[[1]][1:ncol(all_draws_post_prop) * 3 - 2]
  adj_norm_prop_obs_coverage_length_m[i,] <-
    adj_prop_coverage_info[[2]][1:ncol(all_draws_post_prop) * 3 - 2]
  adj_norm_prop_obs_bias_m[i,] <-
    adj_prop_coverage_info[[3]][1:ncol(all_draws_post_prop) * 3 - 2]
  adj_norm_prop_post_coverage_m[i,] <-
     adj_prop_coverage_info[[1]][1:ncol(all_draws_post_prop) * 3 - 1]
  adj_norm_prop_post_coverage_length_m[i,] <-
     adj_prop_coverage_info[[2]][1:ncol(all_draws_post_prop) * 3 - 1]
  adj_norm_prop_post_bias_m[i,] <-
     adj_prop_coverage_info[[3]][1:ncol(all_draws_post_prop) * 3 - 1]

  adj_prop_coverage_info <-
    post_draws_get_ac_prop_coverage_adjusted_multiple_non_parallel(draws_df, all_draws_post_prop[1,], sim_reduced_covid_data)
  adj_ac_prop_coverage_m[i,] <-
    adj_prop_coverage_info[[1]][1:ncol(all_draws_post_prop) * 3]
  adj_ac_prop_coverage_length_m[i,] <-
    adj_prop_coverage_info[[2]][1:ncol(all_draws_post_prop) * 3]
  adj_ac_prop_bias_m[i,] <-
    adj_prop_coverage_info[[3]][1:ncol(all_draws_post_prop) * 3]
  adj_ac_prop_obs_coverage_m[i,] <-
    adj_prop_coverage_info[[1]][1:ncol(all_draws_post_prop) * 3 - 2]
  adj_ac_prop_obs_coverage_length_m[i,] <-
    adj_prop_coverage_info[[2]][1:ncol(all_draws_post_prop) * 3 - 2]
  adj_ac_prop_obs_bias_m[i,] <-
    adj_prop_coverage_info[[3]][1:ncol(all_draws_post_prop) * 3 - 2]
  adj_ac_prop_post_coverage_m[i,] <-
     adj_prop_coverage_info[[1]][1:ncol(all_draws_post_prop) * 3 - 1]
  adj_ac_prop_post_coverage_length_m[i,] <-
     adj_prop_coverage_info[[2]][1:ncol(all_draws_post_prop) * 3 - 1]
  adj_ac_prop_post_bias_m[i,] <-
     adj_prop_coverage_info[[3]][1:ncol(all_draws_post_prop) * 3 - 1]

  print("Done with cell")

  zip_code_draws <- 
    post_draws_get_zip_code_draws(draws_df, sim_reduced_covid_data)
  prop_coverage_info <-
    post_draws_get_prop_coverage(
      inv.logit(zip_code_draws), county_positive_prop_info$pos_prop)
  zip_coverage_m[i,] <- prop_coverage_info[[1]]
  zip_coverage_length_m[i,] <- prop_coverage_info[[2]]
  #prop_coverage_m[i,] <-
  #  post_draws_get_prop_coverage(inv.logit(draws_df), all_draws_post_prop[1,])
  #adj_prop_coverage_m[i,] <-
  #  post_draws_get_prop_coverage_adjusted(draws_df, all_draws_post_prop[1,], sim_reduced_covid_data)
  #adj_zip_coverage_info <-
    #post_draws_get_prop_coverage_adjusted(
      #zip_code_draws, county_positive_prop_info$pos_prop[observed_zip],
      #aggregate_sim_data(sim_reduced_covid_data))
      #mrpN_obs_county_positive_prop_info) 
  #adj_zip_coverage_m[i,] <- adj_zip_coverage_info[[1]]
  #adj_zip_coverage_length_m[i,] <- adj_zip_coverage_info[[2]]

  tmp <- sim_reduced_covid_data
  tmp$total <- tmp$mrpN
  zip_code_draws <- post_draws_get_zip_code_draws(draws_df, tmp)
  prop_coverage_info <-
    post_draws_get_prop_coverage(
      inv.logit(zip_code_draws), mrpN_county_positive_prop_info$pos_prop[observed_zip])
  mrp_zip_coverage_m[i,] <- prop_coverage_info[[1]]
  mrp_zip_coverage_length_m[i,] <- prop_coverage_info[[2]]
  mrp_zip_bias_m[i,] <- colMeans(inv.logit(zip_code_draws)) -
    mrpN_county_positive_prop_info$pos_prop[observed_zip] 
  prop_coverage_info <-
    post_draws_get_direct_post_prop_coverage(zip_code_draws, county_positive_prop_info$zip_total, mrpN_county_positive_prop_info$pos_prop[observed_zip])
  direct_post_zip_coverage_m[i,] <- prop_coverage_info[[1]]
  direct_post_zip_coverage_length_m[i,] <- prop_coverage_info[[2]]
  adj_zip_norm_coverage_info <-
    post_draws_get_norm_prop_coverage_adjusted_multiple_non_parallel(
      #draws_df, mrpN_county_positive_prop_info$pos_prop[observed_zip],
      #sim_reduced_covid_data, reduced_covid_data = reduced_covid_data, 
      #zip_code_db = zip_code_db) 
      zip_code_draws, mrpN_county_positive_prop_info$pos_prop[observed_zip],
      mrpN_obs_county_positive_prop_info) 
      #aggregate_sim_data(sim_reduced_covid_data))
  adj_mrp_zip_norm_coverage_m[i,] <-
    adj_zip_norm_coverage_info[[1]][1:ncol(zip_code_draws) * 3]
  adj_mrp_zip_norm_coverage_length_m[i,] <-
    adj_zip_norm_coverage_info[[2]][1:ncol(zip_code_draws) * 3]
  adj_mrp_zip_norm_bias_m[i,] <-
    adj_zip_norm_coverage_info[[3]][1:ncol(zip_code_draws) * 3]
  adj_mrp_zip_norm_obs_coverage_m[i,] <-
    adj_zip_norm_coverage_info[[1]][1:ncol(zip_code_draws) * 3 - 2]
  adj_mrp_zip_norm_obs_coverage_length_m[i,] <-
    adj_zip_norm_coverage_info[[2]][1:ncol(zip_code_draws) * 3 - 2]
  adj_mrp_zip_norm_obs_bias_m[i,] <-
    adj_zip_norm_coverage_info[[3]][1:ncol(zip_code_draws) * 3 - 2]
  adj_mrp_zip_norm_post_coverage_m[i,] <-
    adj_zip_norm_coverage_info[[1]][1:ncol(zip_code_draws) * 3 - 1]
  adj_mrp_zip_norm_post_coverage_length_m[i,] <-
    adj_zip_norm_coverage_info[[2]][1:ncol(zip_code_draws) * 3 - 1]
  adj_mrp_zip_norm_post_bias_m[i,] <-
    adj_zip_norm_coverage_info[[3]][1:ncol(zip_code_draws) * 3 - 1]
  adj_zip_coverage_info <-
    post_draws_get_ac_prop_coverage_adjusted_multiple_non_parallel(
      #draws_df, mrpN_county_positive_prop_info$pos_prop[observed_zip],
      #sim_reduced_covid_data, reduced_covid_data = reduced_covid_data, 
      #zip_code_db = zip_code_db) 
      zip_code_draws, mrpN_county_positive_prop_info$pos_prop[observed_zip],
      mrpN_obs_county_positive_prop_info) 
      #aggregate_sim_data(sim_reduced_covid_data))
  adj_mrp_zip_ac_coverage_m[i,] <-
    adj_zip_coverage_info[[1]][1:ncol(zip_code_draws) * 3]
  adj_mrp_zip_ac_coverage_length_m[i,] <-
    adj_zip_coverage_info[[2]][1:ncol(zip_code_draws) * 3]
  adj_mrp_zip_ac_bias_m[i,] <-
    adj_zip_coverage_info[[3]][1:ncol(zip_code_draws) * 3]
  adj_mrp_zip_ac_obs_coverage_m[i,] <-
    adj_zip_coverage_info[[1]][1:ncol(zip_code_draws) * 3 - 2]
  adj_mrp_zip_ac_obs_coverage_length_m[i,] <-
    adj_zip_coverage_info[[2]][1:ncol(zip_code_draws) * 3 - 2]
  adj_mrp_zip_ac_obs_bias_m[i,] <-
    adj_zip_coverage_info[[3]][1:ncol(zip_code_draws) * 3 - 2]
  adj_mrp_zip_ac_post_coverage_m[i,] <-
    adj_zip_coverage_info[[1]][1:ncol(zip_code_draws) * 3 - 1]
  adj_mrp_zip_ac_post_coverage_length_m[i,] <-
    adj_zip_norm_coverage_info[[2]][1:ncol(zip_code_draws) * 3 - 1]
  adj_mrp_zip_ac_post_bias_m[i,] <-
    adj_zip_norm_coverage_info[[3]][1:ncol(zip_code_draws) * 3 - 1]
  adj_zip_coverage_info <-
    post_draws_get_prop_coverage_adjusted_multiple_non_parallel(
      #draws_df, mrpN_county_positive_prop_info$pos_prop[observed_zip],
      #sim_reduced_covid_data, reduced_covid_data = reduced_covid_data, 
      #zip_code_db = zip_code_db) 
      zip_code_draws, mrpN_county_positive_prop_info$pos_prop[observed_zip],
      mrpN_obs_county_positive_prop_info) 
      #aggregate_sim_data(sim_reduced_covid_data))
  adj_mrp_zip_coverage_m[i,] <-
    adj_zip_coverage_info[[1]][1:ncol(zip_code_draws) * 3]
  adj_mrp_zip_coverage_length_m[i,] <-
    adj_zip_coverage_info[[2]][1:ncol(zip_code_draws) * 3]
  adj_mrp_zip_bias_m[i,] <-
    adj_zip_coverage_info[[3]][1:ncol(zip_code_draws) * 3]
  adj_mrp_zip_obs_coverage_m[i,] <-
    adj_zip_coverage_info[[1]][1:ncol(zip_code_draws) * 3 - 2]
  adj_mrp_zip_obs_coverage_length_m[i,] <-
    adj_zip_coverage_info[[2]][1:ncol(zip_code_draws) * 3 - 2]
  adj_mrp_zip_obs_bias_m[i,] <-
    adj_zip_coverage_info[[3]][1:ncol(zip_code_draws) * 3 - 2]
  adj_mrp_zip_post_coverage_m[i,] <-
    adj_zip_coverage_info[[1]][1:ncol(zip_code_draws) * 3 - 1]
  adj_mrp_zip_post_coverage_length_m[i,] <-
    adj_zip_norm_coverage_info[[2]][1:ncol(zip_code_draws) * 3 - 1]
  adj_mrp_zip_post_bias_m[i,] <-
    adj_zip_norm_coverage_info[[3]][1:ncol(zip_code_draws) * 3 - 1]

warnings()

#new_gen_prop <- colMeans(inv.logit(draws_df))
#sim_data_100 <- sapply(1:100, function(i) {
#  rbinom(length(new_gen_prop), reduced_covid_data$total, new_gen_prop)      
#})

save(list=ls()[grep("_m$", ls())], file=result_file)

