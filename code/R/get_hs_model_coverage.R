library(tidyverse)
library(brms)

print("HS")

args = commandArgs(trailingOnly=TRUE)
sim_data_file = args[1]
result_file = args[2]

load("data_files/stan_indiana_covid_sim_data.Rdata")
load("data_files/zip_code_db.Rdata")
source("code/R/small_area_estimate_helper_functions.R")
load(sim_data_file)
#sim_data_100 <- spatial_sim_prop_m

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
    summarize(pos_prop = sum(prop * total) / sum(total))

mrpN_county_positive_prop_info <-
  positive_prop_info %>% group_by(zip) %>%
  summarize(pos_prop = sum(prop * total_mrp) / sum(total_mrp))

coverage_m <- matrix(0, nrow = 100, ncol = 30)
prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
prop_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
bias_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
zip_coverage_m <- matrix(0, nrow = 100, ncol = 30)
zip_coverage_length_m <- matrix(0, nrow = 100, ncol = 30)
adj_zip_coverage_m <- matrix(0, nrow = 100, ncol = 30)
adj_zip_coverage_length_m <- matrix(0, nrow = 100, ncol = 30)
mrp_zip_coverage_m <- matrix(0, nrow = 100, ncol = 30)
mrp_zip_coverage_length_m <- matrix(0, nrow = 100, ncol = 30)
adj_mrp_zip_coverage_m <- matrix(0, nrow = 100, ncol = 30)
adj_mrp_zip_coverage_length_m <- matrix(0, nrow = 100, ncol = 30)

new_stan_results <- 
    brm(positive | trials(positive + negative) ~ (1 | sex) + (1 | race) + (1 | age) + (1 | zip),
        data = reduced_covid_data, family = binomial(link="logit"),
        #prior = set_prior(horseshoe(scale_slab = 1e9), class = "sd"),
        control = list(adapt_delta = 0.99), cores = 4, iter = 0)
for (i in 1:100) {
  print(paste0("Iter: ", i))
  sim_reduced_covid_data <- reduced_covid_data
  sim_reduced_covid_data$positive <- sim_data_100[,i]
  sim_reduced_covid_data$negative <- 
    sim_reduced_covid_data$total - sim_reduced_covid_data$positive
  new_stan_results <- 
    update(new_stan_results,
	   positive | trials(positive + negative) ~ (1 | sex) + (1 | race) + (1 | age) + (1 | zip),
        newdata = sim_reduced_covid_data, family = binomial(link="logit"),
        prior = set_prior(horseshoe(scale_slab = 1e9), class = "sd"),
        control = list(adapt_delta = 0.95), cores = 4)

  bias_m[i,] <- all_draws_post_prop[1,] - colMeans(posterior_epred(new_stan_results))  
    
  # fit_summary_info <- new_stan_results$summary(
  #   variables = "loc_param",
  #   posterior::default_summary_measures(),
  #   extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
  # )
  coverage_m[i,] <- 
    rstanarm_results_get_coverage(new_stan_results, sim_reduced_covid_data)
  prop_coverage_info <-
    brms_get_prop_coverage(new_stan_results, all_draws_post_prop[1,])
  prop_coverage_m[i,] <- prop_coverage_info[[1]]
  prop_coverage_length_m[i,] <- prop_coverage_info[[2]]
  #prop_coverage_m[i,] <-
  #  brms_get_prop_coverage(new_stan_results, all_draws_post_prop[1,])
  adj_prop_coverage_info <-
    brms_get_prop_coverage_adjusted(
      new_stan_results, all_draws_post_prop[1,], sim_reduced_covid_data)
  adj_prop_coverage_m[i,] <- adj_prop_coverage_info[[1]]
  adj_prop_coverage_length_m[i,] <- adj_prop_coverage_info[[2]]

  zip_code_draws <- 
    brms_get_zip_code_draws(new_stan_results, sim_reduced_covid_data)
  prop_coverage_info <-
    post_draws_get_prop_coverage(zip_code_draws, county_positive_prop_info$pos_prop)
  zip_coverage_m[i,] <- prop_coverage_info[[1]]
  zip_coverage_length_m[i,] <- prop_coverage_info[[2]]
  #prop_coverage_m[i,] <-
  #  post_draws_get_prop_coverage(inv.logit(draws_df), all_draws_post_prop[1,]) 
  #adj_prop_coverage_m[i,] <-
  #  post_draws_get_prop_coverage_adjusted(draws_df, all_draws_post_prop[1,], sim_reduced_covid_data) 
  adj_zip_coverage_info <-
    post_draws_get_prop_coverage_adjusted(
      logit(zip_code_draws), county_positive_prop_info$pos_prop, 
      aggregate_sim_data(sim_reduced_covid_data))
  adj_zip_coverage_m[i,] <- adj_zip_coverage_info[[1]]
  adj_zip_coverage_length_m[i,] <- adj_zip_coverage_info[[2]]
  #adj_prop_coverage_m[i,] <-
  #  brms_get_prop_coverage_adjusted(
  #    new_stan_results, all_draws_post_prop[1,], sim_reduced_covid_data)


  tmp <- sim_reduced_covid_data
  tmp$total <- tmp$mrpN
  zip_code_draws <-
    brms_get_zip_code_draws(new_stan_results, tmp)
  prop_coverage_info <-
    post_draws_get_prop_coverage(zip_code_draws, mrpN_county_positive_prop_info$pos_prop)
  mrp_zip_coverage_m[i,] <- prop_coverage_info[[1]]
  mrp_zip_coverage_length_m[i,] <- prop_coverage_info[[2]]
  adj_zip_coverage_info <-
    post_draws_get_prop_coverage_adjusted(
      logit(zip_code_draws), mrpN_county_positive_prop_info$pos_prop,
      aggregate_sim_data(sim_reduced_covid_data))
  adj_mrp_zip_coverage_m[i,] <- adj_zip_coverage_info[[1]]
  adj_mrp_zip_coverage_length_m[i,] <- adj_zip_coverage_info[[2]]
}
warnings()
save(coverage_m, prop_coverage_m, prop_coverage_length_m, adj_prop_coverage_m, adj_prop_coverage_length_m, bias_m, zip_coverage_m, zip_coverage_length_m, adj_zip_coverage_m, adj_zip_coverage_length_m, mrp_zip_coverage_m, mrp_zip_coverage_length_m, adj_mrp_zip_coverage_m, adj_mrp_zip_coverage_length_m, file = result_file)
#save(coverage_m, prop_coverage_m, prop_coverage_length_m, adj_prop_coverage_m, adj_prop_coverage_length_m, bias_m, zip_coverage_m, zip_coverage_length_m, adj_zip_coverage_m, adj_zip_coverage_length_m, file = result_file)
