library(tidyverse)
library(brms)

load("data_files/stan_indiana_covid_sim_data.Rdata")
load("data_files/zip_code_db.Rdata")
source("small_area_estimate_helper_functions.R")

county_zip_code_info <-
  sapply(reduced_covid_data$zip, function(zip_code) {
     zip_code_db[zip_code_db$zipcode == zip_code,]$county
  })
reduced_covid_data$county <- county_zip_code_info 
reduced_covid_data <-
  reduced_covid_data %>% group_by(county, sex, race, age) %>% 
  summarize(total = sum(total), positive = sum(positive),
	    negative = sum(total - positive))

coverage_m <- matrix(0, nrow = 100, ncol = 30)
prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
prop_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_coverage_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
adj_prop_coverage_length_m <- matrix(0, nrow = 100, ncol = nrow(reduced_covid_data))
new_stan_results <- 
    brm(positive | trials(positive + negative) ~ (1 | sex) + (1 | race) + (1 | age) + (1 | county),
        data = reduced_covid_data, family = binomial(link="logit"),
        #prior = set_prior(horseshoe(scale_slab = 1e9), class = "sd"),
        control = list(adapt_delta = 0.99), cores = 4, iter = 0)
for (i in 1:100) {
  print(i)
  sim_reduced_covid_data <- reduced_covid_data
  positive_info <- data.frame("cases" = sim_data_100[,i], 
			      "county" = county_zip_code_info)
  county_positive_info <- 
    positive_info %>% group_by(county) %>% 
    summarize(positive = sum(cases)) 
  sim_reduced_covid_data$positive <- county_positive_info$positive 
  sim_reduced_covid_data$negative <- 
    sim_reduced_covid_data$total - sim_reduced_covid_data$positive
  new_stan_results <- 
    update(new_stan_results,
	   positive | trials(positive + negative) ~ (1 | sex) + (1 | race) + (1 | age) + (1 | county),
        newdata = sim_reduced_covid_data, family = binomial(link="logit"),
        #prior = set_prior(horseshoe(scale_slab = 1e9), class = "sd"),
        control = list(adapt_delta = 0.95), cores = 4)
    
  # fit_summary_info <- new_stan_results$summary(
  #   variables = "loc_param",
  #   posterior::default_summary_measures(),
  #   extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
  # )
  positive_prop_info <- data.frame("prop" = all_draws_post_prop[1,],
				   "total" = sim_reduced_covid_data$total,
                                   "county" = county_zip_code_info)
  county_positive_prop_info <-
    positive_prop_info %>% group_by(county) %>%
    summarize(pos_prop = sum(prop * total) / sum(total))
  coverage_m[i,] <- 
    rstanarm_results_get_coverage(new_stan_results, sim_reduced_covid_data)
  prop_coverage_info <-
    brms_get_prop_coverage(new_stan_results, 
			   county_positive_prop_info$pos_prop)
  prop_coverage_m[i,] <- prop_coverage_info[[1]]
  prop_coverage_length_m[i,] <- prop_coverage_info[[2]]
  #prop_coverage_m[i,] <-
  #  brms_get_prop_coverage(new_stan_results, all_draws_post_prop[1,])
  adj_prop_coverage_info <-
    brms_get_prop_coverage_adjusted(
      new_stan_results, county_positive_prop_info$pos_prop, sim_reduced_covid_data)
  adj_prop_coverage_m[i,] <- adj_prop_coverage_info[[1]]
  adj_prop_coverage_length_m[i,] <- adj_prop_coverage_info[[2]]
  #adj_prop_coverage_m[i,] <-
  #  brms_get_prop_coverage_adjusted(
  #    new_stan_results, all_draws_post_prop[1,], sim_reduced_covid_data)
}
warnings()
save(coverage_m, prop_coverage_m, prop_coverage_length_m, adj_prop_coverage_m, adj_prop_coverage_length_m, file = "result_files/base_coverage_m_3_coverage_county.Rdata")

