library(tidyverse)
library(brms)
library(cmdstanr)
library(boot)

args = commandArgs(trailingOnly=TRUE)
sim_data_file = args[1]
result_file = args[2]

load("data_files/stan_indiana_covid_bootstrap_sim_data.Rdata")
load("data_files/zip_code_db.Rdata")
source("code/R/gen_functions/small_area_estimate_helper_functions.R")

load(sim_data_file)
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

new_stan_results_data <- 
    brm(positive | trials(positive + negative) ~ (1 | sex) + (1 | race) + (1 | age) + (1 | zip),
        data = reduced_covid_data, family = binomial(link="logit"),
        prior = set_prior(horseshoe(scale_slab = 1e9), class = "sd"),
        control = list(adapt_delta = 0.99), cores = 4, iter = 0)

mod <- cmdstan_model("code/stan/default_brms_stan_code_normal.stan")

for (i in 1:100) {
  print(paste0("Iter: ", i))
  sim_reduced_covid_data <- reduced_covid_data_list[[i]]
  #sim_reduced_covid_data$positive <- sim_data_100[,i]
  #sim_reduced_covid_data$negative <-
  #  sim_reduced_covid_data$total - sim_reduced_covid_data$positive
  new_stan_results_data <- 
    update(new_stan_results_data,
	   #positive | trials(positive + negative) ~ 0 + Intercept + sex + (1 | race) + (1 | age) + (1 | zip),
	   positive | trials(positive + negative) ~ (1 | sex) + (1 | race) + (1 | age) + (1 | zip),
        newdata = sim_reduced_covid_data, family = binomial(link="logit"),
        #prior = set_prior(horseshoe(scale_slab = 1e9), class = "sd"),
        control = list(adapt_delta = 0.95), cores = 4, iter = 0)
  data_list <- standata(new_stan_results_data)
  data_list$M_4 = max(data_list$J_4)

  data_list$Y <- sim_reduced_covid_data$positive 
  new_stan_results <- mod$sample(
    data = data_list, adapt_delta = 0.95, max_treedepth = 12,
    parallel_chains = 4)

  # fit_summary_info <- new_stan_results$summary(
  #   variables = "post_draws",
  #   posterior::default_summary_measures(),
  #   extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
  # )
  new_stan_results$save_object(file = paste0(result_file, "_", i, ".Rdata"))
}
