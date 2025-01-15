library(tidyverse)
library(cmdstanr)

args = commandArgs(trailingOnly=TRUE)
real_data_file = args[1]
result_file = args[2]
l = as.numeric(args[3])

load(real_data_file)
#source("small_area_estimate_helper_functions.R")

#mod <- cmdstan_model("mrp_model.stan")
mod <- cmdstan_model("mrp_model_sample_sigma.stan")

data_list <- list("N" = nrow(input_data),
		  "y" = input_data$positive,
		  "n_sample" = input_data$total,
		  "K" = 7, 
		  "X" = input_data %>% select(urbanicity:sex),
		  "N_pop" = nrow(new_data),
		  "K_pop" = 7,
		  "X_pop" = new_data %>% select(urbanicity:sex),
		  "N_race" = max(input_data$race),
		  "J_race" = input_data$race,
                  "N_race_pop" = max(new_data$race),
		  "J_race_pop" = new_data$race,
		  "N_age" = max(input_data$age),
		  "J_age" = input_data$age,
                  "N_age_pop" = max(new_data$age),
		  "J_age_pop" = new_data$age,
		  "N_zip" = max(input_data$zip),
		  "J_zip" = input_data$zip,
                  "N_zip_pop" = max(new_data$zip),
		  "J_zip_pop" = new_data$zip,
		  "zip_code_distance" = zip_dist_m,
		  "sens" = 1,
		  "spec" = 1,
		  "l" = l)
		  #"sigma" = 1)

stan_results <- mod$sample(data = data_list, adapt_delta = 0.95, 
			   max_treedepth = 12, parallel_chains = 4)
stan_results$save_object(file = result_file)

#save(stan_results, file = result_file)
