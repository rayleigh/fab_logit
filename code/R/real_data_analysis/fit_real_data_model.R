library(tidyverse)
library(cmdstanr)
library(geosphere)

args = commandArgs(trailingOnly=TRUE)
real_data_file = args[1]
result_file = args[2]
l = as.numeric(args[3])

load(real_data_file)
load("real_data_files/zip_code_db.Rdata")

input_data <- preprocess_data 
input_data <- input_data[input_data$state %in% 17:18,]

num_zip_codes = length(unique(input_data$zip))
zip_dist_m <- matrix(0, nrow = num_zip_codes, ncol = num_zip_codes) 

for (i in 1:(num_zip_codes - 1)) {
  i_ind <- which(zip_code_db$zipcode == unique(input_data$zip)[i])
  for (j in (i + 1):num_zip_codes) {
    j_ind <- which(zip_code_db$zipcode == unique(input_data$zip)[j])
    zip_dist_m[i, j] = distHaversine(
      zip_code_db[c(i_ind, j_ind), c("lng", "lat")], r = 3961)
    zip_dist_m[j, i] = zip_dist_m[i, j]
  }
}

print(input_data)

mod <- cmdstan_model("code/stan/mrp_model_simple_sample_sigma.stan")

data_list <- list("N" = nrow(input_data),
		  "y" = input_data$positive,
		  "K" = 7, 
		  "X" = input_data %>% select(urbanicity:sex_level),
		  "n_sample" = input_data$total,
		  "N_race" = max(input_data$race_level),
		  "J_race" = input_data$race_level,
		  "N_age" = max(input_data$age_level),
		  "J_age" = input_data$age_level,
		  "N_zip" = max(input_data$zip_level),
		  "J_zip" = input_data$zip_level,
		  "zip_code_distance" = zip_dist_m,
		  "sens" = 1,
		  "spec" = 1,
		  "l" = l)

stan_results <- mod$sample(data = data_list, adapt_delta = 0.95, 
			   max_treedepth = 12, parallel_chains = 4)
stan_results$save_object(file = result_file)

