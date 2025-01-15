library(boot)
library(parallel)
library(tidyverse)
library(boot)

source("small_area_estimate_helper_functions.R")

args = commandArgs(trailingOnly=TRUE)
stan_result_file = args[1]
ci_result_file = args[2]

print(stan_result_file)
load("real_data_files/interface_preprocess_peak_data.Rdata")
stan_results <- readRDS(stan_result_file)

draws_df <- stan_results$draws(format = "df")
cell_draws_df <- data.matrix(draws_df[, grep("^p\\[", colnames(draws_df))])

ci_intervals <- apply(cell_draws_df, 2, quantile, probs = c(0.025, 0.975))

obs_prop <- input_data$positive / input_data$total
post_means <- colMeans(logit(cell_draws_df))
post_sd <- apply(logit(cell_draws_df), 2, sd)
adj_ci_intervals <-
  do.call(cbind, mclapply(1:ncol(cell_draws_df), function(i) {
  # sapply(1:ncol(cell_draws_df), function(i) {
    print(i)
    #if (input_data$total[i] < 3) {
    #  obs_prop[i] = inv.logit(post_means[i])
    #}
    derive_fab_intervals_multiple_y(
      c(obs_prop[i], inv.logit(post_means[i])), 0.05, 
      input_data$total[i], post_means[i], post_sd[i])
    # })
}, mc.cores = 4))

direct_ci_intervals <-
  sapply(1:ncol(cell_draws_df), function(i) {
    print(i)
    derive_intervals(obs_prop[i], 0.05, 0.5, input_data$total[i])
  })

#Zip code draws
#pop_draws_df <- data.matrix(draws_df[, grep("^p_pop\\[", colnames(draws_df))])
#zip_code_draws <- post_draws_get_zip_code_draws(logit(pop_draws_df), new_data)
draws_df <- stan_results$draws(format = "df")
cell_draws_df <- data.matrix(draws_df[, grep("^p\\[", colnames(draws_df))])
tmp <- input_data
tmp_dem_data <- new_data %>% unite("dem_label", sex:zip, remove = F)
tmp <- tmp %>% unite("dem_label", sex:zip, remove = F)
tmp <- merge(tmp, tmp_dem_data, by = 'dem_label', all.x = T, all.y= F)
tmp$zip <- tmp$zip.x
tmp$total <- tmp$total.y
#tmp <- tmp %>% arrange(zip)
zip_code_draws <- post_draws_get_zip_code_draws(logit(cell_draws_df), tmp)
zip_code_input_data <- aggregate_sim_data(input_data)

obs_prop <- zip_code_input_data$positive / zip_code_input_data$total
post_means <- colMeans(zip_code_draws)
post_sd <- apply(zip_code_draws, 2, sd)
zip_adj_ci_intervals <-
  do.call(cbind, mclapply(1:ncol(zip_code_draws), function(i) {
  # sapply(1:ncol(zip_code_draws), function(i) {
    print(i)
    #if (zip_code_input_data$total[i] < 3) {
    #  obs_prop[i] = inv.logit(post_means[i])
    #}
    derive_fab_intervals_multiple_y(
      c(obs_prop[i], inv.logit(post_means[i])), 0.05, 
      zip_code_input_data$total[i], post_means[i], post_sd[i])
  }, mc.cores = 4))

zip_direct_ci_intervals <-
  sapply(1:length(obs_prop), function(i) {
    print(i)
    derive_intervals(obs_prop[i], 0.05, 0.5, zip_code_input_data$total[i])
  })

zip_ci_intervals <- apply(inv.logit(zip_code_draws), 2, quantile, probs = c(0.025, 0.975))

save(ci_intervals, adj_ci_intervals, direct_ci_intervals,
     zip_ci_intervals, zip_adj_ci_intervals, zip_direct_ci_intervals,
#save(zip_ci_intervals, zip_adj_ci_intervals, zip_direct_ci_intervals,
     file = ci_result_file)

