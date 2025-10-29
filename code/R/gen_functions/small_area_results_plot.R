library(tidyverse)
library(latex2exp)
library(reshape2)
library(boot)

get_ci_coverage_info <- function(file_list, num_cells, prop, obs,
                                 chunk_size = 1, num_exp = length(file_list)) {

  adj_wilson_ci_m <- matrix(0, nrow = num_exp, ncol = num_cells)
  adj_wald_ci_m <- matrix(0, nrow = num_exp, ncol = num_cells)
  adj_ac_ci_m <- matrix(0, nrow = num_exp, ncol = num_cells)
  for (i in 1:length(file_list)) {
    load(file_list[i])
    iter = (i - 1) * chunk_size + 1:chunk_size
    if (prop) {
      if (obs) {
        adj_wilson_ci_m[iter,] <- adj_prop_obs_coverage_m[iter,]
        adj_wald_ci_m[iter,] <- adj_norm_prop_obs_coverage_m[iter,]
        adj_ac_ci_m[iter,] <- adj_ac_prop_obs_coverage_m[iter,]
      } else {
        adj_wilson_ci_m[iter,] <- adj_prop_post_coverage_m[iter,]
        adj_wald_ci_m[iter,] <- adj_norm_prop_post_coverage_m[iter,]
        adj_ac_ci_m[iter,] <- adj_ac_prop_post_coverage_m[iter,]
      }
    } else {
      if (obs) {
        adj_wilson_ci_m[iter,] <- adj_mrp_zip_obs_coverage_m[iter,]
        adj_wald_ci_m[iter,] <- adj_mrp_zip_norm_obs_coverage_m[iter,]
        adj_ac_ci_m[iter,] <- adj_mrp_zip_ac_obs_coverage_m[iter,]
      } else {
        adj_wilson_ci_m[iter,] <- adj_mrp_zip_post_coverage_m[iter,]
        adj_wald_ci_m[iter,] <- adj_mrp_zip_norm_post_coverage_m[iter,]
        adj_ac_ci_m[iter,] <- adj_mrp_zip_ac_post_coverage_m[iter,]
      }
    }
  }
  adj_ci_post_credible_coverage <-
    cbind(colMeans(adj_wilson_ci_m), colMeans(adj_wald_ci_m),
          colMeans(adj_ac_ci_m))
  colnames(adj_ci_post_credible_coverage) <- c("Wilson", "Wald", "AC")
  return(adj_ci_post_credible_coverage)
}

get_bootstrap_ci_coverage_info <- function(file_list, num_cells, prop, obs, zip_code_list) {
  
  adj_wilson_ci_m <- matrix(NA, nrow = length(file_list), ncol = num_cells)
  adj_wald_ci_m <- matrix(NA, nrow = length(file_list), ncol = num_cells)
  adj_ac_ci_m <- matrix(NA, nrow = length(file_list), ncol = num_cells)
  for (i in 1:length(file_list)) {
    load(file_list[i])
    if (prop) {
      if (obs) {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_prop_obs_coverage_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_norm_prop_obs_coverage_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_ac_prop_obs_coverage_m[i,]
      } else {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_prop_post_coverage_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_norm_prop_post_coverage_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_ac_prop_post_coverage_m[i,]
      }
    } else {
      if (obs) {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_obs_coverage_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_norm_obs_coverage_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_ac_obs_coverage_m[i,]
      } else {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_post_coverage_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_norm_post_coverage_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_ac_post_coverage_m[i,]
      }
    }
  }
  adj_ci_post_credible_coverage <- 
    cbind(colMeans(adj_wilson_ci_m, na.rm = T), colMeans(adj_wald_ci_m, na.rm = T),
          colMeans(adj_ac_ci_m, na.rm = T))
  colnames(adj_ci_post_credible_coverage) <- c("Wilson", "Wald", "AC")
  return(adj_ci_post_credible_coverage)
}

make_coverage_direct_draws_multiple_adj_ci_box_plot <- function(
    file_list, direct_draws, result_plot_file, prop = T) {
  
  ci_credible_coverage <- sapply(file_list, function(file) {
    load(file)
    if (prop) {
      return(colMeans(prop_coverage_m))
    } else {
      return(colMeans(mrp_zip_coverage_m))
    }
  })
  colnames(ci_credible_coverage) <- c("Default", "HS", "GL Lasso", "GP")
  
  adj_ci_post_credible_coverage <- sapply(file_list, function(file) {
    load(file)
    if (prop) {
      # adj_list <- list(adj_prop_coverage_m, 
      #                  adj_prop_obs_coverage_m,
      #                  adj_prop_post_coverage_m)
      # return(colMeans(adj_list[[ind]]))
      return(colMeans(adj_prop_post_coverage_m))
    } else {
      # adj_list <- list(adj_mrp_zip_coverage_m, 
      #                  adj_mrp_zip_obs_coverage_m,
      #                  adj_mrp_zip_post_coverage_m)
      return(colMeans(adj_mrp_zip_post_coverage_m))
    }
  })
  colnames(adj_ci_post_credible_coverage) <- c("Default", "HS", "GL Lasso", "GP")
  
  adj_ci_obs_credible_coverage <- sapply(file_list, function(file) {
    load(file)
    if (prop) {
      # adj_list <- list(adj_prop_coverage_m, 
      #                  adj_prop_obs_coverage_m,
      #                  adj_prop_post_coverage_m)
      return(colMeans(adj_prop_obs_coverage_m))
    } else {
      # adj_list <- list(adj_mrp_zip_coverage_m, 
      #                  adj_mrp_zip_obs_coverage_m,
      #                  adj_mrp_zip_post_coverage_m)
      return(colMeans(adj_mrp_zip_obs_coverage_m))
    }
  })
  colnames(adj_ci_obs_credible_coverage) <- c("Default", "HS", "GL Lasso", "GP")
  
  
  plot_df <- rbind(melt(ci_credible_coverage), 
                   melt(adj_ci_post_credible_coverage),
                   melt(adj_ci_obs_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(direct_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(direct_draws))
  plot_df <- rbind(plot_df, direct_draw_info)
  
  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))), 
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("adj_obs", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(direct_draws)))
  plot_df <- 
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(0:3 * 3.5 + 0.5, each = nrow(ci_credible_coverage)),
                        rep(0:3 * 3.5 + 1.5, each = nrow(adj_ci_post_credible_coverage)),
                        rep(0:3 * 3.5 + 2.5, each = nrow(adj_ci_obs_credible_coverage)),
                        rep(-1, ncol(direct_draws)))
  
  # print(plot_df)
  
  pdf(result_plot_file, height = 9, width = 15)
  print(ggplot(data = plot_df) + 
          geom_hline(yintercept = 0.95, size = 1.25, linetype = "longdash") +
          geom_boxplot(aes(y = value, x = position, 
                           group = label, color = group),
                       size = 2) + xlab("") +
          theme_bw() + ylab("Coverage prob") + ylim(c(0, 1)) +
          scale_x_continuous(breaks = c(-1, 1.5 + 0:3 * 3.5), 
                             labels = c("Direct", "Default", "HS", "GL Lasso", "GP")) + 
          # values=c("black", "grey60", "grey40")) +
          # scale_color_manual(values = c("adj" = "black",
          scale_color_manual(values = c("adj_post" = "grey20",
                                        "adj_obs" = "grey60",
                                        "normal" = "grey80",
                                        "dir" = "grey40"), 
                             guide = "none") +
          theme(axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                axis.title.x = element_text(size = 30),
                axis.title.y = element_text(size = 30)))
  dev.off()
}

make_coverage_credible_intervals_multiple_intervals_adj_ci_box_plot <- function(
    file_list, credible_draws, norm_coverage_m, norm_prop_coverage_m, norm_ac_coverage_m,
    result_plot_file, prop = T, obs = T, 
    bootstrap = F, zip_code_list = NULL) {
  
  ci_credible_coverage <- cbind(colMeans(norm_coverage_m, na.rm = T),
                                colMeans(norm_prop_coverage_m, na.rm = T),
                                colMeans(norm_ac_coverage_m, na.rm = T))
  colnames(ci_credible_coverage) <- c("Wilson", "Wald", "AC")
  
  if (bootstrap) {
    adj_ci_post_credible_coverage <- get_bootstrap_ci_coverage_info(
      file_list, ncol(norm_prop_coverage_m), prop, obs, zip_code_list)
  } else {
    adj_ci_post_credible_coverage <- get_ci_coverage_info(
      file_list, ncol(norm_prop_coverage_m), prop, obs)
  }
  
  plot_df <- rbind(melt(ci_credible_coverage), 
                   melt(adj_ci_post_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws, na.rm = T))
  plot_df <- rbind(plot_df, direct_draw_info)
  
  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))), 
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(credible_draws)))
  plot_df <- 
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(0:2 * 2.5 + 0.5, each = nrow(ci_credible_coverage)),
                        rep(0:2 * 2.5 + 1.5, each = nrow(adj_ci_post_credible_coverage)),
                        rep(-1, ncol(credible_draws)))
  
  # print(plot_df)
  
  pdf(result_plot_file, height = 9, width = 15)
  print(ggplot(data = plot_df) + 
          geom_hline(yintercept = 0.95, size = 1.25, linetype = "longdash") +
          geom_boxplot(aes(y = value, x = position, 
                           group = label, color = group),
                       size = 2) + xlab("") +
          theme_bw() + ylab("Coverage prob") + ylim(c(0, 1)) +
          scale_x_continuous(breaks = c(-1, 1 + 0:2 * 2.5), 
                             labels = c("CI", "Wilson", "Wald", "AC")) + 
          # values=c("black", "grey60", "grey40")) +
          # scale_color_manual(values = c("adj" = "black",
          scale_color_manual(values = c("adj_post" = "black",
                                        # "adj_obs" = "grey60",
                                        "normal" = "grey80",
                                        "dir" = "grey40"), 
                             guide = "none") +
          theme(axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                axis.title.x = element_text(size = 30),
                axis.title.y = element_text(size = 30)))
  dev.off()
}


make_coverage_multiple_intervals_comp_only_adj_ci_box_plot <- function(
    prior_1_file_list, prior_2_file_list,
    credible_draws_prior_1, credible_draws_prior_2,
    result_plot_file, prop = T, obs = T,
    bootstrap = F, zip_code_list = NULL) {

  if (bootstrap) {
    adj_ci_post_credible_coverage <- get_bootstrap_ci_coverage_info(
      prior_1_file_list, ncol(credible_draws_prior_1), prop, obs, zip_code_list)
  } else {
    adj_ci_post_credible_coverage <- get_ci_coverage_info(
      prior_1_file_list, ncol(credible_draws_prior_1), prop, obs)
  }

  if (bootstrap) {
    adj_ci_post_credible_coverage_2 <- get_bootstrap_ci_coverage_info(
      prior_2_file_list, ncol(credible_draws_prior_1), prop, obs, zip_code_list)
  } else {
    adj_ci_post_credible_coverage_2 <- get_ci_coverage_info(
      prior_2_file_list, ncol(credible_draws_prior_1), prop, obs)
  }

  plot_df <- rbind(melt(adj_ci_post_credible_coverage),
                   melt(adj_ci_post_credible_coverage_2))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws_prior_1),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws_prior_1,
						    na.rm = T))
  plot_df <- rbind(plot_df, direct_draw_info)
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws_prior_2),
                                 "Var2" = "Direct 2",
                                 "value" = colMeans(credible_draws_prior_2,
						    na.rm = T))
  plot_df <- rbind(plot_df, direct_draw_info)

  # direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
  #                                "Var2" = "Direct",
  #                                "value" = colMeans(credible_draws))
  # plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(adj_ci_post_credible_coverage))),
                     rep("adj_post", prod(dim(adj_ci_post_credible_coverage))),
                     rep("dir", 2 * ncol(credible_draws_prior_1)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(0:2 * 2.5 + 0.5, each = ncol(credible_draws_prior_1)),
                        rep(0:2 * 2.5 + 1.5, each = ncol(credible_draws_prior_1)),
                        rep(-2, ncol(credible_draws_prior_1)),
                        rep(-1, ncol(credible_draws_prior_2)))

  # print(plot_df)

  pdf(result_plot_file, height = 9, width = 15)
  print(ggplot(data = plot_df) +
          geom_hline(yintercept = 0.95, size = 1.25, linetype = "longdash") +
          geom_boxplot(aes(y = value, x = position,
                           group = label, color = group),
                       size = 2) + xlab("") +
          theme_bw() + ylab("Coverage prob") + ylim(c(0, 1)) +
          scale_x_continuous(breaks = c(-1.5, 1 + 0:2 * 2.5),
                             labels = c("CI", "Wilson", "Wald", "AC")) +
          # values=c("black", "grey60", "grey40")) +
          # scale_color_manual(values = c("adj" = "black",
          scale_color_manual(values = c("adj_post" = "grey60",
                                        # "adj_obs" = "grey60",
                                        "normal" = "black",
                                        "dir" = "grey40"),
                             guide = "none") +
          theme(axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                axis.title.x = element_text(size = 30),
                axis.title.y = element_text(size = 30)))
  dev.off()
}

get_ci_bias_info <- function(
    file_list, num_cells, prop, obs, abs = F,
    chunk_size = 1, num_exp = length(file_list)) {

  adj_wilson_ci_m <- matrix(0, nrow = num_exp, ncol = num_cells)
  adj_wald_ci_m <- matrix(0, nrow = num_exp, ncol = num_cells)
  adj_ac_ci_m <- matrix(0, nrow = num_exp, ncol = num_cells)

  for (i in 1:length(file_list)) {
    # print(i)
    iter = (i - 1) * chunk_size + 1:chunk_size
    load(file_list[i])
    if (prop) {
      if (obs) {
        adj_wilson_ci_m[iter,] <- adj_prop_obs_bias_m[iter,]
        adj_wald_ci_m[iter,] <- adj_norm_prop_obs_bias_m[iter,]
        adj_ac_ci_m[iter,] <- adj_ac_prop_obs_bias_m[iter,]
      } else {
        adj_wilson_ci_m[iter,] <- adj_prop_post_bias_m[iter,]
        adj_wald_ci_m[iter,] <- adj_norm_prop_post_bias_m[iter,]
        adj_ac_ci_m[iter,] <- adj_ac_prop_post_bias_m[iter,]
      }
    } else {
      if (obs) {
        adj_wilson_ci_m[iter,] <- adj_mrp_zip_obs_bias_m[iter,]
        adj_wald_ci_m[iter,] <- adj_mrp_zip_norm_obs_bias_m[iter,]
        adj_ac_ci_m[iter,] <- adj_mrp_zip_ac_obs_bias_m[iter,]
      } else {
        adj_wilson_ci_m[iter,] <- adj_mrp_zip_post_bias_m[iter,]
        adj_wald_ci_m[iter,] <- adj_mrp_zip_norm_post_bias_m[iter,]
        adj_ac_ci_m[iter,] <- adj_mrp_zip_ac_post_bias_m[iter,]
      }
    }
  }
  if (abs) {
    adj_ci_post_credible_coverage <- cbind(colMeans(abs(adj_wilson_ci_m)),
                                           colMeans(abs(adj_wald_ci_m)),
                                           colMeans(abs(adj_ac_ci_m)))
  } else {
    adj_ci_post_credible_coverage <- cbind(colMeans(adj_wilson_ci_m^2),
                                           colMeans(adj_wald_ci_m^2),
                                           colMeans(adj_ac_ci_m^2))
  }
  colnames(adj_ci_post_credible_coverage) <- c("Wilson", "Wald", "AC")
  return(adj_ci_post_credible_coverage)
}

get_bootstrap_ci_bias_info <- function(file_list, num_cells, prop, obs, zip_code_list, abs = F) {

  adj_wilson_ci_m <- matrix(NA, nrow = length(file_list), ncol = num_cells)
  adj_wald_ci_m <- matrix(NA, nrow = length(file_list), ncol = num_cells)
  adj_ac_ci_m <- matrix(NA, nrow = length(file_list), ncol = num_cells)
  for (i in 1:length(file_list)) {
    # print(i)
    load(file_list[i])
    if (prop) {
      if (obs) {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_prop_obs_bias_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_norm_prop_obs_bias_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_ac_prop_obs_bias_m[i,]
      } else {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_prop_post_bias_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_norm_prop_post_bias_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_ac_prop_post_bias_m[i,]
      }
    } else {
      if (obs) {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_obs_bias_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_norm_obs_bias_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_ac_obs_bias_m[i,]
      } else {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_post_bias_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_norm_post_bias_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_ac_post_bias_m[i,]
      }
    }
  }
  if (abs) {
    adj_ci_post_credible_coverage <- cbind(colMeans(abs(adj_wilson_ci_m), na.rm = T),
                                           colMeans(abs(adj_wald_ci_m), na.rm = T),
                                           colMeans(abs(adj_ac_ci_m), na.rm = T))
  } else {
    adj_ci_post_credible_coverage <- cbind(colMeans(adj_wilson_ci_m^2, na.rm = T),
                                           colMeans(adj_wald_ci_m^2, na.rm = T),
                                           colMeans(adj_ac_ci_m^2, na.rm = T))
  }
  colnames(adj_ci_post_credible_coverage) <- c("Wilson", "Wald", "AC")
  return(adj_ci_post_credible_coverage)
}

make_bias_credible_intervals_multiple_intervals_adj_ci_box_plot <- function(
    file_list, credible_draws, norm_coverage_m, norm_prop_coverage_m, norm_ac_coverage_m,
    result_plot_file, prop = T, obs = T, abs = F, ylim_info = c(0, 1),
    bootstrap = F, zip_code_list = NULL) {

  ci_credible_coverage <- cbind(colMeans(norm_coverage_m, na.rm = T),
                                colMeans(norm_prop_coverage_m, na.rm = T),
                                colMeans(norm_ac_coverage_m, na.rm = T))
  colnames(ci_credible_coverage) <- c("Wilson", "Wald", "AC")

  if (bootstrap) {
  adj_ci_post_credible_coverage <-
    get_bootstrap_ci_bias_info(file_list,  ncol(norm_prop_coverage_m), prop, obs, zip_code_list, abs)

  } else {
  adj_ci_post_credible_coverage <-
    get_ci_bias_info(file_list,  ncol(norm_prop_coverage_m), prop, obs, abs)
  }
  # colnames(adj_ci_post_credible_coverage) <- c("Wilson", "Wald", "AC")

  plot_df <- rbind(melt(ci_credible_coverage),
                   melt(adj_ci_post_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws, na.rm = T))
  plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))),
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(credible_draws)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(0:2 * 2.5 + 0.5, each = nrow(ci_credible_coverage)),
                        rep(0:2 * 2.5 + 1.5, each = nrow(adj_ci_post_credible_coverage)),
                        rep(-1, ncol(credible_draws)))

  # print(plot_df)

  if (!abs) {
    plot_df$value <- sqrt(plot_df$value)
    y_label <- c("Root mean squared error")
  } else {
    y_label <- c("Mean absolute bias")
  }

  pdf(result_plot_file, height = 9, width = 15)
  print(ggplot(data = plot_df) +
          geom_boxplot(aes(y = value, x = position,
                           group = label, color = group),
                       size = 2) + xlab("") +
          theme_bw() + ylab(y_label) +
          ylim(ylim_info) +
          scale_x_continuous(breaks = c(-1, 1 + 0:2 * 2.5),
                             labels = c("CI", "Wilson", "Wald", "AC")) +
          # values=c("black", "grey60", "grey40")) +
          # scale_color_manual(values = c("adj" = "black",
          scale_color_manual(values = c("adj_post" = "black",
                                        # "adj_obs" = "grey60",
                                        "normal" = "grey80",
                                        "dir" = "grey40"),
                             guide = "none") +
          theme(axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                axis.title.x = element_text(size = 30),
                axis.title.y = element_text(size = 30)))
  dev.off()
}

make_bias_multiple_intervals_comp_only_adj_ci_box_plot <- function(
    prior_1_file_list, prior_2_file_list,
    credible_draws_prior_1, credible_draws_prior_2,
    result_plot_file, prop = T, obs = T, abs = F, ylim_info = c(0, 1),
    bootstrap = T, zip_code_list = NULL) {


  if (bootstrap) {
  adj_ci_post_credible_coverage <-
    get_bootstrap_ci_bias_info(prior_1_file_list, ncol(credible_draws_prior_1),
                     prop, obs, zip_code_list, abs)
  } else {
  adj_ci_post_credible_coverage <-
    get_ci_bias_info(prior_1_file_list, ncol(credible_draws_prior_1), 
		     prop, obs, abs)
  }

  if (bootstrap) {
    adj_ci_post_credible_coverage_2 <-
    get_bootstrap_ci_bias_info(prior_2_file_list, ncol(credible_draws_prior_1),
                     prop, obs, zip_code_list, abs)
  } else {
  adj_ci_post_credible_coverage_2 <-
    get_ci_bias_info(prior_2_file_list, ncol(credible_draws_prior_1), 
		     prop, obs, abs)
  }

  plot_df <- rbind(melt(adj_ci_post_credible_coverage),
                   melt(adj_ci_post_credible_coverage_2))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws_prior_1),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws_prior_1, 
						    na.rm = T))
  plot_df <- rbind(plot_df, direct_draw_info)
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws_prior_2),
                                 "Var2" = "Direct 2",
                                 "value" = colMeans(credible_draws_prior_2,
						    na.rm = T))
  plot_df <- rbind(plot_df, direct_draw_info)

  # direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
  #                                "Var2" = "Direct",
  #                                "value" = colMeans(credible_draws))
  # plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(adj_ci_post_credible_coverage))),
                     rep("adj_post", prod(dim(adj_ci_post_credible_coverage))),
                     rep("dir", 2 * ncol(credible_draws_prior_1)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(0:2 * 2.5 + 0.5, each = ncol(credible_draws_prior_1)),
                        rep(0:2 * 2.5 + 1.5, each = ncol(credible_draws_prior_1)),
                        rep(-2, ncol(credible_draws_prior_1)),
                        rep(-1, ncol(credible_draws_prior_2)))

  if (!abs) {
    plot_df$value <- sqrt(plot_df$value)
    y_label <- c("Root mean squared error")
  } else {
    y_label <- c("Mean absolute bias")
  }

  # print(plot_df)

  pdf(result_plot_file, height = 9, width = 15)
  print(ggplot(data = plot_df) +
          # geom_hline(yintercept = 0.95, size = 1.25, linetype = "longdash") +
          geom_boxplot(aes(y = value, x = position,
                           group = label, color = group),
                       size = 2) + xlab("") +
          theme_bw() + ylab(y_label) +
          ylim(ylim_info) +
          scale_x_continuous(breaks = c(-1.5, 1 + 0:2 * 2.5),
                             labels = c("CI", "Wilson", "Wald", "AC")) +
          # values=c("black", "grey60", "grey40")) +
          # scale_color_manual(values = c("adj" = "black",
          scale_color_manual(values = c("adj_post" = "grey60",
                                        # "adj_obs" = "grey60",
                                        "normal" = "black",
                                        "dir" = "grey40"),
                             guide = "none") +
          theme(axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                axis.title.x = element_text(size = 30),
                axis.title.y = element_text(size = 30)))
  dev.off()
}

get_ci_coverage_length_info <- function(
    file_list, num_cells, prop, obs,
    chunk_size = 1, num_exp = length(file_list)) {

  adj_wilson_ci_m <- matrix(0, nrow = num_exp, ncol = num_cells)
  adj_wald_ci_m <- matrix(0, nrow = num_exp, ncol = num_cells)
  adj_ac_ci_m <- matrix(0, nrow = num_exp, ncol = num_cells)
  for (i in 1:length(file_list)) {
    # print(i)
    load(file_list[i])
    iter = (i - 1) * chunk_size + 1:chunk_size
    if (prop) {
      if (obs) {
        adj_wilson_ci_m[iter,] <- adj_prop_obs_coverage_length_m[iter,]
        adj_wald_ci_m[iter,] <- adj_norm_prop_obs_coverage_length_m[iter,]
        adj_ac_ci_m[iter,] <- adj_ac_prop_obs_coverage_length_m[iter,]
      } else {
        adj_wilson_ci_m[iter,] <- adj_prop_post_coverage_length_m[iter,]
        adj_wald_ci_m[iter,] <- adj_norm_prop_post_coverage_length_m[iter,]
        adj_ac_ci_m[iter,] <- adj_ac_prop_post_coverage_length_m[iter,]
      }
    } else {
      if (obs) {
        adj_wilson_ci_m[iter,] <- adj_mrp_zip_obs_coverage_length_m[iter,]
        adj_wald_ci_m[iter,] <- adj_mrp_zip_norm_obs_coverage_length_m[iter,]
        adj_ac_ci_m[iter,] <- adj_mrp_zip_ac_obs_coverage_length_m[iter,]
      } else {
        adj_wilson_ci_m[iter,] <- adj_mrp_zip_post_coverage_length_m[iter,]
        adj_wald_ci_m[iter,] <- adj_mrp_zip_norm_post_coverage_length_m[iter,]
        adj_ac_ci_m[iter,] <- adj_mrp_zip_ac_post_coverage_length_m[iter,]
      }
    }
  }
  adj_ci_post_credible_coverage <- cbind(colMeans(adj_wilson_ci_m),
                                         colMeans(adj_wald_ci_m),
                                         colMeans(adj_ac_ci_m))
  colnames(adj_ci_post_credible_coverage) <- c("Wilson", "Wald", "AC")
  return(adj_ci_post_credible_coverage)
}

get_bootstrap_ci_coverage_length_info <- function(file_list, num_cells, prop, obs, zip_code_list) {

  adj_wilson_ci_m <- matrix(NA, nrow = length(file_list), ncol = num_cells)
  adj_wald_ci_m <- matrix(NA, nrow = length(file_list), ncol = num_cells)
  adj_ac_ci_m <- matrix(NA, nrow = length(file_list), ncol = num_cells)
  for (i in 1:length(file_list)) {
    # print(i)
    load(file_list[i])
    if (prop) {
      if (obs) {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_prop_obs_coverage_length_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_norm_prop_obs_coverage_length_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_ac_prop_obs_coverage_length_m[i,]
      } else {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_prop_post_coverage_length_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_norm_prop_post_coverage_length_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_ac_prop_post_coverage_length_m[i,]
      }
    } else {
      if (obs) {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_obs_coverage_length_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_norm_obs_coverage_length_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_ac_obs_coverage_length_m[i,]
      } else {
        adj_wilson_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_post_coverage_length_m[i,]
        adj_wald_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_norm_post_coverage_length_m[i,]
        adj_ac_ci_m[i,zip_code_list[[i]]] <- adj_mrp_zip_ac_post_coverage_length_m[i,]
      }
    }
  }
  adj_ci_post_credible_coverage <- cbind(colMeans(adj_wilson_ci_m, na.rm = T),
                                         colMeans(adj_wald_ci_m, na.rm = T),
                                         colMeans(adj_ac_ci_m, na.rm = T))
  colnames(adj_ci_post_credible_coverage) <- c("Wilson", "Wald", "AC")
  return(adj_ci_post_credible_coverage)
}

make_coverage_credible_intervals_length_multiple_intervals_adj_ci_box_plot <- function(
    file_list, credible_draws, norm_coverage_m, norm_prop_coverage_m, norm_ac_coverage_m,
    result_plot_file, prop = T, obs = T, ylim_info = c(0, 1),
    bootstrap = F, zip_code_list = NULL) {

  ci_credible_coverage <- cbind(colMeans(norm_coverage_m, na.rm = T),
                                colMeans(norm_prop_coverage_m, na.rm = T),
                                colMeans(norm_ac_coverage_m, na.rm = T))
  colnames(ci_credible_coverage) <- c("Wilson", "Wald", "AC")

  if (bootstrap) {
    adj_ci_post_credible_coverage <-
      get_bootstrap_ci_coverage_length_info(
        file_list, ncol(norm_prop_coverage_m), prop, obs, zip_code_list)
  } else {
    adj_ci_post_credible_coverage <-
      get_ci_coverage_length_info(
        file_list, ncol(norm_prop_coverage_m), prop, obs)
  }

  plot_df <- rbind(melt(ci_credible_coverage),
                   melt(adj_ci_post_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws, na.rm = T))
  plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))),
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(credible_draws)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(0:2 * 2.5 + 0.5, each = nrow(ci_credible_coverage)),
                        rep(0:2 * 2.5 + 1.5, each = nrow(adj_ci_post_credible_coverage)),
                        rep(-1, ncol(credible_draws)))

  pdf(result_plot_file, height = 9, width = 15)
  print(ggplot(data = plot_df) +
          geom_boxplot(aes(y = value, x = position,
                           group = label, color = group),
                       size = 2) + xlab("") +
          theme_bw() + ylab("Interval length") + ylim(ylim_info) +
          scale_x_continuous(breaks = c(-1, 1 + 0:2 * 2.5),
                             labels = c("CI", "Wilson", "Wald", "AC")) +
          scale_color_manual(values = c("adj_post" = "black",
                                        # "adj_obs" = "grey60",
                                        "normal" = "grey80",
                                        "dir" = "grey40"),
                             guide = "none") +
          theme(axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                axis.title.x = element_text(size = 30),
                axis.title.y = element_text(size = 30)))
  dev.off()
}

make_coverage_multiple_intervals_length_comp_only_adj_ci_box_plot <- function(
    prior_1_file_list, prior_2_file_list,
    credible_draws_prior_1, credible_draws_prior_2,
    result_plot_file, prop = T, obs = T, ylim_info = c(0, 1),
    bootstrap = F, zip_code_list = NULL) {

  if (bootstrap) {
  adj_ci_post_credible_coverage <-
    get_bootstrap_ci_coverage_length_info(prior_1_file_list, ncol(credible_draws_prior_1), prop, obs, zip_code_list)
  } else {
  adj_ci_post_credible_coverage <-
    get_ci_coverage_length_info(prior_1_file_list, ncol(credible_draws_prior_1), prop, obs)
  }

  if (bootstrap) {
  adj_ci_post_credible_coverage_2 <-
    get_bootstrap_ci_coverage_length_info(prior_2_file_list, ncol(credible_draws_prior_1), prop, obs, zip_code_list)
  } else {
  adj_ci_post_credible_coverage_2 <-
    get_ci_coverage_length_info(prior_2_file_list, ncol(credible_draws_prior_1), prop, obs)
  }

  plot_df <- rbind(melt(adj_ci_post_credible_coverage),
                   melt(adj_ci_post_credible_coverage_2))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws_prior_1),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws_prior_1, 
						    na.rm = T))
  plot_df <- rbind(plot_df, direct_draw_info)
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws_prior_2),
                                 "Var2" = "Direct 2",
                                 "value" = colMeans(credible_draws_prior_2,
						    na.rm = T))
  plot_df <- rbind(plot_df, direct_draw_info)

  # direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
  #                                "Var2" = "Direct",
  #                                "value" = colMeans(credible_draws))
  # plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(adj_ci_post_credible_coverage))),
                     rep("adj_post", prod(dim(adj_ci_post_credible_coverage_2))),
                     rep("dir", 2 * ncol(credible_draws_prior_1)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(0:2 * 2.5 + 0.5, each = ncol(credible_draws_prior_1)),
                        rep(0:2 * 2.5 + 1.5, each = ncol(credible_draws_prior_1)),
                        rep(-2, ncol(credible_draws_prior_1)),
                        rep(-1, ncol(credible_draws_prior_2)))

  # print(plot_df)

  pdf(result_plot_file, height = 9, width = 15)
  print(ggplot(data = plot_df) +
          # geom_hline(yintercept = 0.95, size = 1.25, linetype = "longdash") +
          geom_boxplot(aes(y = value, x = position,
                           group = label, color = group),
                       size = 2) + xlab("") +
          theme_bw() + ylab("Interval length") + ylim(ylim_info) +
          scale_x_continuous(breaks = c(-1.5, 1 + 0:2 * 2.5),
                             labels = c("CI", "Wilson", "Wald", "AC")) +
          # values=c("black", "grey60", "grey40")) +
          # scale_color_manual(values = c("adj" = "black",
          scale_color_manual(values = c("adj_post" = "grey60",
                                        # "adj_obs" = "grey60",
                                        "normal" = "black",
                                        "dir" = "grey40"),
                             guide = "none") +
          theme(axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                axis.title.x = element_text(size = 30),
                axis.title.y = element_text(size = 30)))
  dev.off()
}


plot_coverage_info_shape <- function(
    plot_df, result_plot_file, y_label, ylim_info, nominal_line, interval_type = NULL) {

  if (!is.null(interval_type)) {
    plot_df <- plot_df[plot_df$Var2 == interval_type,]
  }
  g <- ggplot(data = plot_df) +
    # geom_hline(yintercept = 0.95, size = 1.25, linetype = "longdash") +
    geom_point(aes(y = value, x = position,
                    shape = Var2, color = group),
                 size = 8) + xlab("") +
    theme_bw() + ylab(y_label) + ylim(ylim_info) +
    scale_shape_manual(values = c("Wilson" = 16,
                                  # "adj_obs" = "grey60",
                                  "Wald" = 6,
                                  "AC" = 7,
                                  "Direct" = 10),
                       guide = "none") +
    # scale_x_continuous(breaks = c(-1, 1 + 0:2 * 2.5),
    #                    labels = c("CI", "Wilson", "Wald", "AC")) +
    # values=c("black", "grey60", "grey40")) +
    # scale_color_manual(values = c("adj" = "black",
    scale_color_manual(values = c("adj_post" = "black",
                                  # "adj_obs" = "grey60",
                                  "normal" = "grey80",
                                  "dir" = "grey40"),
                       guide = "none") +
    theme(axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          axis.title.x = element_text(size = 40),
          axis.title.y = element_text(size = 40))
  if (nominal_line) {
    g <- g + geom_hline(yintercept = 0.95, size = 1.25, linetype = "longdash")
  }
  pdf(result_plot_file, height = 9, width = 15)
  print(g)
  dev.off()

}

make_coverage_credible_intervals_multiple_intervals_info_shape <- function(
    file_list, credible_draws, norm_coverage_m, norm_prop_coverage_m, norm_ac_coverage_m,
    prop_vals, result_plot_file, prop = T, obs = T, bootstrap = F, zip_code_list = NULL,
    chunk_size = 1, exp_num = length(file_list), interval_type = NULL) {

  ci_credible_coverage <- cbind(colMeans(norm_coverage_m, na.rm = T),
                                colMeans(norm_prop_coverage_m, na.rm = T),
                                colMeans(norm_ac_coverage_m, na.rm = T))
  colnames(ci_credible_coverage) <- c("Wilson", "Wald", "AC")

  if (bootstrap) {
    adj_ci_post_credible_coverage <- get_bootstrap_ci_coverage_info(
      file_list, ncol(norm_prop_coverage_m), prop, obs, zip_code_list)
  } else {
    adj_ci_post_credible_coverage <- get_ci_coverage_info(
      file_list, ncol(norm_prop_coverage_m), prop, obs, chunk_size, exp_num)
  }

  plot_df <- rbind(melt(ci_credible_coverage),
                   melt(adj_ci_post_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws))
  plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))),
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(credible_draws)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(prop_vals, 7))
  # plot_df$interval_type <- c(rep(c("Wald"), each = length(prop_vals)),
  #                            rep(1:3, each = length(prop_vals)),
  #                            rep(4, length(prop_vals)))

  plot_coverage_info_shape(
    plot_df, result_plot_file, "Coverage prob", c(0, 1), T, interval_type)

  # print(plot_df)
}

make_coverage_credible_intervals_lengths_multiple_intervals_info_shape <- function(
    file_list, credible_draws, norm_coverage_m, norm_prop_coverage_m, norm_ac_coverage_m,
    prop_vals, result_plot_file, prop = T, obs = T, ylim_info = c(0, 1),
    bootstrap = F, zip_code_list = NULL,
    chunk_size = 1, exp_num = length(file_list), interval_type = NULL) {

  ci_credible_coverage <- cbind(colMeans(norm_coverage_m, na.rm = T),
                                colMeans(norm_prop_coverage_m, na.rm = T),
                                colMeans(norm_ac_coverage_m, na.rm = T))
  colnames(ci_credible_coverage) <- c("Wilson", "Wald", "AC")

  if (bootstrap) {
    adj_ci_post_credible_coverage <-
      get_bootstrap_ci_coverage_length_info(
        file_list, ncol(norm_prop_coverage_m), prop, obs, zip_code_list)
  } else {
    adj_ci_post_credible_coverage <-
      get_ci_coverage_length_info(
        file_list, ncol(norm_prop_coverage_m), prop, obs, chunk_size, exp_num)
  }

  plot_df <- rbind(melt(ci_credible_coverage),
                   melt(adj_ci_post_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws))
  plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))),
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(credible_draws)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(prop_vals, 7))
  # plot_df$interval_type <- c(rep(c("Wald"), each = length(prop_vals)),
  #                            rep(1:3, each = length(prop_vals)),
  #                            rep(4, length(prop_vals)))

  plot_coverage_info_shape(
    plot_df, result_plot_file, "Interval lengths", ylim_info, F, interval_type)
}

make_bias_credible_intervals_multiple_intervals_adj_ci_info_shape <- function(
    file_list, credible_draws, norm_coverage_m, norm_prop_coverage_m, norm_ac_coverage_m,
    prop_vals, result_plot_file, prop = T, obs = T, abs = F, ylim_info = c(0, 1),
    bootstrap = F, zip_code_list = NULL, chunk_size = 1, exp_num = length(file_list), interval_type = NULL) {

  ci_credible_coverage <- cbind(colMeans(norm_coverage_m, na.rm = T),
                                colMeans(norm_prop_coverage_m, na.rm = T),
                                colMeans(norm_ac_coverage_m, na.rm = T))
  colnames(ci_credible_coverage) <- c("Wilson", "Wald", "AC")

  if (bootstrap) {
    adj_ci_post_credible_coverage <-
      get_bootstrap_ci_bias_info(file_list,  ncol(norm_prop_coverage_m), prop, obs, zip_code_list, abs)

  } else {
    adj_ci_post_credible_coverage <-
      get_ci_bias_info(file_list,  ncol(norm_prop_coverage_m), prop, obs, abs,
                       chunk_size, exp_num)
  }

  # colnames(adj_ci_post_credible_coverage) <- c("Wilson", "Wald", "AC")

  plot_df <- rbind(melt(ci_credible_coverage),
                   melt(adj_ci_post_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws))
  plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))),
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(credible_draws)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(prop_vals, 7))

  if (!abs) {
    plot_df$value <- sqrt(plot_df$value)
    y_label <- c("Root mean squared error")
  } else {
    y_label <- c("Mean absolute bias")
  }

  # print(plot_df)

  plot_coverage_info_shape(
    plot_df, result_plot_file, y_label, ylim_info, F, interval_type)
}

plot_coverage_info_shape <- function(
    plot_df, result_plot_file, y_label, ylim_info, nominal_line, interval_type = NULL) {

  if (!is.null(interval_type)) {
    plot_df <- plot_df[plot_df$Var2 == interval_type,]
  }
  g <- ggplot(data = plot_df) +
    # geom_hline(yintercept = 0.95, size = 1.25, linetype = "longdash") +
    geom_point(aes(y = value, x = position,
                    shape = Var2, color = group),
                 size = 8) + xlab("") +
    theme_bw() + ylab(y_label) + ylim(ylim_info) +
    scale_shape_manual(values = c("Wilson" = 16,
                                  # "adj_obs" = "grey60",
                                  "Wald" = 6,
                                  "AC" = 7,
                                  "Direct" = 10),
                       guide = "none") +
    # scale_x_continuous(breaks = c(-1, 1 + 0:2 * 2.5),
    #                    labels = c("CI", "Wilson", "Wald", "AC")) +
    # values=c("black", "grey60", "grey40")) +
    # scale_color_manual(values = c("adj" = "black",
    scale_color_manual(values = c("adj_post" = "black",
                                  # "adj_obs" = "grey60",
                                  "normal" = "grey80",
                                  "dir" = "grey40"),
                       guide = "none") +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  if (nominal_line) {
    g <- g + geom_hline(yintercept = 0.95, size = 1.25, linetype = "longdash")
  }
  pdf(result_plot_file, height = 9, width = 15)
  print(g)
  dev.off()

}

make_coverage_credible_intervals_multiple_intervals_info_shape <- function(
    file_list, credible_draws, norm_coverage_m, norm_prop_coverage_m, norm_ac_coverage_m,
    prop_vals, result_plot_file, prop = T, obs = T, bootstrap = F, zip_code_list = NULL,
    chunk_size = 1, exp_num = length(file_list), interval_type = NULL) {

  ci_credible_coverage <- cbind(colMeans(norm_coverage_m, na.rm = T),
                                colMeans(norm_prop_coverage_m, na.rm = T),
                                colMeans(norm_ac_coverage_m, na.rm = T))
  colnames(ci_credible_coverage) <- c("Wilson", "Wald", "AC")

  if (bootstrap) {
    adj_ci_post_credible_coverage <- get_bootstrap_ci_coverage_info(
      file_list, ncol(norm_prop_coverage_m), prop, obs, zip_code_list)
  } else {
    adj_ci_post_credible_coverage <- get_ci_coverage_info(
      file_list, ncol(norm_prop_coverage_m), prop, obs, chunk_size, exp_num)
  }

  plot_df <- rbind(melt(ci_credible_coverage),
                   melt(adj_ci_post_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws))
  plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))),
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(credible_draws)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(prop_vals, 7))
  # plot_df$interval_type <- c(rep(c("Wald"), each = length(prop_vals)),
  #                            rep(1:3, each = length(prop_vals)),
  #                            rep(4, length(prop_vals)))

  plot_coverage_info_shape(
    plot_df, result_plot_file, "Coverage prob", c(0, 1), T, interval_type)

  # print(plot_df)
}

make_coverage_credible_intervals_lengths_multiple_intervals_info_shape <- function(
    file_list, credible_draws, norm_coverage_m, norm_prop_coverage_m, norm_ac_coverage_m,
    prop_vals, result_plot_file, prop = T, obs = T, ylim_info = c(0, 1),
    bootstrap = F, zip_code_list = NULL,
    chunk_size = 1, exp_num = length(file_list), interval_type = NULL) {

  ci_credible_coverage <- cbind(colMeans(norm_coverage_m, na.rm = T),
                                colMeans(norm_prop_coverage_m, na.rm = T),
                                colMeans(norm_ac_coverage_m, na.rm = T))
  colnames(ci_credible_coverage) <- c("Wilson", "Wald", "AC")

  if (bootstrap) {
    adj_ci_post_credible_coverage <-
      get_bootstrap_ci_coverage_length_info(
        file_list, ncol(norm_prop_coverage_m), prop, obs, zip_code_list)
  } else {
    adj_ci_post_credible_coverage <-
      get_ci_coverage_length_info(
        file_list, ncol(norm_prop_coverage_m), prop, obs, chunk_size, exp_num)
  }

  plot_df <- rbind(melt(ci_credible_coverage),
                   melt(adj_ci_post_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws))
  plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))),
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(credible_draws)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(prop_vals, 7))
  # plot_df$interval_type <- c(rep(c("Wald"), each = length(prop_vals)),
  #                            rep(1:3, each = length(prop_vals)),
  #                            rep(4, length(prop_vals)))

  plot_coverage_info_shape(
    plot_df, result_plot_file, "Interval lengths", ylim_info, F, interval_type)
}

make_bias_credible_intervals_multiple_intervals_adj_ci_info_shape <- function(
    file_list, credible_draws, norm_coverage_m, norm_prop_coverage_m, norm_ac_coverage_m,
    prop_vals, result_plot_file, prop = T, obs = T, abs = F, ylim_info = c(0, 1),
    bootstrap = F, zip_code_list = NULL, chunk_size = 1, exp_num = length(file_list), interval_type = NULL) {

  ci_credible_coverage <- cbind(colMeans(norm_coverage_m, na.rm = T),
                                colMeans(norm_prop_coverage_m, na.rm = T),
                                colMeans(norm_ac_coverage_m, na.rm = T))
  colnames(ci_credible_coverage) <- c("Wilson", "Wald", "AC")

  if (bootstrap) {
    adj_ci_post_credible_coverage <-
      get_bootstrap_ci_bias_info(file_list,  ncol(norm_prop_coverage_m), prop, obs, zip_code_list, abs)

  } else {
    adj_ci_post_credible_coverage <-
      get_ci_bias_info(file_list,  ncol(norm_prop_coverage_m), prop, obs, abs,
                       chunk_size, exp_num)
  }

  # colnames(adj_ci_post_credible_coverage) <- c("Wilson", "Wald", "AC")

  plot_df <- rbind(melt(ci_credible_coverage),
                   melt(adj_ci_post_credible_coverage))
  direct_draw_info <- data.frame("Var1" = 1:ncol(credible_draws),
                                 "Var2" = "Direct",
                                 "value" = colMeans(credible_draws))
  plot_df <- rbind(plot_df, direct_draw_info)

  plot_df$group <- c(rep("normal", prod(dim(ci_credible_coverage))),
                     rep("adj_post", prod(dim(ci_credible_coverage))),
                     rep("dir", ncol(credible_draws)))
  plot_df <-
    plot_df %>% unite(label, Var2, group, remove = F)
  plot_df$position <- c(rep(prop_vals, 7))

  if (!abs) {
    plot_df$value <- sqrt(plot_df$value)
    y_label <- c("Root mean squared error")
  } else {
    y_label <- c("Mean absolute bias")
  }

  # print(plot_df)

  plot_coverage_info_shape(
    plot_df, result_plot_file, y_label, ylim_info, F, interval_type)
}

