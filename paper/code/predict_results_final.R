

library(this.path)

setwd(paste0(this.path::here(), '/../../'))
load('paper/results/cluster_gridjanuary.RData')
library(fields)
cluster_mat_full_reduced <- tail(cluster_mat_full, 35280)
library(ggplot2)

seasons <- c('january', 'april', 'july', 'october')

source('paper/code/src/plot_src.R')
source('paper/code/src/pred_at_pressure_levels_src.R')


library(ggplot2)
library(dplyr)

grid_use <- grid_locs
cluster_probs <- array(0, dim = c(360, 98, 5, 4), 
                       dimnames = list(grid_locs[1:360,2], unique(grid_locs[,1]), 1:5, 
                                       c('January', 'April', 'July', 'October')))
for (q in 1:length(seasons)) {
  load(paste0('paper/results/cluster_grid', seasons[q], '.RData'))
  clust_grid <-  tail(cluster_mat_full, 35280)
  cluster_probability <- t(apply(clust_grid, 1, function(x) tabulate(x, nbins = 5),
                                 simplify = T)/ncol(clust_grid))
  
  cluster_probs[,,,q] <- cluster_probability
  
  cluster_assigned <- apply(cluster_mat_full_reduced, 1, function(x) which.max(tabulate(x)))
  cluster_probability <- apply(cluster_mat_full_reduced, 
                               1, function(x) tabulate(x)[which.max(tabulate(x))])/
    ncol(cluster_mat_full_reduced)
}

prob_df <- reshape2::melt(cluster_probs)
colnames(prob_df) <- c('longitude', 'latitude', 'cluster', 'month', 'value')

prob_df_sum <- prob_df %>%
  group_by(longitude, latitude, month) %>%
  arrange(cluster) %>%
  summarize(cluster = which.max(value),
            probability =value[which.max(value)])
prob_df_sum <- merge(prob_df_sum, grid_use, by = c('longitude', 'latitude'))

load('paper/results/oxy_none__5_18_18_30_30-matern_spheretime_warp_seasonality/final_results.RData')

library(fda)
library(Matrix)
parameters <- estimated_model[['parameters']]
Phi_pred <- bdiag(eval.basis(params[['basis_predictors']][[1]], evalarg = 1000),
                  eval.basis(params[['basis_predictors']][[1]], evalarg = 1000))
t_vals <- sapply(1:params[['G']], function(g){ (Phi_pred %*% parameters[['means_predictors']][[g]])[1,1]})
new_order <- order(order(t_vals))


prob_df_sum <- mutate(prob_df_sum, cluster_use = new_order[cluster])



theme_set(theme_bw() +theme(legend.position = 'bottom', 
                            text = element_text(family = 'Arial')))


a <- ggplot()+
  geom_tile(data = prob_df_sum,
            aes(x = ifelse(longitude > 180, longitude - 360, longitude), y = latitude,
                fill = factor(cluster_use),
                colour = NULL,
                height = height + .035, width = width + .035, 
                alpha = probability)) +
  SO_coord + SO_theme +
  fronts_dark + continents + latitude_lines + longitude_lines + 
  facet_wrap(~month, ncol = 2, nrow = 2)+
  labs(alpha = 'Probability', fill = 'Cluster')+
  theme(legend.position = 'bottom')
ggsave(filename = paste0('paper/images/product_clusters_season.png'), 
       plot = a, width = 7, 
       height = 7.4)

a <- ggplot()+
  geom_tile(data = prob_df_sum,
            aes(x = ifelse(longitude > 180, longitude - 360, longitude), y = latitude,
                fill = factor(cluster_use),
                height = height + .035, width = width + .035, 
                alpha = probability)) +
  SO_coord + SO_theme +
  fronts_dark + continents + latitude_lines + longitude_lines + 
  facet_wrap(~month, ncol = 4, nrow = 1)+
  labs(alpha = 'Probability', fill = 'Cluster')+
  theme(legend.position = 'bottom')
ggsave(filename = paste0('paper/images/product_clusters_season_row.png'), 
       plot = a, width = 9, 
       height = 3.1)

a <- ggplot()+
  geom_tile(data = prob_df_sum,
            aes(x = ifelse(longitude > 180, longitude - 360, longitude), y = latitude,
                fill = factor(cluster_use),
                height = height + .04, width = width + .04, 
                alpha = probability)) +
  SO_coord + SO_theme +
  fronts_dark + continents + latitude_lines + longitude_lines + 
  facet_wrap(~month, ncol = 4, nrow = 1)+
  labs(alpha = 'Probability', fill = 'Cluster')+
  theme(legend.position = 'bottom')
ggsave(filename = paste0('paper/images/product_clusters_season_row.png'), 
       plot = a, width = 9, 
       height = 3.1)

#####

# oxygen predictions
return_values <- list()
return_variances <- list()
for (q in 1:length(seasons)) {
  load(paste0('paper/results/1.5_one_chunk_oxy_grid_',
              seasons[q], '_1.RData'))
  pcs_array_response <- pred_score_results$pcs_array_response
  pcs_array_predictors <- pred_score_results$pcs_array_predictors
  wt <- pred_score_results$wt 
  
  params <- pred_set_up$prediction_model$params
  parameters <- pred_set_up$prediction_model$parameters
  pressures <- c(50, 150, 300, 500, 1000, 1500)
  Phi <- eval.basis(params[['basis_response']], pressures)
  Phi_pred <- bdiag(eval.basis(params[['basis_predictors']][[1]], evalarg= pressures),
                    eval.basis(params[['basis_predictors']][[2]], evalarg= pressures))
  cluster_mat_pred_use <- tail(cluster_mat_full, 
                               dim(pcs_array_response)[1])[, seq(10, 100, by = 10)]
  cluster_mat_pred_use <- tail(cluster_mat_full, 
                               dim(pcs_array_response)[1])[, round(1:40 * 2.5)]
  t_vals <- sapply(1:params[['G']], function(g){ (Phi_pred %*% parameters[['means_predictors']][[g]])[1,1]})
  new_order <- order(order(t_vals))
  
  Omega_responseg <- lapply(1:params[['G']], function(g) {
    cbind(Phi %*% parameters[['Omega_response1']][[g]],
          Phi %*% parameters[['Omega_response2']][[g]])
  })
  
  Omega_predictorsg <- lapply(1:params[['G']], function(g) {
    Phi_pred %*% parameters[['Omega_predictors']][[g]]
  })
  
  var_values <- pred_score_results[['variances']]
  n_pcs_predictors <- pred_set_up[['prediction_model']][['params']][['pc_predictors']]
  day_use <- tail(pred_set_up$prediction_model$data_inputs$day, dim(pcs_array_response)[1])
  wt_use <- wt
  
  covariates <- diag(c(1,
                  sin(day_use[1]/365.25 * 2 * pi * 1),
                  cos(day_use[1]/365.25 * 2 * pi * 1),
                  sin(day_use[1]/365.25 * 2 * pi * 2),
                  cos(day_use[1]/365.25 * 2 * pi * 2),
                  sin(day_use[1]/365.25 * 2 * pi * 3),
                  cos(day_use[1]/365.25 * 2 * pi * 3)))
  
  values <- array(dim = c(nrow(pcs_array_predictors), length(wt), length(pressures)))
  for (r in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:length(wt)) { 
      z_i_tau <- cluster_mat_pred_use[r,mc]
      values[r,mc, ] <- 1/ncol(cluster_mat_pred_use) * 
        as.double(Phi %*% (rowSums(parameters[['means_response']][[z_i_tau]] %*% covariates) + 
                             parameters[['Omega_response1']][[z_i_tau]] %*% pcs_array_predictors[r,mc,] + 
                             parameters[['Omega_response2']][[z_i_tau]] %*% pcs_array_response[r,mc,]))
    }
  }
  avg_values <- apply(values,c(1,3), function(x) sum(x))
  library(ncdf4)
  grid_file <- nc_open('paper/data/grid.nc', write = F)
  grid_size <- 1.5
  grid_long <- ncvar_get(grid_file, 'XC')[,1]
  grid_lat <- ncvar_get(grid_file, 'YC')[1,]
  grid_depth <- ncvar_get(grid_file, 'Depth')[seq(1, length(grid_long), by = round(6*grid_size)),
                                              seq(1, length(grid_lat), by = round(6*grid_size))]
  grid_long <- grid_long[seq(1, length(grid_long), by = round(6*grid_size))]
  grid_lat <- grid_lat[seq(1, length(grid_lat), by = round(6*grid_size))]
  
  height <- grid_lat[2:length(grid_lat)] - grid_lat[1:(length(grid_lat) - 1)]
  height <- data.frame("latitude" = grid_lat, height = c(height[1], height))
  
  grid_df <- cbind(expand.grid('longitude' = grid_long, 'latitude' = grid_lat),
                   depth = as.double(grid_depth))
  grid_df[['width']] <- grid_size
  grid_df <- merge(grid_df, height, by = "latitude")
  grid_use <- grid_df
  
  var_values <- pred_score_results[['variances']]
  
  pressure_var_values <- array(dim = c(nrow(pcs_array_predictors),
                                       length(wt), length(pressures)))
  for (r in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:length(wt)) { 
      z_i_tau <- cluster_mat_pred_use[r,mc]
      pressure_var_values[r,mc, ] <- 
        diag(Phi %*% cbind(parameters[['Omega_response1']][[z_i_tau]],
                           parameters[['Omega_response2']][[z_i_tau]]) %*% 
               var_values[[mc]][[r]] %*% 
               t(cbind(parameters[['Omega_response1']][[z_i_tau]],
                       parameters[['Omega_response2']][[z_i_tau]])) %*%
               t(Phi))
    }
  }
  avg_var_values <- apply(pressure_var_values,c(1,3), function(x) mean(x))
  var_avg_values <- apply(values,c(1,3), function(x) var(x))

  return_values[[q]] <- avg_values
  return_variances[[q]] <- list(avg_var_values, var_avg_values)
  gc()
  print(q)
}

oxy_df <- data.frame(grid_use, 
                     dplyr::bind_rows(lapply(1:length(seasons), 
                                             function(x) data.frame(value = return_values[[x]][,3],
                                                                    month = paste(c('January', 'April', 'July', 'October')[x], 2020) ))))
oxy_df$month <- factor(oxy_df$month,
                       levels = c('January 2020', 'April 2020', 'July 2020', 'October 2020'))

a <- ggplot()+
  geom_tile(data = oxy_df %>% filter(latitude > -75),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) + 
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  coord_map('ortho', orientation = c(-90, 0, 0), ylim = c(-90, -30)) +
  SO_theme +
  fronts_light + continents +
  latitude_lines + longitude_lines +
  facet_wrap(~month, nrow = 2, ncol = 2) + 
  labs(color ='Oxygen Prediction\n(μmol/kg)', fill = 'Oxygen Prediction\n(μmol/kg)')+
  theme(panel.grid = element_blank(), legend.key.width = unit(1, "cm"),
        text = element_text(family = 'Arial', size = 16))
ggsave(filename = 'paper/images/product_oxy_pred.png', 
       plot = a, width = 6, 
       height = 7, dpi = 500)

oxy_df <- data.frame(grid_use, 
                     dplyr::bind_rows(lapply(1:length(seasons), 
                                             function(x) data.frame(value = sqrt(return_variances[[x]][[1]][,3] + return_variances[[x]][[2]][,3]),
                                                                    month = paste(c('January', 'April', 'July', 'October')[x], 2020) ))))
oxy_df$month <- factor(oxy_df$month,
                       levels = c('January 2020', 'April 2020', 'July 2020', 'October 2020'))

a <- ggplot()+
  geom_tile(data = oxy_df%>% filter(latitude > -75),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) +
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  coord_map('ortho', orientation = c(-90, 0, 0), ylim = c(-90, -30)) +
  SO_theme +
  fronts_light + continents +
  latitude_lines + longitude_lines +
  facet_wrap(~month, nrow = 2, ncol = 2) +
  labs(color ='Oxygen SD\n(μmol/kg)', fill = 'Oxygen SD\n(μmol/kg)')+
  theme(panel.grid = element_blank(), legend.key.width = unit(1, "cm"),
        text = element_text(family = 'Arial', size = 16))
ggsave(filename = 'paper/images/product_oxy_var.png', 
       plot = a, width = 6, 
       height = 7, dpi = 500)

oxy_df <- data.frame(grid_use, value = sqrt(return_variances[[2]][[1]][,3]))
a <- ggplot()+
  geom_tile(data = oxy_df%>% filter(latitude > -75),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) +
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  coord_map('ortho', orientation = c(-90, 0, 0), ylim = c(-90, -30)) +
  SO_theme +
  fronts_light + continents +
  latitude_lines + longitude_lines +
  labs(color ='Oxygen SD\n(μmol/kg)', fill = 'Oxygen SD\n(μmol/kg)')+
  theme(panel.grid = element_blank(), legend.key.width = unit(1, "cm"),
        text = element_text(family = 'Arial', size = 16))
ggsave(filename = 'paper/images/product_oxy_var_p1.png', 
       plot = a, width = 4.5, height = 5.11)

oxy_df <- data.frame(grid_use, value = sqrt(return_variances[[2]][[2]][,3]))
a <- ggplot()+
  geom_tile(data = oxy_df%>% filter(latitude > -75),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) +
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  coord_map('ortho', orientation = c(-90, 0, 0), ylim = c(-90, -30)) +
  SO_theme +
  fronts_light + continents +
  latitude_lines + longitude_lines +
  labs(color ='Oxygen SD\n(μmol/kg)', fill = 'Oxygen SD\n(μmol/kg)')+
  theme(panel.grid = element_blank(), legend.key.width = unit(1, "cm"),
        text = element_text(family = 'Arial', size = 16))
ggsave(filename = 'paper/images/product_oxy_var_p2.png', 
       plot = a, width = 4.5, height = 5.11)

MC_max_pred <- 40
uncon_values <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))

covariates <- c(1,
                sin(105/365.25 * 2 * pi * 1),
                cos(105/365.25 * 2 * pi * 1),
                sin(105/365.25 * 2 * pi * 2),
                cos(105/365.25 * 2 * pi * 2),
                sin(105/365.25 * 2 * pi * 3),
                cos(105/365.25 * 2 * pi * 3))
for (i in 1:nrow(pcs_array_predictors)) {
  for (mc in 1:MC_max_pred) { 
    z_i_tau <- cluster_mat_pred_use[i,mc]
    uncon_values[i,mc, ] <- 
      as.double(Phi %*% (rowSums(parameters[['means_response']][[z_i_tau]] %*% diag(covariates)) ))
  }
}
uncon_avg_values <- apply(uncon_values,c(1,3), function(x) mean(x))

# first part of variance
uncon_vars_response <- lapply(1:length(Omega_responseg), function(g){
  diag((Omega_responseg[[g]] %*% diag(c(parameters[['variances_predictors']][[g]], 
                                        parameters[['variances_response']][[g]]))) %*%
         t(Omega_responseg[[g]]))
})
uncon_var_values1 <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))
for (i in 1:nrow(pcs_array_predictors)) {
  for (mc in 1:MC_max_pred) { 
    z_i_tau <- cluster_mat_pred_use[i,mc]
    uncon_var_values1[i, mc, ] <- 
      uncon_vars_response[[z_i_tau]]
  }
}
uncon_avg_var_values <- apply(uncon_var_values1, c(1,3), function(x) mean(x))
# second part of variance
uncon_var_avg_values <- apply(uncon_values, c(1,3), function(x) mean((x  -
                                                                        mean(x))^2))

#for (p in 1:length(pressures)) {
p = 3 
a <- ggplot(data = data.frame(grid_use, score = 1- (avg_var_values[,p] + var_avg_values[,p])/
                                (uncon_avg_var_values[,p] + uncon_var_avg_values[,p])) %>%
              mutate(score = ifelse(score < 0, 0, score), 
                     score = ifelse(score > 1, 1, score))%>% filter(latitude > -75))+
  geom_tile(aes(x = longitude, y = latitude, color = score,
                fill = score,
                height = height, width = width)) +
  scale_color_viridis_c(limits = c(0, 1)) + scale_fill_viridis_c(limits = c(0, 1)) +   SO_coord + SO_theme +
  fronts_light + continents + latitude_lines + longitude_lines + 
  labs(fill = substitute(paste(nn, ~R^2), list(nn = 'Oxygen'[1])), 
       color = substitute(paste(nn, ~R^2), list(nn = 'Oxygen'[1]))) + 
  theme(text = element_text(size = 16, family = 'Arial'),
        legend.key.width = unit(1, "cm"))
ggsave(filename = 'paper/images/r2.png', a,
       width = 4.5, height = 5.11)

# temperature and salinity

return_values <- list()
return_variances <- list()
for (q in 1:length(seasons)) {
  load(paste0('paper/results/1.5_one_chunk_oxy_grid_',
              seasons[q], '_1.RData'))
  pcs_array_response <- pred_score_results$pcs_array_response
  pcs_array_predictors <- pred_score_results$pcs_array_predictors
  wt <- pred_score_results$wt 
  
  params <- pred_set_up$prediction_model$params
  parameters <- pred_set_up$prediction_model$parameters
  cluster_mat_pred_use <- tail(cluster_mat_full, 
                               dim(pcs_array_response)[1])[, round(1:40 * 2.5)]

  var_values <- pred_score_results[['variances']]
  day_use <- tail(pred_set_up$prediction_model$data_inputs$day, dim(pcs_array_response)[1])

  covariates <- diag(c(1,
                       sin(day_use[1]/365.25 * 2 * pi * 1),
                       cos(day_use[1]/365.25 * 2 * pi * 1),
                       sin(day_use[1]/365.25 * 2 * pi * 2),
                       cos(day_use[1]/365.25 * 2 * pi * 2),
                       sin(day_use[1]/365.25 * 2 * pi * 3),
                       cos(day_use[1]/365.25 * 2 * pi * 3)))
  
  values <- array(dim = c(nrow(pcs_array_predictors), length(wt), 2*length(pressures)))
  for (r in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:length(wt)) { 
      z_i_tau <- cluster_mat_pred_use[r,mc]
      values[r,mc, ] <- 1/ncol(cluster_mat_pred_use) * 
        as.double(Phi_pred %*% (rowSums(parameters[['means_predictors']][[z_i_tau]] %*% diag(covariates)) + 
                             parameters[['Omega_predictors']][[z_i_tau]] %*% pcs_array_predictors[r,mc,]))
    }
  }
  avg_values <- apply(values,c(1,3), function(x) sum(x))
  
  var_values <- pred_score_results[['variances']]
  
  pressure_var_values <- array(dim = c(nrow(pcs_array_predictors),
                                       length(wt), 2*length(pressures)))
  for (r in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:length(wt)) { 
      z_i_tau <- cluster_mat_pred_use[r,mc]
      pressure_var_values[r,mc, ] <- 
        diag(Phi_pred %*% parameters[['Omega_predictors']][[z_i_tau]] %*% 
               var_values[[mc]][[r]][1:18, 1:18] %*% 
               t(parameters[['Omega_predictors']][[z_i_tau]]) %*%
               t(Phi_pred))
    }
  }
  avg_var_values <- apply(pressure_var_values,c(1,3), function(x) mean(x))
  var_avg_values <- apply(values,c(1,3), function(x) var(x))
  
  return_values[[q]] <- avg_values
  return_variances[[q]] <- list(avg_var_values, var_avg_values)
  gc()
  print(q)
}


t_df <- data.frame(grid_use, 
                     dplyr::bind_rows(lapply(1:length(seasons), 
                                             function(x) data.frame(value = return_values[[x]][,3],
                                                                    month = paste(c('January', 'April', 'July', 'October')[x], 2020) ))))
t_df$month <- factor(t_df$month,
                       levels = c('January 2020', 'April 2020', 'July 2020', 'October 2020'))

a <- ggplot()+
  geom_tile(data = t_df%>% filter(latitude > -75),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) + 
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  coord_map('ortho', orientation = c(-90, 0, 0), ylim = c(-90, -30)) +
  SO_theme +
  fronts_light + continents +
  latitude_lines + longitude_lines +
  facet_wrap(~month, nrow = 2, ncol = 2) + 
  labs(color ='Temperature Prediction\n(°C)', fill = 'Temperature Prediction\n(°C)')+
  theme(panel.grid = element_blank(), legend.key.width = unit(1, "cm"),
        text = element_text(family = 'Arial', size = 16))
ggsave(filename = 'paper/images/product_temp_pred.png', 
       plot = a, width = 6, 
       height = 7, dpi = 500)

s_df <- data.frame(grid_use, 
                   dplyr::bind_rows(lapply(1:length(seasons), 
                                           function(x) data.frame(value = return_values[[x]][,9],
                                                                  month = paste(c('January', 'April', 'July', 'October')[x], 2020) ))))
s_df$month <- factor(s_df$month,
                     levels = c('January 2020', 'April 2020', 'July 2020', 'October 2020'))

a <- ggplot()+
  geom_tile(data = s_df%>% filter(latitude > -75),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) + 
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  coord_map('ortho', orientation = c(-90, 0, 0), ylim = c(-90, -30)) +
  SO_theme +
  fronts_light + continents +
  latitude_lines + longitude_lines +
  facet_wrap(~month, nrow = 2, ncol = 2) + 
  labs(color ='Salinity Prediction\n(PSU)', fill = 'Salinity Prediction\n(PSU)')+
  theme(panel.grid = element_blank(), legend.key.width = unit(1, "cm"),
        text = element_text(family = 'Arial', size = 16))
ggsave(filename = 'paper/images/product_psal_pred.png', 
       plot = a, width = 6, 
       height = 7, dpi = 500)



t_df <- data.frame(grid_use, 
                   dplyr::bind_rows(lapply(1:length(seasons), 
                                           function(x) data.frame(value = sqrt(return_variances[[x]][[1]][,3] + 
                                                                                 return_variances[[x]][[2]][,3]),
                                                                  month = paste(c('January', 'April', 'July', 'October')[x], 2020) ))))
t_df$month <- factor(t_df$month,
                     levels = c('January 2020', 'April 2020', 'July 2020', 'October 2020'))

a <- ggplot()+
  geom_tile(data = t_df%>% filter(latitude > -75),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) + 
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  coord_map('ortho', orientation = c(-90, 0, 0), ylim = c(-90, -30)) +
  SO_theme +
  fronts_light + continents +
  latitude_lines + longitude_lines +
  facet_wrap(~month, nrow = 2, ncol = 2) + 
  labs(color ='Temperature SD\n(°C)', fill = 'Temperature SD\n(°C)')+
  theme(panel.grid = element_blank(), legend.key.width = unit(1, "cm"),
        text = element_text(family = 'Arial', size = 16))
ggsave(filename = 'paper/images/product_temp_var.png', 
       plot = a, width = 6, 
       height = 7, dpi = 500)

s_df <- data.frame(grid_use, 
                   dplyr::bind_rows(lapply(1:length(seasons), 
                                           function(x) data.frame(value = sqrt(return_variances[[x]][[1]][,9] + 
                                                                                 return_variances[[x]][[2]][,9]),
                                                                  month = paste(c('January', 'April', 'July', 'October')[x], 2020) ))))
s_df$month <- factor(s_df$month,
                     levels = c('January 2020', 'April 2020', 'July 2020', 'October 2020'))

a <- ggplot()+
  geom_tile(data = s_df%>% filter(latitude > -75),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) + 
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  coord_map('ortho', orientation = c(-90, 0, 0), ylim = c(-90, -30)) +
  SO_theme +
  fronts_light + continents +
  latitude_lines + longitude_lines +
  facet_wrap(~month, nrow = 2, ncol = 2) + 
  labs(color ='Salinity SD\n(PSU)', fill = 'Salinity SD\n(PSU)')+
  theme(panel.grid = element_blank(), legend.key.width = unit(1, "cm"),
        text = element_text(family = 'Arial', size = 16))
ggsave(filename = 'paper/images/product_psal_var.png', 
       plot = a, width = 6, 
       height = 7, dpi = 500)






# Section plots

############# Section plots

pressures <- c(seq(0, 2000, by = 5))
final_prediction <- list()
longitude_use <- -90.0833358764648
tol <- 1

indexes_use <- abs(ifelse( grid_use$longitude > 180, 
                           grid_use$longitude - 360, grid_use$longitude) -longitude_use) < tol

grid_use_subset <- grid_use %>%
  mutate(longitude = ifelse( longitude > 180, longitude - 360, longitude)) %>%
  filter(abs(longitude -longitude_use) < tol)

Phi <- eval.basis(params[['basis_response']], pressures)
Omega_responseg_subset <- lapply(1:params[['G']], function(g) {
  cbind(Phi %*% parameters[['Omega_response1']][[g]],
        Phi %*% parameters[['Omega_response2']][[g]])
})

return_values <- list()
for (q in 1:length(seasons)) {
  load(paste0('paper/results/1.5_one_chunk_oxy_grid_',
              seasons[q], '_1.RData'))
  pcs_array_response <- pred_score_results$pcs_array_response
  pcs_array_predictors <- pred_score_results$pcs_array_predictors
  wt <- pred_score_results$wt 
  
  params <- pred_set_up$prediction_model$params
  parameters <- pred_set_up$prediction_model$parameters
  cluster_mat_pred_use <- tail(cluster_mat_full, 
                               dim(pcs_array_response)[1])[, round(1:40 * 2.5)]
  t_vals <- sapply(1:params[['G']], function(g){ (Phi_pred %*% parameters[['means_predictors']][[g]])[1,1]})
  new_order <- order(order(t_vals))
  
  Omega_responseg <- lapply(1:params[['G']], function(g) {
    cbind(Phi %*% parameters[['Omega_response1']][[g]],
          Phi %*% parameters[['Omega_response2']][[g]])
  })
  
  Omega_predictorsg <- lapply(1:params[['G']], function(g) {
    Phi_pred %*% parameters[['Omega_predictors']][[g]]
  })
  
  var_values <- pred_score_results[['variances']]
  var_values_subset <- lapply(var_values, function(x) {x[indexes_use]})
  day_use <- tail(pred_set_up$prediction_model$data_inputs$day, nrow(pcs_array_response))
  pcs_array_predictors_subset <- pcs_array_predictors[indexes_use,,]
  pcs_array_response_subset <- pcs_array_response[indexes_use,,]
  values_reduced <- array(dim = c(nrow(pcs_array_predictors_subset), MC_max_pred, length(pressures)))
  cluster_mat_subset <- cluster_mat_pred_use[indexes_use,]
  var_values1_reduced <- array(dim = c(nrow(pcs_array_predictors_subset), MC_max_pred, length(pressures)))
  values_df <- list()
  
  covariates <- c(1,
                  sin(tail(day_use, 1)/365.25 * 2 * pi * 1),
                  cos(tail(day_use, 1)/365.25 * 2 * pi * 1),
                  sin(tail(day_use, 1)/365.25 * 2 * pi * 2),
                  cos(tail(day_use, 1)/365.25 * 2 * pi * 2),
                  sin(tail(day_use, 1)/365.25 * 2 * pi * 3),
                  cos(tail(day_use, 1)/365.25 * 2 * pi * 3))
  return_values[[q]] <- list()
  for (i in 1:nrow(pcs_array_predictors_subset)) {
    for (mc in 1:MC_max_pred) { 
      z_i_tau <- cluster_mat_subset[i,mc]
      values_reduced[i,mc, ] <- 
        as.double(Phi %*% (rowSums(parameters[['means_response']][[z_i_tau]] %*% diag(covariates)) + 
                             parameters[['Omega_response1']][[z_i_tau]] %*% pcs_array_predictors[i,mc,] + 
                             parameters[['Omega_response2']][[z_i_tau]] %*% pcs_array_response[i,mc,]))
      var_values1_reduced[i, mc, ] <- 
        rowSums((Omega_responseg_subset[[z_i_tau]] %*% var_values_subset[[mc]][[i]])*
                  Omega_responseg_subset[[z_i_tau]])
      # first part of variance
    }
    return_values[[q]][[i]] <- data.frame('oxygen' = apply(values_reduced[i,,],2, function(x) mean(x)),
                                 'v_p1' = apply(var_values1_reduced[i,,],2, function(x) mean(x)),
                                 'v_p2' = apply(values_reduced[i,,], 2, function(x) mean((x -
                                                                                            mean(x))^2)),
                                 'oxygen_var' = apply(var_values1_reduced[i,,],2, function(x) mean(x )) +
                                   apply(values_reduced[i,,], 2, function(x) mean((x -
                                                                                     mean(x ))^2)),
                                 'pressure' = pressures,
                                 'latitude' = grid_use_subset[i,'latitude'], 
                                 'month' = c('January', 'April', 'July', 'October')[q])
    
  }
  gc()
  print(q)
}


values_df <- bind_rows(return_values) %>% left_join(grid_use_subset)


latitudes_use <- unique(values_df[['latitude']])[seq(1, length(unique(values_df[['latitude']])), 3)]

ggplot(data = values_df %>%left_join(data.frame(month = month.name[c(1, 4, 7, 10)], 
                                                month_label = factor(paste(month.name[c(1, 4, 7, 10)], 2020),
                                                                     levels= c('January 2020', 'April 2020', 'July 2020',
                                                                               'October 2020')))) %>%
         filter(latitude %in% latitudes_use),
       aes(color = latitude, y = pressure, x = oxygen,
                                                                     group = latitude))+
  geom_path() + 
  scale_y_continuous(trans = 'reverse') + 
  labs(x = 'Oxygen Prediction\n(μmol/kg)', y = 'Pressure (decibars)', color = 'Latitude')+ 
  facet_wrap(~month_label, ncol = 4) + 
  
  theme(text = element_text(size = 16), 
        legend.key.width = unit(1, "cm"))+
  scale_color_viridis_c()
ggsave(filename = 'paper/images/section_oxy_lines.png',  
       width = 9, 
       height = 7)

ggplot(data = values_df %>%left_join(data.frame(month = month.name[c(1, 4, 7, 10)], 
                                                month_label = factor(paste(month.name[c(1, 4, 7, 10)], 2020),
                                                                     levels= c('January 2020', 'April 2020', 'July 2020',
                                                                               'October 2020')))) %>%
         filter(latitude %in% latitudes_use), 
       aes(color = latitude, y = pressure, x = sqrt(oxygen_var),
                                                                     group = latitude))+
  geom_path() + 
  facet_wrap(~month_label, ncol = 4) + 
  scale_y_continuous(trans = 'reverse') + 
  labs(x = 'Oxygen SD\n(μmol/kg)', y = 'Pressure (decibars)', color = 'Latitude')+ 
  theme(text = element_text(size = 16), 
        legend.key.width = unit(1, "cm"))+
  scale_color_viridis_c()
ggsave(filename = 'paper/images/section_oxy_var.png',
       width = 9, 
       height = 7)