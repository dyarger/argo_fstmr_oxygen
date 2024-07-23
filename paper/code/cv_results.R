library(this.path)

###############################################################
#### Basic Error Calculation with TS
###############################################################
setwd(paste0(this.path::here(), '/../results/'))

leave_out_types = c('profiles', 'floats')
for(leave_out_type in leave_out_types){
  error_t = c()
  error_145155 = c()
  
  nfolds = ifelse(leave_out_type == 'floats', 10, 5)
  
  for (i in c(1:nfolds)){
    load(paste0('oxy_', leave_out_type, '_', as.character(i), '_5_18_18_30_30-matern_spheretime_warp_seasonality/TS_prediction.RData'))
    
    gdf = split(pred_at_p_levels$pred_df, pred_at_p_levels$pred_df$mc)
    prediction = rep(0, nrow(gdf[[1]]))
    for(j in 1:length(gdf)){
      prediction = prediction + gdf[[j]]$pred_full * pred_score_results$wt[j]
    }
    prediction = prediction / sum(pred_score_results$wt) 
    
    error = prediction - held_out_profs$oxy
    
    held_out_profs$prediction = prediction
    check_gdf = split(held_out_profs, held_out_profs$profile_unique)
    
    error_t = c(error_t, error)
    
    # Collect pressure subset
    ind = held_out_profs$pressure < 155 & held_out_profs$pressure > 145
    error_145155 = c(error_145155, error[ind])
  }

  print(paste0('-- Errors for ', leave_out_type, ' ------------------------'))
  print(paste0("MAD total: ", as.character(median(abs(error_t)))))
  print(paste0("MAD 145155: ", as.character(median(abs(error_145155)))))
  print(paste0("RMSE total: ", as.character((sqrt(mean(error_t^2))))))
  print(paste0("RMSE 145155: ", as.character(sqrt(mean(error_145155^2)))))
}


####################################################################
#### Compute difference between Core Argo nearby vs not nearby no TS
####################################################################
rm(list = ls())

setwd(paste0(this.path::here(), '/../results/'))
leave_out_types = c('profiles', 'floats')

# Load Core data
load('../data/close_core.RData')

one_row_core_data = close_core_df %>%
  group_by(profile_unique) %>%
  slice(1)

one_row_core_data$daynew = julian(as.Date(one_row_core_data$date, format = '%Y-%m-%d'), origin = as.Date('2000-01-01')) 

for(leave_out_type in leave_out_types){
  errs = c()
  errs_145155 = c()
  errs_no_core = c()
  errs_no_core_145155 = c()
  
  nfolds = ifelse(leave_out_type == 'floats', 10, 5)
  
  for (fold in c(1:nfolds)){
    load(paste0('oxy_', leave_out_type, '_', as.character(fold), '_5_18_18_30_30-matern_spheretime_warp_seasonality/prediction.RData'))
    
    one_row_soccom_data = held_out_profs %>%
      group_by(profile_unique) %>%
      slice(1)
    
    pred_at_p_levels$pred_df$pressure = held_out_profs$pressure
    pred_at_p_levels$pred_df$oxy = held_out_profs$oxy
    gdf = split(pred_at_p_levels$pred_df, pred_at_p_levels$pred_df$mc)
    prediction = rep(0, nrow(gdf[[1]]))
    for(j in 1:length(gdf)){
      prediction = prediction + gdf[[j]]$pred_full * pred_score_results$wt[j]
    } 
    prediction = prediction / sum(pred_score_results$wt)
    final_df = gdf[[1]]
    final_df[['prediction_final']] = prediction
    
    unique_prof = unique(held_out_profs$profile_unique)
    for (i in 1:length(unique_prof)){
      
      day_soccom = julian(as.Date(one_row_soccom_data$day[i], format = '%Y-%m-%d'), origin = as.Date('2000-01-01'))
      lon_soccom = one_row_soccom_data$longitude[i]
      lat_soccom = one_row_soccom_data$latitude[i]
      
      ind = which(abs(one_row_core_data$daynew - day_soccom) < 30)
      candidates = one_row_core_data[ind,] %>%
        filter(abs(longitude - lon_soccom) <= 2) %>%
        filter(abs(latitude - lat_soccom) <= 2)
      
      prof = unique_prof[i]
      filter_df = final_df %>%
        filter(profile_unique == prof)
      pind = filter_df$pressure < 155 & filter_df$pressure > 145
      
      if(nrow(candidates) > 0){
        errs_145155 = c(errs_145155, filter_df$prediction[pind] - filter_df$oxy[pind])
        errs = c(errs, filter_df$prediction - filter_df$oxy)
      } else {
        errs_no_core_145155 = c(errs_no_core_145155, filter_df$prediction[pind] - filter_df$oxy[pind])
        errs_no_core = c(errs_no_core, filter_df$prediction - filter_df$oxy)
      }
    }
  }
  
  print(paste0('-- Errors for ', leave_out_type, ' ------------------------'))
  errs_total = c(errs, errs_no_core)
  errs_total_145155 = c(errs_145155, errs_no_core_145155)
  print(paste0('MAD total, core nearby: ', as.character(median(abs(errs)))))
  print(paste0('MAD 145155, core nearby: ', as.character(median(abs(errs_145155)))))
  print(paste0('MAD total, no core nearby: ', as.character(median(abs(errs_no_core)))))
  print(paste0('MAD 145155, no core nearby: ', as.character(median(abs(errs_no_core_145155)))))
  print(paste0('MAD total, total: ', as.character(median(abs(errs_total)))))
  print(paste0('MAD 145155, total: ', as.character(median(abs(errs_total_145155)))))
  
  
  print(paste0('RMSE total, core nearby: ', as.character(sqrt(mean(errs^2)))))
  print(paste0('RMSE 145155, core nearby: ', as.character(sqrt(mean(errs_145155^2)))))
  print(paste0('RMSE total, no core nearby: ', as.character(sqrt(mean(errs_no_core^2)))))
  print(paste0('RMSE 145155, no core nearby: ', as.character(sqrt(mean(errs_no_core_145155^2)))))
  print(paste0('RMSE total, total: ', as.character(sqrt(mean(errs_total^2)))))
  print(paste0('RMSE 145155, total: ', as.character(sqrt(mean(errs_total_145155^2)))))
  
  length(errs) / length(errs_total) 
}


library(dplyr)
library(ggplot2)
library(Matrix)
library(patchwork)
library(fda)
# new leave out results
df_list <- list()
for (fold in 1:5) {
  load(paste0('oxy_profiles_', as.character(fold), '_5_18_18_30_30-matern_spheretime_warp_seasonality/TS_prediction.RData'))
  
  pred_df <- pred_at_p_levels$pred_df
  df_list[[fold]] <- pred_df %>%
    group_by(index) %>%
    summarise(avg_full = mean(pred_full*weight),
              avg_TS = mean(pred_TS*weight),
              avg_mean = mean(pred_mean*weight),
              avg_var = mean(var_at_pressure*weight),
              var_avg = mean((pred_full - mean(pred_full*weight))^2 * weight),
              profile_unique = profile_unique[1], 
              type = 'profile') %>%
    left_join(cbind(held_out_profs, index = 1:nrow(held_out_profs)))
}

df_list_floats <- list()
for (fold in 1:10) {
  load(paste0('oxy_floats_', as.character(fold), '_5_18_18_30_30-matern_spheretime_warp_seasonality/TS_prediction.RData'))
  pred_df <- pred_at_p_levels$pred_df
  df_list_floats[[fold]] <- pred_df %>%
    group_by(index) %>%
    summarise(avg_full = mean(pred_full*weight),
              avg_TS = mean(pred_TS*weight),
              avg_mean = mean(pred_mean*weight),
              avg_var = mean(var_at_pressure*weight),
              var_avg = mean((pred_full - mean(pred_full*weight))^2 * weight),
              profile_unique = profile_unique[1], 
              type = 'float', clusters = 5) %>%
    left_join(cbind(held_out_profs, index = 1:nrow(held_out_profs)))
}
df_list_floats <- lapply(df_list_floats, function(x) x[,c(1:8, 10:19)])
df_list <- lapply(df_list, function(x) x[,c(1:18)])

mean_prediction <- rbind(dplyr::bind_rows(df_list), bind_rows(df_list_floats))
pressure_summary <- mean_prediction %>%
  mutate(pgroup = (seq(0, 2000, by = 100) + 50)[findInterval(pressure, seq(0, 2000, by = 100))]) %>%
  group_by(pgroup, type) %>%
  summarise(med_abs = median(abs(oxy - avg_full)),
            rmse = sqrt(mean((oxy - avg_full)^2)),
            med_abs_mean = median(abs(oxy - avg_mean)),
            rmse_mean = sqrt(mean((oxy - avg_mean)^2)),
            med_abs_TS = median(abs(oxy - avg_TS)),
            rmse_TS = sqrt(mean((oxy - avg_TS)^2)),
            cov_1 = mean(abs(oxy - avg_full) < sqrt(avg_var + var_avg)),
            cov_2 = mean(abs(oxy - avg_full) < 2*sqrt(avg_var + var_avg)),
            cov_3 = mean(abs(oxy - avg_full) < 3*sqrt(avg_var + var_avg)),
            med_length = median(sqrt(avg_var + var_avg)))


df_label <- data.frame('Type' = c('Mean', 'Full\nprediction', 'Temperature\nand salinity'), 
                       'type_var' = c('med_abs_mean', 'med_abs', 'med_abs_TS', 'rmse_mean', 'rmse', 'rmse_TS'))
df_model_label <- data.frame('ModelType' = c('Profiles', 'Floats'), 
                             'type' = c('profile', 'float'))
df_cov_label <- data.frame('type_var' = paste0('cov_', 1:3), 
                           'type_label' = paste0('Within ', 1:3, ' SD', c('', 's', 's')))

theme_set(theme_bw())


ggplot(data = tidyr::pivot_longer(pressure_summary, cols = starts_with('med_abs'),
                                       names_to = 'type_var', values_to = 'value') %>%
              left_join(df_label) %>% left_join(df_model_label) ,
            aes(x = pgroup, shape = Type,color = Type, y = value,
                linetype = Type)) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits = c(0, 20)) + 
  labs(x = 'Pressure (in 100 dbar bins)', y = 'Oxygen Median Absolute Error\n(μmol/kg)',
       linetype = 'Prediction\nType', shape = 'Prediction\nType', 
       color = 'Prediction\nType')+
       #title = 'Leave out point prediction performance'
       #)+
  guides(linetype = guide_legend(),
         color =guide_legend(), shape = guide_legend()) + 
  theme(text = element_text(size = 14, family = 'Arial')) + 
  facet_wrap(~ModelType, ncol = 2)
ggsave('../images/cv_error_pressure.png', height = 4,
       width = 12)

ggplot(data = tidyr::pivot_longer(pressure_summary, cols = starts_with('med_length'),
                                       names_to = 'type_var', values_to = 'value') %>%
              left_join(df_cov_label)%>% left_join(df_model_label) %>% left_join(df_label),
            aes(x = pgroup, y = value, group = factor(ModelType), color = factor(ModelType), 
                linetype =  factor(ModelType), 
                shape =  factor(ModelType))) + 
  geom_point(size = .7) + 
  geom_line() + 
  scale_y_continuous(limits = c(0, 12.5), breaks = seq(0, 12.5, by = 2.5)) + 
  theme(text = element_text(family = 'Arial')) + 
  labs(x = 'Pressure (in 100 dbar bins)', 
       y = 'Median Prediction Standard\nDeviation (μmol/kg)',
       title = 'Size of predicted uncertainties',
       color = 'Leave out\ntype',
       linetype = 'Leave out\ntype',
       shape = 'Leave out\ntype')
ggsave('../images/cv_interval_lengths_by_pressure.png', 
       height = 3*1.5, width = 5.8*2)



head(mean_prediction)
mean_prediction$upper <- 
  mean_prediction$avg_full +
  qnorm(0.975) * sqrt(mean_prediction$var_avg + mean_prediction$avg_var)

mean_prediction$lower <- 
  mean_prediction$avg_full -
  qnorm(0.975) * sqrt(mean_prediction$var_avg + mean_prediction$avg_var)

mean_prediction %>%
  group_by(type) %>%
  summarize(cov = mean(oxy > lower & oxy < upper))

mean_prediction %>%
  mutate(pgroup = (seq(0, 2000, by = 500) + 250)[findInterval(pressure, seq(0, 2000, by = 500))]) %>%
  group_by(type, pgroup) %>%
  summarize(cov = mean(oxy > lower & oxy < upper)) %>%
  tidyr::pivot_wider(names_from = 'pgroup', values_from ='cov')


# ##### BANDS
# find quantile of W_theta
library(CompQuadForm)
######################################
# Util Functions band generation
######################################

compute_normalizing_constants <- function(parameters, G, inner_prod_mat){
  lapply(1:G, function(g){
    normal_constants <- diag(t(parameters[['Omega_response1']][[g]]) %*%
                               crossprod(inner_prod_mat) %*% 
                               parameters[['Omega_response1']][[g]])
    
    Omega_response1_norm <- parameters[['Omega_response1']][[g]] %*%
      diag(1/sqrt(normal_constants))
    return(list(norm_constants = normal_constants, Omega_response1_norm = Omega_response1_norm))
  })
}

compute_bands_input <- function(cond_mean_oxy, cond_var_oxy,
                                con_mean_pred, cond_var_pred,
                                membership, norm_const, quantile_nom){
  
  pred_var_scores_renorm <- cond_var_pred * norm_const[[1]]
  
  var_scores_use <- c(diag(pred_var_scores_renorm), diag(cond_var_oxy))
  
  
  e_vals <- var_scores_use
  c_vec <- sqrt(sqrt(e_vals))
  
  cond_var_norm <- diag(var_scores_use) * 1/sqrt(e_vals)
  cond_var_norm <- t(cond_var_norm) * 1/sqrt(e_vals)
  lambda_in <- eigen(diag(e_vals/c_vec^2) %*% cond_var_norm)$values
  
  quantile <- 1 - (1-quantile_nom)/2
  max_val <- 10000
  
  # Find quantile of W_theta
  quantile_approx <- optimise(f = function(x) {
    y <- imhof(q = x,
               lambda = lambda_in)$Qq
    abs(y - (1-quantile))
  }, lower = 0, upper = max_val)$minimum
  
  # Repeat in case the optimisation doesn't work
  # Set accepted accuracy
  eps = 0.1 * (1 - quantile_nom)/2
  while ((imhof(quantile_approx, lambda_in)$Qq < (1-quantile)-eps | 
          imhof(quantile_approx, lambda_in)$Qq > (1-quantile)+eps) &
         max_val > .2) {
    max_val <- max_val/5
    quantile_approx <- optimise(f = function(x) {
      y <- imhof(q = x,
                 lambda = lambda_in)$Qq
      abs(y - (1-quantile))
    }, lower = 0, upper = max_val)$minimum
  }
  return(list(quantile_approx = quantile_approx, c_vec = c_vec))
}

return_lower_upper <- function(p, g, bands_input, basis, parameters,
                               Omega_norm, pred, quantile = 0.5, use_noise=T){
  n_out <- length(p)
  quantile_noise <- qnorm(p = (1-quantile)/n_out/2 )
  #quantile_noise <- qnorm(p = (1-quantile)/2 )
  
  Phi = fda::eval.basis(basis, evalarg = p)
  pc_values_oxy <- Phi %*% parameters[['Omega_response2']][[g]]
  pc_values_preds <- Phi %*% Omega_norm
  pc_values_all <- cbind(pc_values_preds, pc_values_oxy)
  
  c_vec = bands_input[['c_vec']]
  if(use_noise){
    band_values = sqrt(bands_input[['quantile_approx']] * 
                         rowSums(sapply(1:length(c_vec), function(x){
                           c_vec[x]^2 * pc_values_all[,x]^2})) + abs(quantile_noise) 
                       * parameters$measurement_error_response)
  } else {
    band_values = sqrt( quantile_approx* rowSums(sapply(1:length(c_vec), function(x)
    { c_vec[x]^2 * pc_values_all[,x]^2})))
  }
  return(list(lower = pred - band_values, upper = pred + band_values))
}



######################################
# Floats
######################################

setwd(paste0(this.path::here(), '/../results/'))

return_list <- list()
for (fold in as.character(1:10)) {
  load(paste0('oxy_floats_', fold, '_5_18_18_30_30-matern_spheretime_warp_seasonality/final_results.RData'))
  load(paste0('oxy_floats_', fold, '_5_18_18_30_30-matern_spheretime_warp_seasonality/TS_prediction.RData'))
  MC_max = 10
  ind = 1:((nrow(pred_at_p_levels$pred_df)/MC_max))
  
  pred_at_p_levels$pred_df$pressure = held_out_profs$pressure
  pred_at_p_levels$pred_df$oxy = held_out_profs$oxy
  
  gdf = split(pred_at_p_levels$pred_df[ind,], pred_at_p_levels$pred_df$profile_unique[ind])
  
  norm_consts = compute_normalizing_constants(estimated_model$parameters,
                                              estimated_model$params$G,
                                              estimated_model$data_inputs$inner_prod_mat_response)
  
  clusters = apply(cluster_mat, 1, function(x){
    as.numeric(names(table(x))[which.max(table(x))])
  })
  
  pred_index = pred_set_up$is_prediction
  
  n_pcs_r = estimated_model$params$pc_response2
  n_pcs_p = estimated_model$params$pc_predictors
  n_pcs_total = n_pcs_p + n_pcs_r
  
  cond_means_oxy = pred_score_results$pcs_array_response[pred_index, 1, ]
  cond_means_pred = pred_score_results$pcs_array_predictors[pred_index, 1, ]
  cond_var_oxy = lapply(pred_score_results$variances[[1]], function(x){
    return(x[(n_pcs_p+1):n_pcs_total, (n_pcs_p+1):n_pcs_total])
  })
  cond_var_pred = lapply(pred_score_results$variances[[1]], function(x){
    return(x[1:n_pcs_p,1:n_pcs_p])
  })
  
  lower_band = c()
  upper_band = c()
  for( i in 1:length(gdf) ){
    if( i %% 100 == 0){
      print(i)
    }
    g = clusters[i]
    df = gdf[[i]]
    
    nc = norm_consts[[g]]
    bands_input = compute_bands_input(cond_means_oxy[i], cond_var_oxy[[i]],
                                      cond_means_pred[i], cond_var_pred[[i]],
                                      g, nc, 0.95)
    
    lower_upper = return_lower_upper(df$pressure, g, bands_input, 
                                     estimated_model$params$basis_response,
                                     estimated_model$parameters,
                                     nc[['Omega_response1_norm']], df$pred_full, use_noise=T)
    
    lower_band = c(lower_band, lower_upper[['lower']])
    upper_band = c(upper_band, lower_upper[['upper']])
  }
  
  pred_at_p_levels$pred_df$lower_band = lower_band
  pred_at_p_levels$pred_df$upper_band = upper_band
  
  return_list[[fold]] <- pred_at_p_levels$pred_df[, c('oxy', 'lower_band', 'upper_band', 'profile_unique')]
  print(fold)
}

dplyr::bind_rows(return_list) %>%
  group_by(profile_unique) %>%
  summarize(cov = sum((oxy > upper_band | oxy < lower_band)) == 0) %>%
  summarize(mean(cov))

######################################
# Profiles
######################################

setwd(paste0(this.path::here(), '/../results/'))

return_list_profiles <- list()
for (fold in as.character(1:5)) {
  load(paste0('oxy_profiles_', fold, '_5_18_18_30_30-matern_spheretime_warp_seasonality/final_results.RData'))
  load(paste0('oxy_profiles_', fold, '_5_18_18_30_30-matern_spheretime_warp_seasonality/TS_prediction.RData'))
  MC_max = 10
  ind = 1:((nrow(pred_at_p_levels$pred_df)/MC_max))
  lower = pred_at_p_levels$pred_df$pred_full - sqrt(pred_at_p_levels$pred_df$var_at_pressure) * qnorm(0.975)
  upper = pred_at_p_levels$pred_df$pred_full + sqrt(pred_at_p_levels$pred_df$var_at_pressure) * qnorm(0.975)
  
  pred_at_p_levels$pred_df$upper = upper
  pred_at_p_levels$pred_df$lower = lower
  pred_at_p_levels$pred_df$pressure = held_out_profs$pressure
  pred_at_p_levels$pred_df$oxy = held_out_profs$oxy
  mean(held_out_profs$oxy > lower & held_out_profs$oxy < upper)
  
  gdf = split(pred_at_p_levels$pred_df[ind,], pred_at_p_levels$pred_df$profile_unique[ind])
  
  norm_consts = compute_normalizing_constants(estimated_model$parameters,
                                              estimated_model$params$G,
                                              estimated_model$data_inputs$inner_prod_mat_response)
  
  clusters = apply(cluster_mat, 1, function(x){
    as.numeric(names(table(x))[which.max(table(x))])
  })
  
  pred_index = pred_set_up$is_prediction
  
  n_pcs_r = estimated_model$params$pc_response2
  n_pcs_p = estimated_model$params$pc_predictors
  n_pcs_total = n_pcs_p + n_pcs_r
  
  cond_means_oxy = pred_score_results$pcs_array_response[pred_index, 1, ]
  cond_means_pred = pred_score_results$pcs_array_predictors[pred_index, 1, ]
  cond_var_oxy = lapply(pred_score_results$variances[[1]], function(x){
    return(x[(n_pcs_p+1):n_pcs_total, (n_pcs_p+1):n_pcs_total])
  })
  cond_var_pred = lapply(pred_score_results$variances[[1]], function(x){
    return(x[1:n_pcs_p,1:n_pcs_p])
  })
  
  lower_band = c()
  upper_band = c()
  for( i in 1:length(gdf) ){
    if( i %% 100 == 0){
      print(i)
    }
    g = clusters[i]
    df = gdf[[i]]
    
    nc = norm_consts[[g]]
    bands_input = compute_bands_input(cond_means_oxy[i], cond_var_oxy[[i]],
                                      cond_means_pred[i], cond_var_pred[[i]],
                                      g, nc, 0.95)
    
    lower_upper = return_lower_upper(df$pressure, g, bands_input, 
                                     estimated_model$params$basis_response,
                                     estimated_model$parameters,
                                     nc[['Omega_response1_norm']], df$pred_full, use_noise=T)
    
    lower_band = c(lower_band, lower_upper[['lower']])
    upper_band = c(upper_band, lower_upper[['upper']])
  }
  
  pred_at_p_levels$pred_df$lower_band = lower_band
  pred_at_p_levels$pred_df$upper_band = upper_band
  
  return_list_profiles[[fold]] <- pred_at_p_levels$pred_df[, c('oxy', 'lower_band', 'upper_band', 'profile_unique')]
  print(fold)
}

dplyr::bind_rows(return_list_profiles) %>%
  group_by(profile_unique) %>%
  summarize(cov = sum((oxy > upper_band | oxy < lower_band)) == 0) %>%
  summarize(mean(cov))
