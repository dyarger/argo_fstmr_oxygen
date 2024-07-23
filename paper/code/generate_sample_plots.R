library(fstmr)
library(fda)
library(dplyr)
library(ggplot2)
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
fold = as.character(1)
load(paste0('oxy_floats_', fold, '_5_18_18_30_30-matern_spheretime_warp_seasonality/final_results.RData'))
load(paste0('oxy_floats_', fold, '_5_18_18_30_30-matern_spheretime_warp_seasonality/TS_prediction.RData'))
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

df = pred_at_p_levels$pred_df %>%
  filter(profile_unique == "5906249_13") %>%
  group_by(pressure) %>%
  summarise(pred_full = mean(pred_full), lower = mean(lower),
            upper = mean(upper), lower_band = mean(lower_band),
            upper_band = mean(upper_band), oxy = mean(oxy))

a = ggplot(df, aes(x = pressure)) +  ylim(160, 350) + 
  geom_point(aes(y=oxy, color = 'True Oxygen'), size = 3) +
  geom_line(aes(y=pred_full, color = 'Prediction'), linewidth = 3) +
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.3) + 
  geom_ribbon(aes(ymin=lower_band, ymax = upper_band), alpha=0.2) + 
  theme_bw() + scale_x_reverse() + coord_flip() +
  xlab("Pressure") + ylab("Oxygen (μmol/kg)") + scale_color_manual(values = c('goldenrod1','deepskyblue4')) +
  theme(legend.title = element_text(size = 0), legend.position = "bottom", legend.text = element_text(size = 30)) + 
  theme(axis.text=element_text(size=30), axis.title=element_text(size=30)) + 
  theme(text = element_text(family ="Arial"))
a
w = 7
ggsave(plot = a, file = "../images/float_sample_full.png", height = 1.9*w, width = w)

a = ggplot(df %>% filter(pressure <= 500), aes(x = pressure)) +  ylim(160, 350) + 
  geom_point(aes(y=oxy, color = 'True Oxygen'), size = 3) +
  geom_line(aes(y=pred_full, color = 'Prediction'), linewidth = 3) +
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.3) + 
  geom_ribbon(aes(ymin=lower_band, ymax = upper_band), alpha=0.2) + 
  theme_bw() + scale_x_reverse() + coord_flip() +
  xlab("Pressure") + ylab("Oxygen (μmol/kg)") + scale_color_manual(values = c('goldenrod1','deepskyblue4')) +
  theme(legend.title = element_text(size = 0), legend.position = "bottom", legend.text = element_text(size = 30)) + 
  theme(axis.text=element_text(size=30), axis.title=element_text(size=30)) + 
  theme(text = element_text(family ="Arial"))
a
w = 7
ggsave(plot = a, file = "../images/float_sample_zoom.png", height = 1.9*w, width = w)

######################################
# Profiles
######################################

setwd(paste0(this.path::here(), '/../results/'))
fold = as.character(1)
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

df = pred_at_p_levels$pred_df %>%
  filter(profile_unique == "5905072_1") %>%
  group_by(pressure) %>%
  summarise(pred_full = mean(pred_full), lower = mean(lower),
            upper = mean(upper), lower_band = mean(lower_band),
            upper_band = mean(upper_band), oxy = mean(oxy))


a = ggplot(df, aes(x = pressure)) +  ylim(160, 350) + 
  geom_point(aes(y=oxy, color = 'True Oxygen'), size = 3) +
  geom_line(aes(y=pred_full, color = 'Prediction'), linewidth = 3) +
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.3) + 
  geom_ribbon(aes(ymin=lower_band, ymax = upper_band), alpha=0.2) + 
  theme_bw() + scale_x_reverse() + coord_flip() +
  xlab("Pressure") + ylab("Oxygen (μmol/kg") + scale_color_manual(values = c('goldenrod1','deepskyblue4')) +
  theme(legend.title = element_text(size = 0), legend.position = "bottom", legend.text = element_text(size = 30)) + 
  theme(axis.text=element_text(size=30), axis.title=element_text(size=30)) + 
  theme(text = element_text(family ="Arial"))
a
w = 7
ggsave(plot = a, file = "../images/profile_sample_full.png", height = 1.9*w, width = w)

a = ggplot(df %>% filter(pressure <= 500), aes(x = pressure)) +  ylim(160, 350) + 
  geom_point(aes(y=oxy, color = 'True Oxygen'), size = 3) +
  geom_line(aes(y=pred_full, color = 'Prediction'), linewidth = 3) +
  geom_ribbon(aes(ymin = lower,ymax = upper), alpha = 0.3) + 
  geom_ribbon(aes(ymin=lower_band, ymax = upper_band), alpha=0.2) + 
  theme_bw() + scale_x_reverse() + coord_flip() +
  xlab("Pressure") + ylab("Oxygen (μmol/kg)") + scale_color_manual(values = c('goldenrod1','deepskyblue4')) +
  theme(legend.title = element_text(size = 0), legend.position = "bottom", legend.text = element_text(size = 30)) + 
  theme(axis.text=element_text(size=30), axis.title=element_text(size=30)) +
  theme(text = element_text(family ="Arial"))
a
w = 7
ggsave(plot = a, file = "../images/profile_sample_zoom.png", height = 1.9*w, width = w)

