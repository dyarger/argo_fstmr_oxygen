library(fstmr)
library(dplyr)
library(fda)
library(Matrix)
library(doParallel)
library(foreach)
library(this.path)

# Should be adapted to the specific environment and controls how many folds are run in parallel
# Note that part of the prediction uses up to processing units in parallel
registerDoParallel(3)

# variables and basic setup
set.seed(56)

setwd(paste0(this.path::here(), '/../..'))

params = list()
params[['leave_out_type']] = 'profiles'
params[['pc_response1']]  = params[['pc_predictors']] <- 18
params[['pc_response2']] <- 18
params[['leave_out_prediction_core']] = F

# 1: Leave out type
# 2: Number of principal components
# 3: Whether to use T/S at prediction location
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 3){
  params[['leave_out_type']] = args[1]
  params[['pc_response1']]  = params[['pc_predictors']] <- as.numeric(args[2])
  params[['pc_response2']] <- as.numeric(args[2])
  params[['leave_out_prediction_core']] = predicted_with_TS = as.logical(args[3])
} else {
  stop("The arguments were specified incorrectly (neither length 0 nor 3).")
}

# type of leave out data
if (params[['leave_out_type']] == 'floats'){
  params[['n_folds']] = 10
} else if (params[['leave_out_type']] == 'profiles'){
  params[['n_folds']] = 5
} else if (params[['leave_out_type']] == 'none'){
  params[['n_folds']] = 1
} else {
  stop("Specified invalid leave_out_type, has to be 'floats', 'profiles' or 'none'")
}

params[['results_folder']] <- paste0(getwd(), '/results/')

# 4 spatial covariance parameters
params[['m_pred_predictors']] <- 60
params[['m_pred_response']] <- 60

# number of MCEM and EM iterations
params[['MC_max']] =  10
params[['use_lanczos_prediction']] = T

params[['leave_out_prediction_core']] = F

par_results <- foreach(fold = c(1:params[['n_folds']])) %dopar% {
  
  params[['folds']] = fold
  
  ##### Load and subset data #####
  load('paper/data/soccom_processed_oxy_05_05_2024.RData')
  source('paper/code/misc/prepare_argo_data_src.R')
  
  data_leave_out <- leave_out_data(df_list, params[['leave_out_type']],
                                   params[['folds']], params[['n_folds']], params[['leave_out_prediction_core']],
                                   0, params[['core_sampling_strategy']],
                                   reduce_data = F
  )
  df_list <- data_leave_out[[1]]; df_list_pred <- data_leave_out[[2]]
  held_out_profs <- df_list_pred[[1]]
  rm(data_leave_out)

  params[['dir_name']] <- paste0(params[['results_folder']], 'oxy_',
                                 params[['leave_out_type']], '_', params[['folds']], '_5_',
                                 params[['pc_predictors']], '_', params[['pc_response2']], '_30_30-matern_spheretime_warp_seasonality/')


  load(paste0(params[['dir_name']], '/final_results.RData'))
  save_location = params[['dir_name']]
  
  if (!predicted_with_TS){
    load('paper/data/close_core.RData')
    close_core_df[['day']] <- close_core_df[['date']]
    close_core_df[['dayofyear']] <- julian(as.Date(close_core_df[['date']]), origin = as.Date('2000-01-01')) %% 365.25
  
    one_row_core_data = close_core_df %>%
      group_by(profile_unique) %>%
      slice(1)
  
    one_row_soccom_data = held_out_profs %>%
      group_by(profile_unique) %>%
      slice(1)
  
    day = one_row_core_data$day
    close_profs = c()
    counter = 1
    for (i in 1:nrow(one_row_soccom_data)){
      day_soccom = one_row_soccom_data$day[i]
      lon_soccom = one_row_soccom_data$longitude[i]
      lat_soccom = one_row_soccom_data$latitude[i]
  
      ind = which(abs(day - day_soccom) < 30)
  
      candidates = one_row_core_data[ind,] %>%
        filter(abs(longitude - lon_soccom) <= 2) %>%
        filter(abs(latitude - lat_soccom) <= 2)
  
      if(nrow(candidates) > 0){
        close_profs = unique(c(close_profs, candidates$profile_unique))
      }
  
    }

    close_core_df = close_core_df %>%
      filter(profile_unique %in% close_profs)
    
    # Prediction for leave out
    grid = held_out_profs[,c('longitude', 'latitude', 'dayofyear', 'day')] %>%
      distinct()
    grid[['time']] = julian(as.Date(grid$day, format = '%Y-%m-%d'), origin = as.Date('2000-01-01'))
  }
  
  if(predicted_with_TS){
    estimated_model$params[['m_variance_computation']] <- 5
    params[['m_variance_computation']] <- 5
    params[['compute_variances']] = TRUE
    pred_set_up <- fstmr:::set_up_prediction(estimated_model,
                                             params,
                                             new_data = held_out_profs)
  } else {
    params[['compute_variances']] = FALSE
    pred_set_up <- fstmr:::set_up_prediction(estimated_model,
                                             params,
                                             core_data = close_core_df,
                                             grid = grid)    
  }
  
  cluster_mat <- fstmr:::rPotts_prediction(params[['G']], pred_set_up[['neighbors_all']], pred_set_up[['has_data']],
                                           pred_set_up[['p_mat']], pred_set_up[['like_mat']],
                                           estimated_model[['parameters']][['theta']],
                                           n_samples = params[['MC_max']], skip = 4, dists = pred_set_up[['dist_all']])
  
  print(summary(apply(cluster_mat, 1, var)))

  pred_score_results <- fstmr:::pred_E_step(pred_set_up, params, cluster_mat)
  pred_at_p_levels <- fstmr:::pressure_level_prediction(held_out_profs, params, pred_set_up, pred_score_results,
                                                        cluster_mat, seasonal_mean = T)

  pred_set_up[['prediction_model']][['data_inputs']] <- NULL
  
  prediction_file_name = ifelse(predicted_with_TS, 'TS_prediction.RData', 'prediction.RData')
  
  path = paste0('results/oxy_', params[['leave_out_type']], '_', as.character(fold), '_5_18_18_30_30-matern_spheretime_warp_seasonality/')
  print(paste0(path,prediction_file_name))
  save(pred_set_up, pred_at_p_levels, pred_score_results, held_out_profs, cluster_mat,
       file = paste0(path,prediction_file_name))
  1
}