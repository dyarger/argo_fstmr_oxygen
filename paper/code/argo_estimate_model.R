library(fstmr)
library(dplyr)
library(fda)
library(Matrix)
library(doParallel)
library(foreach)
library(this.path)

# variables and basic setup
set.seed(56)

# Should be adapted to the specific environment and controls how many folds are run in parallel
registerDoParallel(1)

setwd(paste0(this.path::here(), '/../..'))

params = list()

# 1: Leave out type
# 2: Number of principal components
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  params[['leave_out_type']] = 'none'
  params[['pc_response1']]  = params[['pc_predictors']] <- 18
  params[['pc_response2']] <- 18
} else if(length(args) == 2){
  params[['leave_out_type']] = args[1]
  params[['pc_response1']]  = params[['pc_predictors']] <- as.numeric(args[2])
  params[['pc_response2']] <- as.numeric(args[2])
} else {
  stop("The arguments were specified incorrectly (neither length 0 nor 2).")
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

params[['results_folder']] <- paste0(getwd(), '/paper/results/')

# set up data variables
params[['leave_out_prediction_core']] = T
params[['n_core_profiles']] <- 0

params[['id']] = 'profile_unique'
params[['ind_name']] = 'pressure'
params[['loc_names']] = c('longitude', 'latitude')
params[['var']] = params[['Y']] = 'oxy'
params[['X']] = c('temp', 'psal')
params[['domain_Y']] = c(0,2000)
params[['domain_X']] <- list(c(0,2000), c(0,2000))

# initialization
params[['init_strategy']] = 'kmeans'
params[['levels']] = list('nitrate' = seq(0,2000,length.out = 15),
                          'oxy' = seq(0,2000,length.out = 15),
                          'temp'= seq(0,2000,length.out = 15),
                          'psal'= seq(0,2000,length.out = 15))
# MRF parameters
params[['nn_strategy']] = 'nn'
params[['nn_type']] = 'spacewarp_day'
params[['nn_range_time']] = NA
params[['nn_time_param']] = 12
params[['nn_lat_param']] = 3
params[['nn']] = 15
params[['remove_between_land']] <- T

# basis parameters
params[['knots']] = c(seq(0,100,5), seq(110, 290, by = 10), seq(300,1000,50), 1100, 1200, 1300, 1400, 1500, 1600, 2000)
params[['basis_response']] = create.bspline.basis(params[['knots']])
params[['basis_predictors']] = list(create.bspline.basis(params[['knots']]),
                                    create.bspline.basis(params[['knots']])) 

# 4 spatial covariance parameters
params[['m_train_predictors']] <- 30
params[['m_train_response']] <- 30
params[['m_AIC_predictors']] <- 15
params[['m_AIC_response']] <- 15
params[['m_pred_predictors']] <- 60
params[['m_pred_response']] <- 60
params[['maxit']] = 20
params[['time']] <- 'day'
params[['dayofyear']] <- 'dayofyear'
params[['covariance_function']] <- 'matern_spheretime_warp'

# number of clusters
params[['G']] = 5

# smoothing parameters
params[['lambda_mean_response']] <- 10^(3:5)
params[['lambda_mean_predictors']] <- 10^(3:5)
params[['lambda_pcs_response']] <- 10^(4:6)
params[['lambda_lt']] <- 10^(3:6)
params[['lambda_pcs_predictors']] <- 10^(3:5)
params[['cv_skip']] <- 10


# number of MCEM and EM iterations
params[['MC_max']] =  10
params[['MC_max_AIC']] = 10
params[['EM_iter']] = 20
params[['EM_iter_init']] <-  3

params[['use_lanczos_training']] = T
params[['use_lanczos_prediction']] = T

params[['seasonal_mean']] <- T
params[['num_seasonal_basis']] <- 3

params[['use_nugget']] = F
params[['use_dist_mrf']] = T
params[['skip_range_update_init']] = T

log_dir = paste0(getwd(),
                 "/logs/",as.character(params[['pc_response2']]),
                 "_", as.character(params[['pc_response1']]),
                 "_", as.character(params[['m_train_predictors']]),
                 "_", as.character(params[['m_pred_predictors']]),
                 "_", params[['covariance_function']], 
                 "_", params[['leave_out_type']], "_seaonality_test/")

system(paste('mkdir', log_dir))

printf <- function(model, diagnostics, parameters, data_inputs, params, iter) {
  print(diagnostics[['likelihood']])
  write(paste0('MeasurementError is in Iteration' , as.character(iter), ": ",
               as.character(round(parameters[['measurement_error_response']], 3))),
        file = paste0(log_dir, as.character(params[['folds']]), "log.txt"),
        append = T)
}

par_results <- foreach(fold = c(1:params[['n_folds']])) %dopar% {
  
  if (params[['leave_out_type']] != 'none'){
    params[['folds']] <- fold
  } else {
    params[['n_folds']] = NA
  }
  # things to print each MC iteration

  ##### Load and subset data #####
  load('paper/data/soccom_processed_oxy_05_05_2024.RData')
  source('paper/code/misc/prepare_argo_data_src.R')
  data_leave_out <- leave_out_data(df_list, params[['leave_out_type']],
                                   params[['folds']], params[['n_folds']], params[['leave_out_prediction_core']],
                                   params[['n_core_profiles']], params[['core_sampling_strategy']],
                                   reduce_data = F
  )
  df_list <- data_leave_out[[1]]; df_list_pred <- data_leave_out[[2]]
  held_out_profs <- df_list_pred[[1]]
  rm(data_leave_out)

  params[['dir_name']] <- paste0(params[['results_folder']], params[['Y']], '_',
                                 params[['leave_out_type']], '_', params[['folds']], '_', params[['G']], '_',
                                 params[['pc_predictors']], '_', params[['pc_response2']], '_',
                                 params[['m_train_predictors']], '_', params[['m_train_response']], '-',
                                 params[['covariance_function']],'_seasonality','/')
  
  # If not exists, create directory to store files
  system(paste('mkdir', params[['dir_name']]))

  # estimate model
  estimated_model <- fstmr(data = df_list, params,
                                  compute_diagnostics = ifelse(params[['leave_out_type']] == 'none', T, F),
                                  verbose = 50, printf = printf)
  
  save(params, estimated_model, diagnostics,
       file = paste0(params[['dir_name']], 'final_results.RData'))
  1
}
