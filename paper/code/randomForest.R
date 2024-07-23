library(fstmr)
library(dplyr)
library(fda)
library(randomForest)
library(rfinterval)
# 1 variables and basic setup
set.seed(56)
params = list()
params[['Y']] = 'oxy'
# 2 type of leave out data
params[['leave_out_type']] = 'profiles'
params[['n_folds']] = 5
errs = c()
lower = c()
upper = c()
obs = c()
preds = c()
method = "quantreg"
for (n in 1:params[['n_folds']]){#1:params[['n_folds']]){
  params[['folds']] <- n
  params[['leave_out_prediction_core']] = T
  
  ##### Load and subset data #####
  load(paste0('paper/datasoccom_processed_', params[['Y']], '_05_05_2024.RData'))
  
  source('paper/code/misc/prepare_argo_data_src.R')
  data_leave_out <- leave_out_data(df_list, params[['leave_out_type']], 
                                   params[['folds']], params[['n_folds']], params[['leave_out_prediction_core']],
                                   0, params[['core_sampling_strategy']],
                                   reduce_data = F)
  
  df_list <- data_leave_out[[1]]; df_list_pred <- data_leave_out[[2]]
  held_out_profs <- df_list_pred[[1]]; prediction_profs <- df_list_pred[[2]]
  profiles_for_prediction <- data_leave_out[[3]]
  rm(data_leave_out)
  df_unique <- df_list[[2]][!duplicated(df_list[[2]][['profile_unique']]),]
  
  # random forest
  p = 6 # T, S, lat, long, year, month
  m = floor(p/3)
  B = 500
  
  df_train_150 <- df_list[[1]] %>%
    filter(pressure < 155, pressure > 145) %>%
    mutate(year = as.double(substr(day, 1, 4)),
           month = as.double(substr(day, 6, 7)))
  df_test_150 <- df_list_pred[[2]] %>%
    filter(pressure < 155, pressure > 145) %>%
    mutate(year = as.double(substr(day, 1, 4)),
           month = as.double(substr(day, 6, 7)))
  df_test_150 <- df_test_150[!is.na(df_test_150[[params[['Y']]]]),]
  
  params[['X']] = c('temp', 'psal', 'longitude', 'latitude', 'year', 'month')
  
  X <- df_train_150[, params[['X']]]
  data_train = df_train_150[, c(params[['Y']], params[['X']])]
  X_out <- df_test_150[, params[['X']]]
  data_test <- df_test_150[, c(params[['Y']], params[['X']])]
  rf_limit_quantreg <- rfinterval('oxy ~ temp+psal+longitude+latitude+year+month',
                                  train_data = as.data.frame(data_train),
                                  test_data = as.data.frame(data_test),
                                  params_ranger = list(mtry = m),
                                  method = method,
                                  alpha = 0.05)
  if(method == 'quantreg'){
    lower = c(lower, rf_limit_quantreg$quantreg_interval[,1])
    upper = c(upper, rf_limit_quantreg$quantreg_interval[,2])
  } else if (method == 'oob'){
    lower = c(lower, rf_limit_quantreg$oob_interval[,1])
    upper = c(upper, rf_limit_quantreg$oob_interval[,2])
  }
  preds = c(preds, rf_limit_quantreg$testPred)
  errs = c(errs, rf_limit_quantreg$testPred - df_test_150$oxy)
  obs = c(obs, df_test_150$oxy)
  print(median(abs(rf_limit_quantreg$testPred - df_test_150$oxy)))
  print(sqrt(mean((rf_limit_quantreg$testPred - df_test_150$oxy)^2)))
  print(mean(rf_limit_quantreg$quantreg_interval[,1] < df_test_150$oxy &  df_test_150$oxy < rf_limit_quantreg$quantreg_interval[,2]) )
}
print(median(abs(errs)))
print(sqrt(mean(abs(errs^2))))
print(mean(lower < obs &  obs < upper) )
mean(abs(upper - lower))


