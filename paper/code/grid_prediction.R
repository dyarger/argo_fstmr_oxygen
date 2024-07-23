library(fstmr)
library(dplyr)
library(fda)
library(Matrix)
library(doParallel)
library(foreach)
library(this.path)

setwd(paste0(this.path::here(), '/../results/'))
load('oxy_none__5_18_18_30_30-matern_spheretime_warp_seasonality/final_results.RData')
source('../code/misc/pred_at_pressure_levels_src.R')

registerDoParallel(4)

grid_locs_total <- get_grid_info(1.5)

# Get near core data

load('../data/core_processed_06_21.RData')

core_data[['day']] <- julian(as.Date(core_data[['date']]), origin = as.Date('2000-01-01'))
core_data[['dayofyear']] <- julian(as.Date(core_data[['date']]), origin = as.Date('2000-01-01')) %% 365.25

one_row_core_data = core_data %>%
  group_by(profile_unique) %>%
  slice(1)

dates_to_predict = c('2020-01-15','2020-04-15','2020-07-15','2020-10-15')

dates_names = c('january', 'april', 'july', 'october')
grid_chunks = 1
chunk_length = nrow(grid_locs_total) / grid_chunks


par_results <- foreach(j = c(1:4)) %dopar% {
  
  for (l in 1:grid_chunks){
    
    grid_ind = ((l-1)*chunk_length+1):(l*chunk_length)
    
    grid_locs <- cbind(grid_locs_total[grid_ind,], profile_unique = grid_ind,
                       date = dates_to_predict[j],
                       time = as.vector(julian(as.Date(dates_to_predict[j], format = "%Y-%m-%d"), origin = as.Date("2000-01-01"))))
    
    grid_locs[['dayofyear']] <- grid_locs[['time']] %% 365.25
    grid_locs[['day']] <- grid_locs[['time']]  
    one_row_grid_data = grid_locs %>%
      group_by(profile_unique) %>%
      slice(1) %>%
      mutate(longitude = ifelse(longitude > 180, longitude - 360, longitude))
    
    day = julian(as.Date(one_row_core_data$date,format = '%Y-%m-%d'),
                 origin = as.Date('2000-01-01'))%% 365.25
    
    day_grid = one_row_grid_data$time %% 365.25
    close_profs = c()
    counter = 1
    close_counter = 0
    
    
    time_dist1 <- abs(day - day_grid[1])
    time_dist2 <- abs(day - day_grid[1] + 365.25)
    time_dist3 <- abs(day - day_grid[1] - 365.25)
    time_dist <- time_dist3
    b1 <- time_dist1 < time_dist2 & time_dist1 < time_dist3
    b2 <- time_dist2 < time_dist1 & time_dist2 < time_dist3
    time_dist[b1] <- time_dist1[b1]
    time_dist[b2] <- time_dist2[b2]
    ind <- which(abs(time_dist) < 30)
    one_row_core_data_ind = one_row_core_data[ind,]
    
    
    for (i in 1:nrow(one_row_grid_data)){
      
      lon_grid = one_row_grid_data$longitude[i]
      lat_grid = one_row_grid_data$latitude[i]
      
      ind = ( abs(one_row_core_data_ind$latitude - lat_grid) <= 2 ) & ( abs(one_row_core_data_ind$longitude - lon_grid) <= 2 ) 
      
      candidates = one_row_core_data_ind$profile_unique[ind] 
      
      if(length(candidates) > 0){
        close_profs = unique(c(close_profs, candidates))
        close_counter = close_counter + 1
      }
      
      counter = counter + 1
      if(counter %% 100 == 0){
        print(counter)
      }
    }
    
    close_core_df = core_data %>%
      filter(profile_unique %in% close_profs)
    
    estimated_model$params[['m_pred_predictors']] <- 20
    estimated_model$params[['m_pred_response']] <- 20
    estimated_model$params[['m_variance_computation']] <- 4
    estimated_model$params[['m_variance_computation']] <- 4
    estimated_model$params[['time']] <- 'date'
    pred_set_up <- fstmr:::set_up_prediction(estimated_model, estimated_model$params, core_data = close_core_df, grid = grid_locs)
    
    MC_max_cluster <- max(c(100, estimated_model$params$MC_max))
    
    cluster_mat_full <- fstmr:::rPotts_prediction(params[['G']], pred_set_up[['neighbors_all']], pred_set_up[['has_data']],
                                                  pred_set_up[['p_mat']], pred_set_up[['like_mat']],
                                                  estimated_model[['parameters']][['theta']], 
                                                  n_samples = MC_max_cluster, skip = 4, init_pred_cluster = pred_set_up[['memberships']],
                                                  dists = pred_set_up[['dist_all']])
    
    mc_ind = round(1:40  * 2.5)
    cluster_mat_use = cluster_mat_full[,mc_ind]
    estimated_model$params$MC_max = length(mc_ind)
    
    a = Sys.time()
    pred_score_results <- fstmr:::pred_E_step(pred_set_up, estimated_model$params, cluster_mat_use,  prediction_only = T)
    print(Sys.time() - a)
    
    path = paste0(getwd(), '/')
    file_name = paste0('1.5_one_chunk_oxy_grid_', dates_names[j], '_', as.character(l), '.RData')
    save(pred_score_results, cluster_mat_full, pred_set_up,
         file = paste0(path, file_name))
    
  }
}




