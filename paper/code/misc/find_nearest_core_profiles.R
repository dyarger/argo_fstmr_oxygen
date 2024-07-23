library(dplyr)
library(fields)
setwd(paste0(this.path::here(), '/../../..'))

load('paper/data/core_processed_06_21.RData')
load('paper/data/soccom_processed_oxy_05_05_2024.RData')

data = df_list[[1]]

unique_soccom_profiles = unique(data$profile_unique)
unique_core_profiles = unique(core_data$profile_unique)

# Gather data + loc per profile unique
loc_df = data.frame(profile_unique = unique_core_profiles[1],
                    day = core_data[1,'day'],
                    longitude = core_data[1,'longitude'],
                    latitude = core_data[1,'latitude'])

one_row_core_data = core_data %>% 
  group_by(profile_unique) %>% 
  slice(1)
one_row_soccom_data = data %>% 
  group_by(profile_unique) %>% 
  slice(1)

day_soccom = one_row_soccom_data$day[1]
day_core = one_row_core_data$date

close_profs = c()
counter = 1
has_close_prof_counter = 0
has_cn = c()
for (i in 1:nrow(one_row_soccom_data)){
  
  day_soccom = one_row_soccom_data$day[i]
  lon_soccom = one_row_soccom_data$longitude[i]
  lat_soccom = one_row_soccom_data$latitude[i]
  
  ind = which(abs(day_core - day_soccom) < 30)
  
  candidates = one_row_core_data[ind,] %>%
    filter(abs(longitude - lon_soccom) <= 2) %>%
    filter(abs(latitude - lat_soccom) <= 2)
  
  if(nrow(candidates) > 0){
    close_profs = unique(c(close_profs, candidates$profile_unique))
    has_close_prof_counter = has_close_prof_counter + 1
    has_cn = c(has_cn, i)
  }
  
  counter = counter + 1
  if(counter %% 1000 == 0){
    print(counter)
  }
}
print(has_close_prof_counter)


close_core_df = core_data %>%
  filter(profile_unique %in% close_profs)

soccom_unique <- unique(df_list[[1]]$profile_unique)

close_core_df_unique <- close_core_df[!duplicated(close_core_df$profile_unique),]

profiles_to_remove <- close_core_df_unique$profile_unique[close_core_df_unique$profile_unique %in% soccom_unique]

close_core_df_prelim <- close_core_df  %>%
  filter(!(profile_unique %in% profiles_to_remove))

close_core_df <- close_core_df_prelim
length(unique(close_core_df$profile_unique))
length(unique(close_core_df_prelim$profile_unique))


save(close_core_df, file = 'paper/data/close_core.RData')

