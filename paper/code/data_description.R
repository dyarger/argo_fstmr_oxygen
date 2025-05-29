load('paper/data/soccom_processed_oxy_05_05_2024.RData')

df_soccom <- df_list[[1]]
df_core <- df_list[[2]]

df_all <- as.data.frame(rbind(cbind(df_soccom[, 'pressure'], 'Variable' = 'Oxygen'),
                              cbind(df_core[, 'pressure'], 'Variable' = 'Temperature/salinity')))
colnames(df_all)[1] <- 'Pressure (decibars)' 
head(df_all)
df_all$`Pressure (decibars)` <- as.numeric(df_all$`Pressure (decibars)`)

source('paper/code/src/plot_src.R')
ggplot(data = df_all) + 
  geom_histogram(aes(x = `Pressure (decibars)`), color = 'black', 
                 fill = 'white') +
  facet_wrap(~Variable) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse') + 
  labs(y = 'Number of measurements') + 
  theme_bw()
ggsave('paper/images/pressure_hist.png', width = 6, height = 6)





df_unique <- df_soccom[!duplicated(df_soccom$profile_unique),]

# create grid
library(ncdf4)
grid_size <- 1
grid_file <- nc_open('paper/data/grid.nc', write = F)
grid_long <- ncvar_get(grid_file, 'XC')[,1]
grid_lat <- ncvar_get(grid_file, 'YC')[1,]
grid_depth <- ncvar_get(grid_file, 'Depth')[seq(1, length(grid_long), by = round(6*grid_size)),
                                            seq(1, length(grid_lat), by = round(6*grid_size))]
grid_long <- grid_long[seq(1, length(grid_long), by = round(6*grid_size))]
grid_lat <- grid_lat[seq(1, length(grid_lat), by = round(6*grid_size))]

height <- grid_lat[2:length(grid_lat)] - grid_lat[1:(length(grid_lat) - 1)]
height <- data.frame("latitude" = grid_lat, height = c(height[1], height))

grid_df <- cbind(expand.grid('longitude' = grid_long, 'latitude' = grid_lat),
                 depth = as.double(grid_depth), index = 1:length(grid_long),
                 index_all = 1:length(grid_depth))
grid_df[['width']] <- grid_size
grid_df <- merge(grid_df, height, by = "latitude")

closest_grid <- rep(0, nrow(df_unique))
library(fields)
for (i in 1:length(closest_grid)) {
  distances <- rdist.earth.vec(grid_df[,c('longitude', 'latitude')],
                               as.matrix(df_unique[i, c('longitude', 'latitude')]))
  closest_grid[i] <- which.min(distances)
}
table_df <- data.frame(index_all = as.integer(names(table(closest_grid))),
                       value = as.integer(as.vector(table(closest_grid))))
grid_tally <- dplyr::left_join(dplyr::rename(grid_df, id = index), table_df, by = 'index_all', keep = T)
library(dplyr)
ggplot(data =  grid_tally %>% filter(!is.na(value)))+
  geom_tile(aes(x = longitude, y = latitude, fill = value,
                height = height + .022, width = width +.022)) +
  SO_coord + SO_theme +
  fronts_dark + continents +
  latitude_lines + longitude_lines +
  labs(color ='Probability', fill = 'Probability', 
       title = 'Probability of the most likely cluster')+
  theme(panel.grid = element_blank()) +
  scale_color_viridis_c()+  scale_fill_viridis_c()

