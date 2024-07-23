library(dplyr)
library(ggplot2)
library(patchwork)
library(fda)
library(Matrix)
theme_set(theme_bw() + theme(legend.position = 'bottom'))
source('paper/code/misc/plot_src.R')
G <- 5
library(this.path)

setwd(paste0(this.path::here(), '/../results/'))
load('oxy_none__5_18_18_30_30-matern_spheretime_warp_seasonality/final_results.RData')

df <- data.frame(iter = 0:20, 
                 Oxygen= 
                   estimated_model$diagnostics$measurement_error_response,
                 Salinity = estimated_model$diagnostics$measurement_error_predictors[,1],
                 Temperature = estimated_model$diagnostics$measurement_error_predictors[,2])
ggplot(data = df %>% tidyr::pivot_longer(cols = c('Oxygen', 'Salinity', 'Temperature')),
       aes(x = iter, y = value)) +
  geom_point() +
  geom_line() + 
  facet_grid(name~., scales = 'free_y') + 
  scale_y_continuous(limits = c(0, NA)) + 
  labs(x = 'EM iteration', y = 'Measurement Error Variance') +
  theme(text = element_text(family = 'Arial'))
ggsave('../images/meas_errors.png', height = 6, width = 6)

n_pcs <- c(12, 14, 16, 18, 20)
AIC_values <- rep(0, length(n_pcs))

for (i in 1:length(n_pcs)) {
  load(paste0('oxy_none__5_', n_pcs[i], '_', n_pcs[i], '_30_30-matern_spheretime_warp_seasonality/final_results.RData'))
  AIC_values[i] <- tail(estimated_model$AIC, 1)
}

ggplot(data = data.frame(n_pcs = n_pcs, AIC = AIC_values),
       aes(x = n_pcs, y = AIC)) +
  geom_point() + 
  labs(x = 'Number of predictor/response\nprincipal components') + 
  theme(text = element_text(family = 'Arial'))
ggsave('../images/AIC.png', height = 4, width = 6)


parameters <- estimated_model[['parameters']]
data_inputs <- estimated_model[['data_inputs']]
params <- estimated_model[['params']]
# reorder clusters based on temperature at 1000 dbar
Phi_pred <- bdiag(eval.basis(params[['basis_predictors']][[1]],evalarg = 1000),
                  eval.basis(params[['basis_predictors']][[2]],evalarg = 1000))
p_seq <- seq(0, 2000, length.out = 1000)
Phi_test <- eval.basis(p_seq, params$basis_response)
test_seq <- seq(15, 365, by = 91.5)

results <- list()
for (g in 1:params[['G']]) {
  results[[g]] <- list()
  
  covariates <- c(1, 0, 0, 0, 0, 0, 0)
  results[[g]][[1]] <- rbind(data.frame(pressure = p_seq,
                                        mean_value = Phi_test %*% rowSums(parameters[['means_predictors']][[g]][1:64,] %*% diag(covariates)),
                                        cluster = g, month = factor('avg', levels = c(month.name, 'avg')),
                                        variable = 'Temperature (°C)'),
                             data.frame(pressure = p_seq,
                                        mean_value = Phi_test %*% rowSums(parameters[['means_predictors']][[g]][-c(1:64), ] %*% diag(covariates)),
                                        cluster = g, month = factor('avg', levels = c(month.name, 'avg')),
                                        variable = 'Salinity (PSU)'),
                             data.frame(pressure = p_seq,
                                        mean_value = Phi_test %*% rowSums(parameters[['means_response']][[g]] %*% diag(covariates)),
                                        cluster = g, month = factor('avg', levels = c(month.name, 'avg')),
                                        variable = 'Oxygen (μmol/kg)'))
  
  for (i in 1:length(test_seq)) {
    covariates <- c(1,
                    sin(test_seq[i]/365.25 * 2 * pi * 1),
                    cos(test_seq[i]/365.25 * 2 * pi * 1),
                    sin(test_seq[i]/365.25 * 2 * pi * 2),
                    cos(test_seq[i]/365.25 * 2 * pi * 2),
                    sin(test_seq[i]/365.25 * 2 * pi * 3),
                    cos(test_seq[i]/365.25 * 2 * pi * 3))
    results[[g]][[i+1]] <- rbind(data.frame(pressure = p_seq,
                                            mean_value = Phi_test %*% rowSums(parameters[['means_predictors']][[g]][1:64,] %*% diag(covariates)),
                                            cluster = g, month = factor(c(month.name[c(1, 4, 7, 10)][i]), levels = c(month.name, 'avg')),
                                            variable = 'Temperature (°C)'),
                                 data.frame(pressure = p_seq,
                                            mean_value = Phi_test %*% rowSums(parameters[['means_predictors']][[g]][-c(1:64), ] %*% diag(covariates)),
                                            cluster = g, month = factor(c(month.name[c(1, 4, 7, 10)][i]), levels = c(month.name, 'avg')),
                                            variable = 'Salinity (PSU)'),
                                 data.frame(pressure = p_seq,
                                            mean_value = Phi_test %*% rowSums(parameters[['means_response']][[g]] %*% diag(covariates)),
                                            cluster = g, month = factor(c(month.name[c(1, 4, 7, 10)][i]), levels = c(month.name, 'avg')),
                                            variable = 'Oxygen (μmol/kg)'))
  }
}

t_vals <- sapply(1:params[['G']], function(g){ (Phi_test %*% parameters[['means_predictors']][[g]][1:64,])[p_seq == 2000,1]})
new_order <- order(order(t_vals))

df_means <- dplyr::bind_rows(results) %>%
  mutate(cluster_use = new_order[cluster])

theme_set(theme_bw() + theme(text = element_text(family = 'Arial')))
ggplot(data = df_means %>% mutate(cluster_final = paste0('Cluster: ', cluster_use)) %>%
         filter(month != 'avg'), #%>% filter(variable == 'Oxygen (μmol/kg)'),
       aes(x = mean_value, y = pressure, color = month, linetype = month))+
  geom_path() +
  facet_grid(cluster_final~variable, scales = 'free_x') +
  scale_y_continuous(trans = 'reverse') +
  labs(x = 'Mean value', y = 'Pressure (decibars)', color = 'Month', linetype = 'Month') +
  theme(legend.position = 'right', axis.text.x = element_text(angle = 45, vjust = .8),
        axis.title = element_blank())
ggsave('../images/season_means_all.png',
       height = 8.86/1.55, width =8.3/1.55)


ggplot(data = df_means %>% filter(month == 'avg')) + 
  geom_path(size = 1.1,aes(x = mean_value, y = pressure, group = factor(cluster_use), 
                           color = factor(cluster_use), linetype = factor(cluster_use))) + 
  facet_grid(~variable, scales = 'free_x') +
  #coord_flip() + 
  scale_y_continuous(trans = 'reverse')+
  labs(y = 'Pressure (decibars)', x = 'Time-averaged cluster mean functions',
       linetype = 'Cluster', color = 'Cluster')+
  theme(legend.position = 'right', axis.title.x = element_blank(), text = element_text(size = 18))
ggsave('../images/means_all.png',
       height = 6.5/1.22, width = 11/1.22)

library(tidyr)
df_all <- df_means %>%
  filter(month == 'avg') %>%
  pivot_wider(values_from  = 'mean_value', names_from = 'variable') %>%
  mutate(abs_psal = gsw::gsw_SA_from_SP(SP = `Salinity (PSU)`, p = pressure, longitude = 0, latitude = -50),
         cons_temp = gsw::gsw_pt_from_t(SA= abs_psal, t = `Temperature (°C)`, p = pressure, p_ref = 0)) %>%
  arrange(g, pressure, month)

range_SA <- range(df_all$abs_psal)
range_CT <- range(df_all$cons_temp)
grid_TS <- expand.grid(seq(range_SA[1], range_SA[2], length.out = 400),
                       seq(range_CT[1], range_CT[2], length.out = 400))
density <- rep(NA, nrow(grid_TS))

for (i in 1:nrow(grid_TS)) {
  density[i] <- gsw::gsw_rho_t_exact(grid_TS[i,1], t = grid_TS[i,2], p = 0)
}

ggplot() +
  geom_contour(data = data.frame(grid_TS, density), aes(x = Var1, y = Var2, z = density),
               color = 'gray75', breaks = seq(1023, 1029, by = .5)) + 
  geom_text(data = data.frame(salinity = 35.835, pot_t = c(1.5, 4.3, 6.35, 8.1, 9.6, 11.1, 12.4, 13.65, 14.8, 15.85, 17, 18, 19)[seq(1, 13, by = 2)], 
                              label = seq(1028.5, 1025.5, by = -.5) - 1000), 
            aes(x = salinity, y = pot_t, label = format(label, digits = 4)), size = 3,
            color = 'gray75', family = 'Arial') + 
  geom_path(data = df_all ,
            aes(x = abs_psal, y = cons_temp, color = `Oxygen (μmol/kg)`, group = paste(cluster_use, month))) + 
  geom_text(data = df_all %>% filter(pressure == min(pressure)),
            aes(x = abs_psal, y = cons_temp, label = cluster_use), nudge_x = -.015, family = 'Arial') +
  geom_point(data = df_all %>% filter(pressure %in% c(500, 1000, 1500) ),
             aes(x = abs_psal, y = cons_temp), size = .5) +
  scale_color_viridis_c() +
  theme(panel.grid = element_blank(), text = element_text(family = 'Arial')) +
  labs(x = 'Absolute Salinity', y = 'Potential Temperature (°C)', 
       color = 'Oxygen\n(μmol/kg)')
ggsave('../images/TS_plot_cluster.png', width = 6, height = 3.3)


ggplot(data = df_means %>% mutate(cluster_final = paste0('Cluster: ', cluster_use)) %>% filter(pressure < 400, month != 'avg'), #%>% filter(variable == 'Oxygen (μmol/kg)'),
       aes(x = mean_value, y = pressure, color = month, linetype = month))+
  geom_path() +
  facet_grid(cluster_final~variable, scales = 'free_x') +
  scale_y_continuous(trans = 'reverse') +
  labs(x = 'Mean value', y = 'Pressure (decibars)', color = 'Month', linetype = 'Month') +
  theme(legend.position = 'right', axis.text.x = element_text(angle = 45, vjust = .8),
        axis.title = element_blank())
ggsave('../images/season_means_all_reduced.png',
       height = 8.86/1.55, width =8.3/1.55)

ggplot(data = data.frame(day = data_inputs[['day']], 
                         date = as.Date(data_inputs[['day']], origin = as.Date('2020-01-01')),
                         mems = paste0('Cluster: ', new_order[parameters[['cluster_membership']]]), 
                         type=  data_inputs[['BGC']])) +
  geom_histogram(aes(x = date), size = .02, alpha = .6, breaks = as.Date(seq(0, 365.25, by = 14), origin = as.Date('2020-01-01'))) +
  facet_grid(~mems)  +
  labs(y = 'Number of profiles', x= 'Day of year') + 
  scale_x_date(date_labels = '%B', breaks = as.Date(seq(0, 366.25, length.out = 5)[1:4], origin = as.Date('2020-01-01'))) + 
  theme(axis.text.x = element_text(angle = 45, vjust = .85, hjust = .9))
ggsave('../images/dayofyear_hist.png',height = 2.2, width = 7)


label_use <- paste('Oxygen function (μmol/kg)')
df_pcs <- data.frame('pressure' = p_seq, 
                     'pc' = as.double(sapply((1:params[['G']]), function(g){ Phi_test %*% parameters[['Omega_response1']][[g]][,1]})),
                     'g' = factor(rep(new_order[(1:params[['G']])], each = length(p_seq))))
ggplot(data = df_pcs, aes(x = pressure, y = pc, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = label_use,
       linetype = 'Cluster', color = 'Cluster',
       title = 'First transformation functions')+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave('../images/var_pcs_response1.png',
       height = 8.86, width = 5.11)

df_pcs <- data.frame('pressure' = p_seq, 
                     'pc' = as.double(sapply((1:params[['G']]), function(g){ Phi_test %*% parameters[['Omega_response2']][[g]][,1]})),
                     'g' = factor(rep(new_order[(1:params[['G']])], each = length(p_seq))))
ggplot(data = df_pcs, aes(x = pressure, y = pc, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = "Residual function (μmol/kg)",
       linetype = 'Cluster', color = 'Cluster',
       title = paste('First residual PC functions'))+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave('../images/var_pcs_response2.png',
       height = 8.86, width = 5.11)

nbasis_pred <- 64
df_pcs <- data.frame('pressure' = p_seq, 
                     'pc' = as.double(sapply((1:params[['G']]), function(g){ Phi_test %*% parameters[['Omega_predictors']][[g]][1:nbasis_pred,1]})),
                     'g' = factor(rep(new_order[(1:params[['G']])], each = length(p_seq))))
ggplot(data = df_pcs, aes(x = pressure, y = pc, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = 'Temperature PC (°C)',
       linetype = 'Cluster', color = 'Cluster',
       title = 'First temperature PC functions')+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave('../images/var_pcs_predictors1.png',
       height = 8.86, width = 5.11)

df_pcs <- data.frame('pressure' = p_seq, 
                     'pc' = as.double(sapply((1:params[['G']]), function(g){ Phi_test %*% 
                         parameters[['Omega_predictors']][[g]][(nbasis_pred + 1):(2*nbasis_pred),1]})),
                     'g' = factor(rep(new_order[(1:params[['G']])], each = length(p_seq))))
ggplot(data = df_pcs, aes(x = pressure, y = pc, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = 'Salinity PC (PSU)',
       linetype = 'Cluster', color = 'Cluster',
       title = 'First salinity PC functions')+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave('../images/var_pcs_predictors2.png',
       height = 8.86, width = 5.11)

G <- params[['G']]
library(Matrix)
p_seq <- seq(0, 2000, length.out = 400)
Phi_test <- eval.basis(p_seq, params$basis_response)
cov_df <- bind_rows(lapply(1:G, function(g) {
  cov_mat <- Phi_test %*% parameters$Omega_response1[[g]] %*% diag(parameters$variances_predictors[[g]]) %*%
    t(parameters$Omega_predictors[[g]]) %*% t(bdiag(Phi_test, Phi_test))
  
  cov_mat_var_oxy <- Phi_test %*% parameters$Omega_response1[[g]] %*% diag(parameters$variances_predictors[[g]]) %*%
    t(Phi_test %*% parameters$Omega_response1[[g]]) + 
    Phi_test %*% parameters$Omega_response2[[g]] %*% diag(parameters$variances_response[[g]]) %*%
    t(Phi_test %*% parameters$Omega_response2[[g]])
  cov_mat_var_oxy_diag <- diag(1/sqrt(diag(cov_mat_var_oxy)))
  cov_mat_var_pred <- bdiag(Phi_test, Phi_test) %*% parameters$Omega_predictors[[g]] %*% diag(parameters$variances_predictors[[g]]) %*%
    t(parameters$Omega_predictors[[g]]) %*% t(bdiag(Phi_test, Phi_test))
  cov_mat_var_pred_diag <- diag(1/sqrt(diag(cov_mat_var_pred)))
  
  cov_mat_pred_norm <-  cov_mat_var_pred_diag %*% cov_mat_var_pred %*% cov_mat_var_pred_diag
  cov_mat_norm <-  cov_mat_var_oxy_diag %*% cov_mat %*% cov_mat_var_pred_diag
  
  
  data.frame(p_oxy = p_seq, p_pred = rep(p_seq, each = length(p_seq)), 
             var1 = rep(c('Oxygen', 'Oxygen', 'Temperature', 
                          'Oxygen', 'Temperature', 'Salinity'), each = length(p_seq)^2),
             predictor = rep(c('Temperature', 'Salinity', 'Salinity',
                               'Oxygen', 'Temperature', 'Salinity'), each = length(p_seq)^2), 
             cluster = new_order[g], 
             value_norm = c(as.vector(cov_mat_norm),
                            as.vector(cov_mat_pred_norm[1:(nrow(Phi_test)), - c(1:nrow(Phi_test))]),
                            as.vector(cov_mat_var_oxy_diag %*% cov_mat_var_oxy %*% cov_mat_var_oxy_diag),
                            as.vector(cov_mat_pred_norm[1:(nrow(Phi_test)), 1:(nrow(Phi_test))]),
                            as.vector(cov_mat_pred_norm[-c(1:(nrow(Phi_test))), -c(1:nrow(Phi_test))])))
}))

ggplot() +
  facet_grid(`Var 1` + `Var 2`~Cluster, labeller = label_both) + 
  scale_fill_gradient2(name = 'Correlation', 
                       low = scales::muted("blue"),
                       mid = "white",
                       high = scales::muted("red")) + 
  scale_y_continuous(trans = 'reverse') +
  coord_equal() +
  geom_raster(data = filter(cov_df, var1 != predictor) %>%
                rename(`Var 1`=var1, `Var 2` = predictor, Cluster = cluster), 
              aes(x= p_pred, y= p_oxy,
                  fill = value_norm)) + 
  labs(x = 'Pressure of Var 2 (decibars)', y = 'Pressure of Var 1 (decibars)')+
  theme(legend.position = 'right')
ggsave('../images/cor_all.png', height = 5.2, width = 10)

ggplot() +
  facet_grid(`Variable`~Cluster, labeller = label_both) + 
  scale_fill_gradient2(name = 'Correlation', 
                       low = scales::muted("blue"),
                       mid = "white",
                       high = scales::muted("red")) + 
  scale_y_continuous(trans = 'reverse') +
  coord_equal() +
  geom_raster(data = filter(cov_df, var1 == predictor) %>%
                rename(`Variable`=var1,  Cluster = cluster), 
              aes(x= p_pred, y= p_oxy,
                  fill = value_norm)) + 
  labs(x = 'Pressure (decibars)', y = 'Pressure (decibars)')+
  theme(legend.position = 'right')
ggsave('../images/cor_marg_all.png', height = 5.2, width = 10)
