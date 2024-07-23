compute_u_matrices <- function(Phi_x_Phi_resp, weights, cmi, cmi_BGC, n_e_step_samples, G,
                               pcs_r, Phi_x_Phi_pred = NULL, pcs_p1=NULL, pcs_p2=NULL, ind = F, 
                               seasonal_mean = F, num_seasonal_basis = 3, days = NULL, 
                               is_BGC = rep(T, nrow(cmi))){
  u_mats <- list()
  if (seasonal_mean) {
    covariates <- list(rep(1, length(days)), 
                       sin(days/365.25 * 2 * pi * 1),
                       cos(days/365.25 * 2 * pi * 1), 
                       sin(days/365.25 * 2 * pi * 2), 
                       cos(days/365.25 * 2 * pi * 2), 
                       sin(days/365.25 * 2 * pi * 3), 
                       cos(days/365.25 * 2 * pi * 3), 
                       sin(days/365.25 * 2 * pi * 4), 
                       cos(days/365.25 * 2 * pi * 4), 
                       sin(days/365.25 * 2 * pi * 5), 
                       cos(days/365.25 * 2 * pi * 5), 
                       sin(days/365.25 * 2 * pi * 6), 
                       cos(days/365.25 * 2 * pi * 6))
    total_basis <- 1 + 2 * num_seasonal_basis
    
    u_mats[['mean_resp_season']] <- list()
    for (k in 1:total_basis) {
      u_mats[['mean_resp_season']][[k]] <- list()
      for (k_star in 1:k) {
        u_mats[['mean_resp_season']][[k]][[k_star]] <-
          lapply(fstmr:::c_create_summed_U_matrix_pcs_sparse_for_season(Phi_x_Phi_resp, cmi_BGC, as.numeric(weights), 
                                                                as.numeric(covariates[[k]][is_BGC]), 
                                                                as.numeric(covariates[[k_star]][is_BGC]),
                                                                as.integer(G),
                                                                n_e_step_samples), as, 'sparseMatrix')
        if (k_star != k) {
          u_mats[['mean_resp_season']][[k_star]][[k]] <- u_mats[['mean_resp_season']][[k]][[k_star]]
        }
      }
    }
    u_mats[['mean_resp']] <- list()
    for (g in 1:G) {
      u_mats[['mean_resp']][[g]] <- do.call(rbind, 
                                            lapply(u_mats[['mean_resp_season']], 
                                                   function(x) do.call(cbind, lapply(x, function(y) y[[g]]))))
    }
    u_mats[['pc_resp']] = lapply(1:dim(pcs_r)[3], function(q){
      pc_weights_r = pcs_r[,,q]
      Matrix::bdiag(
        fstmr:::c_create_summed_U_matrix_pcs_sparse(Phi_x_Phi_resp, cmi_BGC, as.numeric(weights), pc_weights_r, as.integer(G),
                                            n_e_step_samples)
      )
    })
    if(!is.null(Phi_x_Phi_pred)){
      u_mats[['mean_pred_season']] <- list()
      for (k in 1:total_basis) {
        u_mats[['mean_pred_season']][[k]] <- list()
        for (k_star in 1:k) {
          u_mats[['mean_pred_season']][[k]][[k_star]] <-
            lapply(c_create_summed_U_matrix_pcs_sparse_for_season(Phi_x_Phi_pred, cmi, as.numeric(weights), 
                                                                  as.numeric(covariates[[k]]), as.numeric(covariates[[k_star]]),
                                                                  as.integer(G),
                                                                  n_e_step_samples), as, 'sparseMatrix')
          if (k_star != k) {
            u_mats[['mean_pred_season']][[k_star]][[k]] <- u_mats[['mean_pred_season']][[k]][[k_star]]
          }
        }
      }
      
      u_mats[['mean_pred']] <- list()
      for (g in 1:G) {
        u_mats[['mean_pred']][[g]] <- do.call(rbind, 
                                              lapply(u_mats[['mean_pred_season']], 
                                                     function(x) do.call(cbind, lapply(x, function(y) y[[g]]))))
      }
      u_mats[['pc_pred']] = lapply(1:dim(pcs_p1)[3], function(q){
        pc_weights_p = pcs_p1[,,q]
        Matrix::bdiag(
          c_create_summed_U_matrix_pcs_sparse(Phi_x_Phi_pred, cmi, as.numeric(weights), pc_weights_p, as.integer(G),
                                              n_e_step_samples)
        )
      })
      if (!ind){
        u_mats[['lt']] = lapply(1:dim(pcs_p1)[3], function(q){
          pc_weights_p = pcs_p2[,,q]
          Matrix::bdiag(
            c_create_summed_U_matrix_pcs_sparse(Phi_x_Phi_resp, cmi_BGC, as.numeric(weights), pc_weights_p, as.integer(G),
                                                n_e_step_samples)
          )})
      }
    }
    
  } else {
    u_mats[['mean_resp']] = Matrix::bdiag(
      c_create_summed_U_matrix_sparse(Phi_x_Phi_resp, cmi_BGC, as.numeric(weights), as.integer(G),
                                      n_e_step_samples)
    )
    u_mats[['pc_resp']] = lapply(1:dim(pcs_r)[3], function(q){
      pc_weights_r = pcs_r[,,q]
      Matrix::bdiag(
        c_create_summed_U_matrix_pcs_sparse(Phi_x_Phi_resp, cmi_BGC, as.numeric(weights), pc_weights_r, as.integer(G),
                                            n_e_step_samples)
      )
    })
    if(!is.null(Phi_x_Phi_pred)){
      u_mats[['mean_pred']] = Matrix::bdiag(
        c_create_summed_U_matrix_sparse(Phi_x_Phi_pred, cmi, as.numeric(weights), as.integer(G),
                                        n_e_step_samples)
      )
      
      u_mats[['pc_pred']] = lapply(1:dim(pcs_p1)[3], function(q){
        pc_weights_p = pcs_p1[,,q]
        Matrix::bdiag(
          c_create_summed_U_matrix_pcs_sparse(Phi_x_Phi_pred, cmi, as.numeric(weights), pc_weights_p, as.integer(G),
                                              n_e_step_samples)
        )
      })
      if (!ind){
        u_mats[['lt']] = lapply(1:dim(pcs_p1)[3], function(q){
          pc_weights_p = pcs_p2[,,q]
          Matrix::bdiag(
            c_create_summed_U_matrix_pcs_sparse(Phi_x_Phi_resp, cmi_BGC, as.numeric(weights), pc_weights_p, as.integer(G),
                                                n_e_step_samples)
          )
        })
      }
    }
  }
  return(u_mats)
}


update_means = function(U, Phi_prof, profs, cmi, pcs_array1, Omegas1, weights,
                        pen_mat, G=1, Omegas2=NULL, pcs_array2=NULL, lambda=10^-5, 
                        seasonal_mean = F, num_seasonal_basis= NA, days = NULL){
  
  n_e_step_samples = as.integer(ncol(cmi))

  pcs1 = lapply(1:length(profs), function(x) matrix(pcs_array1[x,,], nrow = n_e_step_samples))
  
  pen_mat_all = Matrix::bdiag(replicate(G, pen_mat))
  
  if (is.null(Omegas2)){
    
    if (seasonal_mean) { 
      covariates <- list(rep(1, length(days)), 
                         sin(days/365.25 * 2 * pi * 1),
                         cos(days/365.25 * 2 * pi * 1), 
                         sin(days/365.25 * 2 * pi * 2), 
                         cos(days/365.25 * 2 * pi * 2), 
                         sin(days/365.25 * 2 * pi * 3), 
                         cos(days/365.25 * 2 * pi * 3), 
                         sin(days/365.25 * 2 * pi * 4), 
                         cos(days/365.25 * 2 * pi * 4), 
                         sin(days/365.25 * 2 * pi * 5), 
                         cos(days/365.25 * 2 * pi * 5), 
                         sin(days/365.25 * 2 * pi * 6), 
                         cos(days/365.25 * 2 * pi * 6))
      total_basis <- 1 + 2 * num_seasonal_basis
      
      v_mats <- list()
      for (k in 1:total_basis) {
        v_mats[[k]] <- c_create_summed_V_matrix_sparse_single_for_season(Phi_prof, profs, as.integer(G), 
                                                                         Omegas1, cmi, pcs1,
                                                                         as.numeric(covariates[[k]]), n_e_step_samples,
                                                                         as.numeric(weights))
      }
      res <- list()
      pen_mat_all = Matrix::bdiag(replicate(total_basis, pen_mat))
      for (g in 1:G) {
        V_g <- as.double(sapply(v_mats, function(x) x[,g]))
        res[[g]] <- matrix(solve(U[[g]] + lambda * pen_mat_all, V_g), nrow = dim(U[[g]])/total_basis, ncol = total_basis)
      }
      
      return(res)
    } else {
      v_mat = c_create_summed_V_matrix_sparse_single(Phi_prof, profs, as.integer(G), Omegas1, cmi,
                                                     pcs1, n_e_step_samples, as.numeric(weights))
      
      V <- as.double(v_mat)
      
      res = as.double(Matrix::solve(U + lambda * pen_mat_all, V))
      
      return(split(res, rep(1:G, each=length(res)/G)))
    }
    
  } else {
    if (seasonal_mean) {
      
      covariates <- list(rep(1, nrow(cmi)), 
                         sin(days/365.25 * 2 * pi * 1),
                         cos(days/365.25 * 2 * pi * 1), 
                         sin(days/365.25 * 2 * pi * 2), 
                         cos(days/365.25 * 2 * pi * 2), 
                         sin(days/365.25 * 2 * pi * 3), 
                         cos(days/365.25 * 2 * pi * 3), 
                         sin(days/365.25 * 2 * pi * 4), 
                         cos(days/365.25 * 2 * pi * 4), 
                         sin(days/365.25 * 2 * pi * 5), 
                         cos(days/365.25 * 2 * pi * 5), 
                         sin(days/365.25 * 2 * pi * 6), 
                         cos(days/365.25 * 2 * pi * 6))
      total_basis <- 1 + 2 * num_seasonal_basis
      pcs2 = lapply(1:length(profs), function(x) matrix(pcs_array2[x,,], nrow = n_e_step_samples))
      
      v_mats <- list()
      for (k in 1:total_basis) {
        v_mats[[k]] <- c_create_summed_V_matrix_sparse_for_season(Phi_prof, profs, as.integer(G), 
                                                                  Omegas1, Omegas2, cmi, pcs1, pcs2,
                                                                  as.numeric(covariates[[k]]), n_e_step_samples,
                                                                  as.numeric(weights))
      }
      res <- list()
      pen_mat_all = Matrix::bdiag(replicate(total_basis, pen_mat))
      for (g in 1:G) {
        V_g <- as.double(sapply(v_mats, function(x) x[,g]))
        res[[g]] <- matrix(solve(U[[g]] + lambda * pen_mat_all, V_g), nrow = dim(U[[g]])/total_basis, ncol = total_basis)
      }
      
      
      return(res)
    } else {
      pcs2 = lapply(1:length(profs), function(x) matrix(pcs_array2[x,,], nrow = n_e_step_samples))
      
      v_mat = c_create_summed_V_matrix_sparse(Phi_prof, profs, as.integer(G), Omegas1, Omegas2,
                                              cmi, pcs1, pcs2, n_e_step_samples, as.numeric(weights))
      
      V <- as.double(v_mat)
      
      res = as.double(Matrix::solve(U + lambda * pen_mat_all, V))
      return(split(res, rep(1:G, each=length(res)/G)))
    }
  }
}

update_principal_components <- function(u_mats, cmi, pcs_r1, profs, Phi_prof, Omegas1,
                                        means, weights, pen_mat, G=1,  pcs_r2 = NULL,
                                        Omegas2=NULL, lambda=10^-5,
                                        seasonal_mean = FALSE, days = NULL){
  
  # The sqeuence pcs/Omegas needs to match, i.e. if one updates the 
  # parameters$Omega_r1, then pcs_r1 = betas (Gamma %*% pcs_p)
  # Its always Omegas1 that is updated
  
  Q = dim(pcs_r1)[3]
  n_e_step_samples = as.integer(ncol(cmi))
  means_mat = as.matrix(do.call(cbind, means))
  
  pcs1 = lapply(1:length(profs), function(x) matrix(pcs_r1[x,,], nrow = n_e_step_samples))
  
  if (!is.null(Omegas2)){
    pcs2 = lapply(1:length(profs), function(x) matrix(pcs_r2[x,,], nrow = n_e_step_samples))
  }
  
  for (q in 1:Q) {
    
    pc_weights = pcs_r1[,,q]
    if (seasonal_mean) {
      if(is.null(Omegas2)){
        v_mat = c_create_summed_V_matrix_pcs_sparse_single_for_season(Phi_prof, profs, means, as.integer(G), Omegas1,
                                                                      cmi, pcs1, n_e_step_samples, weights, as.integer(q-1), as.numeric(days))
      } else {
        v_mat = c_create_summed_V_matrix_pcs_sparse_for_season(Phi_prof, profs, means, as.integer(G), Omegas1, Omegas2,
                                                               cmi, pcs1, pcs2, n_e_step_samples,
                                                               weights, as.integer(q-1), days)      
      }
    } else {
      if(is.null(Omegas2)){
        v_mat = c_create_summed_V_matrix_pcs_sparse_single(Phi_prof, profs, means_mat, as.integer(G), Omegas1,
                                                           cmi, pcs1, n_e_step_samples, weights, as.integer(q-1))
      } else {
        v_mat = c_create_summed_V_matrix_pcs_sparse(Phi_prof, profs, means_mat, as.integer(G), Omegas1, Omegas2,
                                                    cmi, pcs1, pcs2, n_e_step_samples,
                                                    weights, as.integer(q-1))      
      }
    }
    
    U <- u_mats[[q]]
    V <- as.double(v_mat)
    
    pen_mat_all = Matrix::bdiag(replicate(G, pen_mat))
    
    res = as.double(Matrix::solve(U + lambda * pen_mat_all, V))

    for (g in 1:G){
      Omegas1[[g]][,q] = split(res, rep(1:G, each = length(res)/G))[[g]]
    } 
  }
  return(Omegas1)
}


update_measurement_error <- function(cluster_mat, pcs_array1, profs, phi_prof, means, Omegas1,
                                     n_profiles, Omegas2=NULL, pcs_array2=NULL, n_vars = 1,
                                     weights = rep(1, ncol(cluster_mat)), G=1, prof_ns = NULL, 
                                     seasonal_mean = F, days = NULL) {
  
  cluster_mat_i = cluster_mat - 1
  mode(cluster_mat_i) = 'integer'
  means_mat = as.matrix(do.call(cbind, means))
  
  if(!is.null(Omegas2)){

    pcs1 = lapply(1:n_profiles, function(x) {
      matrix(pcs_array1[x,,], nrow = ncol(cluster_mat))
    })
    
    pcs2 = lapply(1:n_profiles, function(x) {
      matrix(pcs_array2[x,,], nrow = ncol(cluster_mat))
    })
    
    if (seasonal_mean) {
      c_update_measurement_error_for_season(cluster_mat_i, phi_prof, pcs1, pcs2, profs, means, Omegas1, Omegas2, n_profiles, weights, G, as.double(days))
    } else {
      c_update_measurement_error(cluster_mat_i, phi_prof, pcs1, pcs2, profs, means_mat, Omegas1, Omegas2, n_profiles, weights, G)
    }
    
  } else {
    
    errors = rep(0, n_vars)
    MC <- ncol(cluster_mat)
    all_prof_lengths_equal <- all(apply(prof_ns, 1, function(x) {diff(range(x)) < .Machine$double.eps ^ 0.5}))
    for (x in 1:n_profiles){
      if (seasonal_mean) {
        diff = c_compute_squared_sparse_for_season(cluster_mat_i[x,,drop= F], matrix(pcs_array1[x, ,], ncol = dim(pcs_array1)[3]),
                                                   profs[[x]], as.integer(G), phi_prof[[x]], means, Omegas1, weights, 
                                                   as.double(days[x]))
      } else {
        diff = c_compute_squared_sparse(cluster_mat_i[x,,drop= F], matrix(pcs_array1[x, ,], ncol = dim(pcs_array1)[3]),
                                        profs[[x]], as.integer(G), phi_prof[[x]], means_mat, Omegas1, weights)
      }
      if (all_prof_lengths_equal) {
        errors <- errors + colSums(matrix(ncol = n_vars, diff))
      } else {
        vec_indexes <- as.double(
          sapply(1:length(prof_ns[x,]), function(i) rep(i, times=prof_ns[x,i]))
        )
        errors <- errors + tapply(diff, vec_indexes, sum)
      }
    }
    return(errors/(sum(weights)*colSums(prof_ns)))
  }
}

update_cluster_membership <- function(cluster_mat, n_profiles, G){
  conditional_probs = matrix(0, nrow=n_profiles, ncol=G)
  for(g in 1:G){
    conditional_probs[,g] = rowMeans(cluster_mat==g)
  }	
  cluster_membership = apply(conditional_probs, 1, which.max)
  return(list(cluster_membership, conditional_probs))
}

# Gauss Markov Random Field
# Optimize functions for MRF (pretty much taken from Liang et al)s
update_theta <- function(interval, nn_list, cluster_mat, n_profiles, G=1, theta, weights= rep(1, ncol(cluster_mat)), dists = NULL){
  tryCatch(expr = {uniroot(gradient_mrf, interval, nn_list = nn_list, weights = weights,
                           cluster_mat = cluster_mat, tol = .Machine$double.eps^0.15, n = n_profiles, G = G, dists = dists)$root },
           error = function(x) {print('MRF optim failed');theta})
}

gradient_mrf <- function(theta, cluster_mat, nn_list, n_profiles, G, weights = rep(1, ncol(cluster_mat)), dists = NULL) {
  MC <- ncol(cluster_mat)
  tmp <- matrix(0, n_profiles, G)
  for(x in 1:n_profiles){
    nn <- length(nn_list[[x]])
    for (mc in 1:MC) {
      if(!is.null(dists)){
        n_each = sapply(1:G, function(g) sum(1 / (dists[[x]][cluster_mat[nn_list[[x]],mc]==g])))
      } else {
        n_each = sapply(1:G, function(g) sum(cluster_mat[nn_list[[x]],mc]==g))
      }
      denominator <-  exp(n_each * theta)
      numerator <- denominator * n_each
      if (sum(numerator) == Inf) {
        next
      }
      for(g in 1:G){
        tmp[x,g] <-tmp[x,g] + 
          weights[mc] * (cluster_mat[x,mc] == g)*(n_each[g] -sum(numerator)/sum(denominator))
      }
    }
  }
  return(sum(tmp))
}


update_variances <- function(cluster_mat, pcs_array, range_params,
                             G, n_pcs, MC, cov_fun, independent_model = F, 
                             Linv, variances, vecchia_info) {
  lapply(1:G, function(g) {
    sapply(1:n_pcs, function(q) {
      mean(sapply(1:MC, function(mc) {
        indexes <- cluster_mat[,mc] == g
        pcs_use <- as.double(pcs_array[indexes,mc,q])
        if (independent_model) {
          return(sum(pcs_use^2)/sum(indexes))
        } else {
          return(sum(as.double((Linv[[mc]][[g]][[q]] *sqrt(variances[[g]][[q]])) %*% pcs_use[vecchia_info[[mc]][['order']][[g]]])^2)/sum(indexes))
        }
      }))
    })
  })
}

# spatial range parameters
update_range_params <- function(G, n_pcs, variances, cluster_mat, pcs_array,
                                vecchia_info,
                                mc_weights = rep(1, ncol(cluster_mat)),
                                covfun_name = 'exponential_spheretime',
                                maxit = 10, start_params_previous = NULL) {
  if (covfun_name == 'exponential_spheretime_warp'){
    start_params <- log(c(.2, 25, exp(c(0,0,0,0,0))))
    lower_params <- log(c(.001, .5, exp(c(-8,-8,-8,-8,-8))))
    upper_params <- log(c(4, 800, exp(c(8,8,8,8,8))))
  } else {
    start_params <- log(c(.2, 25))
    lower_params <- log(c(.000001, .1))
    upper_params <-  log(c(4, 1200))
  }
  no_cluster_variation <- all(apply(cluster_mat, 1, function(x) var(x) == 0))
  # for each cluster and principal component, use optim
  vals <- lapply(1:G, function(g) {
    sapply(1:n_pcs, function(q) {
      if (!is.null(start_params_previous)) {
        if (covfun_name == 'exponential_spheretime_warp') {
          start_params <- c(log(start_params_previous[[g]][1:2,q]),
                            start_params_previous[[g]][-c(1:2),q])
        } else {
          start_params <- log(start_params_previous[[g]][,q])
        }
      }
      # theta, var_score, cluster_mat, pcs_array, q, g2, MC_keep,
      # mc_weights, vecchia_info, covfun_name = 'exponential_spheretime',
      # no_cluster_variation
      return_vals <- optim(method = 'L-BFGS-B', range_optim_vecchia, par = start_params,
                           var_score = variances[[g]][q],
                           cluster_mat = cluster_mat, pcs_array = pcs_array, control = list(factr = 1e7, maxit = maxit),
                           mc_weights = mc_weights, 
                           q = q, g2 = g, MC_keep = dim(cluster_mat)[2], lower = lower_params, upper = upper_params,
                           covfun_name = covfun_name, 
                           vecchia_info = vecchia_info, no_cluster_variation = no_cluster_variation)$par
      return_vals[1:2] <- exp(return_vals[1:2])
      return_vals
    })
  })
  vals
}

# spatial range parameter for 1 pc and 1 cluster
range_optim_vecchia <- function(theta, var_score, cluster_mat, pcs_array, q, g2, MC_keep,
                                mc_weights, vecchia_info, covfun_name = 'exponential_spheretime',
                                no_cluster_variation) {
  scale <- exp(theta[1:2])
  if (covfun_name == 'exponential_spheretime') {
    cov_params_use <-  c(var_score, scale[1], scale[2], .00001)
  } else if (covfun_name == 'exponential_spheretime_warp') {
    cov_params_use <-  c(var_score, scale[1], scale[2], .00001, theta[3:7])
  } else if (covfun_name == 'matern_spheretime') {
    cov_params_use <-  c(var_score, scale[1], scale[2], .5, .00001)
  } else if (covfun_name == 'exponential_isotropic'){
    cov_params_use <-  c(var_score, scale[1], .00001)
  }
  
  like_comp <- tryCatch({
    
    val <- rep(0, MC_keep)
    for (mc in 1:MC_keep) {
      if (mc == 1) {
        Linv_temp <- GpGp::vecchia_Linv(cov_params_use, covfun_name = covfun_name,
                                        locs = vecchia_info[[mc]][['locs_list']][[g2]], 
                                        NNarray = vecchia_info[[mc]][['NN_list']][[g2]])
      } else if (no_cluster_variation) {
        
      } else if (!isTRUE(all.equal(cluster_mat[,mc], cluster_mat[,mc-1]))){
        Linv_temp <- GpGp::vecchia_Linv(cov_params_use, covfun_name = covfun_name,
                                        locs =vecchia_info[[mc]][['locs_list']][[g2]], 
                                        NNarray = vecchia_info[[mc]][['NN_list']][[g2]])
      }
      order_use <- vecchia_info[[mc]][['order']][[g2]]
      vec_to_sq <- GpGp::Linv_mult(Linv_temp, pcs_array[cluster_mat[,mc] == g2,mc,q][order_use],
                                   vecchia_info[[mc]][['NN_list']][[g2]])
      val[mc] <- mc_weights[mc]*(.5 *  sum(vec_to_sq^2) -
                                   sum(log(Linv_temp[,1])) + .5*nrow(Linv_temp) * log(2 * pi))
    }
    val}, 
    error = function(z) {
      message(z)
      NA
    }) 
  (as.double(0.5 * sum(like_comp)))
}

update_range_params_fisher <- function(G, n_pcs, variances, cluster_mat, pcs_array,
                                       vecchia_info,
                                       mc_weights = rep(1, ncol(cluster_mat)),
                                       covfun_name = 'exponential_spheretime',
                                       start_params_previous = NULL,
                                       maxit = 15,
                                       use_box_constraints = F,
                                       nugget = F) {

  no_cluster_variation <- all(apply(cluster_mat, 1, function(x) var(x) == 0))
  
  mc = which.max(mc_weights)

  vals <- lapply(1:G, function(g) {

    order_use <- vecchia_info[[mc]][['order']][[g]]
    locs_list <- vecchia_info[[mc]][['locs_list']][[g]]
    NN_list <- vecchia_info[[mc]][['NN_list']][[g]]

    # Transform lat/lon if necessary
    # Should be correct already from the E step
    # locs_list = get_locs(covfun_name, locs_list)
    
    sapply(1:n_pcs, function(q) {
      
      # Reset
      
      if(is.null(start_params_previous)){
        cov_params_use = init_cov_params(covfun_name, variances[[g]][q], start_params_previous, nugget)
      } else {
        cov_params_use = init_cov_params(covfun_name, variances[[g]][q], start_params_previous[[g]][,q], nugget)
      }
      
      likobj = GpGp::vecchia_meanzero_loglik_info(cov_params_use,
                                                  covfun_name,
                                                  pcs_array[cluster_mat[,mc] == g,mc,q][order_use],
                                                  locs_list,
                                                  NN_list)
      old_lik = likobj$loglik
      conv_cnt = 0
      
      for (i in 1:maxit){
        
        update_ind = get_update_ind(covfun_name, nugget)
        
        result = get_newton_raphson_step(covfun_name, cov_params_use, likobj, update_ind, use_box_constraints)
        step = result[[1]]
        update_ind = result[[2]]
        
        # Make sure we stay with permissible parameters, move along step (tends to work better than gradient)
        res = update_after_check_cov_parms(covfun_name, cov_params_use, step, use_box_constraints, nugget)
        newparms = res[[1]]
        step = res[[2]]
        
        likobj = GpGp::vecchia_meanzero_loglik_info(newparms,
                                                    covfun_name,
                                                    pcs_array[cluster_mat[,mc] == g,mc,q][order_use],
                                                    locs_list,
                                                    NN_list)
        
        if(likobj$loglik < old_lik){
          use_gradient = TRUE
        } else {
          use_gradient = FALSE
        }
        
        if( use_gradient ){
          # Somehow moving along the newton raphson step often works much better.
          result = linesearch(cov_params_use, old_lik, step,
                              pcs_array[cluster_mat[,mc] == g,mc,q][order_use],
                              covfun_name, locs_list, NN_list, use_box_constraints, nugget)
          
          newparms = result[[1]]
          step = result[[2]]
          
          # If Newton-Raphson step still doesn't work, try gradient
          if(min(newparms[update_ind] == cov_params_use[update_ind]) > 0){
            result = prepare_grad(covfun_name, cov_params_use, likobj$grad, update_ind, use_box_constraints, nugget)
            step = result[[1]]
            update_ind = result[[2]]
            
            result = linesearch(cov_params_use, old_lik, step,
                                pcs_array[cluster_mat[,mc] == g,mc,q][order_use],
                                covfun_name, locs_list, NN_list, use_box_constraints, nugget)
            
            newparms = result[[1]]
            step = result[[2]]
            
          }
          
          if(min(newparms == cov_params_use) > 0){
            break
          }
          
          likobj = GpGp::vecchia_meanzero_loglik_info(newparms,
                                                      covfun_name,
                                                      pcs_array[cluster_mat[,mc] == g,mc,q][order_use],
                                                      locs_list,
                                                      NN_list)
        } 
        
        if(max(abs((newparms[update_ind] - cov_params_use[update_ind]) / cov_params_use[update_ind])) < .01){
          conv_cnt = conv_cnt + 1
        }
        # if (abs(old_lik - likobj$loglik) < 0.5){
        #   conv_cnt = conv_cnt + 1
        # }
        if(conv_cnt > 2){
          break
        }
        cov_params_use = newparms
        old_lik = likobj$loglik
        
        # If parameter values are almost equal to box constraints, set them equal
        # to avoid small step sizes
        if(use_box_constraints){
          cov_params_use = set_equal_to_box(covfun_name, cov_params_use)
        }
        if(i == maxit){
          print("Reached maxit during range parameter estimation.")
        }
      }
      # Remove the variance which is estimated separately
      cov_params_use[2:length(cov_params_use)]
    })
  })
  vals
}

set_equal_to_box <- function(covfun_name, params){
  box_ind = get_box_ind(covfun_name)
  res = init_lower_upper_box(covfun_name, nugget)
  lower = res[[1]]
  upper = res[[2]]
  
  for(i in 1:length(box_ind)){
    if(abs(params[box_ind[i]] - upper[i]) < 1e-6){
      params[box_ind[i]] = upper[i]
    }
    if(abs(params[box_ind[i]] - lower[i]) < 1e-6){
      params[box_ind[i]] = lower[i]
    }
  }
  return(params)
}

get_locs <- function(covfun_name, locs, convert_to_xyz = F){
  if(convert_to_xyz ||
     covfun_name == 'exponential_scaledim' ||
     covfun_name == 'matern_scaledim'){
    x = cos( (locs[,1] + 180 ) * 2 * pi / 360 ) * sin( (locs[,2] + 90) * pi / 180 )
    y = sin( (locs[,1] + 180 ) * 2 * pi / 360 ) * sin( (locs[,2] + 90) * pi / 180 )
    z = cos( (locs[,2] + 90) * pi / 180 )
    if(ncol(locs) > 2){
      return(cbind(x,y,z,locs[,ncol(locs)]))
    } else {
      return(cbind(x,y,z))
    }
  }
  if(covfun_name == "matern_isotropic" || 
     covfun_name == "exponential_isotropic" ||
     covfun_name == 'matern_sphere' ||
     covfun_name == 'matern_spheretime' ||
     covfun_name == 'matern_spheretime_warp' ||
     covfun_name == 'exponential_spheretime' ||
     covfun_name == 'exponential_spheretime_warp'){
    return(locs)
  }
  stop("Invalid covariance function specified.")
}

get_newton_raphson_step <- function(covfun_name, params, likobj, update_ind, use_box_constraints){
  
  step = rep(0, length(params))
  
  # Following the GpGp package
  if (fstmr:::condition_number(likobj$info[update_ind, update_ind]) > 1 / 1e-5){
    ee <- eigen(likobj$info[update_ind, update_ind], symmetric = T)
    ee_ratios <- ee$values/max(ee$values)
    ee_ratios[ ee_ratios < 1e-5 ] <- 1e-5
    ee$values <- max(ee$values)*ee_ratios
    likobj$info[update_ind, update_ind] <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
  }
  
  if(max(abs(likobj$info[update_ind, update_ind])) < 1e-15){
    likobj$info = diag(0.01, ncol = ncol(likobj$info), nrow = nrow(likobj$info))
  }
  # First try
  step[update_ind] = solve(likobj$info[update_ind, update_ind], likobj$grad[update_ind])
  
  # Only consider 
  if(use_box_constraints){
    
    res = init_lower_upper_box(covfun_name, nugget)
    lower = res[[1]]
    upper = res[[2]]
    
    box_ind = get_box_ind(covfun_name)
    
    equal_box_ind_upper = box_ind[params[box_ind]==upper]
    equal_box_ind_lower = box_ind[params[box_ind]==lower]
    
    for(i in 1:length(update_ind)){
      if (update_ind[i] %in% equal_box_ind_upper && step[update_ind[i]] > 0){
        update_ind = update_ind[-i]
      }
      if (update_ind[i] %in% equal_box_ind_lower && step[update_ind[i]] < 0){
        update_ind = update_ind[-i]
      }
    }
    
    step = rep(0, length(params))
    
    if (fstmr:::condition_number(likobj$info[update_ind, update_ind]) > 1 / 1e-5){
      ee <- eigen(likobj$info[update_ind, update_ind], symmetric = T)
      ee_ratios <- ee$values/max(ee$values)
      ee_ratios[ ee_ratios < 1e-5 ] <- 1e-5
      ee$values <- max(ee$values)*ee_ratios
      likobj$info[update_ind, update_ind] <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
    }
    step[update_ind] = solve(likobj$info[update_ind, update_ind], likobj$grad[update_ind])
  }
  
  return(list(step, update_ind))
}

prepare_grad <- function(covfun_name, params, grad, update_ind, use_box_constraints, nugget){
  
  grad[-update_ind] = 0
  
  if(use_box_constraints){
    
    result = init_lower_upper_box(covfun_name, nugget)
    lower = result[[1]]
    upper = result[[2]]
    
    box_ind = get_box_ind(covfun_name)
    
    equal_box_ind_upper = box_ind[params[box_ind]==upper]
    equal_box_ind_lower = box_ind[params[box_ind]==lower]
    
    for(i in 1:length(update_ind)){
      if (update_ind[i] %in% equal_box_ind_upper && grad[update_ind[i]] > 0){
        update_ind = update_ind[-i]
      }
      if (update_ind[i] %in% equal_box_ind_lower && grad[update_ind[i]] < 0){
        update_ind = update_ind[-i]
      }
    }
  }
  grad[-update_ind] = 0
  return(list(grad, update_ind))
}

linesearch <- function(params, lik, step, y, covfun_name, locs_list, NN_list, use_box_constraints, nugget, maxiter = 30){  
  
  # Update but ensure that constraints are satisfied
  newparms = params
  check_result = update_after_check_cov_parms(covfun_name, params, step, use_box_constraints, nugget)
  newparms = check_result[[1]]
  step = check_result[[2]]
  
  new_lik = -1/0
  cnt = 1
  reached_maxit = FALSE
  
  while(new_lik < lik & cnt < maxiter){
    step = 0.2 * step  
    check_result = update_after_check_cov_parms(covfun_name, params, step, use_box_constraints, nugget)
    newparms = check_result[[1]]
    step = check_result[[2]]
    
    new_lik = GpGp::vecchia_meanzero_loglik(newparms,
                                            covfun_name,
                                            y,
                                            locs_list,
                                            NN_list)$loglik
    
    cnt = cnt + 1
  }
  if( cnt == maxiter ){
    print("Reached max iter in line search.")
    reached_maxit = TRUE
  }
  list(newparms, step, new_lik, reached_maxit)
}

condition_number <- function(info){
  # assumes that information matrix has finite numbers in it
  
  retval = tryCatch({if(max(diag(info))/min(diag(info)) > 1e6){
    return( max(diag(info))/min(diag(info)) )
  } else {
    ee <- eigen(info, symmetric = T)
    return( max(ee$values)/min(ee$values) )
  }},
  error = function(z){ print(info); 1/0 })
  return(retval)
}   

update_after_check_cov_parms <- function(covfun_name, params, step, use_box_constraints, nugget){

  newparms = params + step
  
  c_counter = 0
  c_counter = c_counter + basecheck_cov_params(covfun_name, newparms, nugget)
  c_counter = c_counter + additional_check_cov_params(covfun_name, newparms, nugget)
  if(c_counter == 2){
    gtg = T
  }
  if (use_box_constraints){
    c_counter = c_counter + boxcheck_range_cov_params(covfun_name, newparms, nugget)
    if(c_counter != 3){
      gtg = F
    }
  }
  
  err_counter = 1
  
  while( !gtg ){
    step <- 0.8 * step
    
    newparms = params + step
    
    c_counter = 0
    c_counter = c_counter + basecheck_cov_params(covfun_name, newparms, nugget)
    c_counter = c_counter + additional_check_cov_params(covfun_name, newparms, nugget)
    if(c_counter == 2){
      gtg = T
    } 
    
    if(use_box_constraints){
      c_counter = c_counter + boxcheck_range_cov_params(covfun_name, newparms, nugget)
      if(c_counter != 3){
        gtg = F
      }
    } 
    
    if(err_counter > 3000){
      stop("Problem with new params: update_after_check")
    }
    err_counter = err_counter + 1
  }
 return(list(newparms, step))
}

get_update_ind <- function(covfun_name, nugget){
  if(!nugget){
    if (covfun_name == 'matern_isotropic'){
      return(2:3)
    }
    if (covfun_name == 'exponential_isotropic'){
      return(2)
    }
    if(covfun_name == 'matern_scaledim'){
      return(2:6)
    }
    if(covfun_name == 'matern_spheretime'){
      return(2:3)
    }
    if (covfun_name == 'matern_spheretime_warp'){
      return(c(2:4, 6:10))
    }
    if (covfun_name == 'exponential_spheretime_warp'){
      return(c(2:3, 5:9))
    }
  } else {
    if (covfun_name == 'matern_isotropic'){
      return(2:4)
    }
    if (covfun_name == 'exponential_isotropic'){
      return(2:3)
    }
    if(covfun_name == 'matern_scaledim'){
      return(2:7)
    }
    if(covfun_name == 'matern_spheretime'){
      return(2:4)
    }
    if (covfun_name == 'matern_spheretime_warp'){
      return(2:10)
    }
    if (covfun_name == 'exponential_spheretime_warp'){
      return(2:9)
    }
  }
  stop("Invalid covfun_name: get_update_ind")
}
get_box_ind <- function(covfun_name){
  if (covfun_name == 'matern_isotropic' ||
      covfun_name == 'exponential_isotropic'){
    return(2)
  }
  if(covfun_name == 'matern_scaledim'){
    return(2:5)
  }
  if(covfun_name == 'matern_spheretime'){
    return(2:3)
  }
  if (covfun_name == 'matern_spheretime_warp'){
    return(c(2:3, 6:10))
  }
  if (covfun_name == 'exponential_spheretime_warp'){
    return(c(2:3, 5:9))
  }
}

# Need to allow parameters to be on the boundary
# Currently only checks range
boxcheck_range_cov_params <- function(covfun_name, params, nugget){
  box_ind = get_box_ind(covfun_name)
  res = init_lower_upper_box(covfun_name, nugget)
  lower = res[[1]]
  upper = res[[2]]
  return(sum(params[box_ind] < lower) + sum(params[box_ind] > upper) == 0)
}

init_cov_params <- function(covfun_name, variance, params=NULL, nugget=F){
  
  if (!is.null(params)){
    return(c(variance, params))
  } else {
    if (!nugget){
      if (covfun_name == 'matern_isotropic'){
        return(c(variance, .5, .5, 0))
      }
      if(covfun_name == 'exponential_isotropic'){
        return(c(variance, .5, 0))
      }
      if(covfun_name == 'matern_scaledim'){
        return(c(variance, .5, .5, .5, 25, 1.5, 0))
      }
      if(covfun_name == 'matern_spheretime'){
        return(c(variance, .5, 25, 0.5, 0))
      }
      if (covfun_name == 'exponential_spheretime_warp'){
        return(c(variance, .05, 25, 0, rep(0, 5)))
      }  
      if (covfun_name == 'matern_spheretime_warp'){
        return(c(variance, .05, 25, .5, 0, rep(0, 5)))
      }
    } else {
      if (covfun_name == 'matern_isotropic'){
        return(c(variance, .5, .5, 0.5*variance))
      }
      if(covfun_name == 'exponential_isotropic'){
        return(c(variance, .5, 0.5*variance))
      }
      if(covfun_name == 'matern_scaledim'){
        return(c(variance, .5, .5, .5, 25, 0.5, 0.5*variance))
      }
      if(covfun_name == 'matern_spheretime'){
        return(c(variance, .5, 25, 0.5, 0.5*variance))
      }
      if (covfun_name == 'exponential_spheretime_warp'){
        return(c(variance, .05, 25, 0, 0.5*variance, rep(0, 5)))
      }  
      if (covfun_name == 'matern_spheretime_warp'){
        return(c(variance, .05, 15, .5, 0.5*variance, rep(0, 5)))
      }
    }
  } 
  stop("Invalid covfun_name: init")
}

basecheck_cov_params <- function(covfun_name, params, nugget=F){
  
  if(nugget){
    if (covfun_name == 'matern_isotropic' ||
        covfun_name == 'exponential_isotropic' ||
        covfun_name == 'matern_scaledim' ||
        covfun_name == 'matern_spheretime'){
      return(sum(params > 0) == length(params))
    }
    if (covfun_name == 'exponential_spheretime_warp'){
      return(sum(params[1:4] > 0) == 4)
    }
    if (covfun_name == 'matern_spheretime_warp'){
      return(sum(params[1:5] > 0) == 5)
    }
  } else {
    if (covfun_name == 'exponential_isotropic'){
      return(sum(params[1:2] > 0) == 2)
    }
    if (covfun_name == 'matern_isotropic'){
      return(sum(params[1:3] > 0) == 3)
    }
    if(covfun_name == 'matern_scaledim'){
      return(sum(params[1:6] > 0) == 6)
    }
    if(covfun_name == 'matern_spheretime'){
      return(sum(params[1:4] > 0) == 4)
    }
    if (covfun_name == 'exponential_spheretime_warp'){
      return(sum(params[1:3] > 0) == 3)
    }
    if (covfun_name == 'matern_spheretime_warp'){
      return(sum(params[1:4] > 0) == 4)
    }
  }

  stop("Invalid covfun_name: basecheck")
}

additional_check_cov_params <- function(covfun_name, params, nugget = F){
  if (covfun_name == 'matern_isotropic' ||
      covfun_name == 'matern_scaledim'){
    # This needs to be changed -  we need to know the number of dimensions
    # This assumes d = 3
    return(params[6] < 7)
  }
  if (covfun_name == 'matern_spheretime' || 
      covfun_name == 'matern_spheretime_warp'){
    return(params[4] < 7)
  }
  TRUE
}

init_lower_upper_box <- function(covfun_name, nugget=F){
  
  if (covfun_name == 'matern_isotropic' ||
      covfun_name == 'exponential_isotropic'){
    lower = c(1e-8)
    upper = 365
    return(list(lower, upper))
  }
  if(covfun_name == 'matern_scaledim'){
    lower = c(1e-8, 1e-8, 1e-8, 0.5)
    upper = c(1, 1, 1,  365)
    return(list(lower, upper))
  }
  if(covfun_name == 'matern_spheretime'){
    lower = c(1e-8, 0.5)
    upper = c(1, 365)
    return(list(lower, upper))
  }
  if (covfun_name == 'matern_spheretime_warp' || 
      covfun_name == 'exponential_spheretime_warp'){
    lower = c(1e-8, .5, rep(-3, 5))
    upper = c(1, 545, rep(3, 5))
    return(list(lower, upper))
  }
  
  stop("Invalid covfun_name: boxcheck")
}

update_Lambda <- function(cluster_mat, response_pcs_array, predictor_pcs_array, G=1){
  # Transform to Matrix so that indexing is easier
  nrows = dim(response_pcs_array)[1] * dim(response_pcs_array)[2]
  response_pcs_array = matrix(response_pcs_array, nrow = nrows, ncol = dim(response_pcs_array)[3])
  nrows = dim(predictor_pcs_array)[1] * dim(predictor_pcs_array)[2]
  predictor_pcs_array = matrix(predictor_pcs_array, nrow = nrows, ncol = dim(predictor_pcs_array)[3])
  lapply(1:G, function(g){
    ind <- as.vector(cluster_mat == g)
    cross = crossprod(response_pcs_array[ind,], predictor_pcs_array[ind,])
    predictor_matrix = crossprod(predictor_pcs_array[ind,])
    tryCatch({ t(solve(predictor_matrix, t(cross)))},
             error = function(x) {
               t(solve(predictor_matrix + diag(nrow(predictor_matrix)), t(cross)))})
  })
}

orthogonalize <- function(Omegas1, variances1, inner_prod_mat1, G,
                          Omegas2=NULL, variances2=NULL, inner_prod_mat2=NULL,
                          Lambdas=NULL, cluster_mat=NULL, pcs1=NULL, pcs2=NULL){
  if(!is.null(cluster_mat)){
    ret = list()
    variances_response = list()
    variances_predictors = list()
    Sigma_eta_inv = list()
    # Transform to Matrix so that indexing is easier
    n_pcs_response = ncol(Omegas2[[1]])
    n_pcs_predictors = ncol(Omegas1[[1]])
    nrows = dim(pcs1)[1] * dim(pcs1)[2]
    response_pcs_array = matrix(pcs1, nrow = nrows, ncol = n_pcs_response)
    nrows = dim(pcs2)[1] * dim(pcs2)[2]
    predictors_pcs_array = matrix(pcs2, nrow = nrows, ncol = n_pcs_predictors)
    
    for (g in 1:G){
      ind <- as.vector(cluster_mat == g)
      #Should the mean be estimated as well?
      cov_estimate_response = 1/sum(ind)*crossprod(response_pcs_array[ind,])
      A = inner_prod_mat1 %*% Omegas1[[g]] %*%
        tcrossprod(cov_estimate_response, Omegas1[[g]])  %*% t(inner_prod_mat1)
      cov_estimate_predictors = 1/sum(ind)*crossprod(predictors_pcs_array[ind,])
      B = inner_prod_mat2 %*% Omegas2[[g]] %*%
        tcrossprod(cov_estimate_predictors, Omegas2[[g]])%*% 
        t(inner_prod_mat2)
      #eig_A = eigen(A, symmetric = T)
      #eig_B = eigen(B, symmetric = T)
      
      svd_A = svd(A, nu = n_pcs_response, nv = n_pcs_response)
      svd_B = svd(B, nu = n_pcs_predictors, nv = n_pcs_predictors)
      
      Lambdas[[g]] = crossprod(make_first_entry_positive(svd_A$u),
                               inner_prod_mat1 %*% Omegas1[[g]]) %*% Lambdas[[g]] %*% 
        solve(crossprod(make_first_entry_positive(svd_B$u), inner_prod_mat2 %*% Omegas2[[g]]))
      
      Omegas1[[g]] = solve(inner_prod_mat1) %*% make_first_entry_positive(svd_A$u)
      Omegas2[[g]] = solve(inner_prod_mat2) %*%  make_first_entry_positive(svd_B$u)
      variances_response[[g]] = svd_A$d[1:n_pcs_response]
      variances_predictors[[g]] = svd_B$d[1:n_pcs_predictors]
      Sigma_eta_inv[[g]] = solve(diag(variances_response[[g]]) - Lambdas[[g]] %*% tcrossprod(diag(variances_predictors[[g]]), Lambdas[[g]]), diag(1, nrow=length(variances_response[[g]]), ncol=length(variances_response[[g]])))
    }
    ret[[1]]=Omegas1
    ret[[2]]=Omegas2
    ret[[3]]=variances_response
    ret[[4]]=variances_predictors
    ret[[5]]=Lambdas
    ret[[6]]=Sigma_eta_inv
    return(ret)
  } else if (!is.null(Omegas2)){
    
    lapply(1:G, function(g) {
      n_pcs_response = ncol(Omegas2[[g]])
      n_pcs_predictors = ncol(Omegas1[[g]])
      A = tcrossprod(inner_prod_mat2 %*% Omegas2[[g]] %*% 
                       diag(variances2[[g]]),
                     Omegas2[[g]])  %*%  t(inner_prod_mat2)
      
      #eig_A = eigen(A, symmetric = T)
      svd_A = svd(A, nu = n_pcs_response, nv = n_pcs_response)
      
      # Make sure matrices aren't too ill-conditioned early
      ee_ratios = (svd_A$d[1:n_pcs_response] / max(svd_A$d))
      ee_ratios[ee_ratios < 1e-4] = 1e-4
      new_variances_r = max(svd_A$d) * ee_ratios 
      
      new_Omega2 = solve(inner_prod_mat2) %*% make_first_entry_positive(svd_A$u)
      
      A = tcrossprod(inner_prod_mat1 %*% Omegas1[[g]] %*%
                       diag(variances1[[g]]),
                     Omegas1[[g]])  %*%
        t(inner_prod_mat1)
      
      #eig_A = eigen(A, symmetric = T)
      svd_A = svd(A, nu = n_pcs_predictors, nv = n_pcs_predictors)
      
      ee_ratios = (svd_A$d[1:n_pcs_predictors] / max(svd_A$d))
      ee_ratios[ee_ratios < 1e-5] = 1e-5
      new_variances_p = max(svd_A$d) * ee_ratios
      
      O_alpha <- make_first_entry_positive(svd_A$u)
      new_Omega1 <- solve(inner_prod_mat1) %*% O_alpha
      new_lambda <- Lambdas[[g]] %*% 
        solve(t(O_alpha) %*% inner_prod_mat1 %*% Omegas1[[g]])
      list(new_variances_p, new_lambda, new_Omega1,
           new_variances_r, new_Omega2)
    })
  } else {
    lapply(1:G, function(g) {
      n_pcs_response = ncol(Omegas1[[g]])
      A = tcrossprod(inner_prod_mat1 %*% Omegas1[[g]] %*% 
                       diag(variances1[[g]]),
                     Omegas1[[g]])  %*%  t(inner_prod_mat1)
      #eig_A = eigen(A, symmetric = T)
      svd_A = svd(A, nu = n_pcs_response, nv = n_pcs_response)
      new_variances_r <- svd_A$d[1:n_pcs_response] #eig_A[['values']][1:ncol(Omegas1[[g]])]
      new_Omegas <- solve(inner_prod_mat1) %*% make_first_entry_positive(svd_A$u)
      list(new_variances_r, new_Omegas)
    })
  }
}

compute_marginal_probabilities <- function(cluster_membership, neighbors, theta, n_profiles, G, dists = NULL) {
  marginal_membership <- matrix(NA, length(neighbors), G)
  cluster_mat <- matrix(0, length(cluster_membership), G)
  cluster_mat[cbind(1:length(cluster_membership), 
                    cluster_membership)] <- 1
  
  if(is.null(dists)){
    dists = lapply(neighbors, function(x) rep(1, length(x)))
  }
  
  for(x in 1:length(neighbors)){
    tmp = apply(cluster_mat[neighbors[[x]],,drop=F], 2,
                function(t){exp(theta*sum(t / dists[[x]], na.rm = T))})
    if (sum(tmp) == 0) {
      marginal_membership[x,] <- rep(1/G, G)
    } else {
      marginal_membership[x,] <- tmp/sum(tmp)
    }
  } 
  marginal_membership
}

match_new_pc_col_index <- function(old_pcs, new_pcs, range_params){
  range_params_copy = range_params
  errs = matrix(0, ncol(old_pcs[[1]]), 2)
  for(g in 1:length(old_pcs)){
    old_mat = old_pcs[[g]]
    new_mat = new_pcs[[g]]
    map = rep(0,ncol(old_mat))
    for(q in 1:ncol(old_mat)){
      norm = sqrt(sum(old_mat[,q]^2))
      for(p in 1:ncol(old_mat)){
        errs[p,1] = (sum(((old_mat[,q] - new_mat[,p]) /norm)^2))
        errs[p,2] = (sum(((old_mat[,q] + new_mat[,p]) /norm)^2))
      }
      index = which.min(errs)  
      if(index > ncol(old_mat)){
        map[q] = index - ncol(old_mat)
      } else {
        map[q] = index
      }
    }
    range_params[[g]] = range_params_copy[[g]][,map]
  }
  return(range_params)
}


