# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

splineMatrixC <- function(order, x, knots) {
    .Call(`_fstmr_splineMatrixC`, order, x, knots)
}

c_compute_UTU <- function(phi_x_phi_r, phi_x_phi_p, Omegas_r1, Omegas_r2, Omegas_p, me_r, me_p, n_basis_p, clust_mem, n_samples, n_samples_TS, G, BGC) {
    .Call(`_fstmr_c_compute_UTU`, phi_x_phi_r, phi_x_phi_p, Omegas_r1, Omegas_r2, Omegas_p, me_r, me_p, n_basis_p, clust_mem, n_samples, n_samples_TS, G, BGC)
}

c_compute_UTX <- function(basis_evals_r, basis_evals_p, profs_r, profs_p, Omegas_r1, Omegas_r2, Omegas_p, means_r, means_p, clust_mem, me_r, me_p, n_basis_p, BGC) {
    .Call(`_fstmr_c_compute_UTX`, basis_evals_r, basis_evals_p, profs_r, profs_p, Omegas_r1, Omegas_r2, Omegas_p, means_r, means_p, clust_mem, me_r, me_p, n_basis_p, BGC)
}

c_compute_centered_obs <- function(basis_evals_r, basis_evals_p, profs_r, profs_p, means_r, means_p, clust_mem, n_basis_p, BGC) {
    .Call(`_fstmr_c_compute_centered_obs`, basis_evals_r, basis_evals_p, profs_r, profs_p, means_r, means_p, clust_mem, n_basis_p, BGC)
}

c_compute_UTU_single <- function(phi_x_phi, Omegas, me, clust_mem, n_samples, G) {
    .Call(`_fstmr_c_compute_UTU_single`, phi_x_phi, Omegas, me, clust_mem, n_samples, G)
}

c_compute_UTX_single <- function(basis_evals, profs, Omegas, means, clust_mem, me) {
    .Call(`_fstmr_c_compute_UTX_single`, basis_evals, profs, Omegas, means, clust_mem, me)
}

c_compute_centered_obs_single <- function(basis_evals, profs, clust_mem, means, me) {
    .Call(`_fstmr_c_compute_centered_obs_single`, basis_evals, profs, clust_mem, means, me)
}

c_compute_E_step_likelihoods <- function(profs_resp, profs_pred, is_bgc, n_profiles, basis_evals_r, basis_evals_p, Omegas_r1, Omegas_r2, Omegas_p, means_resp, means_pred, variances, vars_r, vars_p, profile_lengths_p, G, n_preds) {
    .Call(`_fstmr_c_compute_E_step_likelihoods`, profs_resp, profs_pred, is_bgc, n_profiles, basis_evals_r, basis_evals_p, Omegas_r1, Omegas_r2, Omegas_p, means_resp, means_pred, variances, vars_r, vars_p, profile_lengths_p, G, n_preds)
}

c_compute_E_step_likelihoods_ind <- function(profs_resp, profs_pred, is_bgc, n_profiles, basis_evals_resp, basis_evals_pred, Gammas, Omegas_resp, Omegas_preds, means_resp, means_pred, variances, vars_pred, profile_lengths_p, G, n_preds) {
    .Call(`_fstmr_c_compute_E_step_likelihoods_ind`, profs_resp, profs_pred, is_bgc, n_profiles, basis_evals_resp, basis_evals_pred, Gammas, Omegas_resp, Omegas_preds, means_resp, means_pred, variances, vars_pred, profile_lengths_p, G, n_preds)
}

c_compute_E_step_likelihoods_single <- function(profs, n_profiles, basis_evals, Omegas, means, me, vars, G) {
    .Call(`_fstmr_c_compute_E_step_likelihoods_single`, profs, n_profiles, basis_evals, Omegas, means, me, vars, G)
}

c_create_summed_U_matrix_sparse <- function(phi_x_phi, clust_mem, weights, G, reps) {
    .Call(`_fstmr_c_create_summed_U_matrix_sparse`, phi_x_phi, clust_mem, weights, G, reps)
}

c_create_summed_V_matrix_sparse <- function(basis_evals, profiles, G, Omegas1, Omegas2, clust_mem, pcs1, pcs2, reps, weights) {
    .Call(`_fstmr_c_create_summed_V_matrix_sparse`, basis_evals, profiles, G, Omegas1, Omegas2, clust_mem, pcs1, pcs2, reps, weights)
}

c_create_summed_V_matrix_sparse_single <- function(basis_evals, profiles, G, Omegas, clust_mem, pcs, reps, weights) {
    .Call(`_fstmr_c_create_summed_V_matrix_sparse_single`, basis_evals, profiles, G, Omegas, clust_mem, pcs, reps, weights)
}

c_create_summed_V_matrix_pcs_sparse <- function(basis_evals, profiles, means, G, Omegas1, Omegas2, clust_mem, pcs1, pcs2, reps, weights, q) {
    .Call(`_fstmr_c_create_summed_V_matrix_pcs_sparse`, basis_evals, profiles, means, G, Omegas1, Omegas2, clust_mem, pcs1, pcs2, reps, weights, q)
}

c_create_summed_V_matrix_pcs_sparse_single <- function(basis_evals, profiles, means, G, Omegas, clust_mem, pcs, reps, weights, q) {
    .Call(`_fstmr_c_create_summed_V_matrix_pcs_sparse_single`, basis_evals, profiles, means, G, Omegas, clust_mem, pcs, reps, weights, q)
}

c_create_summed_U_matrix_pcs_sparse <- function(phi_x_phi, clust_mem, weights, pc_weights, G, reps) {
    .Call(`_fstmr_c_create_summed_U_matrix_pcs_sparse`, phi_x_phi, clust_mem, weights, pc_weights, G, reps)
}

c_create_summed_V_matrix_gamma_sparse_r <- function(basis_evals, profiles, means, G, Omegas1, Omegas2, clust_mem, etas, alphas, reps, weights, q) {
    .Call(`_fstmr_c_create_summed_V_matrix_gamma_sparse_r`, basis_evals, profiles, means, G, Omegas1, Omegas2, clust_mem, etas, alphas, reps, weights, q)
}

c_create_summed_U_matrix_gamma_sparse <- function(phi_x_phi, Omegas1, clust_mem, weights, pc_weights, G, reps) {
    .Call(`_fstmr_c_create_summed_U_matrix_gamma_sparse`, phi_x_phi, Omegas1, clust_mem, weights, pc_weights, G, reps)
}

c_update_measurement_error <- function(cluster_mat, basis_evals, pcs1, pcs2, profiles, means_mat, Omegas1, Omegas2, n_profiles, weights, G) {
    .Call(`_fstmr_c_update_measurement_error`, cluster_mat, basis_evals, pcs1, pcs2, profiles, means_mat, Omegas1, Omegas2, n_profiles, weights, G)
}

c_compute_squared_sparse <- function(cluster_mat_i, pcs_mat, profile, G, basis_eval, means_mat, Omegas, weights) {
    .Call(`_fstmr_c_compute_squared_sparse`, cluster_mat_i, pcs_mat, profile, G, basis_eval, means_mat, Omegas, weights)
}

c_compute_squared_sparse_response <- function(cluster_mat_i, pcs_mat1, pcs_mat2, profile, G, basis_eval, means_mat, Lambda1, Lambda2, weights) {
    .Call(`_fstmr_c_compute_squared_sparse_response`, cluster_mat_i, pcs_mat1, pcs_mat2, profile, G, basis_eval, means_mat, Lambda1, Lambda2, weights)
}

c_update_measurement_error_p_space <- function(cluster_mat, basis_evals, pcs, profiles, means_mat, Omegas, n_profiles, G) {
    .Call(`_fstmr_c_update_measurement_error_p_space`, cluster_mat, basis_evals, pcs, profiles, means_mat, Omegas, n_profiles, G)
}

c_compute_conditional_distribution <- function(profs_resp, profs_pred, basis_evals_resp, basis_evals_pred, phi_x_phi_resp, phi_x_phi_pred, n_samples, means_resp, means_pred, Omegas_resp, Omegas_pred, Lambdas, Sigma_eta_inv, me_resp, me_pred, vars_resp, vars_pred, basis_lengths_pred, cond_probs, is_bgc, conditional_distributions) {
    invisible(.Call(`_fstmr_c_compute_conditional_distribution`, profs_resp, profs_pred, basis_evals_resp, basis_evals_pred, phi_x_phi_resp, phi_x_phi_pred, n_samples, means_resp, means_pred, Omegas_resp, Omegas_pred, Lambdas, Sigma_eta_inv, me_resp, me_pred, vars_resp, vars_pred, basis_lengths_pred, cond_probs, is_bgc, conditional_distributions))
}

c_lik_eigen_sherman_pred <- function(x, mean, variances, U, W) {
    .Call(`_fstmr_c_lik_eigen_sherman_pred`, x, mean, variances, U, W)
}

stl_sort <- function(x) {
    .Call(`_fstmr_stl_sort`, x)
}

c_CG <- function(A, w) {
    .Call(`_fstmr_c_CG`, A, w)
}

c_quad_form_log <- function(A, w, maxiter = 1000L, tol = 0.000001, stepsize = 10L, precondition = TRUE) {
    .Call(`_fstmr_c_quad_form_log`, A, w, maxiter, tol, stepsize, precondition)
}

c_reorder <- function(A, ind) {
    .Call(`_fstmr_c_reorder`, A, ind)
}

c_sample_lanczos <- function(A, w, stepsize = 10L, maxiter = 1500L, tol = 1e-3, precondition = TRUE) {
    .Call(`_fstmr_c_sample_lanczos`, A, w, stepsize, maxiter, tol, precondition)
}

c_compute_E_step_likelihoods_for_season <- function(profs_resp, profs_pred, is_bgc, n_profiles, basis_evals_r, basis_evals_p, Omegas_r1, Omegas_r2, Omegas_p, means_resp, means_pred, variances, vars_r, vars_p, profile_lengths_p, G, n_preds, days) {
    .Call(`_fstmr_c_compute_E_step_likelihoods_for_season`, profs_resp, profs_pred, is_bgc, n_profiles, basis_evals_r, basis_evals_p, Omegas_r1, Omegas_r2, Omegas_p, means_resp, means_pred, variances, vars_r, vars_p, profile_lengths_p, G, n_preds, days)
}

c_compute_UTX_for_season <- function(basis_evals_r, basis_evals_p, profs_r, profs_p, Omegas_r1, Omegas_r2, Omegas_p, means_r, means_p, clust_mem, me_r, me_p, n_basis_p, BGC, days) {
    .Call(`_fstmr_c_compute_UTX_for_season`, basis_evals_r, basis_evals_p, profs_r, profs_p, Omegas_r1, Omegas_r2, Omegas_p, means_r, means_p, clust_mem, me_r, me_p, n_basis_p, BGC, days)
}

c_compute_centered_obs_for_season <- function(basis_evals_r, basis_evals_p, profs_r, profs_p, means_r, means_p, clust_mem, n_basis_p, BGC, days) {
    .Call(`_fstmr_c_compute_centered_obs_for_season`, basis_evals_r, basis_evals_p, profs_r, profs_p, means_r, means_p, clust_mem, n_basis_p, BGC, days)
}

c_create_summed_U_matrix_pcs_sparse_for_season <- function(phi_x_phi, clust_mem, weights, pc_weights, pc_weights2, G, reps) {
    .Call(`_fstmr_c_create_summed_U_matrix_pcs_sparse_for_season`, phi_x_phi, clust_mem, weights, pc_weights, pc_weights2, G, reps)
}

c_create_summed_V_matrix_sparse_for_season <- function(basis_evals, profiles, G, Omegas1, Omegas2, clust_mem, pcs1, pcs2, additional_weights, reps, weights) {
    .Call(`_fstmr_c_create_summed_V_matrix_sparse_for_season`, basis_evals, profiles, G, Omegas1, Omegas2, clust_mem, pcs1, pcs2, additional_weights, reps, weights)
}

c_create_summed_V_matrix_sparse_single_for_season <- function(basis_evals, profiles, G, Omegas, clust_mem, pcs, additional_weights, reps, weights) {
    .Call(`_fstmr_c_create_summed_V_matrix_sparse_single_for_season`, basis_evals, profiles, G, Omegas, clust_mem, pcs, additional_weights, reps, weights)
}

c_create_summed_V_matrix_pcs_sparse_single_for_season <- function(basis_evals, profiles, means, G, Omegas, clust_mem, pcs, reps, weights, q, days) {
    .Call(`_fstmr_c_create_summed_V_matrix_pcs_sparse_single_for_season`, basis_evals, profiles, means, G, Omegas, clust_mem, pcs, reps, weights, q, days)
}

c_create_summed_V_matrix_pcs_sparse_for_season <- function(basis_evals, profiles, means, G, Omegas1, Omegas2, clust_mem, pcs1, pcs2, reps, weights, q, days) {
    .Call(`_fstmr_c_create_summed_V_matrix_pcs_sparse_for_season`, basis_evals, profiles, means, G, Omegas1, Omegas2, clust_mem, pcs1, pcs2, reps, weights, q, days)
}

c_compute_squared_sparse_for_season <- function(cluster_mat_i, pcs_mat, profile, G, basis_eval, means_mat, Omegas, weights, day) {
    .Call(`_fstmr_c_compute_squared_sparse_for_season`, cluster_mat_i, pcs_mat, profile, G, basis_eval, means_mat, Omegas, weights, day)
}

c_update_measurement_error_for_season <- function(cluster_mat, basis_evals, pcs1, pcs2, profiles, means_mat, Omegas1, Omegas2, n_profiles, weights, G, days) {
    .Call(`_fstmr_c_update_measurement_error_for_season`, cluster_mat, basis_evals, pcs1, pcs2, profiles, means_mat, Omegas1, Omegas2, n_profiles, weights, G, days)
}

c_compute_conditional_distribution_for_season <- function(profs_resp, profs_pred, basis_evals_resp, basis_evals_pred, phi_x_phi_resp, phi_x_phi_pred, n_samples, means_resp, means_pred, Omegas_resp, Omegas_pred, Lambdas, Sigma_eta_inv, me_resp, me_pred, vars_resp, vars_pred, basis_lengths_pred, cond_probs, is_bgc, days, conditional_distributions) {
    invisible(.Call(`_fstmr_c_compute_conditional_distribution_for_season`, profs_resp, profs_pred, basis_evals_resp, basis_evals_pred, phi_x_phi_resp, phi_x_phi_pred, n_samples, means_resp, means_pred, Omegas_resp, Omegas_pred, Lambdas, Sigma_eta_inv, me_resp, me_pred, vars_resp, vars_pred, basis_lengths_pred, cond_probs, is_bgc, days, conditional_distributions))
}

c_compute_E_step_likelihoods_ind_for_season <- function(profs_resp, profs_pred, is_bgc, n_profiles, basis_evals_resp, basis_evals_pred, Gammas, Omegas_resp, Omegas_preds, means_resp, means_pred, variances, vars_pred, profile_lengths_p, days, G, n_preds) {
    .Call(`_fstmr_c_compute_E_step_likelihoods_ind_for_season`, profs_resp, profs_pred, is_bgc, n_profiles, basis_evals_resp, basis_evals_pred, Gammas, Omegas_resp, Omegas_preds, means_resp, means_pred, variances, vars_pred, profile_lengths_p, days, G, n_preds)
}

c_compute_E_step_likelihoods_single_for_season <- function(profs, n_profiles, basis_evals, Omegas, means, me, vars, G, days) {
    .Call(`_fstmr_c_compute_E_step_likelihoods_single_for_season`, profs, n_profiles, basis_evals, Omegas, means, me, vars, G, days)
}

c_compute_UTX_single_for_season <- function(basis_evals, profs, Omegas, means, clust_mem, me, days) {
    .Call(`_fstmr_c_compute_UTX_single_for_season`, basis_evals, profs, Omegas, means, clust_mem, me, days)
}

c_compute_centered_obs_single_for_season <- function(basis_evals, profs, clust_mem, means, me, G, days) {
    .Call(`_fstmr_c_compute_centered_obs_single_for_season`, basis_evals, profs, clust_mem, means, me, G, days)
}

