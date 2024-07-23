// [[Rcpp::depends("RcppEigen")]]
// [[Rcpp::depends("RcppGSL")]]
// [[Rcpp::depends("Rcpp")]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif 

using Eigen::Map;  
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;

using namespace Eigen;

typedef Eigen::MappedSparseMatrix<double> MSpMatC;
typedef Eigen::Triplet<double> Trip;

static double const log2pi = std::log(2.0 * M_PI);


// [[Rcpp::export]]
Rcpp::NumericMatrix splineMatrixC(const int order,
                                  const Rcpp::NumericVector x, const Rcpp::NumericVector knots){
  // The gsl library creates one object that stores 
  // all necessary b-spline basis attributes. This is this
  // workspace. 
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  
  //Initialize...
  const int n = x.size();
  const int nbreaks = knots.size();
  const int nsplines = nbreaks + order - 2; 
  
  bw = gsl_bspline_alloc(order, nbreaks);
  B = gsl_vector_alloc(nsplines);
  // Could potentially already be a sparse matrix
  // However, if other basis produce dense matrices,
  // then it would be problematic
  Rcpp::NumericMatrix X = Rcpp::NumericMatrix(n, nsplines);
  
  gsl_vector *knots_gsl;
  knots_gsl = gsl_vector_alloc(nbreaks);
  for (int i = 0; i < nbreaks; ++i){
    gsl_vector_set(knots_gsl, i, knots(i));
  }
  
  gsl_bspline_knots(knots_gsl, bw);
  
  /* construct the fit matrix X */
  for (int i = 0; i < n; ++i){
    
    /* compute B_j(xi) for all j */
    gsl_bspline_eval(x[i], B, bw);
    
    /* fill in row i of X */
    for (int j = 0; j < nsplines; ++j){
      X(i,j) = gsl_vector_get(B, j);
    }
  }
  return X;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double, Eigen::ColMajor> c_compute_UTU(const Rcpp::List phi_x_phi_r,
                                                           const Rcpp::List phi_x_phi_p,
                                                           const Rcpp::List Omegas_r1,
                                                           const Rcpp::List Omegas_r2,
                                                           const Rcpp::List Omegas_p,
                                                           const double me_r,
                                                           const Eigen::Map<Eigen::VectorXd> me_p,
                                                           const Eigen::Map<Eigen::VectorXi> n_basis_p,
                                                           const Eigen::Map<Eigen::VectorXi> clust_mem,
                                                           const int n_samples,
                                                           const int n_samples_TS,
                                                           const int G,
                                                           const Eigen::Map<Eigen::VectorXi> BGC){
  
  int n_pcs_p = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas_p[0]).cols();
  int n_pcs_r2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas_r2[0]).cols();
  int n_basis = Rcpp::as<MSpMatC>(phi_x_phi_r[0]).cols();
  int n_pred = me_p.size();
  
  int dim = n_samples *n_pcs_r2 + n_samples_TS * n_pcs_p;
  int nzero = 0.5 * n_samples_TS * (n_pcs_p+1) * n_pcs_p + 
    + 0.5 * n_samples * n_pcs_r2 * (n_pcs_r2 + 1) + 
    + n_samples * n_pcs_p * n_pcs_r2;
  int g_index, q, l, start;
  
  int offset = n_samples_TS * n_pcs_p;
  //int offset2 = 0.5 * n_samples * (n_pcs_p * (2 * n_pcs_r2 + (n_pcs_p+1)));
  //int offset2 = 0.5 * n_samples * (2 * n_pcs_r2 * n_pcs_p  +  n_pcs_p*(n_pcs_p+ 1)));
  int offset2 = n_samples * n_pcs_r2 * n_pcs_p  +  0.5*n_samples_TS * n_pcs_p*(n_pcs_p+ 1);
  
  int* p = new int[dim+1]();
  int* i = new int[nzero]();
  double* v = new double[nzero]();
  
  std::vector<Eigen::MatrixXd > omegas_r1, omegas_r2, omegas_p;
  for (size_t g = 0; g < G; g++){
    omegas_r1.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas_r1[g]));
    omegas_r1[g] *= (1 / std::sqrt(me_r));
    omegas_r2.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas_r2[g]));
    omegas_r2[g] *= (1 / std::sqrt(me_r));
    omegas_p.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas_p[g]));
    start = 0;
    for (int bs = 0; bs < n_pred; bs++){
      omegas_p[g].block(start, 0, n_basis_p(bs), omegas_p[g].cols()) *= (1 / std::sqrt(me_p(bs)));
      start += n_basis_p(bs);
    }
  }
  
  Eigen::MatrixXd temp(n_pcs_p, n_pcs_p);
  Eigen::MatrixXd temp2(n_pcs_p, n_pcs_r2);
  Eigen::MatrixXd temp3(n_pcs_r2, n_pcs_r2);
  Eigen::MatrixXd temp4(n_pcs_p, n_pcs_p);
  Eigen::MatrixXd temp5(n_pcs_p, n_basis);
  
  int vc = 0;
  int pc = 0;
  int pc2 = offset2;
  int vc2 = offset2;
  
  p[0] = 0;
  int BGC_index = 0;
  for (int x = 0; x < n_samples_TS; x++){
    g_index = clust_mem(x);
    MSpMatC psi = Rcpp::as<MSpMatC>(phi_x_phi_p[x]);
    temp.noalias() = omegas_p[g_index].adjoint() * psi * omegas_p[g_index];
    if (BGC(x)) {
      MSpMatC phi = Rcpp::as<MSpMatC>(phi_x_phi_r[BGC_index]);
      temp5 = omegas_r1[g_index].adjoint() * phi;
      temp2.noalias() = temp5 * omegas_r2[g_index];
      temp3.noalias() = omegas_r2[g_index].adjoint() * phi * omegas_r2[g_index];
      temp4.noalias() = temp5 * omegas_r1[g_index]; 

      for (q = 0; q < n_pcs_p; q++){
        for (l = q; l < n_pcs_p; l++){
          v[vc] = temp(l,q) + temp4(l,q);
          i[vc] = x*n_pcs_p + l;
          vc++;
        }
        for (l = 0; l < n_pcs_r2; l++){
          v[vc] = temp2(q,l);
          i[vc] = BGC_index*n_pcs_r2 + offset + l;
          vc++;
        }
        pc += n_pcs_p +  n_pcs_r2 - q;
        p[x*n_pcs_p+q+1] = pc;
      }
      
      for (q = 0; q < n_pcs_r2; q++){
        for (l = q; l < n_pcs_r2; l++){
          v[vc2] = temp3(l,q);
          i[vc2] = BGC_index*n_pcs_r2 + offset + l;
          vc2++;
        }
        pc2 += n_pcs_r2 - q;
        p[BGC_index*n_pcs_r2 + q  + offset+1] = pc2;
      }
      BGC_index += 1;
    } else {
      for (q = 0; q < n_pcs_p; q++){
        for (l = q; l < n_pcs_p; l++){
          v[vc] = temp(l,q);
          i[vc] = x*n_pcs_p + l;
          vc++;
        }
        pc += n_pcs_p  - q;
        p[x*n_pcs_p+q+1] = pc;
      }
    }
  }
  
  p[dim] = nzero;
  
  Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor>> spMap(dim, dim, nzero, p, i, v, 0);
  Eigen::SparseMatrix<double, Eigen::ColMajor> res= spMap.eval();
  
  return res;
}



// [[Rcpp::export]]
Eigen::VectorXd c_compute_UTX(const Rcpp::List basis_evals_r,
                       const Rcpp::List basis_evals_p,
                       const Rcpp::List profs_r,
                       const Rcpp::List profs_p,
                       const Rcpp::List Omegas_r1,
                       const Rcpp::List Omegas_r2,
                       const Rcpp::List Omegas_p,
                       const Eigen::Map<Eigen::MatrixXd> means_r,
                       const Eigen::Map<Eigen::MatrixXd> means_p,
                       const Eigen::Map<Eigen::VectorXi> clust_mem,
                       const double me_r,
                       const Eigen::Map<Eigen::VectorXd> me_p,
                       const Eigen::Map<Eigen::VectorXi> n_basis_p,
                       const Eigen::Map<Eigen::VectorXi> BGC){
  
  const int n_pred = n_basis_p.size();
  const int n_samples = basis_evals_r.size();
  const int n_samples_TS = basis_evals_p.size();
  const int n_pcs_r2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas_r2[0]).cols();
  const int n_pcs_p = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas_p[0]).cols();
  const int n_basis_r = Rcpp::as<MSpMatC>(basis_evals_r[0]).cols();
  const int n_basis_pt = Rcpp::as<MSpMatC>(basis_evals_p[0]).cols();
  const int offset = n_pcs_p * n_samples_TS;
  const int G = Omegas_r1.size();
  int start, g_index;
  
  
  Eigen::VectorXd res = Eigen::VectorXd::Zero(n_samples * n_pcs_r2 + n_samples_TS * n_pcs_p);
  
  std::vector<Eigen::MatrixXd > omegas_r1, omegas_r2, omegas_p;
  for (size_t g = 0; g < G; g++){
    omegas_r1.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas_r1[g]));
    omegas_r1[g] *= (1 / me_r);
    omegas_r2.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas_r2[g]));
    omegas_r2[g] *= (1 / me_r);
    omegas_p.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas_p[g]));
    start = 0;
    for (int bs = 0; bs < n_pred; bs++){
      omegas_p[g].block(start, 0, n_basis_p(bs), omegas_p[g].cols()) *= (1 / me_p(bs));
      start += n_basis_p(bs);
    }
  }
  
  VectorXd temp_r(n_basis_r);
  VectorXd temp_p(n_basis_pt);
  
  int start_r = 0;
  int start_p = 0;
  int BGC_index = 0;
  for (int x = 0; x < n_samples_TS; x++){
    g_index = clust_mem(x);
    MSpMatC psi = Rcpp::as<MSpMatC>(basis_evals_p[x]);
    Eigen::Map<Eigen::VectorXd> X = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_p[x]);
    temp_p = psi.adjoint() * (X - psi * means_p.col(g_index));
    
    if (BGC(x)) {
      MSpMatC phi = Rcpp::as<MSpMatC>(basis_evals_r[BGC_index]);
      Eigen::Map<Eigen::VectorXd> Y = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_r[BGC_index]);
      temp_r = phi.adjoint() * (Y - phi * means_r.col(g_index));
      res.segment(start_p, n_pcs_p) = omegas_p[g_index].adjoint() * temp_p + omegas_r1[g_index].adjoint() * temp_r;
      res.segment(offset + start_r, n_pcs_r2) = omegas_r2[g_index].adjoint() * temp_r;
      start_r += n_pcs_r2;
      start_p += n_pcs_p;
      BGC_index += 1;
    } else {
      res.segment(start_p, n_pcs_p) = omegas_p[g_index].adjoint() * temp_p;
      start_p += n_pcs_p;
    }
  }
  
  return res;
}

// [[Rcpp::export]]
double c_compute_centered_obs(const Rcpp::List basis_evals_r,
                                const Rcpp::List basis_evals_p,
                                const Rcpp::List profs_r,
                                const Rcpp::List profs_p,
                                const Eigen::Map<Eigen::MatrixXd> means_r,
                                const Eigen::Map<Eigen::MatrixXd> means_p,
                                const Eigen::Map<Eigen::VectorXi> clust_mem,
                                const Eigen::Map<Eigen::VectorXi> n_basis_p,
                                const Eigen::Map<Eigen::VectorXi> BGC){
  
  const int n_samples_TS = basis_evals_p.size();
  int g_index;
  
  
  double res = 0;
  int BGC_index = 0;
  Eigen::VectorXd me_preds = Rcpp::as<Eigen::VectorXd>(profs_p[0]);
  for (int x = 0; x < n_samples_TS; x++){
    g_index = clust_mem(x);
    MSpMatC psi = Rcpp::as<MSpMatC>(basis_evals_p[x]);
    Eigen::Map<Eigen::VectorXd> X = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_p[x]);
    res = res +  (X - psi * means_p.col(g_index)).array().square().sum();
    
    if (BGC(x)) {
      MSpMatC phi = Rcpp::as<MSpMatC>(basis_evals_r[BGC_index]);
      Eigen::Map<Eigen::VectorXd> Y = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_r[BGC_index]);
      res = res + ((Y - phi * means_r.col(g_index))).array().square().sum();
      BGC_index += 1;
    } 
  }
  
  return res;
}


// [[Rcpp::export]]
Eigen::SparseMatrix<double, Eigen::ColMajor> c_compute_UTU_single(const Rcpp::List phi_x_phi,
                                                                  const Rcpp::List Omegas,
                                                                  const double me,
                                                                  const Eigen::Map<Eigen::VectorXi> clust_mem,
                                                                  const int n_samples,
                                                                  const int G){
  
  int n_pcs = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[0]).cols();
  int g_index;
  Eigen::MatrixXd temp;
  
  int dim = n_samples * n_pcs;
  int nzero = n_samples * n_pcs * n_pcs;
  
  int* p = new int[dim+1]();
  int* i = new int[nzero]();
  double* v = new double[nzero]();
  
  std::vector<Eigen::MatrixXd > omegas;
  for (size_t g = 0; g < G; g++){
    omegas.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas[g]));
    omegas[g] *= (1 / std::sqrt(me));
  }
  
  int vc = 0;
  int pc = 0;
  
  p[0] = 0;
  for (int x = 0; x < n_samples; x++){
    g_index = clust_mem(x);
    MSpMatC phi = Rcpp::as<MSpMatC>(phi_x_phi[x]);
    temp.noalias() = omegas[g_index].adjoint() * phi * omegas[g_index];
    
    for (int q = 0; q < n_pcs; q++){
      for (int p = 0; p < n_pcs; p++){
        v[vc] = temp(p,q);
        i[vc] = x * n_pcs + p;
        vc++;
      }
      pc++;
      p[pc] = p[pc-1] + n_pcs; 
    }
  }
  
  p[dim] = nzero;
  
  Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor>> spMap(dim, dim, nzero, p, i, v, 0);
  Eigen::SparseMatrix<double, Eigen::ColMajor> res= spMap.eval();
  
  return res;
}



// [[Rcpp::export]]
Eigen::VectorXd c_compute_UTX_single(const Rcpp::List basis_evals,
                                     const Rcpp::List profs,
                                     const Rcpp::List Omegas,
                                     const Eigen::Map<Eigen::MatrixXd> means,
                                     const Eigen::Map<Eigen::VectorXi> clust_mem,
                                     const double me){
  
  const int n_samples = basis_evals.size();
  const int n_pcs = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[0]).cols();
  const int n_basis = Rcpp::as<MSpMatC>(basis_evals[0]).cols();
  const int G = Omegas.size();
  int g_index;
  
  
  Eigen::VectorXd res = Eigen::VectorXd::Zero(n_samples * n_pcs);
  
  std::vector<Eigen::MatrixXd > omegas;
  for (size_t g = 0; g < G; g++){
    omegas.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas[g]));
    omegas[g] *= (1 / me);
  }
  
  Eigen::VectorXd temp(n_basis);
  
  int start = 0;
  
  for (int x = 0; x < n_samples; x++){
    g_index = clust_mem(x);
    MSpMatC phi = Rcpp::as<MSpMatC>(basis_evals[x]);
    Eigen::Map<Eigen::VectorXd> X = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs[x]);
    temp = phi.adjoint() * (X - phi * means.col(g_index));
    res.segment(start, n_pcs) = omegas[g_index].adjoint() * temp;
    start += n_pcs;
  }
  
  return res;
}

// [[Rcpp::export]]
double c_compute_centered_obs_single(const Rcpp::List basis_evals,
                                     const Rcpp::List profs,
                                     const Eigen::Map<Eigen::VectorXi> clust_mem,
                                     const Eigen::Map<Eigen::MatrixXd> means,
                                     const double me){
  
  const int n_samples = basis_evals.size();
  int g_index;
  
  double res = 0;
  for (int x = 0; x < n_samples; x++){
    g_index = clust_mem(x);
    MSpMatC phi = Rcpp::as<MSpMatC>(basis_evals[x]);
    Eigen::Map<Eigen::VectorXd> X = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs[x]);
    res = res + ((X - phi * means.col(g_index))/std::sqrt(me)).array().square().sum();
  }
  return res;
}

/*
 * This is for a covariance of form diag(variances) + UU^{\top}
 * where variances has a block structure
 */

double lik_eigen_sherman(const Eigen::VectorXd& x,
                         const Eigen::VectorXd& mean,
                         const Eigen::VectorXd& variances,
                         const Eigen::VectorXi& block_sizes,
                         Eigen::MatrixXd& U){
  const int dim = x.size();
  const int n = U.cols();
  double out, log_det;
  int start = 0;
  Eigen::VectorXd centered = x-mean;
  Eigen::VectorXd v_temp(dim);
  
  for (size_t s = 0; s < block_sizes.rows(); s++){
    U.block(start, 0, block_sizes[s], n) *= 1 / std::sqrt(variances[s]);
    v_temp.segment(start, block_sizes[s]) = 1 / std::sqrt(variances[s]) * centered.segment(start, block_sizes[s]);
    start += block_sizes[s];
  }
  
  // LLT seems slightly faster but this depends on the matrix size and LDLT is more stable
  Eigen::LDLT<Eigen::MatrixXd> ldlt;
  
  Eigen::MatrixXd temp = Eigen::MatrixXd::Identity(n, n);
  temp.selfadjointView<Lower>().rankUpdate(U.adjoint());
  ldlt.compute(temp.selfadjointView<Lower>());
  
  out = v_temp.adjoint() * (v_temp - (U * ldlt.solve(U.adjoint() * v_temp)));
  
  Eigen::VectorXd d = ldlt.vectorD();
  log_det = std::log(d.prod()) + (block_sizes.cast<double>().array() * variances.array().log()).sum();
  return - 0.5 * (dim * log2pi + log_det + out);
}


// [[Rcpp::export]]
Eigen::MatrixXd c_compute_E_step_likelihoods(const Rcpp::List profs_resp,
                                             const Rcpp::List profs_pred,
                                             const Eigen::Map<Eigen::VectorXi> is_bgc,
                                             const int n_profiles,
                                             const Rcpp::List basis_evals_r,
                                             const Rcpp::List basis_evals_p,
                                             const Rcpp::List Omegas_r1,
                                             const Rcpp::List Omegas_r2,
                                             const Rcpp::List Omegas_p,
                                             const Eigen::Map<Eigen::MatrixXd> means_resp,
                                             const Eigen::Map<Eigen::MatrixXd> means_pred,
                                             const Eigen::Map<Eigen::VectorXd> variances,
                                             const Eigen::Map<Eigen::MatrixXd> vars_r,
                                             const Eigen::Map<Eigen::MatrixXd> vars_p,
                                             const Eigen::Map<Eigen::MatrixXi> profile_lengths_p,
                                             const int G,
                                             const int n_preds){
  
  int n_total;
  int bgc_index = 0;
  const int n_pcs_r2 = Rcpp::as<Eigen::MatrixXd> (Omegas_r2[0]).cols();
  const int n_pcs_p = Rcpp::as<Eigen::MatrixXd> (Omegas_p[0]).cols();
  const int n_pcs = n_pcs_p + n_pcs_r2;

  Eigen::VectorXd x, mean;
  Eigen::VectorXi profile_lengths(n_preds+1);
  Eigen::MatrixXd liks(n_profiles, G);
  Eigen::MatrixXd L;
  
  std::vector<Eigen::MatrixXd> L1, L2, L3;
  for (size_t g = 0; g < G; g++){
    // Need the square root here
    L1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas_p[g]) * vars_p.col(g).array().sqrt().matrix().asDiagonal());
    L2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas_r1[g]) * vars_p.col(g).array().sqrt().matrix().asDiagonal());
    L3.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas_r2[g]) * vars_r.col(g).array().sqrt().matrix().asDiagonal());
  }
  
  for (size_t x = 0; x < n_profiles; x++){
    
    Eigen::MappedSparseMatrix<double> phi_p = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals_p[x]);
    
    if (is_bgc[x]){
      
      Eigen::MappedSparseMatrix<double> phi_r = Rcpp::as<Eigen::MappedSparseMatrix<double>  >(basis_evals_r[bgc_index]);
      
      n_total = phi_r.rows() + phi_p.rows();
      
      Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(n_total, n_pcs);
      Eigen::VectorXd x_temp(n_total);
      Eigen::VectorXd mean(n_total);
      
      x_temp << (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs_pred[x])), Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs_resp[bgc_index]);
      profile_lengths.head(n_preds) = profile_lengths_p.row(x);
      profile_lengths(profile_lengths.size()-1) = phi_r.rows();

      for (size_t g = 0; g < G; g++){
        
        mean << phi_p * means_pred.col(g), phi_r * means_resp.col(g);
        temp.topLeftCorner(phi_p.rows(), n_pcs_p) = phi_p * L1[g];
        temp.bottomLeftCorner(phi_r.rows(), n_pcs_p) = phi_r * L2[g];
        temp.bottomRightCorner(phi_r.rows(), n_pcs_r2) = phi_r * L3[g];
        
        liks(x, g) = lik_eigen_sherman(x_temp, mean, variances, profile_lengths, temp);
      }
      bgc_index++;
      
    } else {
      
      n_total = phi_p.rows();
      
      Eigen::MatrixXd temp(n_total, n_pcs_p);
      Eigen::VectorXd x_temp(n_total);
      Eigen::VectorXd mean(n_total);
      
      x_temp << Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_pred[x]);
      
      profile_lengths.head(n_preds) = profile_lengths_p.row(x);
      
      for (size_t g = 0; g < G; g++){
        
        mean << phi_p * means_pred.col(g);
        temp = phi_p * L1[g];
        liks(x, g) = lik_eigen_sherman(x_temp, mean, variances.head(n_preds), profile_lengths.head(n_preds), temp);
      }
    }
  }
  return liks;
}

// [[Rcpp::export]]
Eigen::MatrixXd c_compute_E_step_likelihoods_ind(const Rcpp::List profs_resp,
                                                 const Rcpp::List profs_pred,
                                                 const Eigen::Map<Eigen::VectorXi> is_bgc,
                                                 const int n_profiles,
                                                 const Rcpp::List basis_evals_resp,
                                                 const Rcpp::List basis_evals_pred,
                                                 Rcpp::List Gammas,
                                                 const Rcpp::List Omegas_resp,
                                                 const Rcpp::List Omegas_preds,
                                                 const Eigen::Map<Eigen::MatrixXd> means_resp,
                                                 const Eigen::Map<Eigen::MatrixXd> means_pred,
                                                 const Eigen::Map<Eigen::VectorXd> variances,
                                                 const Eigen::Map<Eigen::MatrixXd> vars_pred,
                                                 const Eigen::Map<Eigen::MatrixXi> profile_lengths_p,
                                                 const int G,
                                                 const int n_preds){
  
  int n_total;
  int bgc_index = 0;
  const int n_pcs_resp = Rcpp::as<Eigen::MatrixXd> (Omegas_resp[0]).cols();
  const int n_pcs_pred = Rcpp::as<Eigen::MatrixXd> (Omegas_preds[0]).cols();
  const int n_pcs = n_pcs_resp + n_pcs_pred;
  Eigen::VectorXd x, mean;
  Eigen::VectorXi profile_lengths(n_preds+1);
  Eigen::MatrixXd liks(n_profiles, G);
  Eigen::MatrixXd L;
  
  std::vector<Eigen::MatrixXd> L1, L2, L_pred;
  for (size_t g = 0; g < G; g++){
    // Need the square root here
    L = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Gammas[g]).llt().matrixL();
    L1.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_resp[g]) * L.topRows(n_pcs_resp));
    L2.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_preds[g]) * L.bottomRows(n_pcs_pred));
    L_pred.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_preds[g]) * vars_pred.col(g).array().sqrt().matrix().asDiagonal());
  }
  
  for (size_t x = 0; x < n_profiles; x++){
    
    Eigen::MappedSparseMatrix<double> phi_p = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals_pred[x]);
    
    if (is_bgc[x]){
      
      Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals_resp[bgc_index]);
      
      n_total = phi.rows() + phi_p.rows();
      
      Eigen::MatrixXd temp(n_total, n_pcs);
      Eigen::VectorXd x_temp(n_total);
      Eigen::VectorXd mean(n_total);
      
      x_temp << (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs_resp[bgc_index])), Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs_pred[x]);
      profile_lengths[0] = phi.rows();
      profile_lengths.tail(n_preds) = profile_lengths_p.row(x);
      
      for (size_t g = 0; g < G; g++){
        
        mean << phi * means_resp.col(g), phi_p * means_pred.col(g);
        temp.topRows(phi.rows()) = phi * L1[g];
        temp.bottomRows(phi_p.rows()) = phi_p * L2[g];
        
        liks(x, g) = lik_eigen_sherman(x_temp, mean, variances, profile_lengths, temp);
      }
      bgc_index++;
      
    } else {
      
      n_total = phi_p.rows();
      
      Eigen::MatrixXd temp(n_total, n_pcs_pred);
      Eigen::VectorXd x_temp(n_total);
      Eigen::VectorXd mean(n_total);
      
      x_temp << Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_pred[x]);
      
      profile_lengths.tail(n_preds) = profile_lengths_p.row(x);
      
      for (size_t g = 0; g < G; g++){
        
        mean << phi_p * means_pred.col(g);
        temp = phi_p * L_pred[g];
        
        liks(x, g) = lik_eigen_sherman(x_temp, mean, variances.tail(n_preds), profile_lengths.tail(n_preds), temp);
      }
    }
  }
  return liks;
}


// [[Rcpp::export]]
Eigen::MatrixXd c_compute_E_step_likelihoods_single(const Rcpp::List profs,
                                                    const int n_profiles,
                                                    const Rcpp::List basis_evals,
                                                    const Rcpp::List Omegas,
                                                    const Eigen::Map<Eigen::MatrixXd> means,
                                                    const double me,
                                                    const Eigen::Map<Eigen::MatrixXd> vars,
                                                    const int G){
  Eigen::VectorXd x, mean;
  Eigen::MatrixXd liks(n_profiles, G);
  Eigen::MatrixXd L;
  
  Eigen::VectorXi profile_lengths = Eigen::VectorXi::Zero(1);
  Eigen::VectorXd variances = Eigen::VectorXd::Zero(1);
  variances(0) = me;
  
  std::vector<Eigen::MatrixXd> l;
  for (size_t g = 0; g < G; g++){
    // Need the square root here
    l.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas[g]) * vars.col(g).array().sqrt().matrix().asDiagonal());
  }
  
  for (size_t x = 0; x < n_profiles; x++){
    
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    
    Eigen::VectorXd x_temp = (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs[x]));
    
    for (size_t g = 0; g < G; g++){
      
      Eigen::VectorXd mean = phi * means.col(g);
      Eigen::MatrixXd temp = phi * l[g];
      
      profile_lengths(0) = phi.rows();
      
      liks(x, g) = lik_eigen_sherman(x_temp, mean, variances, profile_lengths, temp);
    }
  }
  return liks;
}

// [[Rcpp::export]]
Rcpp::List c_create_summed_U_matrix_sparse(const Rcpp::List phi_x_phi,
                                           const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                           const Eigen::Map<Eigen::VectorXd> weights,
                                           const int G,
                                           const int reps){
  
  long g_index;
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(phi_x_phi[0]).cols();
  const int n_samples = phi_x_phi.size();
  Eigen::VectorXd summed_weights = Eigen::VectorXd::Zero(G);
  
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n_basis, G);
  
  std::vector<Eigen::SparseMatrix<double> > return_matrices;
  for (size_t i = 0; i < G; i++){
    return_matrices.push_back(Eigen::SparseMatrix<double>(n_basis, n_basis));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double>>(phi_x_phi[x]);
    
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major storage 
      g_index = clust_mem(x,t);
      summed_weights(g_index) += weights(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){
        continue;
      }
      // Not sure if this can't be further optimized. 
      // Hopefully eigen recognizes that it only needs to add the lower half of 
      // the matrix. This should be tested a bit
      return_matrices[g] += summed_weights(g) * phi;
    }
  }
  // This is of course not very elegant but making typed Rcpp lists
  // doesnt work and returning the std::vector neither
  Rcpp::List ret;
  for (size_t g = 0; g < G; g++){
    ret.push_back(return_matrices[g]);
  }
  return(ret);
}

// [[Rcpp::export]]
Eigen::MatrixXd c_create_summed_V_matrix_sparse(const Rcpp::List basis_evals,
                                                const Rcpp::List profiles,
                                                const int G,
                                                const Rcpp::List Omegas1,
                                                const Rcpp::List Omegas2,
                                                const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                const Rcpp::List pcs1,
                                                const Rcpp::List pcs2,
                                                const int reps,
                                                const Eigen::Map<Eigen::VectorXd> weights){
  
  long g;
  const int n_samples = basis_evals.size();
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[0]).cols();
  const int n_components1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[0]).cols();
  const int n_components2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[0]).cols();
  Eigen::MatrixXd summed_components1 = Eigen::MatrixXd::Zero(n_components1, G);
  Eigen::MatrixXd summed_components2 = Eigen::MatrixXd::Zero(n_components2, G);
  Eigen::VectorXd summed_weights = Eigen::VectorXd::Zero(G);
  Eigen::MatrixXd V_matrix = Eigen::MatrixXd::Zero(n_basis, G);
  
  // This is of course terrible but the list won't work directly
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas1, omegas2;
  for (size_t g = 0; g < G; g++){
    omegas1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[g]));
    omegas2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[g]));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    summed_components1.setZero();
    summed_components2.setZero();
    Eigen::MappedSparseMatrix<double> phi = (Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]));
    Eigen::Map<Eigen::MatrixXd > pc_mat1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs1[x]);
    Eigen::Map<Eigen::MatrixXd > pc_mat2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs2[x]);
    
    // Sum scores belonging to the same cluster.
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major order 
      g = clust_mem(x,t);
      summed_weights(g) += weights(t);
      summed_components1.col(g) += weights(t) * pc_mat1.row(t).adjoint();
      summed_components2.col(g) += weights(t) * pc_mat2.row(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){
        continue;
      }
      // Probably want to check if this aumatically optimizes with the parantheses
      V_matrix.col(g) += phi.adjoint() * (summed_weights(g) * (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profiles[x])) - (phi * (omegas1[g] * summed_components1.col(g) + omegas2[g] * summed_components2.col(g))));
    }
  }
  return V_matrix;
}


// [[Rcpp::export]]
Eigen::MatrixXd c_create_summed_V_matrix_sparse_single(const Rcpp::List basis_evals,
                                                       const Rcpp::List profiles,
                                                       const int G,
                                                       const Rcpp::List Omegas,
                                                       const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                       const Rcpp::List pcs,
                                                       const int reps,
                                                       const Eigen::Map<Eigen::VectorXd> weights){
  
  long g;
  const int n_samples = basis_evals.size();
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[0]).cols();
  const int n_components = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[0]).cols();
  Eigen::MatrixXd summed_components = MatrixXd::Zero(n_components, G);
  Eigen::VectorXd summed_weights = VectorXd::Zero(G);
  Eigen::MatrixXd V_matrix = MatrixXd::Zero(n_basis, G);
  
  // This is of course terrible but the list won't work directly
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas;
  for (size_t i = 0; i < G; i++){
    omegas.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[i]));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    summed_components.setZero();
    Eigen::MappedSparseMatrix<double> phi = (Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]));
    Eigen::Map<Eigen::MatrixXd > pc_mat = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs[x]);
    
    // Sum scores belonging to the same cluster.
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major order 
      g = clust_mem(x,t);
      summed_weights(g) += weights(t);
      summed_components.col(g) += weights(t) * pc_mat.row(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){
        continue;
      }
      
      // Probably want to check if this aumatically optimizes with the parantheses
      V_matrix.col(g) += phi.adjoint() * (summed_weights(g) * (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profiles[x])) - (phi * (omegas[g] * summed_components.col(g))));
    }
  }
  return V_matrix;
}

// [[Rcpp::export]]
Eigen::MatrixXd c_create_summed_V_matrix_pcs_sparse(const Rcpp::List basis_evals,
                                                    const Rcpp::List profiles,
                                                    const Eigen::Map<Eigen::MatrixXd> means,
                                                    const int G,
                                                    const Rcpp::List Omegas1,
                                                    const Rcpp::List Omegas2,
                                                    const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                    const Rcpp::List pcs1,
                                                    const Rcpp::List pcs2,
                                                    const int reps,
                                                    const Eigen::Map<Eigen::VectorXd> weights,
                                                    const int q){
  long g_index;
  const int n_samples = basis_evals.size();
  const int n_basis = Rcpp::as<MappedSparseMatrix<double> >(basis_evals[0]).cols();
  const int n_components1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[0]).cols();
  const int n_components2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[0]).cols();
  Eigen::VectorXd summed_weights = Eigen::VectorXd::Zero(G);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(n_basis);
  Eigen::MatrixXd summed_components1 = Eigen::MatrixXd::Zero(n_components1, G);
  Eigen::MatrixXd summed_components2 = Eigen::MatrixXd::Zero(n_components2, G);
  Eigen::MatrixXd V_matrix = Eigen::MatrixXd::Zero(n_basis, G);
  
  // This is of course terrible but the list won't work directly
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas1, omegas2;
  for (size_t g = 0; g < G; g++){
    omegas1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[g]));
    omegas2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[g]));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    summed_components1.setZero();
    summed_components2.setZero();
    
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    Eigen::Map<Eigen::MatrixXd > pc_mat1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs1[x]);
    Eigen::Map<Eigen::MatrixXd > pc_mat2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs2[x]);
    
    // Sum scores belonging to the same cluster.
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major order 
      g_index = clust_mem(x,t);
      summed_weights(g_index) += weights(t) * pc_mat1(t, q);
      summed_components1.col(g_index) += pc_mat1(t, q) * weights(t) * pc_mat1.row(t);
      summed_components2.col(g_index) += pc_mat1(t, q) * weights(t) * pc_mat2.row(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){ continue; }
      temp.setZero();
      for (size_t l = 0; l < n_components1; l++){
        if (l == q){ continue; }
        temp += summed_components1(l, g) * omegas1[g].col(l);
      }
      V_matrix.col(g) += phi.adjoint() * (summed_weights(g) * Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profiles[x]) - (phi * (summed_weights(g) * means.col(g) + omegas2[g] * summed_components2.col(g) + temp)));
    }
  }
  return V_matrix;
}

// [[Rcpp::export]]
Eigen::MatrixXd c_create_summed_V_matrix_pcs_sparse_single(const Rcpp::List basis_evals,
                                                           const Rcpp::List profiles,
                                                           const Eigen::Map<Eigen::MatrixXd> means,
                                                           const int G,
                                                           const Rcpp::List Omegas,
                                                           const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                           const Rcpp::List pcs,
                                                           const int reps,
                                                           const Eigen::Map<Eigen::VectorXd> weights,
                                                           const int q){
  long g_index;
  const int n_samples = basis_evals.size();
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[0]).cols();
  const int n_components = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[0]).cols();
  Eigen::VectorXd summed_weights = Eigen::VectorXd::Zero(G);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(n_basis);
  Eigen::MatrixXd summed_components = Eigen::MatrixXd::Zero(n_components, G);
  Eigen::MatrixXd V_matrix = Eigen::MatrixXd::Zero(n_basis, G);
  
  // This is of course terrible but the list won't work directly
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas;
  for (size_t i = 0; i < G; i++){
    omegas.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[i]));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    summed_components.setZero();
    
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    Eigen::Map<Eigen::MatrixXd > pc_mat = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs[x]);
    
    // Sum scores belonging to the same cluster.
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major order 
      g_index = clust_mem(x,t);
      summed_weights(g_index) += weights(t) * pc_mat(t, q);
      summed_components.col(g_index) += pc_mat(t, q) * weights(t) * pc_mat.row(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){ continue; }
      temp.setZero();
      for (size_t l = 0; l < n_components; l++){
        if (l == q){ continue; }
        temp += summed_components(l, g) * omegas[g].col(l);
      }
      V_matrix.col(g) += phi.adjoint() * (summed_weights(g) * Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profiles[x]) - (phi * (summed_weights(g) * means.col(g) + temp)));
    }
  }
  return V_matrix;
}

// [[Rcpp::export]]
Rcpp::List c_create_summed_U_matrix_pcs_sparse(const Rcpp::List phi_x_phi,
                                               const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                               const Eigen::Map<Eigen::VectorXd> weights,
                                               const Eigen::Map<Eigen::MatrixXd> pc_weights,
                                               const int G,
                                               const int reps){
  
  long g_index;
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(phi_x_phi[0]).cols();
  const int n_samples = phi_x_phi.size();
  Eigen::VectorXd summed_weights = VectorXd::Zero(G);
  
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n_basis, G);
  
  std::vector<Eigen::MatrixXd > return_matrices;
  for (size_t i = 0; i < G; i++){
    return_matrices.push_back(Eigen::SparseMatrix<double>(n_basis, n_basis));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(phi_x_phi[x]);
    
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major storage 
      g_index = clust_mem(x,t);
      summed_weights(g_index) += weights(t) * pc_weights(x, t) * pc_weights(x, t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){
        continue;
      }
      // Not sure if this can't be further optimized. 
      // Hopefully eigen recognizes that it only needs to add the lower half of 
      // the matrix. This should be tested a bit
      return_matrices[g] += summed_weights(g) * phi;
    }
  }
  // This is of course not very elegant but making typed Rcpp lists
  // doesnt work and returning the std::vector neither
  Rcpp::List ret;
  for (size_t g = 0; g < G; g++){
    ret.push_back(return_matrices[g]);
  }
  return(ret);
}


// [[Rcpp::export]]
Eigen::MatrixXd c_create_summed_V_matrix_gamma_sparse_r(const Rcpp::List basis_evals,
                                                        const Rcpp::List profiles,
                                                        const Eigen::Map<Eigen::MatrixXd> means,
                                                        const int G,
                                                        const Rcpp::List Omegas1,
                                                        const Rcpp::List Omegas2,
                                                        const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                        const Rcpp::List etas,
                                                        const Rcpp::List alphas,
                                                        const int reps,
                                                        const Eigen::Map<Eigen::VectorXd> weights,
                                                        const int q){
  long g_index;
  const int n_samples = basis_evals.size();
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[0]).cols();
  const int n_components1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[0]).cols();
  const int n_components2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[0]).cols();
  Eigen::VectorXd summed_weights = Eigen::VectorXd::Zero(G);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(n_basis);
  Eigen::MatrixXd summed_components1 = Eigen::MatrixXd::Zero(n_components1, G);
  Eigen::MatrixXd summed_components2 = Eigen::MatrixXd::Zero(n_components2, G);
  Eigen::MatrixXd V_matrix = Eigen::MatrixXd::Zero(n_components1, G);
  
  // This is of course terrible but the list won't work directly
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas1, omegas2;
  std::vector<Eigen::MatrixXd> L;
  for (size_t g = 0; g < G; g++){
    omegas1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[g]));
    omegas2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[g]));
    L.push_back(omegas1[g]);
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    summed_components1.setZero();
    summed_components2.setZero();
    
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    Eigen::Map<Eigen::MatrixXd > Alpha = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(alphas[x]);
    Eigen::Map<Eigen::MatrixXd > Eta = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(etas[x]);
    
    // Sum scores belonging to the same cluster.
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major order 
      g_index = clust_mem(x,t);
      summed_weights(g_index) += weights(t) * Alpha(t, q);
      summed_components1.col(g_index) += Alpha(t, q) * weights(t) * Alpha.row(t);
      summed_components2.col(g_index) += Alpha(t, q) * weights(t) * Eta.row(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){ continue; }
      temp.setZero();
      for (size_t l = 0; l < n_components2; l++){
        if (l == q){ continue; }
        temp += summed_components1(l, g) * L[g].col(l);
      }
      V_matrix.col(g) += (phi * omegas1[g]).adjoint() * (summed_weights(g) * Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profiles[x]) - (phi * (summed_weights(g) * means.col(g) + omegas2[g] * summed_components2.col(g) + temp)));
    }
  } 
  return V_matrix;
}


// [[Rcpp::export]]
Rcpp::List c_create_summed_U_matrix_gamma_sparse(const Rcpp::List phi_x_phi,
                                                 const Rcpp::List Omegas1,
                                                 const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                 const Eigen::Map<Eigen::VectorXd> weights,
                                                 const Eigen::Map<Eigen::MatrixXd> pc_weights,
                                                 const int G,
                                                 const int reps){
  
  long g_index;
  const int n_basis = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[0]).cols();
  const int n_samples = phi_x_phi.size();
  Eigen::VectorXd summed_weights = VectorXd::Zero(G);
  
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas1;
  for (size_t g = 0; g < G; g++){
    omegas1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[g]));
  }
  
  std::vector<Eigen::MatrixXd > return_matrices;
  for (size_t i = 0; i < G; i++){
    return_matrices.push_back(Eigen::MatrixXd::Zero(n_basis, n_basis));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(phi_x_phi[x]);
    
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major storage 
      g_index = clust_mem(x,t);
      summed_weights(g_index) += weights(t) * pc_weights(x, t) * pc_weights(x, t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){
        continue;
      }
      // Not sure if this can't be further optimized. 
      // Hopefully eigen recognizes that it only needs to add the lower half of 
      // the matrix. This should be tested a bit
      return_matrices[g] += summed_weights(g) * omegas1[g].adjoint() * phi * omegas1[g];
    }
  }
  // This is of course not very elegant but making typed Rcpp lists
  // doesnt work and returning the std::vector neither
  Rcpp::List ret;
  for (size_t g = 0; g < G; g++){
    ret.push_back(return_matrices[g]);
  }
  return(ret);
}

// [[Rcpp::export]]
double c_update_measurement_error(const Eigen::Map<Eigen::MatrixXi> cluster_mat,
                                  const Rcpp::List basis_evals,
                                  const Rcpp::List pcs1,
                                  const Rcpp::List pcs2,
                                  const Rcpp::List profiles,
                                  const Eigen::Map<Eigen::MatrixXd> means_mat,
                                  const Rcpp::List Omegas1,
                                  const Rcpp::List Omegas2,
                                  const int n_profiles,
                                  const Eigen::Map<Eigen::VectorXd> weights,
                                  const int G){
  
  double res = 0;
  int n = 0;
  int g_index;
  
  std::vector<Eigen::Map<Eigen::MatrixXd> >omegas1, omegas2;
  for (size_t g = 0; g < Omegas1.size(); g++){
    omegas1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[g]));
    omegas2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[g]));
  }
  
  for (size_t x = 0; x < n_profiles; x++){
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    Eigen::Map<Eigen::MatrixXd> pcs_mat1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(pcs1[x]);
    Eigen::Map<Eigen::MatrixXd> pcs_mat2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(pcs2[x]);
    Eigen::Map<Eigen::VectorXd> prof = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profiles[x]);
    VectorXd temp = VectorXd::Zero(phi.rows());
    n += phi.rows();
    for (size_t mc = 0; mc < cluster_mat.cols(); mc++){
      g_index = cluster_mat(x, mc);
      temp = weights(mc) * (prof - phi * (means_mat.col(g_index) + omegas1[g_index] * pcs_mat1.row(mc).adjoint() + omegas2[g_index] * pcs_mat2.row(mc).adjoint())).array().square();
      res += temp.sum();
    }
  }
  
  return res / (weights.sum() * n);
}

// [[Rcpp::export]]
Eigen::VectorXd c_compute_squared_sparse(const Eigen::Map<Eigen::VectorXi> cluster_mat_i,
                                         const Eigen::Map<Eigen::MatrixXd> pcs_mat,
                                         const Eigen::Map<Eigen::VectorXd> profile,
                                         const int G,
                                         const Eigen::MappedSparseMatrix<double> basis_eval,
                                         const Eigen::Map<Eigen::MatrixXd> means_mat,
                                         const Rcpp::List Omegas,
                                         const Eigen::Map<Eigen::VectorXd> weights){
  
  Eigen::VectorXd res = Eigen::VectorXd::Zero(basis_eval.rows());
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(basis_eval.rows());
  int g;
  
  std::vector<Eigen::Map<Eigen::MatrixXd> >omegas;
  for (size_t g = 0; g < Omegas.size(); g++){
    omegas.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[g]));
  }
  
  for (size_t mc = 0; mc < cluster_mat_i.size(); mc++){
    g = cluster_mat_i(mc);
    temp = weights(mc) * (profile - basis_eval * (means_mat.col(g) + omegas[g] * pcs_mat.row(mc).adjoint())).array().square();
    res += temp;
  }
  
  return res;
}


// [[Rcpp::export]]
Eigen::VectorXd c_compute_squared_sparse_response(const Eigen::Map<Eigen::VectorXi> cluster_mat_i,
                                                  const Eigen::Map<Eigen::MatrixXd> pcs_mat1,
                                                  const Eigen::Map<Eigen::MatrixXd> pcs_mat2,
                                                  const Eigen::Map<Eigen::VectorXd> profile,
                                                  const int G,
                                                  const Eigen::MappedSparseMatrix<double> basis_eval,
                                                  const Eigen::Map<Eigen::MatrixXd> means_mat,
                                                  const Rcpp::List Lambda1,
                                                  const Rcpp::List Lambda2,
                                                  const Eigen::Map<Eigen::VectorXd> weights){
  
  Eigen::VectorXd res = Eigen::VectorXd::Zero(basis_eval.rows());
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(basis_eval.rows());
  int g;
  
  std::vector<Eigen::Map<Eigen::MatrixXd> >lambda1;
  std::vector<Eigen::Map<Eigen::MatrixXd> >lambda2;
  for (size_t g = 0; g < Lambda1.size(); g++){
    lambda1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Lambda1[g]));
    lambda2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Lambda2[g]));
  }
  
  for (size_t mc = 0; mc < cluster_mat_i.size(); mc++){
    g = cluster_mat_i(mc);
    temp = (profile - basis_eval * (means_mat.col(g) + lambda1[g] * pcs_mat1.row(mc).adjoint() + 
      lambda2[g] * pcs_mat2.row(mc).adjoint())).array().square();
    res += temp;
  }
  
  return res;
}

// [[Rcpp::export]]
double c_update_measurement_error_p_space(const Eigen::Map<Eigen::MatrixXi> cluster_mat,
                                          const Rcpp::List basis_evals,
                                          const Rcpp::List pcs,
                                          const Rcpp::List profiles,
                                          const Eigen::Map<Eigen::MatrixXd> means_mat,
                                          const Rcpp::List Omegas,
                                          const int n_profiles,
                                          const int G){
  
  double res = 0;
  int T = cluster_mat.cols();
  int n = 0;
  int g_index;
  
  std::vector<Eigen::Map<Eigen::MatrixXd> >omegas;
  for (size_t g = 0; g < Omegas.size(); g++){
    omegas.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[g]));
  }
  
  for (size_t x = 0; x < n_profiles; x++){
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    Eigen::Map<Eigen::MatrixXd> pcs_mat = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(pcs[x]);
    Eigen::Map<Eigen::VectorXd> prof = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profiles[x]);
    Eigen::VectorXd temp = Eigen::VectorXd::Zero(phi.rows());
    n += phi.rows();
    for (size_t mc = 0; mc < cluster_mat.cols(); mc++){
      g_index = cluster_mat(x, mc);
      temp = (prof - phi * (means_mat.col(g_index) + omegas[g_index] * pcs_mat.row(mc).adjoint())).array().square();
      res += temp.sum();
    }
  }
  
  return res / (T * n);
}

// [[Rcpp::export]]
void c_compute_conditional_distribution(Rcpp::List profs_resp,
                                        Rcpp::List profs_pred,
                                        Rcpp::List basis_evals_resp,
                                        Rcpp::List basis_evals_pred,
                                        Rcpp::List phi_x_phi_resp,
                                        Rcpp::List phi_x_phi_pred,
                                        const int n_samples,
                                        const Eigen::Map<Eigen::MatrixXd> means_resp,
                                        const Eigen::Map<Eigen::MatrixXd> means_pred,
                                        Rcpp::List Omegas_resp,
                                        Rcpp::List Omegas_pred,
                                        Rcpp::List Lambdas,
                                        Rcpp::List Sigma_eta_inv, 
                                        const double me_resp,
                                        const Eigen::Map<Eigen::VectorXd> me_pred,
                                        const Eigen::Map<Eigen::MatrixXd> vars_resp,
                                        const Eigen::Map<Eigen::MatrixXd> vars_pred,
                                        const Eigen::Map<Eigen::VectorXi> basis_lengths_pred,
                                        const Eigen::Map<Eigen::MatrixXd> cond_probs,
                                        const Eigen::Map<Eigen::VectorXi> is_bgc,
                                        Rcpp::List conditional_distributions){
  Eigen::MatrixXd step;
  const size_t n_resp = vars_resp.rows();
  const size_t n_pred = vars_pred.rows();
  const size_t n_total = n_resp + n_pred;
  const size_t G = Omegas_resp.size();
  int bgc_index = 0;
  const Eigen::MatrixXd id = Eigen::MatrixXd::Identity(n_total, n_total);
  
  Eigen::MatrixXd S_i_aa = Eigen::MatrixXd::Zero(n_resp, n_resp);
  Eigen::MatrixXd S_i_bb = Eigen::MatrixXd::Zero(n_resp, n_pred);
  
  std::vector<Eigen::MatrixXd> omegas_resp, omegas_pred, omegas_pred_sqrt, S_i_ba, S_i_bb_temp, S_i_bb_temp_p;
  std::vector<Eigen::Map<Eigen::MatrixXd> > sigma_eta_inv;
  for (size_t g = 0; g < G; g++){
    omegas_resp.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_resp[g]));
    omegas_pred.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_pred[g]));
    omegas_pred_sqrt.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_pred[g]));
    // Multiplying by the measurement errors here because it makes things easier downstream
    omegas_pred_sqrt[g].topRows(basis_lengths_pred[0]) *= 1 / std::sqrt(me_pred[0]);
    omegas_pred_sqrt[g].bottomRows(basis_lengths_pred[1]) *= 1 / std::sqrt(me_pred[1]);
    omegas_pred[g].topRows(basis_lengths_pred[0]) *= 1 / me_pred[0];
    omegas_pred[g].bottomRows(basis_lengths_pred[1]) *= 1 / me_pred[1];
    sigma_eta_inv.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Sigma_eta_inv[g]));
    S_i_ba.push_back(-Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Lambdas[g]).adjoint() * sigma_eta_inv[g]);
    step = vars_pred.col(g).array().inverse().matrix().asDiagonal();
    S_i_bb_temp_p.push_back(step);
    S_i_bb_temp.push_back(step - S_i_ba[g] * Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Lambdas[g]));
  }
  
  Eigen::LLT<Eigen::MatrixXd> llt;
  Eigen::VectorXd res_resp(n_resp);
  Eigen::VectorXd res_pred(n_pred);
  Eigen::MatrixXd chol(n_total, n_total);
  for (size_t x = 0; x < n_samples; x++){
    
    Eigen::MappedSparseMatrix<double> phi_p = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals_pred[x]);
    Eigen::MappedSparseMatrix<double> phi_x_phi_p = Rcpp::as<Eigen::MappedSparseMatrix<double> >(phi_x_phi_pred[x]);
    Eigen::Map<Eigen::VectorXd> prof_pred = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_pred[x]);
    
    Rcpp::List temp = conditional_distributions[x];
    if (is_bgc[x]){
      Eigen::Map<Eigen::MatrixXd> phi = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(basis_evals_resp[bgc_index]);
      Eigen::Map<Eigen::MatrixXd> phi_x_phi = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(phi_x_phi_resp[bgc_index]);
      Eigen::Map<Eigen::VectorXd> prof_resp = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_resp[bgc_index]);
      bgc_index += 1;
      for (size_t g = 0; g < G; g++){
        
        if (cond_probs(x, g) == 0) { continue; }
        
        Rcpp::List temp2 = temp[g];
        S_i_bb = S_i_bb_temp[g] + omegas_pred_sqrt[g].adjoint() * (phi_x_phi_p * omegas_pred_sqrt[g]);
        S_i_aa = sigma_eta_inv[g] + 1/me_resp * omegas_resp[g].adjoint() * (phi_x_phi * omegas_resp[g]);
        Eigen::MatrixXd new_inv(n_total, n_total); 
        new_inv << S_i_aa, S_i_ba[g].adjoint(), S_i_ba[g], S_i_bb;
        temp2[0] = new_inv;
        Eigen::MatrixXd U = llt.compute(new_inv).matrixU();
        // Somehow householder is faster than the triangular view method
        // TODO: Shouldn't this be ldlt???
        temp2[1] = U.householderQr().solve(id);
        
        // Multiplying only the sqrt with me since the gammas were already multiplied
        res_resp = (1 / me_resp) * omegas_resp[g].adjoint() * (phi.adjoint() * (prof_resp - (phi * means_resp.col(g))));
        res_pred = omegas_pred[g].adjoint() * (phi_p.adjoint() * (prof_pred - (phi_p * means_pred.col(g))));
        Eigen::VectorXd res(n_total);
        res << res_resp, res_pred;
        // Alternatively, this could be computed via the Cholesky of the inverse, speed wise almost equal
        temp2[2] = llt.solve(res);
      }
    } else {
      for (size_t g = 0; g < G; g++){
        
        if (cond_probs(x, g) == 0) { continue; }
        
        Rcpp::List temp2 = temp[g];
        
        S_i_bb = S_i_bb_temp_p[g] + omegas_pred_sqrt[g].adjoint() * (phi_x_phi_p * omegas_pred_sqrt[g]);
        
        res_pred = omegas_pred[g].adjoint() * (phi_p.adjoint() * (prof_pred - (phi_p * means_pred.col(g))));
        
        Eigen::MatrixXd U = llt.compute(S_i_bb).matrixU(); 
        
        temp2[0] = S_i_bb;
        temp2[1] = U.householderQr().solve(Eigen::MatrixXd::Identity(n_pred, n_pred));
        temp2[2] = llt.solve(res_pred);
      }
    }
  }
}

// [[Rcpp::export]]
double c_lik_eigen_sherman_pred(const Eigen::VectorXd& x,
                                const Eigen::VectorXd& mean,
                                const Eigen::VectorXd& variances,
                                const Eigen::Map<Eigen::MatrixXd>& U,
                                const Eigen::VectorXd& W){
  int dim = x.size();
  double out, log_det;
  
  Eigen::VectorXd inv_vars = variances.array().inverse();  
  Eigen::MatrixXd utau = U.adjoint() * (inv_vars.asDiagonal() * U);
  Eigen::MatrixXd W_mat = W.asDiagonal();
  Eigen::VectorXd W_inv = W.array().inverse();

  Eigen::MatrixXd M = utau * W_mat + Eigen::MatrixXd::Identity(utau.cols(), utau.cols());

  Eigen::VectorXd diff, rss;
  diff = x - mean; 
  rss = inv_vars.array() * diff.array();
  rss = U.adjoint() * rss;
  rss = M.colPivHouseholderQr().solve(rss);
  rss = U * (W_mat * rss);
  rss = inv_vars.array() * (diff - rss).array();
  out = diff.dot(rss);
  utau += W_inv.asDiagonal();
  log_det = log(utau.determinant()) + W.array().log().sum() + variances.array().log().sum(); 
  
  return -(double)dim/2.0 * log2pi - 0.5*log_det - 0.5*out;
}

// [[Rcpp::export]]
Rcpp::NumericVector stl_sort(Rcpp::NumericVector x) {
  Rcpp::NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
Eigen::VectorXd c_CG(Eigen::MappedSparseMatrix<double> A, Eigen::VectorXd& w){
  ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
  Eigen::setNbThreads(6);
  cg.setMaxIterations(3000);
  cg.setTolerance(1 / (100000000000 * w.lpNorm<2>()));
  cg.compute(A);
  Eigen::VectorXd result = cg.solve(w);
  Rcpp::Rcout << cg.error() << "\n";
  Rcpp::Rcout << w.lpNorm<2>() << "\n";
  Rcpp::Rcout << cg.iterations() << "\n";
  return result;
}

// [[Rcpp::export]]
double c_quad_form_log(Eigen::MappedSparseMatrix<double> A, Eigen::VectorXd& w, int maxiter = 1000, double tol = 0.000001, int stepsize = 10, bool precondition = true){
  MatrixXd T(maxiter+1, 2);
  VectorXd v(w.size());
  VectorXd t(w.size());
  double b, res_t;
  SelfAdjointEigenSolver<MatrixXd> es;
  
  // Quadratic forms are always positive, and so is the output of the algorithm
  double res = -10000;
  bool converged = false;
  double offset = 0;
  double diff = 0;
  
  SparseMatrix<double> Anew = A;
  if (precondition){
    VectorXd D = Anew.diagonal().array().cwiseInverse().sqrt().matrix();
    Anew = D.asDiagonal() * (Anew * D.asDiagonal());
    offset = 2 * w.dot( (D.array().log() * w.array()).matrix() );  
  } 
  
  double w_norm = w.norm();
  
  w = w / w_norm;
  v = Anew * w;
  T(0,0) = w.dot(v);
  
  
  for (size_t i = 0; i < maxiter; i++){
    v = v - T(i, 0) * w;
    b = v.norm();
    T(i,1) = b;
    
    t = w;
    w = v / b;
    v = (- b) * t;
    v = v + Anew * w;
    T(i+1,0) = w.dot( v );
    
    if ( remainder(i, stepsize) == 0 && i > 1 ){
      
      es.computeFromTridiagonal(T.col(0).head(i+1), T.col(1).head(i));
      res_t = es.eigenvectors().row(0).array().square().matrix().dot(es.eigenvalues().array().log().matrix());
      
      diff = std::abs(res_t - res);
      
      if (diff < tol){
        res = res_t;
        converged = true;
        break;
      }
      res = res_t;
    }
  }
  if (!converged){
    Rcpp::Rcout << diff << "\n";
    Rcpp::Rcout << "Did not converge: quad form\n";
  }
  return res * w_norm * w_norm - offset;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> c_reorder(Eigen::MappedSparseMatrix<double> A, Eigen::VectorXi ind){
  std::vector<Trip> trps;
  Eigen::SparseMatrix<double> Anew = A;
  
  for (int k=0; k < A.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(Anew,k); it; ++it){
      trps.push_back(Trip(ind(it.row())-1,ind(it.col())-1,it.value()));
    }
  }
  
  Anew.setFromTriplets(trps.begin(), trps.end());
  return Anew;
}


// [[Rcpp::export]]
Eigen::VectorXd c_sample_lanczos(Eigen::MappedSparseMatrix<double> A, Eigen::VectorXd w, int stepsize = 10, int maxiter = 1500, double tol = 1e-6, bool precondition=true){
  SelfAdjointEigenSolver<MatrixXd> es;
  int i = 0;
  bool converged = false;
  double diff = 0;
  
  MatrixXd T = MatrixXd::Zero(maxiter+1, 2);
  MatrixXd W = MatrixXd::Zero(w.size(), maxiter+1);
  VectorXd v(w.size());
  VectorXd t(w.size());
  VectorXd D; 
  
  VectorXd res = VectorXd::Ones(w.size());
  VectorXd res_old = VectorXd::Ones(w.size()) * 10000;
  double b;
  
  SparseMatrix<double> Anew = A;
  
  // We for any symmetric D > 0, ( D^(-1/2) (D^(-1/2)AD^(-1/2))^(-1/2) )^2 = A^(-1)
  if (precondition){
    D = Anew.diagonal().array().cwiseInverse().sqrt().matrix();
    Anew = D.asDiagonal() * (Anew * D.asDiagonal());
  } 
  
  double w_norm = w.norm();
  
  
  w = w / w_norm;
  v = Anew * w;
  
  T(0,0) = w.dot( v );
  W.col(0) = w;
  
  // Do at least  steps
  for (size_t i = 0; i < maxiter; i++){
    
    v = v - T(i, 0) * w;
    b = v.norm();
    T(i,1) = b;
    
    t = w;
    w = v / b;
    v = (- b) * t;
    v = v + Anew * w;
    T(i+1,0) = w.dot( v );
    
    W.col(i+1) = w;
    
    // Do at least a couple of iterations
    if ( remainder(i, stepsize) == 0 && i > 0){
      es.computeFromTridiagonal(T.col(0).head(i+1), T.col(1).head(i));
      res = W.leftCols(i+1) * (w_norm * (es.eigenvectors() * (es.eigenvalues().array().cwiseInverse().sqrt() * es.eigenvectors().row(0).transpose().array()).matrix()));
      if(precondition){
        res = (res.array() * D.array()).matrix();
      }
      diff = (res_old - res).lpNorm<Infinity>();
      if(diff < tol){
        converged = true;
        break;
      }
      res_old = res;
    }
  }
  if (!converged){
    Rcpp::Rcout << diff << "\n";
    Rcpp::Rcout << "Did not converge: sample.\n";
  }
  return res;
}
 /*
Eigen::VectorXd c_sample_lanczos_block(Eigen::MappedSparseMatrix<double> A, Eigen::MatrixXd w, int stepsize = 10, int maxiter = 1500, double tol = 1e-3, bool precondition=true){
  SelfAdjointEigenSolver<MatrixXd> es;
  int i = 0;
  bool converged = false;
  double diff = 0;
  
  MatrixXd T = MatrixXd::Zero(maxiter+1, 2);
  MatrixXd W = MatrixXd::Zero(w.size(), maxiter+1);
  MatrixXd v(w.size());
  VectorXd t(w.size());
  VectorXd D; 
  
  VectorXd res = VectorXd::Ones(w.size());
  VectorXd res_old = VectorXd::Ones(w.size()) * 10000;
  double b;
  
  SparseMatrix<double> Anew = A;
  
  if (precondition){
    D = Anew.diagonal().array().cwiseInverse().sqrt().matrix();
    Anew = D.asDiagonal() * (Anew * D.asDiagonal());
  } 
  
  HouseholderQR<MatrixXd> qr(w).compute(). householderQ();
  
  
  v = Anew * w;
  
  T(0,0) = w.dot( v );
  W.col(0) = w;
  
  
}
*/

// [[Rcpp::export]]
Eigen::MatrixXd c_compute_E_step_likelihoods_for_season(const Rcpp::List profs_resp,
                                                        const Rcpp::List profs_pred,
                                                        const Eigen::Map<Eigen::VectorXi> is_bgc,
                                                        const int n_profiles,
                                                        const Rcpp::List basis_evals_r,
                                                        const Rcpp::List basis_evals_p,
                                                        const Rcpp::List Omegas_r1,
                                                        const Rcpp::List Omegas_r2,
                                                        const Rcpp::List Omegas_p,
                                                        const Rcpp::List means_resp,
                                                        const Rcpp::List means_pred,
                                                        const Eigen::Map<Eigen::VectorXd> variances,
                                                        const Eigen::Map<Eigen::MatrixXd> vars_r,
                                                        const Eigen::Map<Eigen::MatrixXd> vars_p,
                                                        const Eigen::Map<Eigen::MatrixXi> profile_lengths_p,
                                                        const int G,
                                                        const int n_preds,
                                                        const Eigen::Map<Eigen::VectorXd> days){
  
  int n_total;
  int bgc_index = 0;
  const int n_pcs_r2 = Rcpp::as<Eigen::MatrixXd> (Omegas_r2[0]).cols();
  const int n_pcs_p = Rcpp::as<Eigen::MatrixXd> (Omegas_p[0]).cols();
  Eigen::VectorXd temp_season_pred = Eigen::VectorXd::Zero(Rcpp::as<Eigen::MatrixXd> (Omegas_p[0]).rows());
  Eigen::VectorXd temp_season_resp = Eigen::VectorXd::Zero(Rcpp::as<Eigen::MatrixXd> (Omegas_r2[0]).rows());
  int g;
  const int n_pcs = n_pcs_p + n_pcs_r2;
  
  Eigen::VectorXd x, mean;
  Eigen::VectorXi profile_lengths(n_preds+1);
  Eigen::MatrixXd liks(n_profiles, G);
  Eigen::MatrixXd L;
  
  std::vector<Eigen::MatrixXd> L1, L2, L3;
  std::vector<Eigen::Map<Eigen::MatrixXd> >means_resp_list;
  std::vector<Eigen::Map<Eigen::MatrixXd> >means_pred_list;
  for (size_t g = 0; g < G; g++){
    // Need the square root here
    means_resp_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means_resp[g]));
    means_pred_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means_pred[g]));
    L1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas_p[g]) * vars_p.col(g).array().sqrt().matrix().asDiagonal());
    L2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas_r1[g]) * vars_p.col(g).array().sqrt().matrix().asDiagonal());
    L3.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas_r2[g]) * vars_r.col(g).array().sqrt().matrix().asDiagonal());
  }
  
  for (size_t x = 0; x < n_profiles; x++){
    
    Eigen::MappedSparseMatrix<double> phi_p = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals_p[x]);
    
    if (is_bgc[x]){
      
      Eigen::MappedSparseMatrix<double> phi_r = Rcpp::as<Eigen::MappedSparseMatrix<double>  >(basis_evals_r[bgc_index]);
      
      n_total = phi_r.rows() + phi_p.rows();
      
      Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(n_total, n_pcs);
      Eigen::VectorXd x_temp(n_total);
      Eigen::VectorXd mean(n_total);
      
      x_temp << (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs_pred[x])), Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs_resp[bgc_index]);
      profile_lengths.head(n_preds) = profile_lengths_p.row(x);
      profile_lengths(profile_lengths.size()-1) = phi_r.rows();
      double day = days[x];
      for (size_t g = 0; g < G; g++){
        
        temp_season_pred.setZero();
        temp_season_pred += means_pred_list[g].col(0);
        temp_season_resp.setZero();
        temp_season_resp += means_resp_list[g].col(0);
        for (size_t k =1; k < means_pred_list[g].cols(); k++) {
          if (k == 1) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
          } 
          if (k == 2) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
          } 
          if (k == 3) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
          } 
          if (k == 4) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
          } 
          if (k == 5) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
          } 
          if (k == 6) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
          } 
          if (k == 7) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
          } 
          if (k == 8) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
          } 
          if (k == 9) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
          } 
          if (k == 10) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
          } 
          if (k == 11) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
          } 
          if (k == 12) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
            temp_season_resp += means_resp_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
          } 
        }
        
        mean << phi_p * temp_season_pred, phi_r * temp_season_resp;
        temp.topLeftCorner(phi_p.rows(), n_pcs_p) = phi_p * L1[g];
        temp.bottomLeftCorner(phi_r.rows(), n_pcs_p) = phi_r * L2[g];
        temp.bottomRightCorner(phi_r.rows(), n_pcs_r2) = phi_r * L3[g];
        
        liks(x, g) = lik_eigen_sherman(x_temp, mean, variances, profile_lengths, temp);
      }
      bgc_index++;
      
    } else {
      
      n_total = phi_p.rows();
      
      Eigen::MatrixXd temp(n_total, n_pcs_p);
      Eigen::VectorXd x_temp(n_total);
      Eigen::VectorXd mean(n_total);
      
      x_temp << Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_pred[x]);
      
      profile_lengths.head(n_preds) = profile_lengths_p.row(x);
      double day = days[x];
      
      for (size_t g = 0; g < G; g++){
        
        temp_season_pred.setZero();
        temp_season_pred += means_pred_list[g].col(0);
        for (size_t k =1; k < means_pred_list[g].cols(); k++) {
          if (k == 1) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
          } 
          if (k == 2) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
          } 
          if (k == 3) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
          } 
          if (k == 4) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
          } 
          if (k == 5) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
          } 
          if (k == 6) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
          } 
          if (k == 7) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
          } 
          if (k == 8) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
          } 
          if (k == 9) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
          } 
          if (k == 10) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
          } 
          if (k == 11) {
            temp_season_pred += means_pred_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
          } 
          if (k == 12) {
            temp_season_pred += means_pred_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
          } 
        }
        
        mean << phi_p * temp_season_pred;
        temp = phi_p * L1[g];
        liks(x, g) = lik_eigen_sherman(x_temp, mean, variances.head(n_preds), profile_lengths.head(n_preds), temp);
      }
    }
  }
  return liks;
}


// [[Rcpp::export]]
Eigen::VectorXd c_compute_UTX_for_season(const Rcpp::List basis_evals_r,
                              const Rcpp::List basis_evals_p,
                              const Rcpp::List profs_r,
                              const Rcpp::List profs_p,
                              const Rcpp::List Omegas_r1,
                              const Rcpp::List Omegas_r2,
                              const Rcpp::List Omegas_p,
                              const Rcpp::List means_r,
                              const Rcpp::List means_p,
                              const Eigen::Map<Eigen::VectorXi> clust_mem,
                              const double me_r,
                              const Eigen::Map<Eigen::VectorXd> me_p,
                              const Eigen::Map<Eigen::VectorXi> n_basis_p,
                              const Eigen::Map<Eigen::VectorXi> BGC,
                              const Eigen::Map<Eigen::VectorXd> days){
  
  const int n_pred = n_basis_p.size();
  const int n_samples = basis_evals_r.size();
  const int n_samples_TS = basis_evals_p.size();
  const int n_pcs_r2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas_r2[0]).cols();
  const int n_pcs_p = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas_p[0]).cols();
  const int n_basis_r = Rcpp::as<MSpMatC>(basis_evals_r[0]).cols();
  const int n_basis_pt = Rcpp::as<MSpMatC>(basis_evals_p[0]).cols();
  const int offset = n_pcs_p * n_samples_TS;
  const int G = Omegas_r1.size();
  int start, g_index;
  Eigen::VectorXd temp_season_pred = Eigen::VectorXd::Zero(Rcpp::as<Eigen::MatrixXd> (Omegas_p[0]).rows());
  Eigen::VectorXd temp_season_resp = Eigen::VectorXd::Zero(Rcpp::as<Eigen::MatrixXd> (Omegas_r2[0]).rows());
  
  
  Eigen::VectorXd res = Eigen::VectorXd::Zero(n_samples * n_pcs_r2 + n_samples_TS * n_pcs_p);
  std::vector<Eigen::Map<Eigen::MatrixXd> >means_resp_list;
  std::vector<Eigen::Map<Eigen::MatrixXd> >means_pred_list;
  std::vector<Eigen::MatrixXd > omegas_r1, omegas_r2, omegas_p;
  for (size_t g = 0; g < G; g++){
    means_resp_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means_r[g]));
    means_pred_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means_p[g]));
    omegas_r1.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas_r1[g]));
    omegas_r1[g] *= (1 / me_r);
    omegas_r2.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas_r2[g]));
    omegas_r2[g] *= (1 / me_r);
    omegas_p.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas_p[g]));
    start = 0;
    for (int bs = 0; bs < n_pred; bs++){
      omegas_p[g].block(start, 0, n_basis_p(bs), omegas_p[g].cols()) *= (1 / me_p(bs));
      start += n_basis_p(bs);
    }
  }
  
  VectorXd temp_r(n_basis_r);
  VectorXd temp_p(n_basis_pt);
  
  int start_r = 0;
  int start_p = 0;
  int BGC_index = 0;
  for (int x = 0; x < n_samples_TS; x++){
    g_index = clust_mem(x);
    MSpMatC psi = Rcpp::as<MSpMatC>(basis_evals_p[x]);
    Eigen::Map<Eigen::VectorXd> X = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_p[x]);
    
    temp_season_pred.setZero();
    temp_season_pred += means_pred_list[g_index].col(0);
    temp_season_resp.setZero();
    temp_season_resp += means_resp_list[g_index].col(0);
    double day = days[x];
    for (size_t k =1; k < means_pred_list[g_index].cols(); k++) {
      if (k == 1) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(2 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(2 * M_PI / 365.25 * day);
      }
      if (k == 2) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(2 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(2 * M_PI / 365.25 * day);
      }
      if (k == 3) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(4 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(4 * M_PI / 365.25 * day);
      }
      if (k == 4) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(4 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(4 * M_PI / 365.25 * day);
      }
      if (k == 5) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(6 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(6 * M_PI / 365.25 * day);
      }
      if (k == 6) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(6 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(6 * M_PI / 365.25 * day);
      }
      if (k == 7) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(8 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(8 * M_PI / 365.25 * day);
      } 
      if (k == 8) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(8 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(8 * M_PI / 365.25 * day);
      } 
      if (k == 9) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(10 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(10 * M_PI / 365.25 * day);
      } 
      if (k == 10) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(10 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(10 * M_PI / 365.25 * day);
      } 
      if (k == 11) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(12 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(12 * M_PI / 365.25 * day);
      } 
      if (k == 12) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(12 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(12 * M_PI / 365.25 * day);
      } 
    }
    
    
    temp_p = psi.adjoint() * (X - psi * temp_season_pred);
    if (BGC(x)) {
      MSpMatC phi = Rcpp::as<MSpMatC>(basis_evals_r[BGC_index]);
      Eigen::Map<Eigen::VectorXd> Y = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_r[BGC_index]);
      temp_r = phi.adjoint() * (Y - phi *temp_season_resp);
      res.segment(start_p, n_pcs_p) = omegas_p[g_index].adjoint() * temp_p + omegas_r1[g_index].adjoint() * temp_r;
      res.segment(offset + start_r, n_pcs_r2) = omegas_r2[g_index].adjoint() * temp_r;
      start_r += n_pcs_r2;
      start_p += n_pcs_p;
      BGC_index += 1;
    } else {
      res.segment(start_p, n_pcs_p) = omegas_p[g_index].adjoint() * temp_p;
      start_p += n_pcs_p;
    }
  }
  
  return res;
}

// [[Rcpp::export]]
double c_compute_centered_obs_for_season(const Rcpp::List basis_evals_r,
                              const Rcpp::List basis_evals_p,
                              const Rcpp::List profs_r,
                              const Rcpp::List profs_p,
                              const Rcpp::List means_r,
                              const Rcpp::List means_p,
                              const Eigen::Map<Eigen::VectorXi> clust_mem,
                              const Eigen::Map<Eigen::VectorXi> n_basis_p,
                              const Eigen::Map<Eigen::VectorXi> BGC,
                              const Eigen::Map<Eigen::VectorXd> days){
  
  const int n_samples_TS = basis_evals_p.size();
  int g_index;
  Eigen::VectorXd temp_season_pred = Eigen::VectorXd::Zero(n_basis_p.sum());
  Eigen::VectorXd temp_season_resp = Eigen::VectorXd::Zero(Rcpp::as<Eigen::MatrixXd> (means_r[0]).rows());
  
  std::vector<Eigen::Map<Eigen::MatrixXd> >means_resp_list;
  std::vector<Eigen::Map<Eigen::MatrixXd> >means_pred_list;
  for (size_t g = 0; g < means_r.size(); g++){
    means_resp_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means_r[g]));
    means_pred_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means_p[g]));
  }
  
  double res = 0;
  int BGC_index = 0;
  Eigen::VectorXd me_preds = Rcpp::as<Eigen::VectorXd>(profs_p[0]);
  for (int x = 0; x < n_samples_TS; x++){
    g_index = clust_mem(x);
    MSpMatC psi = Rcpp::as<MSpMatC>(basis_evals_p[x]);
    Eigen::Map<Eigen::VectorXd> X = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_p[x]);
    
    temp_season_pred.setZero();
    temp_season_pred += means_pred_list[g_index].col(0);
    temp_season_resp.setZero();
    temp_season_resp += means_resp_list[g_index].col(0);
    double day = days[x];
    for (size_t k =1; k < means_pred_list[g_index].cols(); k++) {
      if (k == 1) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(2 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(2 * M_PI / 365.25 * day);
      }
      if (k == 2) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(2 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(2 * M_PI / 365.25 * day);
      }
      if (k == 3) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(4 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(4 * M_PI / 365.25 * day);
      }
      if (k == 4) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(4 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(4 * M_PI / 365.25 * day);
      }
      if (k == 5) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(6 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(6 * M_PI / 365.25 * day);
      }
      if (k == 6) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(6 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(6 * M_PI / 365.25 * day);
      }
      if (k == 7) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(8 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(8 * M_PI / 365.25 * day);
      } 
      if (k == 8) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(8 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(8 * M_PI / 365.25 * day);
      } 
      if (k == 9) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(10 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(10 * M_PI / 365.25 * day);
      } 
      if (k == 10) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(10 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(10 * M_PI / 365.25 * day);
      } 
      if (k == 11) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::sin(12 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::sin(12 * M_PI / 365.25 * day);
      } 
      if (k == 12) {
        temp_season_pred += means_pred_list[g_index].col(k) * std::cos(12 * M_PI / 365.25 * day);
        temp_season_resp += means_resp_list[g_index].col(k) * std::cos(12 * M_PI / 365.25 * day);
      } 
    }
    
    res = res +  (X - psi * temp_season_pred).array().square().sum();
    
    if (BGC(x)) {
      MSpMatC phi = Rcpp::as<MSpMatC>(basis_evals_r[BGC_index]);
      Eigen::Map<Eigen::VectorXd> Y = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_r[BGC_index]);
      res = res + ((Y - phi * temp_season_resp)).array().square().sum();
      BGC_index += 1;
    } 
  }
  
  return res;
}

// [[Rcpp::export]]
Rcpp::List c_create_summed_U_matrix_pcs_sparse_for_season(const Rcpp::List phi_x_phi,
                                                          const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                          const Eigen::Map<Eigen::VectorXd> weights,
                                                          const Eigen::Map<Eigen::VectorXd> pc_weights,
                                                          const Eigen::Map<Eigen::VectorXd> pc_weights2,
                                                          const int G,
                                                          const int reps){
  
  long g_index;
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(phi_x_phi[0]).cols();
  const int n_samples = phi_x_phi.size();
  Eigen::VectorXd summed_weights = VectorXd::Zero(G);
  
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n_basis, G);
  
  std::vector<Eigen::MatrixXd > return_matrices;
  for (size_t i = 0; i < G; i++){
    return_matrices.push_back(Eigen::SparseMatrix<double>(n_basis, n_basis));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(phi_x_phi[x]);
    
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major storage 
      g_index = clust_mem(x,t);
      summed_weights(g_index) += weights(t) * pc_weights(x) * pc_weights2(x);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){
        continue;
      }
      // Not sure if this can't be further optimized. 
      // Hopefully eigen recognizes that it only needs to add the lower half of 
      // the matrix. This should be tested a bit
      return_matrices[g] += summed_weights(g) * phi;
    }
  }
  // This is of course not very elegant but making typed Rcpp lists
  // doesnt work and returning the std::vector neither
  Rcpp::List ret;
  for (size_t g = 0; g < G; g++){
    ret.push_back(return_matrices[g]);
  }
  return(ret);
}

// [[Rcpp::export]]
Eigen::MatrixXd c_create_summed_V_matrix_sparse_for_season(const Rcpp::List basis_evals,
                                                           const Rcpp::List profiles,
                                                           const int G,
                                                           const Rcpp::List Omegas1,
                                                           const Rcpp::List Omegas2,
                                                           const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                           const Rcpp::List pcs1,
                                                           const Rcpp::List pcs2,
                                                           const Eigen::Map<Eigen::VectorXd>  additional_weights,
                                                           const int reps,
                                                           const Eigen::Map<Eigen::VectorXd> weights){
  
  long g;
  const int n_samples = basis_evals.size();
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[0]).cols();
  const int n_components1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[0]).cols();
  const int n_components2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[0]).cols();
  Eigen::MatrixXd summed_components1 = Eigen::MatrixXd::Zero(n_components1, G);
  Eigen::MatrixXd summed_components2 = Eigen::MatrixXd::Zero(n_components2, G);
  Eigen::VectorXd summed_weights = Eigen::VectorXd::Zero(G);
  Eigen::MatrixXd V_matrix = Eigen::MatrixXd::Zero(n_basis, G);
  
  // This is of course terrible but the list won't work directly
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas1, omegas2;
  for (size_t g = 0; g < G; g++){
    omegas1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[g]));
    omegas2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[g]));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    summed_components1.setZero();
    summed_components2.setZero();
    Eigen::MappedSparseMatrix<double> phi = (Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]));
    Eigen::Map<Eigen::MatrixXd > pc_mat1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs1[x]);
    Eigen::Map<Eigen::MatrixXd > pc_mat2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs2[x]);
    
    // Sum scores belonging to the same cluster.
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major order 
      g = clust_mem(x,t);
      summed_weights(g) += weights(t) * additional_weights(x);
      summed_components1.col(g) += additional_weights(x) * weights(t) * pc_mat1.row(t);
      summed_components2.col(g) += additional_weights(x) * weights(t) * pc_mat2.row(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){
        continue;
      }
      // Probably want to check if this aumatically optimizes with the parantheses
      V_matrix.col(g) += phi.adjoint() * (summed_weights(g) * (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profiles[x])) - (phi * (omegas1[g] * summed_components1.col(g) + omegas2[g] * summed_components2.col(g))));
    }
  }
  return V_matrix;
}

// [[Rcpp::export]]
Eigen::MatrixXd c_create_summed_V_matrix_sparse_single_for_season(const Rcpp::List basis_evals,
                                                                  const Rcpp::List profiles,
                                                                  const int G,
                                                                  const Rcpp::List Omegas,
                                                                  const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                                  const Rcpp::List pcs,
                                                                  const Eigen::Map<Eigen::VectorXd> additional_weights,
                                                                  const int reps,
                                                                  const Eigen::Map<Eigen::VectorXd> weights){
  
  long g;
  const int n_samples = basis_evals.size();
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[0]).cols();
  const int n_components = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[0]).cols();
  Eigen::MatrixXd summed_components = MatrixXd::Zero(n_components, G);
  Eigen::VectorXd summed_weights = VectorXd::Zero(G);
  Eigen::MatrixXd V_matrix = MatrixXd::Zero(n_basis, G);
  
  // This is of course terrible but the list won't work directly
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas;
  for (size_t i = 0; i < G; i++){
    omegas.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[i]));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    
    summed_weights.setZero();
    summed_components.setZero();
    Eigen::MappedSparseMatrix<double> phi = (Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]));
    Eigen::Map<Eigen::MatrixXd > pc_mat = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs[x]);
    
    // Sum scores belonging to the same cluster.
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major order 
      g = clust_mem(x,t);
      summed_weights(g) += weights(t) * additional_weights(x);
      summed_components.col(g) += additional_weights(x) * weights(t) * pc_mat.row(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){
        continue;
      }
      
      // Probably want to check if this aumatically optimizes with the parantheses
      V_matrix.col(g) += phi.adjoint() * (summed_weights(g) * (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profiles[x])) - (phi * (omegas[g] * summed_components.col(g))));
    }
  }
  return V_matrix;
}

// [[Rcpp::export]]
Eigen::MatrixXd c_create_summed_V_matrix_pcs_sparse_single_for_season(const Rcpp::List basis_evals,
                                                           const Rcpp::List profiles,
                                                           const Rcpp::List means,
                                                           const int G,
                                                           const Rcpp::List Omegas,
                                                           const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                           const Rcpp::List pcs,
                                                           const int reps,
                                                           const Eigen::Map<Eigen::VectorXd> weights,
                                                           const int q,
                                                           const Eigen::Map<Eigen::VectorXd> days){
  long g_index;
  const int n_samples = basis_evals.size();
  const int n_basis = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[0]).cols();
  const int n_components = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[0]).cols();
  Eigen::VectorXd summed_weights = Eigen::VectorXd::Zero(G);
  Eigen::MatrixXd summed_mean = Eigen::MatrixXd::Zero(n_basis, G);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(n_basis);
  Eigen::MatrixXd summed_components = Eigen::MatrixXd::Zero(n_components, G);
  Eigen::MatrixXd V_matrix = Eigen::MatrixXd::Zero(n_basis, G);
  
  // This is of course terrible but the list won't work directly
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas;
  std::vector<Eigen::Map<Eigen::MatrixXd> >means_list;
  for (size_t i = 0; i < G; i++){
    omegas.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[i]));
    means_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means[i]));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    double day = days[x];
    summed_weights.setZero();
    summed_mean.setZero();
    summed_components.setZero();
    
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    Eigen::Map<Eigen::MatrixXd > pc_mat = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs[x]);
    
    // Sum scores belonging to the same cluster.
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major order 
      g_index = clust_mem(x,t);
      summed_weights(g_index) += weights(t) * pc_mat(t, q);
      summed_mean.col(g_index) += pc_mat(t, q)* weights(t)* means_list[g_index].col(0);
      for (size_t k = 1; k < means_list[0].cols(); k++) {
        if (k == 1) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(2 * M_PI / 365.25 * day);
        } 
        if (k == 2) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(2 * M_PI / 365.25 * day);
        } 
        if (k == 3) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(4 * M_PI / 365.25 * day);
        } 
        if (k == 4) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(4 * M_PI / 365.25 * day);
        } 
        if (k == 5) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(6 * M_PI / 365.25 * day);
        } 
        if (k == 6) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(6 * M_PI / 365.25 * day);
        } 
        if (k == 7) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(8 * M_PI / 365.25 * day);
        } 
        if (k == 8) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(8 * M_PI / 365.25 * day);
        } 
        if (k == 9) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(10 * M_PI / 365.25 * day);
        } 
        if (k == 10) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(10 * M_PI / 365.25 * day);
        } 
        if (k == 11) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(12 * M_PI / 365.25 * day);
        } 
        if (k == 12) {
          summed_mean.col(g_index) += pc_mat(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(12 * M_PI / 365.25 * day);
        } 
      }
      summed_components.col(g_index) += pc_mat(t, q) * weights(t) * pc_mat.row(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){ continue; }
      temp.setZero();
      for (size_t l = 0; l < n_components; l++){
        if (l == q){ continue; }
        temp += summed_components(l, g) * omegas[g].col(l);
      }
      V_matrix.col(g) += phi.adjoint() * (summed_weights(g) * Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profiles[x]) - (phi * (summed_mean.col(g) + temp)));
    }
  }
  return V_matrix;
}

// [[Rcpp::export]]
Eigen::MatrixXd c_create_summed_V_matrix_pcs_sparse_for_season(const Rcpp::List basis_evals,
                                                    const Rcpp::List profiles,
                                                    const Rcpp::List means,
                                                    const int G,
                                                    const Rcpp::List Omegas1,
                                                    const Rcpp::List Omegas2,
                                                    const Eigen::Map<Eigen::MatrixXi> clust_mem,
                                                    const Rcpp::List pcs1,
                                                    const Rcpp::List pcs2,
                                                    const int reps,
                                                    const Eigen::Map<Eigen::VectorXd> weights,
                                                    const int q,
                                                    const Eigen::Map<Eigen::VectorXd> days){
  long g_index;
  const int n_samples = basis_evals.size();
  const int n_basis = Rcpp::as<MappedSparseMatrix<double> >(basis_evals[0]).cols();
  const int n_components1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[0]).cols();
  const int n_components2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[0]).cols();
  Eigen::VectorXd summed_weights = Eigen::VectorXd::Zero(G);
  Eigen::MatrixXd summed_mean = Eigen::MatrixXd::Zero(n_basis, G);
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(n_basis);
  Eigen::MatrixXd summed_components1 = Eigen::MatrixXd::Zero(n_components1, G);
  Eigen::MatrixXd summed_components2 = Eigen::MatrixXd::Zero(n_components2, G);
  Eigen::MatrixXd V_matrix = Eigen::MatrixXd::Zero(n_basis, G);
  
  // This is of course terrible but the list won't work directly
  std::vector<Eigen::Map<Eigen::MatrixXd> > omegas1, omegas2, means_list;
  for (size_t g = 0; g < G; g++){
    omegas1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[g]));
    omegas2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[g]));
    means_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means[g]));
  }
  
  for (size_t x = 0; x < n_samples; x++){
    double day = days[x];
    summed_weights.setZero();
    summed_components1.setZero();
    summed_components2.setZero();
    summed_mean.setZero();
    
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    Eigen::Map<Eigen::MatrixXd > pc_mat1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs1[x]);
    Eigen::Map<Eigen::MatrixXd > pc_mat2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(pcs2[x]);
    
    // Sum scores belonging to the same cluster.
    for (size_t t = 0; t < reps; t++){
      // TODO: This should be faster in row-major order 
      g_index = clust_mem(x,t);
      summed_weights(g_index) += weights(t) * pc_mat1(t, q);
      summed_mean.col(g_index) += pc_mat1(t, q)* weights(t)* means_list[g_index].col(0);
      for (size_t k = 1; k < means_list[0].cols(); k++) {
        if (k == 1) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(2 * M_PI / 365.25 * day);
        } 
        if (k == 2) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(2 * M_PI / 365.25 * day);
        } 
        if (k == 3) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(4 * M_PI / 365.25 * day);
        } 
        if (k == 4) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(4 * M_PI / 365.25 * day);
        } 
        if (k == 5) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(6 * M_PI / 365.25 * day);
        } 
        if (k == 6) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(6 * M_PI / 365.25 * day);
        } 
        if (k == 7) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(8 * M_PI / 365.25 * day);
        } 
        if (k == 8) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(8 * M_PI / 365.25 * day);
        } 
        if (k == 9) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(10 * M_PI / 365.25 * day);
        } 
        if (k == 10) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(10 * M_PI / 365.25 * day);
        } 
        if (k == 11) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::sin(12 * M_PI / 365.25 * day);
        } 
        if (k == 12) {
          summed_mean.col(g_index) += pc_mat1(t, q) * weights(t)*means_list[g_index].col(k) * std::cos(12 * M_PI / 365.25 * day);
        } 
      }
      summed_components1.col(g_index) += pc_mat1(t, q) * weights(t) * pc_mat1.row(t);
      summed_components2.col(g_index) += pc_mat1(t, q) * weights(t) * pc_mat2.row(t);
    }
    
    for (size_t g = 0; g < G; g++){
      if (summed_weights(g) == 0){ continue; }
      temp.setZero();
      for (size_t l = 0; l < n_components1; l++){
        if (l == q){ continue; }
        temp += summed_components1(l, g) * omegas1[g].col(l);
      }
      V_matrix.col(g) += phi.adjoint() * (summed_weights(g) * Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profiles[x]) - (phi * (summed_mean.col(g) + omegas2[g] * summed_components2.col(g) + temp)));
    }
  }
  return V_matrix;
}

// [[Rcpp::export]]
Eigen::VectorXd c_compute_squared_sparse_for_season(const Eigen::Map<Eigen::VectorXi> cluster_mat_i,
                                                    const Eigen::Map<Eigen::MatrixXd> pcs_mat,
                                                    const Eigen::Map<Eigen::VectorXd> profile,
                                                    const int G,
                                                    const Eigen::MappedSparseMatrix<double> basis_eval,
                                                    const Rcpp::List means_mat,
                                                    const Rcpp::List Omegas,
                                                    const Eigen::Map<Eigen::VectorXd> weights,
                                                    const double day) {
  // std::cout << 'start';
  // std::cout << '\n';
  // unsigned int microsecond = 1000000;
  // usleep(microsecond/4);
  
  Eigen::VectorXd res = Eigen::VectorXd::Zero(basis_eval.rows());
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(basis_eval.rows());
  Eigen::VectorXd temp_season = Eigen::VectorXd::Zero(basis_eval.cols());
  int g;
  
  std::vector<Eigen::Map<Eigen::MatrixXd> >omegas;
  std::vector<Eigen::Map<Eigen::MatrixXd> >means_list;
  for (size_t g = 0; g < Omegas.size(); g++){
    omegas.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[g]));
    means_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means_mat[g]));
  }
  
  for (size_t mc = 0; mc < cluster_mat_i.size(); mc++){
    g = cluster_mat_i(mc);
    temp_season.setZero();
    temp_season += means_list[g].col(0);
    for (size_t k =1; k < means_list[g].cols(); k++) {
      if (k == 1) {
        temp_season += means_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
      } 
      if (k == 2) {
        temp_season += means_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
      } 
      if (k == 3) {
        temp_season += means_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
      } 
      if (k == 4) {
        temp_season += means_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
      } 
      if (k == 5) {
        temp_season += means_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
      } 
      if (k == 6) {
        temp_season += means_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
      } 
      if (k == 7) {
        temp_season += means_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
      } 
      if (k == 8) {
        temp_season += means_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
      } 
      if (k == 9) {
        temp_season += means_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
      } 
      if (k == 10) {
        temp_season += means_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
      } 
      if (k == 11) {
        temp_season += means_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
      } 
      if (k == 12) {
        temp_season += means_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
      } 
      
    }
    
    temp = weights(mc) * (profile - basis_eval * (temp_season + omegas[g] * pcs_mat.row(mc).adjoint())).array().square();
    res += temp;
  }
  return res;
}

// [[Rcpp::export]]
double c_update_measurement_error_for_season(const Eigen::Map<Eigen::MatrixXi> cluster_mat,
                                             const Rcpp::List basis_evals,
                                             const Rcpp::List pcs1,
                                             const Rcpp::List pcs2,
                                             const Rcpp::List profiles,
                                             const Rcpp::List means_mat,
                                             const Rcpp::List Omegas1,
                                             const Rcpp::List Omegas2,
                                             const int n_profiles,
                                             const Eigen::Map<Eigen::VectorXd> weights,
                                             const int G,
                                             const Eigen::Map<Eigen::VectorXd> days){
  double res = 0;
  int n = 0;
  int g_index;
  Eigen::VectorXd temp_season = Eigen::VectorXd::Zero(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (Omegas1[0]).rows());
  
  std::vector<Eigen::Map<Eigen::MatrixXd> >omegas1, omegas2, means_list;
  for (size_t g = 0; g < Omegas1.size(); g++){
    omegas1.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas1[g]));
    omegas2.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas2[g]));
    means_list.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(means_mat[g]));
  }

  for (size_t x = 0; x < n_profiles; x++){
    double day = days(x);
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    Eigen::Map<Eigen::MatrixXd> pcs_mat1 = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(pcs1[x]);
    Eigen::Map<Eigen::MatrixXd> pcs_mat2 = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(pcs2[x]);
    Eigen::Map<Eigen::VectorXd> prof = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profiles[x]);
    VectorXd temp = VectorXd::Zero(phi.rows());
    n += phi.rows();
    for (size_t mc = 0; mc < cluster_mat.cols(); mc++){
      g_index = cluster_mat(x, mc);
      
      temp_season.setZero();
      temp_season += means_list[g_index].col(0);
      for (size_t k =1; k < means_list[g_index].cols(); k++) {
        if (k == 1) {
          temp_season += means_list[g_index].col(k) * std::sin(2 * M_PI / 365.25 * day);
        } 
        if (k == 2) {
          temp_season += means_list[g_index].col(k) * std::cos(2 * M_PI / 365.25 * day);
        } 
        if (k == 3) {
          temp_season += means_list[g_index].col(k) * std::sin(4 * M_PI / 365.25 * day);
        } 
        if (k == 4) {
          temp_season += means_list[g_index].col(k) * std::cos(4 * M_PI / 365.25 * day);
        } 
        if (k == 5) {
          temp_season += means_list[g_index].col(k) * std::sin(6 * M_PI / 365.25 * day);
        } 
        if (k == 6) {
          temp_season += means_list[g_index].col(k) * std::cos(6 * M_PI / 365.25 * day);
        } 
        if (k == 7) {
          temp_season += means_list[g_index].col(k) * std::sin(8 * M_PI / 365.25 * day);
        } 
        if (k == 8) {
          temp_season += means_list[g_index].col(k) * std::cos(8 * M_PI / 365.25 * day);
        } 
        if (k == 9) {
          temp_season += means_list[g_index].col(k) * std::sin(10 * M_PI / 365.25 * day);
        } 
        if (k == 10) {
          temp_season += means_list[g_index].col(k) * std::cos(10 * M_PI / 365.25 * day);
        } 
        if (k == 11) {
          temp_season += means_list[g_index].col(k) * std::sin(12 * M_PI / 365.25 * day);
        } 
        if (k == 12) {
          temp_season += means_list[g_index].col(k) * std::cos(12 * M_PI / 365.25 * day);
        } 
      }
      
      temp = weights(mc) * (prof - phi * (temp_season + omegas1[g_index] * pcs_mat1.row(mc).adjoint() + omegas2[g_index] * pcs_mat2.row(mc).adjoint())).array().square();
      res += temp.sum();
    }
  }
  
  return res / (weights.sum() * n);
}

// [[Rcpp::export]]
void c_compute_conditional_distribution_for_season(Rcpp::List profs_resp,
                                        Rcpp::List profs_pred,
                                        Rcpp::List basis_evals_resp,
                                        Rcpp::List basis_evals_pred,
                                        Rcpp::List phi_x_phi_resp,
                                        Rcpp::List phi_x_phi_pred,
                                        const int n_samples,
                                        Rcpp::List means_resp,
                                        Rcpp::List means_pred,
                                        Rcpp::List Omegas_resp,
                                        Rcpp::List Omegas_pred,
                                        Rcpp::List Lambdas,
                                        Rcpp::List Sigma_eta_inv, 
                                        const double me_resp,
                                        const Eigen::Map<Eigen::VectorXd> me_pred,
                                        const Eigen::Map<Eigen::MatrixXd> vars_resp,
                                        const Eigen::Map<Eigen::MatrixXd> vars_pred,
                                        const Eigen::Map<Eigen::VectorXi> basis_lengths_pred,
                                        const Eigen::Map<Eigen::MatrixXd> cond_probs,
                                        const Eigen::Map<Eigen::VectorXi> is_bgc,
                                        const Eigen::Map<Eigen::VectorXd> days,
                                        Rcpp::List conditional_distributions){
  Eigen::MatrixXd step;
  const size_t n_resp = vars_resp.rows();
  const size_t n_pred = vars_pred.rows();
  const size_t n_total = n_resp + n_pred;
  const size_t G = Omegas_resp.size();
  int bgc_index = 0;
  const Eigen::MatrixXd id = Eigen::MatrixXd::Identity(n_total, n_total);
  
  Eigen::MatrixXd S_i_aa = Eigen::MatrixXd::Zero(n_resp, n_resp);
  Eigen::MatrixXd S_i_bb = Eigen::MatrixXd::Zero(n_resp, n_pred);
  
  std::vector<Eigen::MatrixXd> means_resp_list, means_pred_list, omegas_resp, omegas_pred, omegas_pred_sqrt, S_i_ba, S_i_bb_temp, S_i_bb_temp_p;
  std::vector<Eigen::Map<Eigen::MatrixXd> > sigma_eta_inv;
  for (size_t g = 0; g < G; g++){
    omegas_resp.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_resp[g]));
    omegas_pred.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_pred[g]));
    means_resp_list.push_back(Rcpp::as<Eigen::MatrixXd>(means_resp[g]));
    means_pred_list.push_back(Rcpp::as<Eigen::MatrixXd>(means_pred[g]));
    omegas_pred_sqrt.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_pred[g]));
    // Multiplying by the measurement errors here because it makes things easier downstream
    omegas_pred_sqrt[g].topRows(basis_lengths_pred[0]) *= 1 / std::sqrt(me_pred[0]);
    omegas_pred_sqrt[g].bottomRows(basis_lengths_pred[1]) *= 1 / std::sqrt(me_pred[1]);
    omegas_pred[g].topRows(basis_lengths_pred[0]) *= 1 / me_pred[0];
    omegas_pred[g].bottomRows(basis_lengths_pred[1]) *= 1 / me_pred[1];
    sigma_eta_inv.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Sigma_eta_inv[g]));
    S_i_ba.push_back(-Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Lambdas[g]).adjoint() * sigma_eta_inv[g]);
    step = vars_pred.col(g).array().inverse().matrix().asDiagonal();
    S_i_bb_temp_p.push_back(step);
    S_i_bb_temp.push_back(step - S_i_ba[g] * Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Lambdas[g]));
  }
  Eigen::VectorXd temp_resp_season = Eigen::VectorXd::Zero(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (Omegas_resp[0]).rows());
  Eigen::VectorXd temp_pred_season = Eigen::VectorXd::Zero(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (Omegas_pred[0]).rows());
  
  Eigen::LLT<Eigen::MatrixXd> llt;
  Eigen::VectorXd res_resp(n_resp);
  Eigen::VectorXd res_pred(n_pred);
  Eigen::MatrixXd chol(n_total, n_total);
  for (size_t x = 0; x < n_samples; x++){
    double day = days(x);
    
    Eigen::MappedSparseMatrix<double> phi_p = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals_pred[x]);
    Eigen::MappedSparseMatrix<double> phi_x_phi_p = Rcpp::as<Eigen::MappedSparseMatrix<double> >(phi_x_phi_pred[x]);
    Eigen::Map<Eigen::VectorXd> prof_pred = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_pred[x]);

    
    
    Rcpp::List temp = conditional_distributions[x];
    if (is_bgc[x]){
      Eigen::Map<Eigen::MatrixXd> phi = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(basis_evals_resp[bgc_index]);
      Eigen::Map<Eigen::MatrixXd> phi_x_phi = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(phi_x_phi_resp[bgc_index]);
      Eigen::Map<Eigen::VectorXd> prof_resp = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_resp[bgc_index]);
      bgc_index += 1;
      for (size_t g = 0; g < G; g++){
        
        if (cond_probs(x, g) == 0) { continue; }
        temp_resp_season.setZero();
        temp_pred_season.setZero();
        temp_resp_season += means_resp_list[g].col(0);
        temp_pred_season += means_pred_list[g].col(0);
        for (size_t k =1; k < means_resp_list[g].cols(); k++) {
          if (k == 1) {
            temp_resp_season += means_resp_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
          } 
          if (k == 2) {
            temp_resp_season += means_resp_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
          } 
          if (k == 3) {
            temp_resp_season += means_resp_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
          } 
          if (k == 4) {
            temp_resp_season += means_resp_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
          } 
          if (k == 5) {
            temp_resp_season += means_resp_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
          } 
          if (k == 6) {
            temp_resp_season += means_resp_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
          } 
          if (k == 7) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
            temp_resp_season += means_resp_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
          } 
          if (k == 8) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
            temp_resp_season += means_resp_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
          } 
          if (k == 9) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
            temp_resp_season += means_resp_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
          } 
          if (k == 10) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
            temp_resp_season += means_resp_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
          } 
          if (k == 11) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
            temp_resp_season += means_resp_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
          } 
          if (k == 12) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
            temp_resp_season += means_resp_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
          } 
        }
        
        Rcpp::List temp2 = temp[g];
        S_i_bb = S_i_bb_temp[g] + omegas_pred_sqrt[g].adjoint() * (phi_x_phi_p * omegas_pred_sqrt[g]);
        S_i_aa = sigma_eta_inv[g] + 1/me_resp * omegas_resp[g].adjoint() * (phi_x_phi * omegas_resp[g]);
        Eigen::MatrixXd new_inv(n_total, n_total); 
        new_inv << S_i_aa, S_i_ba[g].adjoint(), S_i_ba[g], S_i_bb;
        temp2[0] = new_inv;
        Eigen::MatrixXd U = llt.compute(new_inv).matrixU();
        // Somehow householder is faster than the triangular view method
        // TODO: Shouldn't this be ldlt???
        temp2[1] = U.householderQr().solve(id);
        
        // Multiplying only the sqrt with me since the gammas were already multiplied
        res_resp = (1 / me_resp) * omegas_resp[g].adjoint() * (phi.adjoint() * (prof_resp - (phi * temp_resp_season)));
        res_pred = omegas_pred[g].adjoint() * (phi_p.adjoint() * (prof_pred - (phi_p * temp_pred_season)));
        Eigen::VectorXd res(n_total);
        res << res_resp, res_pred;
        // Alternatively, this could be computed via the Cholesky of the inverse, speed wise almost equal
        temp2[2] = llt.solve(res);
      }
    } else {
      for (size_t g = 0; g < G; g++){
        
        if (cond_probs(x, g) == 0) { continue; }
        
        temp_pred_season.setZero();
        temp_pred_season += means_pred_list[g].col(0);
        for (size_t k =1; k < means_pred_list[g].cols(); k++) {
          if (k == 1) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
          } 
          if (k == 2) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
          } 
          if (k == 3) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
          } 
          if (k == 4) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
          } 
          if (k == 5) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
          } 
          if (k == 6) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
          } 
          if (k == 7) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
          } 
          if (k == 8) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
          } 
          if (k == 9) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
          } 
          if (k == 10) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
          } 
          if (k == 11) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
          } 
          if (k == 12) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
          } 
        }
        
        Rcpp::List temp2 = temp[g];
        
        S_i_bb = S_i_bb_temp_p[g] + omegas_pred_sqrt[g].adjoint() * (phi_x_phi_p * omegas_pred_sqrt[g]);
        
        res_pred = omegas_pred[g].adjoint() * (phi_p.adjoint() * (prof_pred - (phi_p * temp_pred_season)));
        
        Eigen::MatrixXd U = llt.compute(S_i_bb).matrixU(); 
        
        temp2[0] = S_i_bb;
        temp2[1] = U.householderQr().solve(Eigen::MatrixXd::Identity(n_pred, n_pred));
        temp2[2] = llt.solve(res_pred);
      }
    }
  }
}


// [[Rcpp::export]]
Eigen::MatrixXd c_compute_E_step_likelihoods_ind_for_season(const Rcpp::List profs_resp,
                                                 const Rcpp::List profs_pred,
                                                 const Eigen::Map<Eigen::VectorXi> is_bgc,
                                                 const int n_profiles,
                                                 const Rcpp::List basis_evals_resp,
                                                 const Rcpp::List basis_evals_pred,
                                                 Rcpp::List Gammas,
                                                 const Rcpp::List Omegas_resp,
                                                 const Rcpp::List Omegas_preds,
                                                 const Rcpp::List means_resp,
                                                 const Rcpp::List means_pred,
                                                 const Eigen::Map<Eigen::VectorXd> variances,
                                                 const Eigen::Map<Eigen::MatrixXd> vars_pred,
                                                 const Eigen::Map<Eigen::MatrixXi> profile_lengths_p,
                                                 const Eigen::Map<Eigen::VectorXd> days,
                                                 const int G,
                                                 const int n_preds){
  
  int n_total;
  int bgc_index = 0;
  const int n_pcs_resp = Rcpp::as<Eigen::MatrixXd> (Omegas_resp[0]).cols();
  const int n_pcs_pred = Rcpp::as<Eigen::MatrixXd> (Omegas_preds[0]).cols();
  const int n_pcs = n_pcs_resp + n_pcs_pred;
  Eigen::VectorXd x, mean;
  Eigen::VectorXi profile_lengths(n_preds+1);
  Eigen::MatrixXd liks(n_profiles, G);
  Eigen::MatrixXd L;
  
  std::vector<Eigen::MatrixXd> L1, L2, L_pred;
  for (size_t g = 0; g < G; g++){
    // Need the square root here
    L = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Gammas[g]).llt().matrixL();
    L1.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_resp[g]) * L.topRows(n_pcs_resp));
    L2.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_preds[g]) * L.bottomRows(n_pcs_pred));
    L_pred.push_back(Rcpp::as<Eigen::MatrixXd>(Omegas_preds[g]) * vars_pred.col(g).array().sqrt().matrix().asDiagonal());
  }
  
  std::vector<Eigen::MatrixXd> means_resp_list, means_pred_list;
  std::vector<Eigen::Map<Eigen::MatrixXd> > sigma_eta_inv;
  for (size_t g = 0; g < G; g++){
    means_resp_list.push_back(Rcpp::as<Eigen::MatrixXd>(means_resp[g]));
    means_pred_list.push_back(Rcpp::as<Eigen::MatrixXd>(means_pred[g]));
  }
  Eigen::VectorXd temp_resp_season = Eigen::VectorXd::Zero(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (Omegas_resp[0]).rows());
  Eigen::VectorXd temp_pred_season = Eigen::VectorXd::Zero(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (Omegas_preds[0]).rows());
  
  
  for (size_t x = 0; x < n_profiles; x++){
    double day = days(x);
    Eigen::MappedSparseMatrix<double> phi_p = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals_pred[x]);
    
    if (is_bgc[x]){
      
      Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals_resp[bgc_index]);
      
      n_total = phi.rows() + phi_p.rows();
      
      Eigen::MatrixXd temp(n_total, n_pcs);
      Eigen::VectorXd x_temp(n_total);
      Eigen::VectorXd mean(n_total);
      
      x_temp << (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs_resp[bgc_index])), Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs_pred[x]);
      profile_lengths[0] = phi.rows();
      profile_lengths.tail(n_preds) = profile_lengths_p.row(x);
      
      for (size_t g = 0; g < G; g++){
        temp_resp_season.setZero();
        temp_pred_season.setZero();
        temp_resp_season += means_resp_list[g].col(0);
        temp_pred_season += means_pred_list[g].col(0);
        for (size_t k =1; k < means_resp_list[g].cols(); k++) {
          if (k == 1) {
            temp_resp_season += means_resp_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
          } 
          if (k == 2) {
            temp_resp_season += means_resp_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
          } 
          if (k == 3) {
            temp_resp_season += means_resp_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
          } 
          if (k == 4) {
            temp_resp_season += means_resp_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
          } 
          if (k == 5) {
            temp_resp_season += means_resp_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
          } 
          if (k == 6) {
            temp_resp_season += means_resp_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
          } 
          if (k == 7) {
            temp_resp_season += means_resp_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
          } 
          if (k == 8) {
            temp_resp_season += means_resp_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
          } 
          if (k == 9) {
            temp_resp_season += means_resp_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
          } 
          if (k == 10) {
            temp_resp_season += means_resp_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
          } 
          if (k == 11) {
            temp_resp_season += means_resp_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
          } 
          if (k == 12) {
            temp_resp_season += means_resp_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
            temp_pred_season += means_pred_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
          } 
        }
        mean << phi * temp_resp_season, phi_p * temp_pred_season;
        temp.topRows(phi.rows()) = phi * L1[g];
        temp.bottomRows(phi_p.rows()) = phi_p * L2[g];
        
        liks(x, g) = lik_eigen_sherman(x_temp, mean, variances, profile_lengths, temp);
      }
      bgc_index++;
      
    } else {
      
      n_total = phi_p.rows();
      
      Eigen::MatrixXd temp(n_total, n_pcs_pred);
      Eigen::VectorXd x_temp(n_total);
      Eigen::VectorXd mean(n_total);
      
      x_temp << Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs_pred[x]);
      
      profile_lengths.tail(n_preds) = profile_lengths_p.row(x);
      
      for (size_t g = 0; g < G; g++){
        temp_pred_season.setZero();
        temp_pred_season += means_pred_list[g].col(0);
        for (size_t k =1; k < means_pred_list[g].cols(); k++) {
          if (k == 1) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
          } 
          if (k == 2) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
          } 
          if (k == 3) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
          } 
          if (k == 4) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
          } 
          if (k == 5) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
          } 
          if (k == 6) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
          } 
          if (k == 7) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
          } 
          if (k == 8) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
          } 
          if (k == 9) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
          } 
          if (k == 10) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
          } 
          if (k == 11) {
            temp_pred_season += means_pred_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
          } 
          if (k == 12) {
            temp_pred_season += means_pred_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
          } 
        }
        mean << phi_p * temp_pred_season;
        temp = phi_p * L_pred[g];
        
        liks(x, g) = lik_eigen_sherman(x_temp, mean, variances.tail(n_preds), profile_lengths.tail(n_preds), temp);
      }
    }
  }
  return liks;
}

// [[Rcpp::export]]
Eigen::MatrixXd c_compute_E_step_likelihoods_single_for_season(const Rcpp::List profs,
                                                    const int n_profiles,
                                                    const Rcpp::List basis_evals,
                                                    const Rcpp::List Omegas,
                                                    const Rcpp::List means,
                                                    const double me,
                                                    const Eigen::Map<Eigen::MatrixXd> vars,
                                                    const int G,
                                                    const Eigen::Map<Eigen::VectorXd> days){
  Eigen::VectorXd x, mean;
  Eigen::MatrixXd liks(n_profiles, G);
  Eigen::MatrixXd L;
  
  Eigen::VectorXi profile_lengths = Eigen::VectorXi::Zero(1);
  Eigen::VectorXd variances = Eigen::VectorXd::Zero(1);
  variances(0) = me;
  
  std::vector<Eigen::MatrixXd> l;
  for (size_t g = 0; g < G; g++){
    // Need the square root here
    l.push_back(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Omegas[g]) * vars.col(g).array().sqrt().matrix().asDiagonal());
  }
  
  std::vector<Eigen::MatrixXd> means_list;
  for (size_t g = 0; g < G; g++){
    means_list.push_back(Rcpp::as<Eigen::MatrixXd>(means[g]));
  }
  Eigen::VectorXd temp_season = Eigen::VectorXd::Zero(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (Omegas[0]).rows());

  
  for (size_t x = 0; x < n_profiles; x++){
    double day = days(x);
    Eigen::MappedSparseMatrix<double> phi = Rcpp::as<Eigen::MappedSparseMatrix<double> >(basis_evals[x]);
    
    Eigen::VectorXd x_temp = (Rcpp::as<Eigen::Map<Eigen::VectorXd> >(profs[x]));
    
    for (size_t g = 0; g < G; g++){
      
      temp_season.setZero();
      temp_season += means_list[g].col(0);
      for (size_t k =1; k < means_list[g].cols(); k++) {
        if (k == 1) {
          temp_season += means_list[g].col(k) * std::sin(2 * M_PI / 365.25 * day);
        } 
        if (k == 2) {
          temp_season += means_list[g].col(k) * std::cos(2 * M_PI / 365.25 * day);
        } 
        if (k == 3) {
          temp_season += means_list[g].col(k) * std::sin(4 * M_PI / 365.25 * day);
        } 
        if (k == 4) {
          temp_season += means_list[g].col(k) * std::cos(4 * M_PI / 365.25 * day);
        } 
        if (k == 5) {
          temp_season += means_list[g].col(k) * std::sin(6 * M_PI / 365.25 * day);
        } 
        if (k == 6) {
          temp_season += means_list[g].col(k) * std::cos(6 * M_PI / 365.25 * day);
        } 
        if (k == 7) {
          temp_season += means_list[g].col(k) * std::sin(8 * M_PI / 365.25 * day);
        } 
        if (k == 8) {
          temp_season += means_list[g].col(k) * std::cos(8 * M_PI / 365.25 * day);
        } 
        if (k == 9) {
          temp_season += means_list[g].col(k) * std::sin(10 * M_PI / 365.25 * day);
        } 
        if (k == 10) {
          temp_season += means_list[g].col(k) * std::cos(10 * M_PI / 365.25 * day);
        } 
        if (k == 11) {
          temp_season += means_list[g].col(k) * std::sin(12 * M_PI / 365.25 * day);
        } 
        if (k == 12) {
          temp_season += means_list[g].col(k) * std::cos(12 * M_PI / 365.25 * day);
        } 
      }
      Eigen::VectorXd mean = phi * temp_season;
      Eigen::MatrixXd temp = phi * l[g];
      
      profile_lengths(0) = phi.rows();
      
      liks(x, g) = lik_eigen_sherman(x_temp, mean, variances, profile_lengths, temp);
    }
  }
  return liks;
}

// [[Rcpp::export]]
Eigen::VectorXd c_compute_UTX_single_for_season(const Rcpp::List basis_evals,
                                     const Rcpp::List profs,
                                     const Rcpp::List Omegas,
                                     const Rcpp::List means,
                                     const Eigen::Map<Eigen::VectorXi> clust_mem,
                                     const double me,
                                     const Eigen::Map<Eigen::VectorXd> days){
  
  const int n_samples = basis_evals.size();
  const int n_pcs = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Omegas[0]).cols();
  const int n_basis = Rcpp::as<MSpMatC>(basis_evals[0]).cols();
  const int G = Omegas.size();
  int g_index;
  
  
  Eigen::VectorXd res = Eigen::VectorXd::Zero(n_samples * n_pcs);
  
  std::vector<Eigen::MatrixXd > omegas;
  for (size_t g = 0; g < G; g++){
    omegas.push_back(Rcpp::as<Eigen::MatrixXd >(Omegas[g]));
    omegas[g] *= (1 / me);
  }
  std::vector<Eigen::MatrixXd> means_list;
  for (size_t g = 0; g < G; g++){
    means_list.push_back(Rcpp::as<Eigen::MatrixXd>(means[g]));
  }
  Eigen::VectorXd temp_season = Eigen::VectorXd::Zero(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (Omegas[0]).rows());
  
  
  
  Eigen::VectorXd temp(n_basis);
  
  int start = 0;
  
  for (int x = 0; x < n_samples; x++){
    double day = days(x);
    
    g_index = clust_mem(x);
    MSpMatC phi = Rcpp::as<MSpMatC>(basis_evals[x]);
    Eigen::Map<Eigen::VectorXd> X = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs[x]);
    
    temp_season.setZero();
    temp_season += means_list[g_index].col(0);
    for (size_t k =1; k < means_list[g_index].cols(); k++) {
      if (k == 1) {
        temp_season += means_list[g_index].col(k) * std::sin(2 * M_PI / 365.25 * day);
      } 
      if (k == 2) {
        temp_season += means_list[g_index].col(k) * std::cos(2 * M_PI / 365.25 * day);
      } 
      if (k == 3) {
        temp_season += means_list[g_index].col(k) * std::sin(4 * M_PI / 365.25 * day);
      } 
      if (k == 4) {
        temp_season += means_list[g_index].col(k) * std::cos(4 * M_PI / 365.25 * day);
      } 
      if (k == 5) {
        temp_season += means_list[g_index].col(k) * std::sin(6 * M_PI / 365.25 * day);
      } 
      if (k == 6) {
        temp_season += means_list[g_index].col(k) * std::cos(6 * M_PI / 365.25 * day);
      } 
      if (k == 7) {
        temp_season += means_list[g_index].col(k) * std::sin(8 * M_PI / 365.25 * day);
      } 
      if (k == 8) {
        temp_season += means_list[g_index].col(k) * std::cos(8 * M_PI / 365.25 * day);
      } 
      if (k == 9) {
        temp_season += means_list[g_index].col(k) * std::sin(10 * M_PI / 365.25 * day);
      } 
      if (k == 10) {
        temp_season += means_list[g_index].col(k) * std::cos(10 * M_PI / 365.25 * day);
      } 
      if (k == 11) {
        temp_season += means_list[g_index].col(k) * std::sin(12 * M_PI / 365.25 * day);
      } 
      if (k == 12) {
        temp_season += means_list[g_index].col(k) * std::cos(12 * M_PI / 365.25 * day);
      } 
    }
    
    temp = phi.adjoint() * (X - phi * temp_season);
    res.segment(start, n_pcs) = omegas[g_index].adjoint() * temp;
    start += n_pcs;
  }
  
  return res;
}


// [[Rcpp::export]]
double c_compute_centered_obs_single_for_season(const Rcpp::List basis_evals,
                                     const Rcpp::List profs,
                                     const Eigen::Map<Eigen::VectorXi> clust_mem,
                                     const Rcpp::List means,
                                     const double me, 
                                     const int G,
                                     const Eigen::Map<Eigen::VectorXd> days){
  
  const int n_samples = basis_evals.size();
  int g_index;
  
  std::vector<Eigen::MatrixXd> means_list;
  for (size_t g = 0; g < G; g++){
    means_list.push_back(Rcpp::as<Eigen::MatrixXd>(means[g]));
  }
  Eigen::VectorXd temp_season = Eigen::VectorXd::Zero(Rcpp::as<Eigen::Map<Eigen::MatrixXd> > (means[0]).rows());
  
  
  
  double res = 0;
  for (int x = 0; x < n_samples; x++){
    double day = days(x);
    g_index = clust_mem(x);
    MSpMatC phi = Rcpp::as<MSpMatC>(basis_evals[x]);
    Eigen::Map<Eigen::VectorXd> X = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(profs[x]);
    
    temp_season.setZero();
    temp_season += means_list[g_index].col(0);
    for (size_t k =1; k < means_list[g_index].cols(); k++) {
      if (k == 1) {
        temp_season += means_list[g_index].col(k) * std::sin(2 * M_PI / 365.25 * day);
      } 
      if (k == 2) {
        temp_season += means_list[g_index].col(k) * std::cos(2 * M_PI / 365.25 * day);
      } 
      if (k == 3) {
        temp_season += means_list[g_index].col(k) * std::sin(4 * M_PI / 365.25 * day);
      } 
      if (k == 4) {
        temp_season += means_list[g_index].col(k) * std::cos(4 * M_PI / 365.25 * day);
      } 
      if (k == 5) {
        temp_season += means_list[g_index].col(k) * std::sin(6 * M_PI / 365.25 * day);
      } 
      if (k == 6) {
        temp_season += means_list[g_index].col(k) * std::cos(6 * M_PI / 365.25 * day);
      } 
      if (k == 7) {
        temp_season += means_list[g_index].col(k) * std::sin(8 * M_PI / 365.25 * day);
      } 
      if (k == 8) {
        temp_season += means_list[g_index].col(k) * std::cos(8 * M_PI / 365.25 * day);
      } 
      if (k == 9) {
        temp_season += means_list[g_index].col(k) * std::sin(10 * M_PI / 365.25 * day);
      } 
      if (k == 10) {
        temp_season += means_list[g_index].col(k) * std::cos(10 * M_PI / 365.25 * day);
      } 
      if (k == 11) {
        temp_season += means_list[g_index].col(k) * std::sin(12 * M_PI / 365.25 * day);
      } 
      if (k == 12) {
        temp_season += means_list[g_index].col(k) * std::cos(12 * M_PI / 365.25 * day);
      } 
    }
    
    res = res + ((X - phi * temp_season)/std::sqrt(me)).array().square().sum();
  }
  return res;
}

