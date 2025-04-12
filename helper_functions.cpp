#include <RcppArmadillo.h>
#include <Rcpp.h>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat compute_B(const arma::mat& Lambda, const double& sigma_sq) {
  arma::mat Lambda_outer = Lambda * Lambda.t(); 
  int p = Lambda.n_rows; 
  arma::mat B = arma::zeros<arma::mat>(p, p);
  
  for (int j = 0; j < p - 1; ++j) {
    for (int l = j + 1; l < p; ++l) {
      B(j, l) = sqrt(1 + (Lambda_outer(j, j) * Lambda_outer(l, l) + Lambda_outer(j, l) * Lambda_outer(j, l)) / 
        (sigma_sq * Lambda_outer(l, l) + sigma_sq * Lambda_outer(j, j)));
      B(l, j) = B(j, l); 
    }
    B(j, j) = sqrt(1 + Lambda_outer(j, j) / (2 * sigma_sq));
  }
  B(p - 1, p - 1) = sqrt(1 + Lambda_outer(p - 1, p - 1) / (2 * sigma_sq));
  
  return B;
}

// [[Rcpp::export]]
List posterior_samples(const arma::mat& Y, 
                   const arma::mat& Lambda_C_hat,
                   const arma::mat& Lambda_N_hat,
                   double tau_C,
                   double tau_N,
                   double sigma_sq_hat, 
                   const arma::mat P_C,
                   double v_0 = 1, 
                   double sigma_sq_0 = 1, 
                   int N_mc = 100, 
                   double rho_C = 1, 
                   double rho_N = 1) {
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  int k = Lambda_C_hat.n_cols;

  arma::cube Lambda_samples(p, k, N_mc, arma::fill::zeros);

  arma::mat Lambda_t(p, k, arma::fill::zeros);
  arma::vec sigma_sq_samples(N_mc, arma::fill::zeros);
  double residuals;
  
  arma::mat Q_C = arma::eye(p, p) - P_C;
  
  double SST = arma::accu(arma::square(Y));
  double SSE = (n + 1/tau_C) * arma::accu(arma::square(Lambda_C_hat)) + 
    (n + 1/tau_N) * arma::accu(arma::square(Lambda_N_hat));
  double v_n = v_0 + n*p;
  double gamma_n = v_0 * sigma_sq_0 + SST - SSE;

  //arma::mat D = arma::diagmat((residuals + gamma_0 * delta_0) / (gamma_0 + n - 2)) * k * rho * rho / (n + 1 / tau_2);
  //Lambda_outer_mean = mu_js * mu_js.t() + D;
  
  // posterior sampling
  for (int it = 0; it < N_mc; ++it) {
    sigma_sq_samples(it) = 1.0 / R::rgamma(v_n / 2, 2.0 / gamma_n);
    
    arma::mat eps_C = sqrt(sigma_sq_samples(it)) * rho_C * arma::randn<arma::mat>(p, k) / sqrt(n + 1/tau_C);
    arma::mat eps_C_P = P_C * eps_C;
    arma::mat eps_N = sqrt(sigma_sq_samples(it)) * rho_N * arma::randn<arma::mat>(p, k) / sqrt(n + 1/tau_N);
    arma::mat eps_N_Q = Q_C * eps_N;
    
    Lambda_t = Lambda_C_hat + Lambda_N_hat + eps_C_P + eps_N_Q;
    
    Lambda_samples.slice(it) = Lambda_t;
  }
  
  return List::create(
    Named("Lambda_samples") = Lambda_samples,
    Named("sigma_sq_samples") = sigma_sq_samples
  );
}


// [[Rcpp::export]]
arma::cube sample_Lambda_outer(arma::cube Lambda_samples){
  
  int n_MC = Lambda_samples.n_slices;
  int k = Lambda_samples.n_cols;
  int p = Lambda_samples.n_rows;
  
  arma::mat Lambda(p, k);
  arma::mat Lambda_outer(p, p);
  arma::cube Lambda_outer_samples(p, p, n_MC);
  
  for(int s=0; s<n_MC; ++s) {
    Lambda = Lambda_samples.slice(s);
    Lambda_outer = Lambda * Lambda.t();
    Lambda_outer_samples.slice(s) = Lambda_outer;
  }
  
  return Lambda_outer_samples;
}



