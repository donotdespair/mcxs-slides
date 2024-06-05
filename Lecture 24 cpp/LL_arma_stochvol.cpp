#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 Functions:
 cholesky_tridiagonal
 forward_algorithm
 backward_algorithm
 are taken from R package stochvol by Hosszejni, Kastner (2016, JSS)Dealing with Stochastic Volatility in Time Series Using the R Package stochvol
 Function cholesky_tridiagonal has been modified: output type changed to List
 Functions forward_algorithm and backward_algorithm have not been modified
 */
List cholesky_tridiagonal(
    const vec& omega_diag,
    const double omega_offdiag) {
  const int T = omega_diag.n_elem - 1;
  vec chol_diag(T+1);
  vec chol_offdiag(T+1);
  chol_diag[0] = std::sqrt(omega_diag[0]);
  for (int j = 1; j < T+1; j++) {
    chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
    chol_diag[j] = std::sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
  }
  return List::create(_["chol_diag"]=chol_diag, _["chol_offdiag"]=chol_offdiag);
}

vec forward_algorithm(
    const vec& chol_diag,
    const vec& chol_offdiag,
    const vec& covector) {
  const int T = chol_diag.n_elem - 1;
  vec htmp(T+1);
  htmp[0] = covector[0]/chol_diag[0];
  for (int j = 1; j < T+1; j++) {
    htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
  }
  return htmp;
}

vec backward_algorithm(
    const vec& chol_diag,
    const vec& chol_offdiag,
    const vec& htmp) {
  const int T = chol_diag.size() - 1;
  vec h(T+1);
  h[T] = htmp[T] / chol_diag[T];
  for (int j = T-1; j >= 0; j--) {
    h[j] = (htmp[j] - chol_offdiag[j] * h[j+1]) / chol_diag[j];
  }
  return h;
}

// [[Rcpp::export]]
List LL_arma_stochvol(vec y, 
                      int S = 10,
                      Nullable<List>  starting_values = R_NilValue,
                      NumericVector   Hyper = NumericVector::create(10,1,3)) {
  
  vec hyper = as<vec>(Hyper);                 // read argument with a default value from the input
  int T = y.n_rows;                           // read sample size
  
  RNGScope scp;                               // facilitates reading R's set.seed()
  
  vec    aux_mu(T, fill::zeros);                           // create auxiliary objects with the current MCMC state
  double aux_mu0      = 0;                              // and fill them with their default values
  double aux_sigma2   = 1;
  double aux_sigma2m  = 1;
  
  if (starting_values.isNotNull()){                     // overwrite these values if starting_values are provided as inputs
    List Starting_values(starting_values);              // read Starting_values from the function argument
    aux_mu            = as<vec>(Starting_values["mu"]);
    aux_mu0           = Starting_values["mu0"];
    aux_sigma2        = Starting_values["sigma2"];
    aux_sigma2m       = Starting_values["sigma2"];
  } 
  
  mat posterior_mu(T, S);                     // create output objects containing posterior draws
  vec posterior_mu0(S);
  vec posterior_sigma2(S);
  vec posterior_sigma2m(S);
  
  vec e1(T);
  e1(0) = 1;
  
  mat H(T, T, fill::eye);                     // create the H matrix
  H.diag(-1) += -1;
  mat HH = H.t() * H;                         // ... and its crossprod
  mat IT(T, T, fill::eye);                    // identity matrix
  
  
  for (int s=0; s<S; s++){                    // Gibbs sampler begins here
    // sample mu0
    double vm     = 1/((1/hyper[0]) + (1/aux_sigma2m));                           // variance
    aux_mu0       = R::rnorm((vm * aux_mu[0])/aux_sigma2m, sqrt(vm));             // sample mu0 R::rnorm is used for compatibility with R's set.seed 
    
    // sample sigma2
    double res_ss = sum(pow(y - aux_mu,2));                                       // residual sum of squares
    aux_sigma2    = (hyper[1] + res_ss)/R::rchisq(hyper[2] + T);                  // sample sigma2 R::rchisq is used for compatibility with R's set.seed 
    
    // sample sigma2m
    vec mu0_vec(1, fill::value(aux_mu0));                                         // mu0_vec is of type vec that can be passed to join_cols() in the next line
    double mu_ss  = sum(pow(diff(join_cols(mu0_vec,aux_mu)), 2));                 // state-space equation residual sum of squares
    aux_sigma2m   = (hyper[1] + mu_ss)/R::rchisq(hyper[2] + T);                   // sample sigma2m R::rchisq is used for compatibility with R's set.seed 
    
    // simulation smoother
    mat precision = IT/aux_sigma2 + HH/aux_sigma2m;                               // computes the precision matrix of mu
    vec location  = y/aux_sigma2 + aux_mu0*IT.col(0)/aux_sigma2m;                 // computes the covector of mu
    
    // BEGIN: algorithm for simulation smoother
    // ######################################################
    vec     precision_diag    = precision.diag();   // read the main diagonal of the precision
    double  precision_offdiag = precision(1,0);     // read the offdiagonal element of the precision
    
    List precision_chol = cholesky_tridiagonal(precision_diag, precision_offdiag);    // Cholesky decomposition using a dedicated technique
    vec  aa             = forward_algorithm(precision_chol["chol_diag"],              // this forward substitution can be performed outside of the loop
                                              precision_chol["chol_offdiag"], 
                                                            location);
    vec  epsilon        = rnorm(T);                       // sample normal draws using Rcpp::rnorm for compatibility with R's set.seed()
    aux_mu              = backward_algorithm(precision_chol["chol_diag"],
                                             precision_chol["chol_offdiag"],
                                                           aa + epsilon);     // this has to be done in the loop as function backward_algorithm requires covector to be a vector (not a matrix)
    // FIN: algorithm for simulation smoother
    // ######################################################
    
    posterior_mu.col(s)     = aux_mu;         // put the draws in the MCMC output
    posterior_mu0(s)        = aux_mu0;
    posterior_sigma2(s)     = aux_sigma2;
    posterior_sigma2m(s)    = aux_sigma2m;
    
  }
  List last_draw;                             // create output lists
  last_draw["mu"]       = aux_mu;
  last_draw["mu0"]      = aux_mu0;
  last_draw["sigma2"]   = aux_sigma2;
  last_draw["sigma2m"]  = aux_sigma2m;
                                                   
  Rcpp::List posterior;
  posterior["mu"]       = posterior_mu;
  posterior["mu0"]      = posterior_mu0;
  posterior["sigma2"]   = posterior_sigma2;
  posterior["sigma2m"]  = posterior_sigma2m;
  
  List output;
  output["last.draw"]   = last_draw;
  output["posterior"]   = posterior;

  return output;
}


/*** R
# y = rnorm(5)
# set.seed(123)
# ll.sv = LL_arma_stochvol(y, 100)
# set.seed(123)
# ll.no = LL_arma_solve(y, 100)
# ll.sv$last.draw$mu0
# ll.no$last.draw$mu0
# ll.sv$last.draw$mu
# ll.no$last.draw$mu
*/
