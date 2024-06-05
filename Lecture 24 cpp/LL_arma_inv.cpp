#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List LL_arma_inv(vec y, 
                 int S = 10,
                 Nullable<List> starting_values = R_NilValue,
                 NumericVector  Hyper = NumericVector::create(10,1,3)) {
  
  vec hyper = as<vec>(Hyper);                 // read argument with a default value from the input
  int T = y.n_rows;                           // read sample size
  
  RNGScope scp;                               // facilitates reading R's set.seed()
  
  vec    aux_mu(T, fill::zeros);                        // create auxiliary objects with the current MCMC state
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
    mat precision_chol_inv  = trans(inv(trimatu(chol(precision))));      // t(solve(chol())) applied to an upper-triangular Cholesky decomposition
    vec epsilon             = rnorm(T);                                               // Rcpp:rnorm is used not Armadillo's randn() for compatibility with R's set.seed(); also somehow rnorm(T) cannot simply be put in the place of epsilon in the next line
    aux_mu                  = trans(precision_chol_inv) * (precision_chol_inv * location + epsilon);    // computes the draws from the target mv normal distribution
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
# sv = list(mu=rep(10,5), mu0=0, sigma2=1, sigma2=1)
# ll.sv = LL_arma_inv(y, 2, starting_values=sv)
# ll.no = LL_arma_inv(y, 2)
# ll.sv$posterior$mu0
# ll.no$posterior$mu0
*/
