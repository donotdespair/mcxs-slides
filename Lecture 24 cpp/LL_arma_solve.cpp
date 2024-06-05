#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List LL_arma_solve(
    vec y, 
    int S = 10,
    Nullable<List> starting_values = R_NilValue,
    NumericVector  Hyper = NumericVector::create(10,1,3)) {
  
  vec hyper = as<vec>(Hyper);
  int T = y.n_rows;
  
  vec    aux_mu(T, fill::zeros);
  double aux_mu0      = 0;
  double aux_sigma2   = 1;
  double aux_sigma2m  = 1;
  
  if (starting_values.isNotNull()){
    List Starting_values(starting_values);
    aux_mu            = as<vec>(Starting_values["mu"]);
    aux_mu0           = Starting_values["mu0"];
    aux_sigma2        = Starting_values["sigma2"];
    aux_sigma2m       = Starting_values["sigma2"];
  } 
  
  mat posterior_mu(T, S);
  vec posterior_mu0(S);
  vec posterior_sigma2(S);
  vec posterior_sigma2m(S);
  
  mat H(T, T, fill::eye);
  H.diag(-1) += -1;
  mat HH = H.t() * H;
  mat IT(T, T, fill::eye);
  
  for (int s=0; s<S; s++){
    double vm     = 1/((1/hyper[0]) + (1/aux_sigma2m));    
    aux_mu0       = R::rnorm((vm*aux_mu[0])/aux_sigma2m, sqrt(vm));
    
    double res_ss = sum(pow(y - aux_mu,2));                 
    aux_sigma2    = (hyper[1] + res_ss)/R::rchisq(hyper[2] + T);
    
    vec mu0_vec(1, fill::value(aux_mu0));                      
    double mu_ss  = sum(pow(diff(join_cols(mu0_vec,aux_mu)), 2)); 
    aux_sigma2m   = (hyper[1] + mu_ss)/R::rchisq(hyper[2] + T);  
    
    mat precision = IT/aux_sigma2 + HH/aux_sigma2m;             
    vec location  = y/aux_sigma2 + aux_mu0*IT.col(0)/aux_sigma2m; 
    mat precision_chol  = trimatu(chol(precision));  
    vec epsilon         = rnorm(T);
    aux_mu              = solve(precision_chol, 
                                solve(trans(precision_chol),
                                  location) + epsilon);
    posterior_mu.col(s)     = aux_mu;
    posterior_mu0(s)        = aux_mu0;
    posterior_sigma2(s)     = aux_sigma2;
    posterior_sigma2m(s)    = aux_sigma2m;
  }
  
  List last_draw; 
  last_draw["mu"]       = aux_mu;
  last_draw["mu0"]      = aux_mu0;
  last_draw["sigma2"]   = aux_sigma2;
  last_draw["sigma2m"]  = aux_sigma2m;
                                                   
  List posterior;
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
# ll.sv = LL_arma_solve(y, 2, starting_values=sv)
# ll.no = LL_arma_solve(y, 2)
# ll.sv$posterior$mu0
# ll.no$posterior$mu0
*/
