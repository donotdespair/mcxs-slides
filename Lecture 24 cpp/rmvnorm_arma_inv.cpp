#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat rmvnorm_arma_inv (int n, mat precision, vec location){
  int T = precision.n_rows;
  
  mat epsilon(T, n);
  for (int i=0; i<n; i++){
    epsilon.col(i) = as<vec>(rnorm(T));
  }
  
  mat location_matrix(T, n, fill::zeros);
  location_matrix.each_col() += location;
  mat precision_chol_inv = trans(inv(trimatu(chol(precision))));
  mat draw    = trans(precision_chol_inv) * 
                (precision_chol_inv * location_matrix 
                 + epsilon);
  return draw;
}

/*** R
set.seed(123)
rmvnorm_arma_inv(10, diag(3), rep(0,3))
*/
