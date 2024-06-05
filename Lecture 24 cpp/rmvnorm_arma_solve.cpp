#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat rmvnorm_arma_solve(int n, mat precision, vec location){
  int T = precision.n_rows;
  
  mat epsilon(T, n);
  for (int i=0; i<n; i++){
    epsilon.col(i) = as<vec>(rnorm(T));
  }
  
  mat location_matrix(T, n);
  location_matrix.each_col() = location;
  mat precision_chol  = trimatu(arma::chol(precision));
  mat draw            = solve(precision_chol, 
                          solve(trans(precision_chol), 
                          location_matrix) + epsilon);
  return draw;
}



/*** R
# N   = 500
# pr = 5*diag(N)
# lo = rep(0,N)

# set.seed(123)
# rmvnorm_arma_solve(2, pr, lo)
# set.seed(123)
# rmvnorm.nothing.special(2, pr, lo)

# library(microbenchmark)
# microbenchmark(
#   rcpp <- rmvnorm_arma_set_seed(100, pr, lo),
#   rr   <- rmvnorm.nothing.special(100, pr, lo),
#   check = "equal",
#   setup = set.seed(123)
# )
# plot(density(rcpp))
# lines(density(rr),col="red")
*/
