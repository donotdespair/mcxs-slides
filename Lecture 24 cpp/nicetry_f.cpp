#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec nicetry (mat y, mat x) {
  vec beta_hat = solve(x.t()*x, x.t()*y);
  return beta_hat;
}

/*** R
x = cbind(rep(1,5),1:5)
y = x %*% c(1,2) + rnorm(5)
nicetry(y, x)
solve(crossprod(x), crossprod(x,y))
*/
