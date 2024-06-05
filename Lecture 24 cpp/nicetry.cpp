#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec nicetry (mat y, mat x) {
  vec beta_hat = solve(x.t()*x, x.t()*y);
  return beta_hat;
}