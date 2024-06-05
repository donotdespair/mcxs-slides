############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 10: Forecasting with Large Bayesian VARs
############################################################

rm(list=ls())
library(mvtnorm)
library(microbenchmark)

set.seed(123456)


# Kronecker product inverse
############################################################
N       = 10
p       = 12

Sigma     = rWishart(1,N+2,diag(N))[,,1]
XX        = rWishart(1,p*N+3,diag(1+p*N))[,,1]

microbenchmark(
  regular   = solve(kronecker(Sigma,XX)),
  kronecker = kronecker(solve(Sigma),solve(XX))
)


# Diagonal matrix inverse
############################################################
K       = 1 + p*N

V.inv   = diag(rgamma(K,1,1))
microbenchmark(
  regular   = solve(V.inv),
  diagonal  = diag(1/diag(V.inv))
)


# Triangular matrix inverse
############################################################
A.bar.tmp     = as.matrix(rnorm(K))
V.bar.inv     = XX + diag(1/diag(V.inv))

dedicated     = function(A.bar.tmp,V.bar.inv){
  C = chol(V.bar.inv); 
  return(backsolve(C, forwardsolve(t(C), A.bar.tmp)))
}

microbenchmark(
  regular   = solve(V.bar.inv) %*% A.bar.tmp,
  dedicated = dedicated(A.bar.tmp,V.bar.inv)
)

