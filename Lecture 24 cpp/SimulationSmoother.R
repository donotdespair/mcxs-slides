rm(list=ls())

library(mgcv)
library(Rcpp)
library(microbenchmark)

source("rmvnorm_trichol.R")
source("rmvnorm_bandchol.R")
source("rmvnorm_solve.R")
sourceCpp("rmvnorm_arma_inv.cpp")
sourceCpp("rmvnorm_arma_solve.cpp")
sourceCpp("rmvnorm_arma_stochvol.cpp")

set.seed(12345)
n     = 100
T     = 250
s     = rgamma(1, shape=10, scale=10)
precision = rgamma(1, shape=10, scale=10)*diag(T) + 2*s*diag(T)
sdiag(precision, 1) = -s
sdiag(precision, -1) = -s
location = as.matrix(rnorm(T))

microbenchmark(
  R.solve   = rmvnorm_solve(n, precision, location),
  R.band    = rmvnorm_bandchol(n, precision, location),
  R.tridiag = rmvnorm_trichol(n, precision, location),
  cpp.inv = rmvnorm_arma_inv(n, precision, location),
  cpp.sol = rmvnorm_arma_solve(n, precision, location),
  cpp.sto = rmvnorm_arma_stochvol(n, precision, location),
  check  = "equal",
  setup  = set.seed(123)
)

# Create a positive definite and symmetric matrix - the precision matrix in this example
set.seed(12345)
T     = 720
Omega = rWishart(1, T+10, diag(T))[,,1]
alpha     = as.matrix(rnorm(T))

microbenchmark(usual   = rmvnorm_nothing_special(n=100, precision=Omega, location=alpha),
               cpp.inv = rmvnorm_arma_inv(n=100, precision=Omega, location=alpha),
               cpp.sol = rmvnorm_arma_solve(n=100, precision=Omega, location=alpha),
               check  = "equal",
               setup  = set.seed(123)
)


# library(SparseM)
# 
# set.seed(12345)
# T     = 50
# s     = rgamma(1, shape=10, scale=10)
# 
# Omega = 2*s*diag(T)
# sdiag(Omega, 1) = -s
# sdiag(Omega, -1) = -s
# alpha     = as.matrix(rep(2,T))
# epsilon     = rnorm(T)
# 
# omega = as.matrix.csr(Omega)
# all.equal(
#   as.vector(solve(Omega,alpha)), 
#   solve(omega,alpha)
# )
# 
# omega.chol  = chol(omega)
# Omega.chol  = t(bandchol(Omega))
# 
# all.equal(
#   backsolve(omega.chol, epsilon, twice = TRUE),
#   backsolve(t(Omega.chol), epsilon)
# )
