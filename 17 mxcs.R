
rm(list=ls())

# Two functions to sample from the multivariate normal
rmvnorm.tridiag.precision = function(n, D, b){
  N           = dim(D)[1]
  lead.diag   = diag(D)
  sub.diag    = sdiag(D, -1)
  
  D.chol    = trichol(ld = lead.diag, sd=sub.diag)
  D.L       = diag(D.chol$ld)
  sdiag(D.L,-1) = D.chol$sd
  
  x           = matrix(rnorm(n*N), ncol=n)
  a           = forwardsolve(D.L, b)
  draw        = backsolve(t(D.L), matrix(rep(a,n), ncol=n) + x)
  
  return(draw)
}


rmvnorm.band.precision = function(n, D, b){
  N           = dim(D)[1]
  D.L         = t(bandchol(D))
  
  x           = matrix(rnorm(n*N), ncol=n)
  a           = forwardsolve(D.L, b)
  draw        = backsolve(t(D.L), matrix(rep(a,n), ncol=n) + x)
  
  return(draw)
}


rmvnorm.nothing.special = function(n, D, b){
  N           = dim(D)[1]
  D.chol    = t(chol(D))
  variance.chol = solve(D.chol)
  
  x           = matrix(rnorm(n*N), ncol=n)
  draw        = t(variance.chol) %*% (matrix(rep(variance.chol%*%b,n), ncol=n) + x)
  
  return(draw)
}


# Create a tri-diagonal symmetric matrix - the precision matrix in this example
library(mgcv)
library(microbenchmark)

set.seed(12345)
T     = 240
md    = rgamma(T, shape=10, scale=10)
od    = rgamma(T-1, shape=10, scale=1)
D     = 2*diag(md)
sdiag(D, 1) = -od
sdiag(D, -1) = -od
b     = as.matrix(rnorm(T))

microbenchmark(
  trid  = rmvnorm.tridiag.precision(n=100, D=D, b=b),
  usual = rmvnorm.nothing.special(n=100, D=D, b=b),
  check = "equal", setup=set.seed(123456)
)



set.seed(12345)
T     = 720
md    = rgamma(T, shape=10, scale=10)
od    = rgamma(T-1, shape=10, scale=1)
D     = 2*diag(md)
sdiag(D, 1) = -od
sdiag(D, -1) = -od
b     = as.matrix(rnorm(T))

microbenchmark(
  trid  = rmvnorm.tridiag.precision(n=100, D=D, b=b),
  usual = rmvnorm.nothing.special(n=100, D=D, b=b),
  check = "equal", setup=set.seed(123456)
)
