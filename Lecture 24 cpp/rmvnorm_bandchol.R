library(mgcv)

rmvnorm_bandchol = function(n, precision, location){
  T           = dim(precision)[1]
  precision.L = t(bandchol(precision))
  
  epsilon     = matrix(rnorm(n*T), ncol=n)
  a           = forwardsolve(precision.L, location)
  draw        = backsolve(t(precision.L), 
                          matrix(rep(a,n), ncol=n) + epsilon)
  
  return(draw)
}
