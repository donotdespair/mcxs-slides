library(mgcv)

rmvnorm_trichol = function(n, precision, location){
  T           = dim(precision)[1]
  lead.diag   = diag(precision)
  sub.diag    = sdiag(precision, -1)
  
  precision.chol    = trichol(ld = lead.diag, sd=sub.diag)
  precision.L       = diag(precision.chol$ld)
  sdiag(precision.L,-1) = precision.chol$sd
  
  epsilon     = matrix(rnorm(n*T), ncol=n)
  a           = forwardsolve(precision.L, location)
  draw        = backsolve(t(precision.L), 
                          matrix(rep(a,n), ncol=n) 
                          + epsilon)
  return(draw)
}