rmvnorm_solve = function(n, precision, location){
  T           = dim(precision)[1]
  precision.chol.inv    = solve(t(chol(precision)))
  
  epsilon     = matrix(rnorm(n*T), ncol=n)
  draw        = t(precision.chol.inv) %*% 
                 (matrix(rep(precision.chol.inv%*%
                    location,n), ncol=n) + epsilon)
  return(draw)
}