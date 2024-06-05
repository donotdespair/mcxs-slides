# requires library(mgcv)

# Gibbs sampler for the local-level model
LL_band_precision = function(y, 
                             S = 100, 
                             starting.values=list(mu=rep(0, length(y)), mu0=0, sigma2=1, sigma2m=1), 
                             hyper = c(10, 1, 3)) {

  T         = length(y)
  aux       = starting.values
  
  posterior = list(
    mu      = matrix(NA, T, S),
    mu0     = rep(NA, S),
    sigma2  = rep(NA, S),
    sigma2m = rep(NA, S)
  )
  
  H = IT    = diag(T)
  sdiag(H, -1) = -1
  HH        = crossprod(H)
  
  for (s in 1:S){
    # sample mu0
    vm          = 1/((1/hyper[1]) + (1/aux$sigma2m))
    aux$mu0     = rnorm(1, (vm * aux$mu[1])/aux$sigma2m, sqrt(vm) )
    
    # sample sigma2
    aux$sigma2  = (hyper[2] + sum((y-aux$mu)^2))/rchisq(1, hyper[3]+T)
    
    # sample sigma2m
    aux$sigma2m = (hyper[2] + sum(diff(c(aux$mu0, aux$mu))^2))/rchisq(1, hyper[3]+T)
    
    # simulation smoother
    precision   = IT/aux$sigma2 + HH/aux$sigma2m
    location    = y/aux$sigma2 + aux$mu0*IT[,1]/aux$sigma2m
      
    # BEGIN: algorithm for simulation smoother
    ######################################################
    precision.L = t(bandchol(precision))
    a           = forwardsolve(precision.L, location)
    draw        = backsolve(t(precision.L), (a + rnorm(T)))
    # FIN: algorithm for simulation smoother
    ######################################################
    aux$mu      = as.vector(draw)
    
    posterior$mu[,s]    = aux$mu
    posterior$mu0[s]    = aux$mu0
    posterior$sigma2[s] = aux$sigma2
    posterior$sigma2m[s] = aux$sigma2m
  }
  
  out    = list(
    last.draw   = aux,
    posterior   = posterior
  )
  return(out)
}