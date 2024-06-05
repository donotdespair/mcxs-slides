############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 19: Modeling trend inflation
############################################################

library(mvtnorm)
library(mgcv)
library(Matrix)

# Gibbs sampler for a UC model using simulation smoother
############################################################
UC.Gibbs.sampler    = function(S, starting.values, priors){
  # Estimation of an UC-AR(p) model with sigma_{\eta e}=0
  # and with time-varying drift parameter
  # with hierarchical prior for sigmas
  # as specified in Lecture 19 slides 

  aux     = starting.values
  p       = length(aux$alpha)
  T       = nrow(aux$Y)
  
  posteriors    = list(
    tau     = matrix(NA,T,S),
    epsilon = matrix(NA,T,S),
    mu      = matrix(NA,T,S),
    tau0    = matrix(NA,S),
    mu0     = matrix(NA,S),
    alpha   = matrix(NA,p,S),
    sigma   = matrix(NA,3,S),
    sigma.s = matrix(NA,S),
    non.stationary.iterations = rep(NA,S)
  )
  
  e1      = diag(T)[,1]
  
  for (s in 1:S){
    
    # Sampling tau0
    ###########################
    V.tau0.bar    = 1/(1/aux$sigma[1] + 1/priors$tau0.v)
    tau0.bar      = V.tau0.bar %*% ( (aux$tau[1] - aux$mu[1])/aux$sigma[1] + priors$tau0.m/priors$tau0.v )
    tau0.draw     = rnorm(1,as.vector(tau0.bar),V.tau0.bar)
    aux$tau0      = as.vector(tau0.draw)
    
    # Sampling mu0
    ###########################
    V.mu0.bar     = 1/(1/aux$sigma[3] + 1/priors$mu0.v)
    mu0.bar       = V.tau0.bar %*% (aux$mu[1]/aux$sigma[3] + priors$mu0.m/priors$mu0.v )
    mu0.draw      = rnorm(1,as.vector(mu0.bar),V.mu0.bar)
    aux$mu0       = as.vector(mu0.draw)
    
    # Sampling alpha
    ###########################
    X.epsilon     = matrix(NA,T,0)
    for (i in 1:p){
      X.epsilon   = cbind(X.epsilon, c(rep(0,i),aux$epsilon[1:(T-i),]))
    }
    
    alpha.v.inv   = diag(1/diag(priors$alpha.v))
    V.alpha.bar   = solve((1/aux$sigma[2])*crossprod(X.epsilon) + alpha.v.inv )
    V.alpha.bar   = 0.5*(V.alpha.bar + t(V.alpha.bar))
    alpha.bar     = V.alpha.bar %*% ( (1/aux$sigma[2])*crossprod(X.epsilon, aux$epsilon) + alpha.v.inv%*%priors$alpha.m )

    # non.stationary= TRUE
    # ns.i          = 1
    # while(non.stationary){
      alpha.draw    = rmvnorm(1,as.vector(alpha.bar),V.alpha.bar)
    #   non.stationary= max(Mod(eigen(rbind(alpha.draw,cbind(diag(p-1),rep(0,p-1))))$values))>=1
    #   ns.i          = ns.i + 1
    # }
    aux$alpha     = as.vector(alpha.draw)
    
    H.alpha       = diag(T)
    for (i in 1:p){
      sdiag(H.alpha,-i) =  -aux$alpha[i]
    }
    aux$H.alpha   = H.alpha
    
    # Sampling sigmas
    ###########################
    s.a           = priors$s.a + 1.5*priors$sigma.nu
    s.s           = 1/((1/priors$s.a) + .5*sum(1/aux$sigma))
    sigma.s.draw  = rgamma(1, shape=s.a, scale=s.s)
    aux$sigma.s   = sigma.s.draw
    
    sigma.eta.s   = as.numeric(aux$sigma.s + crossprod( c(aux$tau[1],diff(aux$tau)) - e1*aux$tau0 - aux$mu))
    sigma.eta.nu  = priors$sigma.nu + T
    sigma.eta.draw= sigma.eta.s/rchisq(1,sigma.eta.nu)
    
    sigma.e.s     = as.numeric(aux$sigma.s + crossprod(aux$H.alpha%*%aux$epsilon))
    sigma.e.nu    = priors$sigma.nu + T
    sigma.e.draw  = sigma.e.s/rchisq(1,sigma.e.nu)
    
    sigma.m.s     = as.numeric(aux$sigma.s + crossprod( diff(c(aux$mu0,aux$mu)) ))
    sigma.m.nu    = priors$sigma.nu + T
    sigma.m.draw  = sigma.m.s/rchisq(1,sigma.m.nu)
    
    aux$sigma     = c(sigma.eta.draw,sigma.e.draw,sigma.m.draw)
    
    # Sampling mu
    ###########################
    V.mu.inv      = HH/aux$sigma[3] + diag(T)/aux$sigma[1]
    b.mu          = diff(c(aux$tau0,aux$tau))/aux$sigma[1] + e1*aux$mu0/aux$sigma[3]
    lead.diag     = diag(V.mu.inv)
    sub.diag      = sdiag(V.mu.inv, -1)
    D.chol        = trichol(ld = lead.diag, sd=sub.diag)
    D.L           = diag(D.chol$ld)
    sdiag(D.L,-1) = D.chol$sd
    x             = as.matrix(rnorm(T))
    a             = forwardsolve(D.L, b.mu)
    draw          = backsolve(t(D.L), a + x)
    aux$mu        = as.matrix(draw)
    
    # Sampling tau
    ###########################
    V.tau.inv     = (1/aux$sigma[2])*crossprod(aux$H.alpha) + (1/aux$sigma[1])*HH
    V.tau.inv     = 0.5*(V.tau.inv + t(V.tau.inv))
    b.tau         = (1/aux$sigma[2])*crossprod(aux$H.alpha, aux$H.alpha%*%aux$Y) + (1/aux$sigma[1])*crossprod(H, aux$mu + e1*aux$tau0)
    precision.L   = t(bandchol(V.tau.inv))
    epsilon       = rnorm(T)
    b.tau.tmp     = forwardsolve(precision.L, b.tau)
    tau.draw      = backsolve(t(precision.L), b.tau.tmp + epsilon)
    aux$tau       = tau.draw
    
    # Sampling epsilon
    ###########################
    b.epsilon     = (1/aux$sigma[1])*HH%*%(aux$Y- cumsum(aux$mu) - rep(aux$tau0,T))
    epsilon.n     = rnorm(T)
    b.epsilon.tmp = forwardsolve(precision.L, b.epsilon)
    epsilon.draw  = backsolve(t(precision.L), b.epsilon.tmp + epsilon.n)
    aux$epsilon   = epsilon.draw
    
    # Storage
    ###########################
    posteriors$tau[,s]     = aux$tau
    posteriors$epsilon[,s] = aux$epsilon
    posteriors$mu[,s]      = aux$mu
    posteriors$tau0[s]     = aux$tau0
    posteriors$mu0[s]      = aux$mu0
    posteriors$alpha[,s]   = aux$alpha
    posteriors$sigma[,s]   = aux$sigma
    posteriors$sigma.s[s]   = aux$sigma.s
    # posteriors$non.stationary.iterations[s] = ns.i
    
    if (s%%1000==0){cat(" ",s)}
  }
  
  output      = list(
    posterior = posteriors,
    last.draw = aux
  )
  return(output)
}

