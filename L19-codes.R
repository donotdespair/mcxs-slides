############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 19: Modeling trend inflation
############################################################

library(mvtnorm)
library(mgcv)
library(Matrix)
library(GIGrvg)

# Gibbs sampler for a simple UC model using simulation smoother
############################################################

UC.AR.Gibbs.sampler    = function(S, starting.values, priors){
  # Estimation of a simple UC-AR(p) model with sigma_{\eta e}=0
  # as specified in Lecture 17 slides 
  
  aux     = starting.values
  p       = length(aux$alpha)
  T       = nrow(aux$Y)
  
  posteriors    = list(
    tau     = matrix(NA,T,S),
    epsilon = matrix(NA,T,S),
    beta    = matrix(NA,2,S),
    alpha   = matrix(NA,p,S),
    sigma   = matrix(NA,2,S),
    non.stationary.iterations = rep(NA,S)
  )
  
  for (s in 1:S){
    
    # Sampling beta
    ###########################
    beta.v.inv    = diag(1/diag(priors$beta.v))
    V.beta.bar    = solve((1/aux$sigma[1])*crossprod(priors$X.tau) + beta.v.inv )
    beta.bar      = V.beta.bar %*% ( (1/aux$sigma[1])*crossprod(priors$X.tau, c(aux$tau[1],diff(aux$tau))) + beta.v.inv%*%priors$beta.m )
    beta.draw     = rmvnorm(1,as.vector(beta.bar),V.beta.bar)
    aux$beta      = as.vector(beta.draw)
    
    X.epsilon     = matrix(NA,T,0)
    for (i in 1:p){
      X.epsilon   = cbind(X.epsilon, c(rep(0,i),aux$epsilon[1:(T-i),]))
    }
    
    # Sampling alpha
    ###########################
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
    
    # Sampling sigma
    ###########################
    sigma.eta.s   = as.numeric(priors$sigma.s + crossprod( (c(aux$tau[1],diff(aux$tau)) - X.tau%*%aux$beta), (priors$H%*%aux$tau - X.tau%*%aux$beta)))
    sigma.eta.nu  = priors$sigma.nu + T
    sigma.eta.draw= sigma.eta.s/rchisq(1,sigma.eta.nu)
    
    sigma.e.s     = as.numeric(priors$sigma.s + crossprod(aux$H.alpha%*%aux$epsilon))
    sigma.e.nu    = priors$sigma.nu + T
    sigma.e.draw  = sigma.e.s/rchisq(1,sigma.e.nu)
    aux$sigma     = c(sigma.eta.draw,sigma.e.draw)
    
    # Sampling tau
    ###########################
    V.tau.inv     = (1/aux$sigma[2])*crossprod(aux$H.alpha) + (1/aux$sigma[1])*HH
    V.tau.inv     = 0.5*(V.tau.inv + t(V.tau.inv))
    b.tau         = (1/aux$sigma[2])*crossprod(aux$H.alpha, aux$H.alpha%*%aux$Y) + (1/aux$sigma[1])*crossprod(H, X.tau%*%aux$beta)
    precision.L   = t(bandchol(V.tau.inv))
    epsilon       = rnorm(T)
    b.tau.tmp     = forwardsolve(precision.L, b.tau)
    tau.draw      = backsolve(t(precision.L), b.tau.tmp + epsilon)
    aux$tau       = tau.draw
    
    # Sampling epsilon
    ###########################
    b.epsilon     = (1/aux$sigma[1])*HH%*%(aux$Y- cumsum(X.tau%*%aux$beta))
    epsilon.n     = rnorm(T)
    b.epsilon.tmp = forwardsolve(precision.L, b.epsilon)
    epsilon.draw  = backsolve(t(precision.L), b.epsilon.tmp + epsilon.n)
    aux$epsilon   = epsilon.draw
    
    posteriors$tau[,s]     = aux$tau
    posteriors$epsilon[,s] = aux$epsilon
    posteriors$beta[,s]    = aux$beta
    posteriors$alpha[,s]   = aux$alpha
    posteriors$sigma[,s]   = aux$sigma
    # posteriors$non.stationary.iterations[s] = ns.i
    if (s%%1000==0){cat(" ",s)}
  }
  
  output      = list(
    posterior = posteriors,
    last.draw = aux
  )
  return(output)
}





UC.AR.Gibbs.sampler.sigma    = function(S, starting.values, priors){
  # Estimation of a simple UC-AR(p) model with sigma_{\eta e}=0
  # with hierarchical prior for sigmas
  # this specification allows for the estimation of a deterministic broken trend
  # in Y implemented by providing appropriately specified matrix priors$X.tau
  # as specified in Lecture 19 slides 
  
  aux     = starting.values
  p       = length(aux$alpha)
  T       = nrow(aux$Y)
  
  posteriors    = list(
    tau     = matrix(NA,T,S),
    epsilon = matrix(NA,T,S),
    beta    = matrix(NA,ncol(priors$X.tau),S),
    alpha   = matrix(NA,p,S),
    sigma   = matrix(NA,2,S),
    sigma.s = matrix(NA,2,S),
    non.stationary.iterations = rep(NA,S)
  )
  
  for (s in 1:S){
    
    # Sampling beta
    ###########################
    beta.v.inv    = diag(1/diag(priors$beta.v))
    V.beta.bar    = solve((1/aux$sigma[1])*crossprod(priors$X.tau) + beta.v.inv )
    beta.bar      = V.beta.bar %*% ( (1/aux$sigma[1])*crossprod(priors$X.tau, c(aux$tau[1],diff(aux$tau))) + beta.v.inv%*%priors$beta.m )
    beta.draw     = rmvnorm(1,as.vector(beta.bar),V.beta.bar)
    aux$beta      = as.vector(beta.draw)
    
    X.epsilon     = matrix(NA,T,0)
    for (i in 1:p){
      X.epsilon   = cbind(X.epsilon, c(rep(0,i),aux$epsilon[1:(T-i),]))
    }
    
    # Sampling alpha
    ###########################
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
    s.a           = priors$s.a + priors$sigma.nu
    s.s           = 1/((1/priors$s.s) + .5*sum(1/aux$sigma))
    sigma.s.draw  = rgamma(1, shape=s.a, scale=s.s)
    aux$sigma.s   = sigma.s.draw
    
    sigma.eta.s   = as.numeric(aux$sigma.s + crossprod( (c(aux$tau[1],diff(aux$tau)) - priors$X.tau%*%aux$beta), (priors$H%*%aux$tau - priors$X.tau%*%aux$beta)))
    sigma.eta.nu  = priors$sigma.nu + T
    sigma.eta.draw= sigma.eta.s/rchisq(1,sigma.eta.nu)
    
    sigma.e.s     = as.numeric(aux$sigma.s + crossprod(aux$H.alpha%*%aux$epsilon))
    sigma.e.nu    = priors$sigma.nu + T
    sigma.e.draw  = sigma.e.s/rchisq(1,sigma.e.nu)
    aux$sigma     = c(sigma.eta.draw,sigma.e.draw)
    
    # Sampling tau
    ###########################
    V.tau.inv     = (1/aux$sigma[2])*crossprod(aux$H.alpha) + (1/aux$sigma[1])*HH
    V.tau.inv     = 0.5*(V.tau.inv + t(V.tau.inv))
    b.tau         = (1/aux$sigma[2])*crossprod(aux$H.alpha, aux$H.alpha%*%aux$Y) + (1/aux$sigma[1])*crossprod(H, priors$X.tau%*%aux$beta)
    precision.L   = t(bandchol(V.tau.inv))
    epsilon       = rnorm(T)
    b.tau.tmp     = forwardsolve(precision.L, b.tau)
    tau.draw      = backsolve(t(precision.L), b.tau.tmp + epsilon)
    aux$tau       = tau.draw
    
    # Sampling epsilon
    ###########################
    b.epsilon     = (1/aux$sigma[1])*HH%*%(aux$Y- cumsum(priors$X.tau%*%aux$beta))
    epsilon.n     = rnorm(T)
    b.epsilon.tmp = forwardsolve(precision.L, b.epsilon)
    epsilon.draw  = backsolve(t(precision.L), b.epsilon.tmp + epsilon.n)
    aux$epsilon   = epsilon.draw
    
    posteriors$tau[,s]     = aux$tau
    posteriors$epsilon[,s] = aux$epsilon
    posteriors$beta[,s]    = aux$beta
    posteriors$alpha[,s]   = aux$alpha
    posteriors$sigma[,s]   = aux$sigma
    posteriors$sigma.s[,s]   = aux$sigma.s
    # posteriors$non.stationary.iterations[s] = ns.i
    if (s%%1000==0){cat(" ",s)}
  }
  
  output      = list(
    posterior = posteriors,
    last.draw = aux
  )
  return(output)
}





UC.AR.Gibbs.sampler.sigma.sigmaetae    = function(S, starting.values, priors){
  # Estimation of a simple UC-AR(p) model
  # with estimated sigma_{\eta e}=0 (although in an alternative parameterization)
  # with hierarchical prior for sigmas
  # this specification allows for the estimation of a deterministic broken trend
  # in Y implemented by providing appropriately specified matrix priors$X.tau
  # as specified in Lecture 19 slides 
  
  aux     = starting.values
  p       = length(aux$alpha)
  T       = nrow(aux$Y)
  
  posteriors    = list(
    tau     = matrix(NA,T,S),
    epsilon = matrix(NA,T,S),
    beta    = matrix(NA,ncol(priors$X.tau)+1,S),
    alpha   = matrix(NA,p,S),
    sigma   = matrix(NA,2,S),
    sigma.s = matrix(NA,1,S),
    non.stationary.iterations = rep(NA,S)
  )
  
  for (s in 1:S){
    
    # Sampling beta
    ###########################
    beta.v.inv    = diag(1/diag(priors$beta.v))
    X.tau         = cbind(priors$X.tau, aux$epsilon)
    V.beta.bar    = solve((1/aux$sigma[1])*crossprod(X.tau) + beta.v.inv )
    beta.bar      = V.beta.bar %*% ( (1/aux$sigma[1])*crossprod(X.tau, c(aux$tau[1],diff(aux$tau))) + beta.v.inv%*%priors$beta.m )
    beta.draw     = rmvnorm(1,as.vector(beta.bar),V.beta.bar)
    aux$beta      = as.vector(beta.draw)
    
    X.epsilon     = matrix(NA,T,0)
    for (i in 1:p){
      X.epsilon   = cbind(X.epsilon, c(rep(0,i),aux$epsilon[1:(T-i),]))
    }
    
    # Sampling alpha
    ###########################
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
    s.a           = priors$s.a + priors$sigma.nu
    s.s           = 1/((1/priors$s.s) + .5*sum(1/aux$sigma))
    sigma.s.draw  = rgamma(1, shape=s.a, scale=s.s)
    aux$sigma.s   = sigma.s.draw
    
    sigma.eta.s   = as.numeric(aux$sigma.s + crossprod( (c(aux$tau[1],diff(aux$tau)) - X.tau%*%aux$beta), (priors$H%*%aux$tau - X.tau%*%aux$beta)))
    sigma.eta.nu  = priors$sigma.nu + T
    sigma.eta.draw= sigma.eta.s/rchisq(1,sigma.eta.nu)
    
    sigma.e.s     = as.numeric(aux$sigma.s + crossprod(aux$H.alpha%*%aux$epsilon))
    sigma.e.nu    = priors$sigma.nu + T
    sigma.e.draw  = sigma.e.s/rchisq(1,sigma.e.nu)
    aux$sigma     = c(sigma.eta.draw,sigma.e.draw)
    
    # Sampling tau
    ###########################
    V.tau.inv     = (1/aux$sigma[2])*crossprod(aux$H.alpha) + (1/aux$sigma[1])*HH
    V.tau.inv     = 0.5*(V.tau.inv + t(V.tau.inv))
    b.tau         = (1/aux$sigma[2])*crossprod(aux$H.alpha, aux$H.alpha%*%aux$Y) + (1/aux$sigma[1])*crossprod(H, X.tau%*%aux$beta)
    precision.L   = t(bandchol(V.tau.inv))
    epsilon       = rnorm(T)
    b.tau.tmp     = forwardsolve(precision.L, b.tau)
    tau.draw      = backsolve(t(precision.L), b.tau.tmp + epsilon)
    aux$tau       = tau.draw
    
    # Sampling epsilon
    ###########################
    b.epsilon     = (1/aux$sigma[1])*HH%*%(aux$Y- cumsum(X.tau%*%aux$beta))
    epsilon.n     = rnorm(T)
    b.epsilon.tmp = forwardsolve(precision.L, b.epsilon)
    epsilon.draw  = backsolve(t(precision.L), b.epsilon.tmp + epsilon.n)
    aux$epsilon   = epsilon.draw
    
    posteriors$tau[,s]     = aux$tau
    posteriors$epsilon[,s] = aux$epsilon
    posteriors$beta[,s]    = aux$beta
    posteriors$alpha[,s]   = aux$alpha
    posteriors$sigma[,s]   = aux$sigma
    posteriors$sigma.s[,s]   = aux$sigma.s
    # posteriors$non.stationary.iterations[s] = ns.i
    if (s%%1000==0){cat(" ",s)}
  }
  
  output      = list(
    posterior = posteriors,
    last.draw = aux
  )
  return(output)
}





UC.AR.Gibbs.sampler.gamma    = function(S, starting.values, priors){
  # Estimation of a simple UC-AR(p) model with sigma_{\eta e}=0
  # with hierarchical prior for sigmas
  # and with gamma prior for sigma2_eta

  aux     = starting.values
  p       = length(aux$alpha)
  T       = nrow(aux$Y)
  
  posteriors    = list(
    tau     = matrix(NA,T,S),
    epsilon = matrix(NA,T,S),
    beta    = matrix(NA,ncol(priors$X.tau),S),
    alpha   = matrix(NA,p,S),
    sigma   = matrix(NA,2,S),
    sigma.s = matrix(NA,2,S),
    non.stationary.iterations = rep(NA,S)
  )
  
  for (s in 1:S){
    
    # Sampling beta
    ###########################
    beta.v.inv    = diag(1/diag(priors$beta.v))
    V.beta.bar    = solve((1/aux$sigma[1])*crossprod(priors$X.tau) + beta.v.inv )
    beta.bar      = V.beta.bar %*% ( (1/aux$sigma[1])*crossprod(priors$X.tau, c(aux$tau[1],diff(aux$tau))) + beta.v.inv%*%priors$beta.m )
    beta.draw     = rmvnorm(1,as.vector(beta.bar),V.beta.bar)
    aux$beta      = as.vector(beta.draw)
    
    X.epsilon     = matrix(NA,T,0)
    for (i in 1:p){
      X.epsilon   = cbind(X.epsilon, c(rep(0,i),aux$epsilon[1:(T-i),]))
    }
    
    # Sampling alpha
    ###########################
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
    sigma.s.draw  = rgig(1, priors$s.a - (priors$sigma.nu+3)/2, 2*aux$sigma[1], (1/aux$sigma[2]) + (2/priors$s.s))
    aux$sigma.s   = sigma.s.draw
    
    sigma.eta.chi = as.numeric(crossprod( priors$X.tau %*% aux$beta - H %*% aux$tau))
    sigma.eta.draw= rgig(1, -(T-1)/2, , 2/sigma.s.draw)
    
    sigma.e.s     = as.numeric(aux$sigma.s + crossprod(aux$H.alpha%*%aux$epsilon))
    sigma.e.nu    = priors$sigma.nu + T
    sigma.e.draw  = sigma.e.s/rchisq(1,sigma.e.nu)
    aux$sigma     = c(sigma.eta.draw,sigma.e.draw)
    
    # Sampling tau
    ###########################
    V.tau.inv     = (1/aux$sigma[2])*crossprod(aux$H.alpha) + (1/aux$sigma[1])*HH
    V.tau.inv     = 0.5*(V.tau.inv + t(V.tau.inv))
    b.tau         = (1/aux$sigma[2])*crossprod(aux$H.alpha, aux$H.alpha%*%aux$Y) + (1/aux$sigma[1])*crossprod(H, priors$X.tau%*%aux$beta)
    precision.L   = t(bandchol(V.tau.inv))
    epsilon       = rnorm(T)
    b.tau.tmp     = forwardsolve(precision.L, b.tau)
    tau.draw      = backsolve(t(precision.L), b.tau.tmp + epsilon)
    aux$tau       = tau.draw
    
    # Sampling epsilon
    ###########################
    b.epsilon     = (1/aux$sigma[1])*HH%*%(aux$Y- cumsum(priors$X.tau%*%aux$beta))
    epsilon.n     = rnorm(T)
    b.epsilon.tmp = forwardsolve(precision.L, b.epsilon)
    epsilon.draw  = backsolve(t(precision.L), b.epsilon.tmp + epsilon.n)
    aux$epsilon   = epsilon.draw
    
    posteriors$tau[,s]     = aux$tau
    posteriors$epsilon[,s] = aux$epsilon
    posteriors$beta[,s]    = aux$beta
    posteriors$alpha[,s]   = aux$alpha
    posteriors$sigma[,s]   = aux$sigma
    posteriors$sigma.s[,s]   = aux$sigma.s
    # posteriors$non.stationary.iterations[s] = ns.i
    if (s%%1000==0){cat(" ",s)}
  }
  
  output      = list(
    posterior = posteriors,
    last.draw = aux
  )
  return(output)
}


