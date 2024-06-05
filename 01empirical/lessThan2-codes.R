
library(mvtnorm)
library(tmvtnorm)

lessThan2.Gibbs = function(S, data, starting.values ){
  
  aux         = starting.values
  
  F           = data$F
  XF          = data$XF
  G           = data$G
  XG          = data$XG
  tau         = data$tau
  Xtau        = data$Xtau
  tech        = data$tech
  tech.trend  = data$tech.trend
  
  T           = nrow(F)
  N           = length(aux$sigma2)
  
  posterior = list(
    gamma     = matrix(NA,2,S),
    phi       = matrix(NA,N-1,S),
    betatau   = matrix(NA,N+3,S),
    sigmaf2   = rep(NA,S),
    sigmag2   = matrix(NA,N-1,S),
    sigma2    = matrix(NA,N,S),
    phi.mu    = rep(NA,S),
    phi.sig   = rep(NA,S),
    delta.mu  = rep(NA,S),
    delta.sig = rep(NA,S),
    s         = rep(NA,S),
    s.sig     = rep(NA,S)
  )
  
  epsilong    = matrix(0,T,N)
  
  for (s in 1:S){
    # sample gamma
    V.gamma     = aux$sigmaf2*solve(t(XF)%*%XF)
    gamma.bar   = solve(t(XF)%*%XF)%*%t(XF)%*%F
    aux$gamma   = as.vector(rtmvnorm(1, mean=as.vector(gamma.bar), sigma=V.gamma, lower=c(0,-0.1), upper=c(1,0.1)))
    epsilong[,1]= F - XF%*%aux$gamma
    
    # sample phi
    Sigma.G.inv = kronecker(diag(1/aux$sigmag2), diag(T))
    V.phi       = solve(t(XG)%*%Sigma.G.inv%*%XG + (1/aux$phi.sig)*diag(N-1)) 
    phi.bar     = V.phi%*%( t(XG)%*%Sigma.G.inv%*%G + (aux$phi.mu/aux$phi.sig)*matrix(1,N-1,1))
    aux$phi     = as.vector(rtmvnorm(1, mean=as.vector(phi.bar), sigma=V.phi, lower=rep(0,N-1), upper=rep(1,N-1)))
    epsilong[,2:N] = matrix(G - XG%*%aux$phi,T,N-1)
    
    # sample betatau
    V.beta.prior.inv = diag(c(100,0,rep(1/aux$delta.sig,N),0))
    mu.beta.prior    = as.matrix(c(0.1,0,rep(aux$delta.mu,N),0))
    Tc          = rep(NA,N)
    Sigma.diag  = rep(NA,0)
    epsilongc   = rep(NA,0)
    for (c in 1:N){
      Tc[c]     = length(tech.trend[[c]])
      Sigma.diag= c(Sigma.diag,rep(aux$sigma2[c],Tc[c]))
      if (c==1){
        epsilongc = c(epsilongc,sqrt(aux$sigma2[c]/aux$sigmaf2)*epsilong[tech.trend[[c]]+28,c])
      } else {
        epsilongc = c(epsilongc,sqrt(aux$sigma2[c]/aux$sigmag2[c-1])*epsilong[tech.trend[[c]]+28,c])
      }
    }
    X           = cbind(Xtau,epsilongc)
    V.beta      = solve( t(X)%*%diag(1/Sigma.diag)%*%X + V.beta.prior.inv )
    beta.bar    = V.beta %*% (t(X)%*%diag(1/Sigma.diag)%*%tau + V.beta.prior.inv%*%mu.beta.prior)
    aux$betatau = as.vector(rtmvnorm(1, mean=as.vector(beta.bar), sigma=V.beta, lower=c(-Inf, 0, rep(-Inf,N), -1), upper=c(Inf, 1, rep(Inf,N), 1)))
    epsilon     = tau - X%*%aux$betatau
    
    # sample sigmas
    ete         = diag(crossprod(epsilong))
    aux$sigmaf2 = (0.0003712422+ete[1])/rchisq(1, df=T+3)
    aux$sigmag2 = ete[2:N]/rchisq(N-1, df=T+3)
    for (c in 1:N){
      if (c==1){
        aux$sigma2[c]   = as.numeric(crossprod(epsilon[1:Tc[1],]))/rchisq(1, df=Tc[c]+3)
      } else {
        aux$sigma2[c]   = as.numeric(crossprod(epsilon[(sum(Tc[1:(c-1)])+1):sum(Tc[1:c]),]))/rchisq(1, df=Tc[c]+3)
      }
    }
    
    # sample hyper-parameters
    aux$phi.mu    = as.vector(rtmvnorm(1, mean=(aux$phi.sig/(N-1))*sum(aux$phi), sigma=(aux$phi.sig/(N-1)), lower=0, upper=1))
    draw          = 1.1
    while (draw>1){
      draw        = as.numeric(crossprod(aux$phi-aux$phi.mu))/rchisq(1, df=N+1)
    }
    aux$phi.sig   = draw
    aux$delta.mu  = rnorm(1, mean=sum(aux$betatau[3:(2+N)])/(N*aux$delta.sig + 1), sd=sqrt(1/(N*aux$delta.sig + 1)))
    aux$delta.sig = (1+as.numeric(crossprod(aux$betatau[3:(2+N)] - aux$delta.mu)))/rchisq(1, df=N+1)
    aux$s         = rgamma(1, scale=1/(0.5*sum(1/aux$sigmag2) + 1), shape=1.5*(N-1)+1)
    aux$s.sig     = rgamma(1, scale=1/(0.5*sum(1/aux$sigma2) + 1), shape=1.5*N+1)
    
    posterior$gamma[,s]     = aux$gamma
    posterior$phi[,s]       = aux$phi
    posterior$betatau[,s]   = aux$betatau
    posterior$sigmaf2[s]    = aux$sigmaf2
    posterior$sigmag2[,s]   = aux$sigmag2
    posterior$sigma2[,s]    = aux$sigma2
    posterior$phi.mu[s]     = aux$phi.mu
    posterior$phi.sig[s]    = aux$phi.sig
    posterior$delta.mu[s]   = aux$delta.mu
    posterior$delta.sig[s]  = aux$delta.sig
    posterior$s[s]          = aux$s
    posterior$s.sig[s]      = aux$s.sig
    
    if (s%%100==0){cat(s," ")}
  }
  
  output = list(
    posterior     = posterior,
    last.draw     = aux
  )
  return(output)
}
