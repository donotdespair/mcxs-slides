
############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 16: Modeling effects of monetary policy
############################################################
library(mvtnorm)
library(parallel)

# useful functions
############################################################
orthogonal.complement.matrix.TW = function(x){
  # x is a mxn matrix and m>n
  # the function returns a mx(m-n) matrix, out, that is an orthogonal complement of x, i.e.:
  # t(x)%*%out = 0 and det(cbind(x,out))!=0
  N     = dim(x)
  tmp   = qr.Q(qr(x, tol = 1e-10),complete=TRUE)
  out   = as.matrix(tmp[,(N[2]+1):N[1]])
  return(out)
}

r.conditional.generalized.normal = function(S.inv, nu, Vn, n, B0){
  # A function to sample a random draw from a conditional generalized normal distribution
  # of the unrestricted elements of the n-th row of matrix B0 
  # given the parameters from the remaining rows of B0
  # Depends on package mvtnorm
  # use: library(rmvtnorm)
  
  rn            = nrow(Vn)
  Un            = chol(nu*solve(Vn%*%S.inv%*%t(Vn)))
  w             = t(orthogonal.complement.matrix.TW(t(B0[-n,])))
  w1            = w %*% t(Vn) %*% t(Un) / sqrt(as.numeric(w %*% t(Vn) %*% t(Un) %*% Un %*% Vn %*% t(w)))
  if (rn>1){
    Wn          = cbind(t(w1),orthogonal.complement.matrix.TW(t(w1)))
  } else {
    Wn          = w1
  }
  alpha         = rep(NA,rn)
  u             = rmvnorm(1,rep(0,nu+1),(1/nu)*diag(nu+1))
  alpha[1]      = sqrt(as.numeric(u%*%t(u)))
  if (runif(1)<0.5){
    alpha[1]    = -alpha[1]
  }
  if (rn>1){
    alpha[2:rn] = rmvnorm(1,rep(0,nrow(Vn)-1),(1/nu)*diag(rn-1))
  }
  bn            = alpha %*% Wn %*% Un
  B0n           = bn %*% Vn
  
  output        = list(bn=bn, B0n=B0n)
  return(output)
}




rgn             = function(n,S.inv,nu,V,B0.initial){
  # This function simulates draws for the unrestricted elements 
  # of the conteporaneous relationships matrix of an SVAR model
  # from a generalized-normal distribution according to algorithm 
  # by Waggoner & Zha (2003, JEDC)
  # n     - a positive integer, the number of draws to be sampled
  # S     - an NxN positive definite matrix, a parameter of the generalized-normal distribution
  # nu    - a positive scalar, degrees of freedom parameter
  # V     - an N-element list, with fixed matrices
  # B0.initial - an NxN matrix, of initial values of the parameters
  
  N             = nrow(B0.initial)
  no.draws      = n
  
  B0            = array(NA, c(N,N,no.draws))
  B0.aux        = B0.initial
  
  for (i in 1:no.draws){
    for (n in 1:N){
      rn            = nrow(V[[n]])
      Un            = chol(nu*solve(V[[n]]%*%S.inv%*%t(V[[n]])))
      w             = t(orthogonal.complement.matrix.TW(t(B0.aux[-n,])))
      w1            = w %*% t(V[[n]]) %*% t(Un) / sqrt(as.numeric(w %*% t(V[[n]]) %*% t(Un) %*% Un %*% V[[n]] %*% t(w)))
      if (rn>1){
        Wn          = cbind(t(w1),orthogonal.complement.matrix.TW(t(w1)))
      } else {
        Wn          = w1
      }
      alpha         = rep(NA,rn)
      u             = rmvnorm(1,rep(0,nu+1),(1/nu)*diag(nu+1))
      alpha[1]      = sqrt(as.numeric(u%*%t(u)))
      if (runif(1)<0.5){
        alpha[1]    = -alpha[1]
      }
      if (rn>1){
        alpha[2:rn] = rmvnorm(1,rep(0,nrow(V[[n]])-1),(1/nu)*diag(rn-1))
      }
      bn            = alpha %*% Wn %*% Un
      B0.aux[n,]    = bn %*% V[[n]]
    }
    B0[,,i]         = B0.aux
  }
  
  return(B0)
}


normalization.wz2003  = function(B0,B0.hat.inv, Sigma.inv, diag.signs){
  # This function normalizes a matrix of contemporaneous effects
  # according to the algorithm by Waggoner & Zha (2003, JOE)
  # B0        - an NxN matrix, to be normalized
  # B0.hat    - an NxN matrix, a normalized matrix
  
  N                 = nrow(B0)
  K                 = 2^N
  distance          = rep(NA,K)
  for (k in 1:K){
    B0.tmp.inv      = solve(diag(diag.signs[k,]) %*% B0)
    distance[k]     = sum(
      unlist(
        lapply(1:N,
               function(n){
                 t(B0.tmp.inv - B0.hat.inv)[n,] %*%Sigma.inv %*% t(B0.tmp.inv - B0.hat.inv)[n,]
               }
        )))
  }
  B0.out            = diag(diag.signs[which.min(distance),]) %*% B0
  
  return(B0.out)
}

normalize.Gibbs.output.parallel          = function(B0.posterior,B0.hat){
  # This function normalizes the Gibbs sampler output from function rgn
  # using function normalization.wz2003 
  # B0.posterior  - a list, output from function rgn
  # B0.hat        - an NxN matrix, a normalized matrix
  
  N                 = nrow(B0.hat)
  K                 = 2^N
  
  B0.hat.inv        = solve(B0.hat)
  Sigma.inv         = t(B0.hat)%*%B0.hat
  
  diag.signs        = matrix(NA,2^N,N)
  for (n in 1:N){
    diag.signs[,n]  = kronecker(c(-1,1),rep(1,2^(n-1)))
  }
  
  B0.posterior.n    = mclapply(1:dim(B0.posterior)[3],function(i){
    normalization.wz2003(B0=B0.posterior[,,i],B0.hat.inv, Sigma.inv, diag.signs)
  },mc.cores=4
  )
  B0.posterior.n  = simplify2array(B0.posterior.n)
  
  return(B0.posterior.n)
}




normalize.Gibbs.output          = function(B0.posterior,B0.hat){
  # This function normalizes the Gibbs sampler output from function rgn
  # using function normalization.wz2003 
  # B0.posterior  - a list, output from function rgn
  # B0.hat        - an NxN matrix, a normalized matrix
  
  N                 = nrow(B0.hat)
  K                 = 2^N
  
  B0.hat.inv        = solve(B0.hat)
  Sigma.inv         = solve(B0.hat.inv %*% t(B0.hat.inv))
  
  diag.signs        = matrix(NA,2^N,N)
  for (n in 1:N){
    diag.signs[,n]  = kronecker(c(-1,1),rep(1,2^(n-1)))
  }
  
  for (i in 1:dim(B0.posterior)[3]){
    if (i%%100==0){ cat(i," ")}
    norm.post                   = normalization.wz2003(B0=B0.posterior[,,i],B0.hat.inv, Sigma.inv, diag.signs)
    B0.posterior[,,i]           = norm.post
  }
  return(B0.posterior)
}

rnorm.ngn       = function(B0.posterior,B,Omega){
  # This function simulates draws for the multivariate normal distribution
  # of the autoregressive slope matrix of an SVAR model
  # from a normal-generalized-normal distribution according to algorithm 
  # by Waggoner & Zha (2003, JEDC)
  # B0.posterior  - a list, output from function rgn
  # B             - an NxK matrix, a parameter determining the mean of the multivariate conditionally normal distribution given B0
  # Omega         - a KxK positive definite matrix, a covariance matrix of the multivariate normal distribution
  
  N             = nrow(B)
  K             = ncol(B)
  no.draws      = dim(B0.posterior)[3]
  L             = t(chol(Omega))
  
  Bp.posterior  = lapply(1:no.draws,function(i){
    Bp          = matrix(NA, N, K)
    for (n in 1:N){
      Bp[n,]    = as.vector(t(B0.posterior[n,,i] %*% B) + L%*%rnorm(K))
    }
    return(Bp)
  })
  Bp.posterior  = simplify2array(Bp.posterior)
  return(Bp.posterior)
}
