############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 16: Modeling effects of monetary policy
############################################################

rm(list=ls())
set.seed(1234567)
source("L16 codes.R")

T       = 200
N       = 3
p       = 1
K       = 1 + N * p
S.burnin = 100
S       = 100



y       = apply(matrix(rnorm(T * N), ncol = N), 2, cumsum)


# create Y and X
############################################################
Y       = y[(p+1):nrow(y),]
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[((p+1):nrow(y))-i,])
}
Y       = t(Y)
X       = t(X)

# set the priors
############################################################
kappa1  = .1       # autoregressive slope shrinkage
kappa4  = 1       # VAR prior persistence
kappa2  = 10      # constant term shrinkage
kappa3  = 10      # contemporaneous effects shrinkage

priors  = list(
  B     = cbind(rep(0,N), kappa4*diag(N), matrix(0, N, (p-1)*N)),
  Omega = diag(c(kappa2,kappa1*((1:p)^(-2))%x%rep(1,N))),
  # Omega = diag(c(kappa2,kappa1*rep(1,N*p))),
  S     = kappa3*diag(N),
  nu    = N
)

# compute posterior distribution parameters, matrices V, and starting values
############################################################
Omega.inv   = solve(priors$Omega)
Omega.post.inv = X%*%t(X) + Omega.inv
Omega.post  = solve( Omega.post.inv )
B.post      = (Y%*%t(X) + priors$B%*%Omega.inv) %*% Omega.post
S.post      = Y%*%t(Y) + solve(priors$S) + priors$B%*%Omega.inv%*%t(priors$B) - B.post%*%Omega.post.inv%*%t(B.post) 
nu.post     = ncol(Y) + priors$nu

posteriors = list(
  B        = B.post,
  Omega    = Omega.post,
  S        = S.post,
  nu       = nu.post
)

FF.V          = vector("list",N)
for (n in 1:N){
  FF.V[[n]]   = cbind(diag(n),matrix(0,n,N-n))
}

B0.initial = matrix(0,N,N)
for (n in 1:N){
  unrestricted    = apply(FF.V[[n]],2,sum)==1
  B0.initial[n,unrestricted] = rnorm(sum(unrestricted))
}


# run the posterior simulations and save their output
############################################################
# burn-in run: sampling B0 from the posterior distribution using Gibbs
t0                  = proc.time()
B0.posterior        = rgn(n=S.burnin, S.inv=S.post, nu=nu.post, V=FF.V, B0.initial=B0.initial)
t1                  = proc.time()
(t1-t0)/60
# sampling B0 from the posterior distribution using Gibbs
t0                  = proc.time()
B0.posterior        = rgn(n=S, S.inv=S.post, nu=nu.post, V=FF.V, B0.initial=B0.posterior[,,S.burnin])
t1                  = proc.time()
(t1-t0)/60
# normalisation
B0.hat              = t(chol((nu.post-N)*S.post))                       # normalisation using this B0.hat should work
FF.B0.posterior    = normalize.Gibbs.output.parallel(B0.posterior,B0.hat=B0.hat)
t2                  = proc.time()
(t2-t1)/60


apply(FF.B0.posterior,1:2, mean)






# sample B+ from the normal conditional posterior
t2                  = proc.time()
FF.Bp.posterior    = rnorm.ngn(FF.B0.posterior, B=B.post,Omega=Omega.post)
t3                  = proc.time()
(t3-t2)/60

save(FF.B0.posterior,FF.Bp.posterior, priors, posteriors, file="FF-posterior.RData")





# Impulse response functions
# Forecast Error Variance Decomposition
############################################################

t4          = proc.time()
B.posterior       = array(NA,c(N,N,S))
A.posterior       = array(NA,c(N,K,S))
for (s in 1:S){
  B               = solve(FF.B0.posterior[,,s])
  B.posterior[,,s]= B
  A.posterior[,,s]= B %*% FF.Bp.posterior[,,s]
}

IRF.posterior     = array(NA,c(N,N,h+1,S))
IRF.inf.posterior = array(NA,c(N,N,S))
FEVD.posterior    = array(NA,c(N,N,h+1,S))
J                 = cbind(diag(N),matrix(0,N,N*(p-1)))
for (s in 1:S){
  A.bold          = rbind(A.posterior[,2:(1+N*p),s],cbind(diag(N*(p-1)),matrix(0,N*(p-1),N)))
  IRF.inf.posterior[,,s]          = J %*% solve(diag(N*p)-A.bold) %*% t(J) %*% B.posterior[,,s]
  A.bold.power    = A.bold
  for (i in 1:(h+1)){
    if (i==1){
      IRF.posterior[,,i,s]        = B.posterior[,,s]
    } else {
      IRF.posterior[,,i,s]        = J %*% A.bold.power %*% t(J) %*% IRF.posterior[,,1,s]
      A.bold.power                = A.bold.power %*% A.bold
    }
    for (n in 1:N){
      for (nn in 1:N){
        FEVD.posterior[n,nn,i,s]  = sum(IRF.posterior[n,nn,1:i,s]^2)
      }
    }
    FEVD.posterior[,,i,s]         = diag(1/apply(FEVD.posterior[,,i,s],1,sum))%*%FEVD.posterior[,,i,s]
  }
}
FEVD.posterior    = 100*FEVD.posterior

t5          = proc.time()
(t5-t4)[3]/60 # Time of computations in minutes

save(IRF.posterior,IRF.inf.posterior, FEVD.posterior, file="FF-irf.RData")





# plot IRFs and FEVDs
############################################################
IRF.posterior.mps = IRF.posterior[,4,,]
IRFs.k1           = apply(IRF.posterior.mps,1:2,median)
IRF.posterior.mps = IRF.posterior.mps*(0.25/IRFs.k1[4,1])
IRFs.k1           = apply(IRF.posterior.mps,1:2,median)
IRFs.inf.k1       = apply(IRF.posterior.mps,1,mean)
rownames(IRFs.k1) = colnames(y)

IRFs.k1.hdi    = apply(IRF.posterior.mps,1:2,hdi, credMass=0.68)
hh          = 1:(h+1)

pdf(file="FF-irf-mps.pdf", height=9, width=12)
par(mfrow=c(4,2), mar=c(4,4.5,2,2),cex.axis=1.5, cex.lab=1.5)
for (n in 1:N){
  ylims     = range(IRFs.k1[n,hh],IRFs.k1.hdi[,n,1:4],0)
  plot(hh,IRFs.k1[n,hh], type="l", ylim=ylims, axes=FALSE, xlab="", ylab=rownames(IRFs.k1)[n])
  if (n==5 | n==6){
    axis(1,c(1,2,5,7),c("","1 quarter","1 year","6 quarters"))
  } else {
    axis(1,c(1,2,5,7),c("","","",""))
  }
  axis(2,c(ylims[1],0,ylims[2]),round(c(ylims[1],0,ylims[2]),3))
  polygon(c(hh,(h+1):1), c(IRFs.k1.hdi[1,n,hh],IRFs.k1.hdi[2,n,(h+1):1]), col=mcxs1.shade1,border=mcxs1.shade1)
  abline(h=0)
  lines(hh, IRFs.k1[n,hh],lwd=2,col=mcxs1)
}
dev.off()

t1         = proc.time()
(t1-t0)[3]/60

system("say do not despair")
