############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 15: Modeling effects of monetary policy
############################################################

############################################################
# SVAR model of the Australian economy:
############################################################
rm(list=ls())
library(mvtnorm)
library(plot3D)
library(HDInterval)
set.seed(123456)

# upload the codes
source("L16 codes.R")

# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"
purple = "#b02442"

mcxs1.rgb   = col2rgb(mcxs1)
mcxs1.shade1= rgb(mcxs1.rgb[1],mcxs1.rgb[2],mcxs1.rgb[3], alpha=120, maxColorValue=255)
mcxs2.rgb   = col2rgb(mcxs2)
mcxs2.shade1= rgb(mcxs2.rgb[1],mcxs2.rgb[2],mcxs2.rgb[3], alpha=120, maxColorValue=255)

# setup
############################################################
N       = 12
p       = 4
K       = 1+N*p
h       = 6
S.burnin= 100
S       = 5000

# upload and transform the variables
############################################################

us.rgdp     = as.data.frame(read.csv("./AU-SVAR-data/US-RGDP.csv", header=TRUE))
us.rgdp     = ts(log(us.rgdp[,2]), start=c(1990,3), frequency=4)

us.cpi      = as.data.frame(read.csv("./AU-SVAR-data/US-CPI.csv", header=TRUE))
us.cpi      = ts(log(us.cpi[,2]), start=c(1990,3), frequency=4)

us.ffr      = as.data.frame(read.csv("./AU-SVAR-data/US-FEDFUNDS.csv", header=TRUE))
us.ffr      = ts(matrix(us.ffr[,2],nrow=3)[1,], start=c(1990,3), frequency=4)

us.sp5      = as.data.frame(read.csv("./AU-SVAR-data/US-SP500.csv", header=TRUE))
us.sp5      = ts(log(matrix(us.sp5[,2],nrow=3)[1,]), start=c(1990,3), frequency=4)

au.tot      = as.data.frame(read.csv("./AU-SVAR-data/AU-ToT.csv", header=TRUE))
au.tot      = ts(log(au.tot[,2]), start=c(1990,3), frequency=4)

au.cpi      = as.data.frame(read.csv("./AU-SVAR-data/AU-CPI.csv", header=TRUE))
au.cpi      = ts(log(au.cpi[,2]), start=c(1990,3), frequency=4)

au.rex      = as.data.frame(read.csv("./AU-SVAR-data/AU-EX.csv", header=TRUE))
au.rex      = ts(log(au.rex[,2]), start=c(1990,3), frequency=4)
au.rex      = au.rex - au.cpi

au.aord     = as.data.frame(read.csv("./AU-SVAR-data/AU-AORD.csv", header=TRUE))
au.aord     = ts(log(matrix(au.aord[,2],nrow=3)[1,]), start=c(1990,3), frequency=4)

au.rgne     = as.data.frame(read.csv("./AU-SVAR-data/AU-GNE.csv", header=TRUE))
au.rgne     = ts(log(au.rgne[,2]), start=c(1990,3), frequency=4)
au.rgne     = au.rgne - au.cpi

au.rgdp     = as.data.frame(read.csv("./AU-SVAR-data/AU-GDP.csv", header=TRUE))
au.rgdp     = ts(log(au.rgdp[,2]), start=c(1990,3), frequency=4)
au.rgdp     = au.rgdp - au.cpi

au.cr       = as.data.frame(read.csv("./AU-SVAR-data/AU-CR.csv", header=TRUE))
au.cr       = ts(matrix(au.cr[,2],nrow=3)[1,], start=c(1990,3), frequency=4)

au.rtwi     = as.data.frame(read.csv("./AU-SVAR-data/AU-RTWI.csv", header=TRUE))
au.rtwi     = ts(log(au.rtwi[,2]), start=c(1990,3), frequency=4)

y           = cbind(us.rgdp,
                    us.cpi,
                    us.ffr,
                    us.sp5,
                    au.tot,
                    au.rex,
                    au.rgne,
                    au.rgdp,
                    au.cpi,
                    au.cr,
                    au.rtwi,
                    au.aord
)

# pdf(file="data-foreign.pdf",height=9, width=6)
# plot(y[,1:6], main="", xlab="", nc=1, lwd=2, col=mcxs1)
# dev.off()
# 
# pdf(file="data-domestic.pdf",height=9, width=6)
# plot(y[,7:12], main="", xlab="", nc=1, lwd=2, col=purple)
# dev.off()

# create Y and X
############################################################
Y       = ts(y[5:122,], start=c(1991,3), frequency=4)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[5:122-i,])
}
Y       = t(Y)
X       = t(X)


# set the priors
############################################################
kappa1  = 1       # autoregressive slope shrinkage
kappa4  = 0.95       # VAR prior persistence
kappa2  = 1      # constant term shrinkage
kappa3  = 1      # contemporaneous effects shrinkage

priors  = list(
  B     = cbind(rep(0,N), kappa4*diag(N), matrix(0, N, (p-1)*N)),
  # Omega = diag(c(kappa2,kappa1*rep(1,N*p))),
  Omega = diag(c(kappa2,kappa1*((1:p)^(-2))%x%rep(1,N))),
  S     = kappa3*diag(N),
  nu    = N
)

# compute posterior distribution parameters, matrices V, and starting values
############################################################
Omega.inv   = solve(priors$Omega)
Omega.post.inv = X%*%t(X) + Omega.inv
Omega.post  = solve( Omega.post.inv )
B.post      = (Y%*%t(X) + priors$B%*%Omega.inv) %*% Omega.post
S.post      = solve(Y%*%t(Y) + solve(priors$S) + priors$B%*%Omega.inv%*%t(priors$B) - B.post%*%Omega.post.inv%*%t(B.post) )
nu.post     = ncol(Y) + priors$nu

posteriors = list(
  B        = B.post,
  Omega    = Omega.post,
  S        = S.post,
  nu       = nu.post
)

# V_n matrices
SOE.V      = vector("list",N)
for (n in 1:N){
  SOE.V[[n]]   = cbind(diag(n),matrix(0,n,N-n))
}

B0.initial = matrix(0,N,N)
for (n in 1:N){
  unrestricted    = apply(SOE.V[[n]],2,sum)==1
  B0.initial[n,unrestricted] = rnorm(sum(unrestricted))
}

# run the posterior simulations and save their output
############################################################
# burn-in run: sampling B0 from the posterior distribution using Gibbs
t0                  = proc.time()
B0.posterior        = rgn(n=S.burnin, S=S.post, nu=nu.post, V=SOE.V, B0.initial=B0.initial)
t1                  = proc.time()
(t1-t0)/60
# sampling B0 from the posterior distribution using Gibbs
t0                  = proc.time()
B0.posterior        = rgn(n=S, S=S.post, nu=nu.post, V=SOE.V, B0.initial=B0.posterior[,,S.burnin])
t1                  = proc.time()
(t1-t0)/60
# normalisation
B0.hat              = t(chol((nu.post-N)*S.post))                       # normalisation using this B0.hat should work
SOE.B0.posterior    = normalize.Gibbs.output.parallel(B0.posterior,B0.hat=B0.hat)
t2                  = proc.time()
(t2-t1)/60

# sample B+ from the normal conditional posterior
t2                  = proc.time()
SOE.Bp.posterior    = rnorm.ngn(SOE.B0.posterior, B=B.post,Omega=Omega.post)
t3                  = proc.time()
(t3-t2)/60

save(SOE.B0.posterior,SOE.Bp.posterior, priors, posteriors, file="AU-posterior.RData")



# Impulse response functions
# Forecast Error Variance Decomposition
############################################################

t4          = proc.time()
B.posterior       = array(NA,c(N,N,S))
A.posterior       = array(NA,c(N,K,S))
for (s in 1:S){
  B               = solve(SOE.B0.posterior[,,s])
  B.posterior[,,s]= B
  A.posterior[,,s]= B %*% SOE.Bp.posterior[,,s]
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
      IRF.posterior[,,i,s]        = J %*% A.bold.power %*% t(J) %*% IRF.posterior[,,i-1,s]
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

# save(IRF.posterior,IRF.inf.posterior, FEVD.posterior, file="irf-fevd-k002.RData")
save(IRF.posterior,IRF.inf.posterior, FEVD.posterior, file="AU-irf.RData")


# plot IRFs and FEVDs
############################################################
IRF.posterior.mps = IRF.posterior[7:12,10,,]
IRFs.k1           = apply(IRF.posterior.mps,1:2,median)
IRF.posterior.mps = IRF.posterior.mps*(0.25/IRFs.k1[4,1])
IRFs.k1           = apply(IRF.posterior.mps,1:2,median)
IRFs.inf.k1       = apply(IRF.posterior.mps,1,mean)
rownames(IRFs.k1) = colnames(y)[7:12]

IRFs.k1.hdi    = apply(IRF.posterior.mps,1:2,hdi, credMass=0.68)
hh          = 1:(h+1)

pdf(file="irf-au-mps.pdf", height=9, width=12)
par(mfrow=c(3,2), mar=c(4,4.5,2,2),cex.axis=1.5, cex.lab=1.5)
for (n in 1:6){
  ylims     = range(IRFs.k1[n,hh],IRFs.k1.hdi[,n,1:4],0)
  plot(hh,IRFs.k1[n,hh], type="l", ylim=ylims, axes=FALSE, xlab="", ylab=rownames(IRFs.k1)[n])
  # if (n==5 | n==6){
  #   axis(1,c(1,2,5,9),c("","1 quarter","1 year","2 years"))
  # } else {
  #   axis(1,c(1,2,5,9),c("","","",""))
  # }
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




# Plots of FEVD of Australian rgdp and p
############################################################
hh            = 1:(h+1)
fevd.au.rgdp  = apply(FEVD.posterior[8,,,],1:2,mean)
fevd.au.rgdp  = rbind(rep(0,h+1),apply(fevd.au.rgdp,2,cumsum))

fevd.au.p     = apply(FEVD.posterior[9,,,],1:2,mean)
fevd.au.p     = rbind(rep(0,h+1),apply(fevd.au.p,2,cumsum))

colors = c("deepskyblue1","deepskyblue2","deepskyblue","deepskyblue3","deepskyblue4","dodgerblue",
           "maroon1","maroon","maroon2","magenta","maroon3","maroon4")

pdf(file="fevd-au-gdp.pdf", height=7, width=12)
par(mar=rep(4,4),cex.axis=1, cex.lab=0.8)
plot(hh,fevd.au.rgdp[1,], type="n", ylim=c(0,100), axes=FALSE, xlab="", ylab="")
# axis(1,hh,c("","1 quarter","","","1 year","","","","2 years"))
axis(1,hh,c("","1 quarter","","","1 year","","6 quarters"))
axis(2,c(0,50,100),c("","FEVD[au.rgdp]",""))
for (n in 1:N){
  polygon(c(hh,(h+1):1), c(fevd.au.rgdp[n,hh],fevd.au.rgdp[n+1,(h+1):1]), col=colors[n],border=colors[n])
}
axis(4, (0.5*(fevd.au.rgdp[1:12,7]+fevd.au.rgdp[2:13,7]))[c(3,8,10)], c("us.mps","","au.mps"))
dev.off()


pdf(file="fevd-au-p.pdf", height=7, width=12)
par(mar=rep(4,4),cex.axis=1, cex.lab=0.8)
plot(hh,fevd.au.p[1,], type="n", ylim=c(0,100), axes=FALSE, xlab="", ylab="")
axis(1,hh,c("","1 quarter","","","1 year","","6 quarters"))
axis(2,c(0,50,100),c("","FEVD[au.p]",""))
for (n in 1:N){
  polygon(c(hh,(h+1):1), c(fevd.au.p[n,hh],fevd.au.p[n+1,(h+1):1]), col=colors[n],border=colors[n])
}
axis(4, (0.5*(fevd.au.p[1:12,7]+fevd.au.p[2:13,7]))[c(3,10)], c("us.mps","au.mps"))
dev.off()




# plot IRFs to the US mps
############################################################
IRF.posterior.mps = IRF.posterior[1:12,3,,]
IRFs.k1           = apply(IRF.posterior.mps,1:2,median)
IRF.posterior.mps = IRF.posterior[7:12,3,,]
IRF.posterior.mps = IRF.posterior.mps*(0.25/IRFs.k1[3,1])
IRFs.k1           = apply(IRF.posterior.mps,1:2,median)
IRFs.inf.k1       = apply(IRF.posterior.mps,1,mean)
rownames(IRFs.k1) = colnames(y)[7:12]

IRFs.k1.hdi    = apply(IRF.posterior.mps,1:2,hdi, credMass=0.68)
hh          = 1:(h+1)

pdf(file="irf-us-mps.pdf", height=9, width=12)
par(mfrow=c(3,2), mar=c(4,4.5,2,2),cex.axis=1.5, cex.lab=1.5)
for (n in 1:6){
  ylims     = range(IRFs.k1[n,hh],IRFs.k1.hdi[,n,1:4],0)
  plot(hh,IRFs.k1[n,hh], type="l", ylim=ylims, axes=FALSE, xlab="", ylab=rownames(IRFs.k1)[n])
  # if (n==5 | n==6){
  #   axis(1,c(1,2,5,9),c("","1 quarter","1 year","2 years"))
  # } else {
  #   axis(1,c(1,2,5,9),c("","","",""))
  # }
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


(t5-t0)[3]/60 # Time of computations in minutes
system("say do not despair")

























# Foreign vs. domestic shocks
############################################################
B.posterior.id    = B.posterior
B.posterior       = array(NA,c(N,N,S))

for (s in 1:S){
  # compute the block diagonal rotation matrix
  X           = matrix(rnorm(N^2),6,6)
  QR          = qr(X, tol = 1e-10)
  Q1          = qr.Q(QR,complete=TRUE)
  R1          = qr.R(QR,complete=TRUE)
  Q1          = t(Q1 %*% diag(sign(diag(R1))))
  X           = matrix(rnorm(N^2),6,6)
  QR          = qr(X, tol = 1e-10)
  Q2          = qr.Q(QR,complete=TRUE)
  R2          = qr.R(QR,complete=TRUE)
  Q2          = t(Q2 %*% diag(sign(diag(R2))))
  Q           = rbind(cbind(Q1,matrix(0,6,6)),cbind(matrix(0,6,6),Q2))
  
  B.posterior[,,s]= B.posterior.id[,,s]%*%t(Q)
}

# Forecast Error Variance Decomposition
############################################################
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
      IRF.posterior[,,i,s]        = J %*% A.bold.power %*% t(J) %*% IRF.posterior[,,i-1,s]
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

save(A.posterior,B.posterior,IRF.posterior,IRF.inf.posterior, FEVD.posterior, file="AU-fd-irf.RData")




# Plots of FEVD of Australian rgdp
############################################################

hh            = 1:(h+1)
fevd.au.rgdp.f= apply(FEVD.posterior[8,1:6,,],2:3,sum)
fevd.au.rgdp.d= apply(FEVD.posterior[8,7:12,,],2:3,sum)

fevd.au.rgdp  = rbind(apply(fevd.au.rgdp.f,1,mean), apply(fevd.au.rgdp.d,1,mean))
fevd.au.rgdp  = rbind(rep(0,h+1),apply(fevd.au.rgdp,2,cumsum))

fevd.au.p.f= apply(FEVD.posterior[9,1:6,,],2:3,sum)
fevd.au.p.d= apply(FEVD.posterior[9,7:12,,],2:3,sum)

fevd.au.p  = rbind(apply(fevd.au.p.f,1,mean), apply(fevd.au.p.d,1,mean))
fevd.au.p  = rbind(rep(0,h+1),apply(fevd.au.p,2,cumsum))

colors = c("dodgerblue","maroon")

pdf(file="fevd-au-fd-gdp.pdf", height=7, width=12)
par(mar=rep(4,4),cex.axis=1, cex.lab=0.8)
plot(hh,fevd.au.rgdp[1,], type="n", ylim=c(0,100), axes=FALSE, xlab="", ylab="")
axis(1,c(1,2,5,7),c("","1 quarter","1 year","6 quarters"))
axis(2,c(0,50,100),c("","FEVD[au.rgdp]",""))
for (n in 1:2){
  polygon(c(hh,(h+1):1), c(fevd.au.rgdp[n,hh],fevd.au.rgdp[n+1,(h+1):1]), col=colors[n],border=colors[n])
}
axis(4, c(0.5*(fevd.au.rgdp[1:2,7]+fevd.au.rgdp[2:3,7])), c("foreign","domestic"))
dev.off()


pdf(file="fevd-au-fd-p.pdf", height=7, width=12)
par(mar=rep(4,4),cex.axis=1, cex.lab=0.8)
plot(hh,fevd.au.p[1,], type="n", ylim=c(0,100), axes=FALSE, xlab="", ylab="")
axis(1,c(1,2,5,7),c("","1 quarter","1 year","6 quarters"))
axis(2,c(0,50,100),c("","FEVD[au.p]",""))
for (n in 1:2){
  polygon(c(hh,(h+1):1), c(fevd.au.p[n,hh],fevd.au.p[n+1,(h+1):1]), col=colors[n],border=colors[n])
}
axis(4, c(0.5*(fevd.au.p[1:2,7]+fevd.au.p[2:3,7])), c("foreign","domestic"))
dev.off()





