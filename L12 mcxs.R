############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 12: SVAR Tools
############################################################

############################################################
# SVAR model of the Australian economy:
############################################################
rm(list=ls())
library(mvtnorm)
library(plot3D)
library(HDInterval)
set.seed(123456)

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
S       = 50000
h       = 8

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

pdf(file="data-foreign.pdf",height=9, width=6)
plot(y[,1:6], main="", xlab="", nc=1, lwd=2, col=mcxs1)
dev.off()

pdf(file="data-domestic.pdf",height=9, width=6)
plot(y[,7:12], main="", xlab="", nc=1, lwd=2, col=purple)
dev.off()


# create Y and X
############################################################
Y       = ts(y[5:122,], start=c(1991,3), frequency=4)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[5:122-i,])
}

t0          = proc.time() # read processor time

# MLE
############################################################
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/nrow(Y)
# round(A.hat,3)
# round(Sigma.hat,3)
# round(cov2cor(Sigma.hat),3)

# prior distribution
############################################################
kappa.1     = 1^2
kappa.2     = 100
kappa.3     = 1
A.prior     = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:13,] = kappa.3*diag(N)
V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N)))
S.prior     = diag(diag(Sigma.hat))
nu.prior    = N+1

# normal-inverse Wishard posterior parameters
############################################################
V.bar.inv   = t(X)%*%X + diag(1/diag(V.prior))
V.bar       = solve(V.bar.inv)
A.bar       = V.bar%*%(t(X)%*%Y + diag(1/diag(V.prior))%*%A.prior)
nu.bar      = nrow(Y) + nu.prior
S.bar       = S.prior + t(Y)%*%Y + t(A.prior)%*%diag(1/diag(V.prior))%*%A.prior - t(A.bar)%*%V.bar.inv%*%A.bar
S.bar.inv   = solve(S.bar)

# posterior draws 
############################################################
Sigma.posterior   = rWishart(S, df=nu.bar, Sigma=S.bar.inv)
Sigma.posterior   = apply(Sigma.posterior,3,solve)
Sigma.posterior   = array(Sigma.posterior,c(N,N,S))
A.posterior       = array(rnorm(prod(c(dim(A.bar),S))),c(dim(A.bar),S))
B.posterior       = array(NA,c(N,N,S))
L                 = t(chol(V.bar))
for (s in 1:S){
  cholSigma.s     = chol(Sigma.posterior[,,s])
  B.posterior[,,s]= t(cholSigma.s)
  A.posterior[,,s]= A.bar + L%*%A.posterior[,,s]%*%cholSigma.s
}
# round(apply(A.posterior,1:2,mean),4)
# round(apply(B.posterior,1:2,mean),4)

# Impulse response functions
# Forecast Error Variance Decomposition
############################################################
IRF.posterior     = array(NA,c(N,N,h+1,S))
IRF.inf.posterior = array(NA,c(N,N,S))
FEVD.posterior    = array(NA,c(N,N,h+1,S))
J                 = cbind(diag(N),matrix(0,N,N*(p-1)))
for (s in 1:S){
  A.bold          = rbind(t(A.posterior[2:(1+N*p),,s]),cbind(diag(N*(p-1)),matrix(0,N*(p-1),N)))
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

t1          = proc.time()
(t1-t0)[3]/60 # Time of computations in minutes

# save(IRF.posterior,IRF.inf.posterior, FEVD.posterior, file="irf-fevd-k002.RData")
save(IRF.posterior,IRF.inf.posterior, FEVD.posterior, file="irf-fevd-k1.RData")


# Plots of responses to domestic monetary policy shock
############################################################
load("irf-fevd-k1.RData")
IRFs.k1           = apply(IRF.posterior[7:12,10,,],1:2,mean)
IRFs.inf.k1       = apply(IRF.inf.posterior[7:12,10,],1,mean)
rownames(IRFs.k1) = colnames(Y)[7:12]

IRFs.k1.hdi    = apply(IRF.posterior[7:12,10,,],1:2,hdi, credMass=0.68)
hh          = 1:9

pdf(file="irf-au-mps.pdf", height=9, width=12)
par(mfrow=c(3,2), mar=c(4,4.5,2,2),cex.axis=1.5, cex.lab=1.5)
for (n in 1:6){
  ylims     = range(IRFs.k1[n,hh],IRFs.k1.hdi[,n,hh])
  plot(hh,IRFs.k1[n,hh], type="l", ylim=ylims, axes=FALSE, xlab="", ylab=rownames(IRFs.k1)[n])
  if (n==5 | n==6){
    axis(1,c(1,2,5,9),c("","1 quarter","1 year","2 years"))
  } else {
    axis(1,c(1,2,5,9),c("","","",""))
  }
  axis(2,c(ylims[1],0,ylims[2]),round(c(ylims[1],0,ylims[2]),3))
  polygon(c(hh,(h+1):1), c(IRFs.k1.hdi[1,n,hh],IRFs.k1.hdi[2,n,(h+1):1]), col=mcxs1.shade1,border=mcxs1.shade1)
  abline(h=0)
  lines(hh, IRFs.k1[n,hh],lwd=2,col=mcxs1)
}
dev.off()




load("irf-fevd-k002.RData")
IRFs.k002         = apply(IRF.posterior[7:12,10,,],1:2,mean)
IRFs.inf.k002     = apply(IRF.inf.posterior[7:12,10,],1,mean)
rownames(IRFs.k002) = colnames(Y)[7:12]
IRFs.k002.hdi    = apply(IRF.posterior[7:12,10,,],1:2,hdi, credMass=0.68)



pdf(file="irf-au-mps-both.pdf", height=9, width=12)
par(mfrow=c(3,2), mar=c(4,4.5,2,2),cex.axis=1.5, cex.lab=1.5)
for (n in 1:6){
  ylims     = range(IRFs.k1[n,hh],IRFs.k1.hdi[,n,hh],IRFs.k002[n,hh],IRFs.k002.hdi[,n,hh])
  plot(hh,IRFs.k1[n,hh], type="l", ylim=ylims, axes=FALSE, xlab="", ylab=rownames(IRFs.k1)[n])
  if (n==5 | n==6){
    axis(1,c(1,2,5,9),c("","1 quarter","1 year","2 years"))
  } else {
    axis(1,c(1,2,5,9),c("","","",""))
  }
  axis(2,c(ylims[1],0,ylims[2]),round(c(ylims[1],0,ylims[2]),3))
  polygon(c(hh,(h+1):1), c(IRFs.k1.hdi[1,n,hh],IRFs.k1.hdi[2,n,(h+1):1]), col=mcxs1.shade1,border=mcxs1.shade1)
  polygon(c(hh,(h+1):1), c(IRFs.k002.hdi[1,n,hh],IRFs.k002.hdi[2,n,(h+1):1]), col=mcxs2.shade1,border=mcxs2.shade1)
  abline(h=0)
  lines(hh, IRFs.k1[n,hh],lwd=2,col=mcxs1)
  lines(hh,IRFs.k002[n,hh],lwd=2,col=mcxs2)
}
dev.off()







# Plots of FEVD of Australian rgdp and p
############################################################
load("irf-fevd-k1.RData")
hh            = 1:(h+1)
fevd.au.rgdp  = apply(FEVD.posterior[8,,,],1:2,mean)
fevd.au.rgdp  = rbind(rep(0,h+1),apply(fevd.au.rgdp,2,cumsum))

fevd.au.p     = apply(FEVD.posterior[9,,,],1:2,mean)
fevd.au.p     = rbind(rep(0,h+1),apply(fevd.au.p,2,cumsum))

colors = c("deepskyblue1","deepskyblue2","deepskyblue","deepskyblue3","deepskyblue4","dodgerblue",
           "maroon1","maroon","maroon2","magenta","maroon3","maroon4")

pdf(file="fevd-au-gdp-k1.pdf", height=7, width=12)
par(mar=rep(4,4),cex.axis=1, cex.lab=0.8)
plot(hh,fevd.au.rgdp[1,], type="n", ylim=c(0,100), axes=FALSE, xlab="", ylab="")
axis(1,hh,c("","1 quarter","","","1 year","","","","2 years"))
axis(2,c(0,50,100),c("","FEVD[au.rgdp]",""))
for (n in 1:N){
  polygon(c(hh,(h+1):1), c(fevd.au.rgdp[n,hh],fevd.au.rgdp[n+1,(h+1):1]), col=colors[n],border=colors[n])
}
axis(4, (0.5*(fevd.au.rgdp[1:12,9]+fevd.au.rgdp[2:13,9]))[c(3,8,10)], c("us.mps","","au.mps"))
dev.off()


pdf(file="fevd-au-p-k1.pdf", height=7, width=12)
par(mar=rep(4,4),cex.axis=1, cex.lab=0.8)
plot(hh,fevd.au.p[1,], type="n", ylim=c(0,100), axes=FALSE, xlab="", ylab="")
axis(1,hh,c("","1 quarter","","","1 year","","","","2 years"))
axis(2,c(0,50,100),c("","FEVD[au.p]",""))
for (n in 1:N){
  polygon(c(hh,(h+1):1), c(fevd.au.p[n,hh],fevd.au.p[n+1,(h+1):1]), col=colors[n],border=colors[n])
}
axis(4, (0.5*(fevd.au.p[1:12,9]+fevd.au.p[2:13,9]))[c(3,10)], c("us.mps","au.mps"))
dev.off()






load("irf-fevd-k002.RData")
hh            = 1:(h+1)
fevd.au.rgdp  = apply(FEVD.posterior[8,,,],1:2,mean)
fevd.au.rgdp  = rbind(rep(0,h+1),apply(fevd.au.rgdp,2,cumsum))

fevd.au.p     = apply(FEVD.posterior[9,,,],1:2,mean)
fevd.au.p     = rbind(rep(0,h+1),apply(fevd.au.p,2,cumsum))

colors = c("deepskyblue1","deepskyblue2","deepskyblue","deepskyblue3","deepskyblue4","dodgerblue",
           "maroon1","maroon","maroon2","magenta","maroon3","maroon4")

pdf(file="fevd-au-gdp-k002.pdf", height=7, width=12)
par(mar=rep(4,4),cex.axis=1, cex.lab=0.8)
plot(hh,fevd.au.rgdp[1,], type="n", ylim=c(0,100), axes=FALSE, xlab="", ylab="")
axis(1,hh,c("","1 quarter","","","1 year","","","","2 years"))
axis(2,c(0,50,100),c("","FEVD[au.rgdp]",""))
for (n in 1:N){
  polygon(c(hh,(h+1):1), c(fevd.au.rgdp[n,hh],fevd.au.rgdp[n+1,(h+1):1]), col=colors[n],border=colors[n])
}
axis(4, (0.5*(fevd.au.rgdp[1:12,9]+fevd.au.rgdp[2:13,9]))[c(3,8,10)], c("us.mps","","au.mps"))
dev.off()


pdf(file="fevd-au-p-k002.pdf", height=7, width=12)
par(mar=rep(4,4),cex.axis=1, cex.lab=0.8)
plot(hh,fevd.au.p[1,], type="n", ylim=c(0,100), axes=FALSE, xlab="", ylab="")
axis(1,hh,c("","1 quarter","","","1 year","","","","2 years"))
axis(2,c(0,50,100),c("","FEVD[au.p]",""))
for (n in 1:N){
  polygon(c(hh,(h+1):1), c(fevd.au.p[n,hh],fevd.au.p[n+1,(h+1):1]), col=colors[n],border=colors[n])
}
axis(4, (0.5*(fevd.au.p[1:12,9]+fevd.au.p[2:13,9]))[c(3,10)], c("us.mps","au.mps"))
dev.off()




