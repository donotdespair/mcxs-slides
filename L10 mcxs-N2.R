############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 10: Forecasting with Large Bayesian VARs
############################################################

############################################################
# Forecasting a bivariate series of Australian:
# + real GDP
# + real GDP deflator - inflation
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

mcxs1.rgb   = col2rgb(mcxs1)
mcxs1.shade1= rgb(mcxs1.rgb[1],mcxs1.rgb[2],mcxs1.rgb[3], alpha=50, maxColorValue=255)
mcxs2.rgb   = col2rgb(mcxs2)
mcxs2.shade1= rgb(mcxs2.rgb[1],mcxs2.rgb[2],mcxs2.rgb[3], alpha=50, maxColorValue=255)

# setup
############################################################
N       = 2
p       = 4
K       = 1+p*N
S       = 50000
h       = 8

# upload data and create Y and X
############################################################
y       = as.data.frame(read.csv("ausmacrodata-2016.csv",header=TRUE))
y       = ts(y[,2:118], start=c(1985,2), frequency=4)
y       = y[,c(1,72,2:71,73:117)]
y       = y[,1:2]  

Y       = ts(y[(p+1):nrow(y),], start=c(1986,2), frequency=4)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[(p+1):nrow(y)-i,])
}
T       = nrow(Y)


# prior distribution
############################################################
kappa.1     = 0.02^2
kappa.2     = 100
A.prior     = matrix(0,K,N)
V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N)))
V.prior.inv = diag(1/c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N)))
s.ols       = rep(NA,N)
for (n in 1:N){
  s.ols[n]  = var(ar(x=Y[,n], aic=FALSE, order.max=8, method="ols")$resid[9:T])
}
S.prior     = diag(s.ols)
nu.prior    = N+1

# normal-inverse Wishard posterior parameters
############################################################
V.bar.inv     = t(X)%*%X + V.prior.inv
V.bar.inv.chol= chol(V.bar.inv)
A.bar.tmp     = t(X)%*%Y + diag(1/diag(V.prior))%*%A.prior
A.bar         = backsolve(V.bar.inv.chol, forwardsolve(t(V.bar.inv.chol), A.bar.tmp))

nu.bar        = T + nu.prior
S.bar         = S.prior + t(Y)%*%Y + t(A.prior)%*%diag(1/diag(V.prior))%*%A.prior - t(A.bar)%*%V.bar.inv%*%A.bar
S.bar         = 0.5*(S.bar + t(S.bar))
S.bar.chol    = chol(S.bar)
S.bar.inv     = backsolve(S.bar.chol, forwardsolve(t(S.bar.chol), diag(N)))

# posterior distribution and predictive density draws 
# WARNING! This takes a while
############################################################
L                 = t(solve(chol(V.bar.inv)))
Y.h               = array(NA,c(h,2,S))
Y.h.m             = array(NA,c(h,2))
for (s in 1:S){
  draw.norm       = array(rnorm(prod(N*K)),c(K,N))
  Sigma.posterior = solve(rWishart(1, df=nu.bar, Sigma=S.bar.inv)[,,1])
  A.posterior.draw= A.bar + L%*%draw.norm%*%chol(Sigma.posterior)
  if (p==1){
    x.Ti          = matrix(Y[(nrow(Y)-p+1):nrow(Y),],nrow=1)
    x.Ti.m        = x.Ti
  } else {
    x.Ti          = Y[(nrow(Y)-p+1):nrow(Y),]
    x.Ti          = x.Ti[p:1,]
    x.Ti.m        = x.Ti
  }
  for (i in 1:h){
    x.T           = c(1,as.vector(t(x.Ti)))
    x.T.m         = c(1,as.vector(t(x.Ti.m)))
    Y.f           = rmvnorm(1, mean = x.T%*%A.posterior.draw, sigma=Sigma.posterior)
    Y.f.m         = x.T.m%*%A.bar
    if (p==1){
      x.Ti        = Y.f
      x.Ti.m      = Y.f.m
    } else {
      x.Ti        = rbind(Y.f,x.Ti[1:(p-1),])
      x.Ti.m      = rbind(Y.f.m,x.Ti.m[1:(p-1),])
    }
    Y.h[i,,s]     = Y.f[1:2]
    Y.h.m[i,]     = Y.f.m[1:2]
  }
}

# save(Y.h,Y.h.m,file="forecasts-N2-k1.RData")
save(Y.h,Y.h.m,file="forecasts-N2-k002.RData")

# plots
############################################################
# plot rgdp growth forecasts
load("forecasts-N2-k002.RData")
Yh02        = Y.h
Yhm02       = Y.h.m
load("forecasts-N2-k1.RData")
Yh1         = Y.h
Yhm1        = Y.h.m

Y.recent    = ts(rbind(Y[92:117,],matrix(NA,h,2)),start=c(2009,1),frequency=4)
interval.fs02   = apply(Yh02,1:2,hdi,credMass=0.68)
interval.fs1    = apply(Yh1,1:2,hdi,credMass=0.68)
lims.gdp.s  = range(Y.recent[1:26,1],interval.fs02[,,1],interval.fs1[,,1])
lims.pi.s   = range(Y.recent[1:26,2],interval.fs02[,,2],interval.fs1[,,2])

median.forecasts.s02  = apply(Yh02,1:2,quantile,probs=0.5)
median.forecasts.s1   = apply(Yh1,1:2,quantile,probs=0.5)

pdf(file=paste("forecasts-N2.pdf",sep=""), width=15,height=6)
par(mfrow=c(1,2), mar=rep(2.2,4),cex.axis=1.5, cex=1.5)
plot(1:nrow(Y.recent),Y.recent[,1], type="l", ylim=lims.gdp.s, axes=FALSE, xlab="", ylab="", main=expression(rgdp[t]), lwd=2, col=mcxs1)
axis(1,c(seq(from=1, to=34, by=4),35),c("","2010","","2012","","2014","","2016","",""), col=mcxs1)
axis(2,c(lims.gdp.s[1],0,lims.gdp.s[2]),round(c(lims.gdp.s[1],0,lims.gdp.s[2]),2), col=mcxs1)
polygon(c(26:34,34:26), c(Y.recent[26,1],interval.fs1[1,,1],interval.fs1[2,h:1,1],Y.recent[26,1]), col=mcxs1.shade1,border=mcxs1.shade1)
polygon(c(26:34,34:26), c(Y.recent[26,1],interval.fs02[1,,1],interval.fs02[2,h:1,1],Y.recent[26,1]), col=mcxs2.shade1,border=mcxs2.shade1)
abline(h=0, col=mcxs1)
lines(26:34, c(Y.recent[26,1],median.forecasts.s02[,1]),lwd=2,col=mcxs2)
lines(26:34, c(Y.recent[26,1],median.forecasts.s1[,1]),lwd=2,col=mcxs1)
plot(1:nrow(Y.recent),Y.recent[,2], type="l", ylim=lims.pi.s, axes=FALSE, xlab="", ylab="", main=expression(cpi[t]), col=mcxs2, lwd=2)
axis(1,c(seq(from=1, to=34, by=4),35),c("","2010","","2012","","2014","","2016","",""), col=mcxs1)
axis(2,c(lims.pi.s[1],0,lims.pi.s[2]),round(c(lims.pi.s[1],0,lims.pi.s[2]),2), col=mcxs1)
polygon(c(26:34,34:26), c(Y.recent[26,2],interval.fs1[1,,2],interval.fs1[2,h:1,2],Y.recent[26,2]), col=mcxs1.shade1,border=mcxs1.shade1)
polygon(c(26:34,34:26), c(Y.recent[26,2],interval.fs02[1,,2],interval.fs02[2,h:1,2],Y.recent[26,2]), col=mcxs2.shade1,border=mcxs2.shade1)
abline(h=0, col=mcxs1)
lines(26:34, c(Y.recent[26,2],median.forecasts.s02[,2]),lwd=2,col=mcxs2)
lines(26:34, c(Y.recent[26,2],median.forecasts.s1[,2]),lwd=2,col=mcxs1)
dev.off()
