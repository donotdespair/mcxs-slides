############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 9: Forecasting with Bayesian VARs
############################################################

rm(list=ls())

# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"

mcxs1.rgb   = col2rgb(mcxs1)
mcxs1.shade1= rgb(mcxs1.rgb[1],mcxs1.rgb[2],mcxs1.rgb[3], alpha=120, maxColorValue=255)
mcxs2.rgb   = col2rgb(mcxs2)
mcxs2.shade1= rgb(mcxs2.rgb[1],mcxs2.rgb[2],mcxs2.rgb[3], alpha=120, maxColorValue=255)
  
############################################################
# Forecasting a bivariate series of Australian:
# + real GDP
# + CPI
############################################################
library(mvtnorm)
library(plot3D)
library(MASS)
library(HDInterval)
set.seed(123456)

# setup
############################################################
N       = 2
p       = 4
S       = 50000
h       = 20

# upload data and create Y and X
############################################################
y       = read.csv("data.csv")
y       = ts(log(y[,c(3,2)]), start=c(1959,3), frequency=4, names=c("rdgp","prices"))

Y       = ts(y[5:nrow(y),], start=c(1960,3), frequency=4)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[5:nrow(y)-i,])
}

# time series plots
############################################################
rgdp.range    = range(y[,1])
p.range       = range(y[,2])

pdf(file="data.pdf", width=15,height=6)
par(mfrow=c(1,2), mar=rep(3,4),cex.axis=1.5)
plot(1:length(y[,1]),y[,1], type="l", ylim=rgdp.range, axes=FALSE, xlab="", ylab="", lwd=4, col=mcxs1)
axis(1,c(3,43,83,123,163,203,243),c("1960","1970","1980","1990","2000","2010","2020"), col=mcxs1)
axis(2,c(rgdp.range[1],mean(rgdp.range),rgdp.range[2]),c("","rgdp",""), col=mcxs1)
plot(1:length(y[,2]),y[,2], type="l", ylim=p.range, axes=FALSE, xlab="", ylab="", lwd=4, col=mcxs2)
axis(1,c(3,43,83,123,163,203,243),c("1960","1970","1980","1990","2000","2010","2020"), col=mcxs1)
axis(2,c(p.range[1],mean(p.range),p.range[2]),c("","prices",""), col=mcxs1)
dev.off()









# MLE
############################################################
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/T
round(A.hat,3)
round(Sigma.hat,3)
round(cov2cor(Sigma.hat),3)

# prior distribution
############################################################
kappa.1     = 0.02^2
kappa.2     = 100
A.prior     = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:3,] = diag(2)
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
L                 = t(chol(V.bar))
for (s in 1:S){
  A.posterior[,,s]= A.bar + L%*%A.posterior[,,s]%*%chol(Sigma.posterior[,,s])
}

round(apply(A.posterior,1:2,mean),3)

# report posterior means and sd of parameters
A.E         = apply(A.posterior,1:2,mean)
A.sd        = apply(A.posterior,1:2,sd)
Sigma.E     = apply(Sigma.posterior,1:2,mean)
Sigma.sd    = apply(Sigma.posterior,1:2,sd)

library(xtable)
xtable(rbind(format(t(A.E),digits=1, scientific=FALSE)[1,],
      format(t(A.sd),digits=1, scientific=FALSE)[1,],
      format(t(A.E),digits=1, scientific=FALSE)[2,],
      format(t(A.sd),digits=1, scientific=FALSE)[2,]))

xtable(rbind(format(Sigma.E,digits=1, scientific=FALSE)[1,],
             format(Sigma.sd,digits=1, scientific=FALSE)[1,],
             format(Sigma.E,digits=1, scientific=FALSE)[2,],
             format(Sigma.sd,digits=1, scientific=FALSE)[2,]))


# simulate draws from the predictive density
# WARNING! This takes a while
############################################################
Y.h         = array(NA,c(h,N,S))

for (s in 1:S){
  x.Ti        = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti        = x.Ti[4:1,]
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1, mean = x.T%*%A.posterior[,,s], sigma=Sigma.posterior[,,s])
    x.Ti        = rbind(Y.h[i,,s],x.Ti[1:3,])
  }
}




# joint predictive density 1-period ahead
############################################################
limits.rgdp   = range(Y.h[1,1,])
limits.infl   = range(Y.h[1,2,])

bands       = 100
xlime       = limits.rgdp
ylime       = limits.infl
predictive.one.kernel = kde2d(x=Y.h[1,1,],y=Y.h[1,2,], n =bands, lims=c(xlime,ylime))

marginal.rgdp = apply(predictive.one.kernel$z,1,sum)
marginal.rgdp = max(predictive.one.kernel$z)*marginal.rgdp/max(marginal.rgdp)

marginal.infl = apply(predictive.one.kernel$z,2,sum)
marginal.infl = max(predictive.one.kernel$z)*marginal.infl/max(marginal.infl)

pdf(file=paste("joint-predictive-1ph.pdf",sep=""), width=9,height=7)
f1    = persp3D(x=predictive.one.kernel$x, y=predictive.one.kernel$y, z=predictive.one.kernel$z, phi=25, theta=25, xlab="\nrgdp", ylab="\nprices", zlab="\npredictive density h=1", shade=0, border=NA, ticktype="detailed", nticks=2,cex.lab=1, col="white")
f1.l1 = trans3d(x=predictive.one.kernel$x, y=predictive.one.kernel$y[100], z=marginal.rgdp, pmat=f1)
lines(f1.l1, lwd=2, col=mcxs1)
f1.l2 = trans3d(x=predictive.one.kernel$x[1], y=predictive.one.kernel$y, z=marginal.infl, pmat=f1)
lines(f1.l2, lwd=2, col=mcxs2)
f1    = persp3D(x=predictive.one.kernel$x, y=predictive.one.kernel$y, z=predictive.one.kernel$z, phi=25, theta=35, xlab="\nrgdp", ylab="\nprices", zlab="\npredictive density h=1", shade=.5, border=NA, ticktype="detailed", nticks=2,cex.lab=1, col=mcxs2.shade1, add=TRUE)
dev.off()

cor(t(Y.h[1,,]))

# joint predictive density of rgdp 1 and 2 periods ahead
############################################################
limits.rgdp1   = range(Y.h[1,1,])
limits.rgdp2   = range(Y.h[2,1,])

bands       = 100
xlime       = limits.rgdp1
ylime       = limits.rgdp2
predictive.one.kernel = kde2d(x=Y.h[1,1,],y=Y.h[2,1,], n =bands, lims=c(xlime,ylime))

marginal.rgdp1 = apply(predictive.one.kernel$z,1,sum)
marginal.rgdp1 = max(predictive.one.kernel$z)*marginal.rgdp1/max(marginal.rgdp1)

marginal.rgdp2 = apply(predictive.one.kernel$z,2,sum)
marginal.rgdp2 = max(predictive.one.kernel$z)*marginal.rgdp2/max(marginal.rgdp2)

pdf(file=paste("joint-predictive-rgdp-1and2ph.pdf",sep=""), width=9,height=7)
f2    = persp3D(x=predictive.one.kernel$x, y=predictive.one.kernel$y, z=predictive.one.kernel$z, phi=25, theta=25, xlab="\nrgdp[t+1|t]", ylab="\nrgdp[t+2|t]", zlab="\npredictive density rgdp h=1:2", shade=0, border=NA, ticktype="detailed", nticks=2,cex.lab=1, col="white")
f2.l1 = trans3d(x=predictive.one.kernel$x, y=predictive.one.kernel$y[100], z=marginal.rgdp1, pmat=f2)
lines(f2.l1, lwd=2, col=mcxs1)
f2.l2 = trans3d(x=predictive.one.kernel$x[1], y=predictive.one.kernel$y, z=marginal.rgdp2, pmat=f2)
lines(f2.l2, lwd=2, col=mcxs1)
f2    = persp3D(x=predictive.one.kernel$x, y=predictive.one.kernel$y, z=predictive.one.kernel$z, shade=.5, border=NA, col=mcxs1.shade1, add=TRUE)
dev.off()

cor(t(Y.h[1:2,1,]))


# joint predictive density of prices 1 and 2 periods ahead
############################################################
limits.i1   = range(Y.h[1,2,])
limits.i2   = range(Y.h[2,2,])

bands       = 100
xlime       = limits.i1
ylime       = limits.i2
predictive.one.kernel = kde2d(x=Y.h[1,2,],y=Y.h[2,2,], n =bands, lims=c(xlime,ylime))

marginal.i1 = apply(predictive.one.kernel$z,1,sum)
marginal.i1 = max(predictive.one.kernel$z)*marginal.i1/max(marginal.i1)

marginal.i2 = apply(predictive.one.kernel$z,2,sum)
marginal.i2 = max(predictive.one.kernel$z)*marginal.i2/max(marginal.i2)

pdf(file=paste("joint-predictive-prices-1and2ph.pdf",sep=""), width=9,height=7)
f3    = persp3D(x=predictive.one.kernel$x, y=predictive.one.kernel$y, z=predictive.one.kernel$z, phi=25, theta=25, xlab="\nprices[t+1|t]", ylab="\nprices[t+2|t]", zlab="\npredictive density prices h=1:2", shade=0, border=NA, ticktype="detailed", nticks=2,cex.lab=1, col="white")
f3.l1 = trans3d(x=predictive.one.kernel$x, y=predictive.one.kernel$y[100], z=marginal.i1, pmat=f3)
lines(f3.l1, lwd=2, col=mcxs2)
f3.l2 = trans3d(x=predictive.one.kernel$x[1], y=predictive.one.kernel$y, z=marginal.i2, pmat=f3)
lines(f3.l2, lwd=2, col=mcxs2)
f3    = persp3D(x=predictive.one.kernel$x, y=predictive.one.kernel$y, z=predictive.one.kernel$z, shade=.5, border=NA, col=mcxs2.shade1, add=TRUE)
dev.off()

cor(t(Y.h[1:2,2,]))



# plot of rgdp forecasts
############################################################
limits.1    = range(Y.h[,1,])
point.f     = apply(Y.h[,1,],1,mean)
interval.f  = apply(Y.h[,1,],1,hdi,credMass=0.90)

x           = seq(from=limits.1[1], to=limits.1[2], length.out=100)
z           = matrix(NA,h,99)
for (i in 1:h){
  z[i,]     = hist(Y.h[i,1,], breaks=x, plot=FALSE)$density
}
x           = hist(Y.h[i,1,], breaks=x, plot=FALSE)$mids
yy          = 1:h
z           = t(z)

pdf(file=paste("predictive-rgdp-horizons.pdf",sep=""), width=9,height=7)
theta = 180
phi   = 15.5
f4    = persp3D(x=x, y=yy, z=z, phi=phi, theta=theta, xlab="\nrgdp[t+h|t]", ylab="h", zlab="\npredictive densities of rgdp", shade=NA, border=NA, ticktype="detailed", nticks=3,cex.lab=1, col=NA,plot=FALSE)
perspbox (x=x, y=yy, z=z, bty="f", col.axis="black", phi=phi, theta=theta, xlab="\nrgdp[t+h|t]", ylab="h", zlab="\npredictive densities of rgdp", ticktype="detailed", nticks=3,cex.lab=1, col = NULL, plot = TRUE)
polygon3D(x=c(interval.f[1,],interval.f[2,h:1]), y=c(1:h,h:1), z=rep(0,2*h), col = mcxs1.shade1, NAcol = "white", border = NA, add = TRUE, plot = TRUE)
for (i in 1:h){
  f4.l = trans3d(x=x, y=yy[i], z=z[,i], pmat=f4)
  lines(f4.l, lwd=0.5, col="black")
}
f4.l1 = trans3d(x=point.f, y=yy, z=0, pmat=f4)
lines(f4.l1, lwd=2, col=mcxs1)
dev.off()




# # for canvas banner
# jpeg(file="predictive-rgdp-horizons.jpg", width=900,height=600)
# theta = 180
# phi   = 15.5
# f4    = persp3D(x=x, y=yy, z=z, phi=phi, theta=theta, xlab="", ylab="", zlab="", shade=NA, border=NA, ticktype="simple", nticks=3,cex.lab=1, col=NA,plot=FALSE)
# perspbox (x=x, y=yy, z=z, bty="f", col.axis="black", phi=phi, theta=theta, xlab="", ylab="", zlab="", ticktype="simple", nticks=3,cex.lab=1, col = NULL, plot = TRUE)
# polygon3D(x=c(interval.f[1,],interval.f[2,h:1]), y=c(1:h,h:1), z=rep(0,2*h), col = "maroon1", NAcol = "white", border = NA, add = TRUE, plot = TRUE)
# for (i in 1:h){
#   f4.l = trans3d(x=x, y=yy[i], z=z[,i], pmat=f4)
#   lines(f4.l, lwd=2, col="black")
# }
# f4.l1 = trans3d(x=point.f, y=yy, z=0, pmat=f4)
# lines(f4.l1, lwd=2)
# dev.off()




# plot of prices forecasts
############################################################
limits.1    = range(Y.h[,2,])
point.f     = apply(Y.h[,2,],1,mean)
interval.f  = apply(Y.h[,2,],1,hdi,credMass=0.90)

x           = seq(from=limits.1[1], to=limits.1[2], length.out=100)
z           = matrix(NA,h,99)
for (i in 1:h){
  z[i,]     = hist(Y.h[i,2,], breaks=x, plot=FALSE)$density
}
x           = hist(Y.h[i,2,], breaks=x, plot=FALSE)$mids
yy          = 1:h
z           = t(z)

pdf(file=paste("predictive-prices-horizons.pdf",sep=""), width=9,height=7)
theta = 180
phi   = 15.5
f4    = persp3D(x=x, y=yy, z=z, phi=phi, theta=theta, xlab="\nprices[t+h|t]", ylab="h", zlab="\npredictive densities of prices", shade=NA, border=NA, ticktype="detailed", nticks=3,cex.lab=1, col=NA,plot=FALSE)
perspbox (x=x, y=yy, z=z, bty="f", col.axis="black", phi=phi, theta=theta, xlab="\nprices[t+h|t]", ylab="h", zlab="\npredictive densities of prices", ticktype="detailed", nticks=3,cex.lab=1, col = NULL, plot = TRUE)
polygon3D(x=c(interval.f[1,],interval.f[2,h:1]), y=c(1:h,h:1), z=rep(0,2*h), col = mcxs2.shade1, NAcol = "white", border = NA, add = TRUE, plot = TRUE)
for (i in 1:h){
  f4.l = trans3d(x=x, y=yy[i], z=z[,i], pmat=f4)
  lines(f4.l, lwd=.5, col="black")
}
f4.l1 = trans3d(x=point.f, y=yy, z=0, pmat=f4)
lines(f4.l1, lwd=2,col=mcxs2)
dev.off()




# 2D plots of forecasts
############################################################
rgdp.point.f    = apply(Y.h[,1,],1,mean)
rgdp.interval.f = apply(Y.h[,1,],1,hdi,credMass=0.90)
rgdp.range      = range(y[,1],rgdp.interval.f)

p.point.f       = apply(Y.h[,2,],1,mean)
p.interval.f    = apply(Y.h[,2,],1,hdi,credMass=0.90)
p.range         = range(y[,2],p.interval.f)

pdf(file="forecasts.pdf", width=15,height=6)
par(mfrow=c(1,2), mar=rep(3,4),cex.axis=1.5)
plot(1:(length(y[,1])+h),c(y[,1],rgdp.point.f), type="l", ylim=rgdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=mcxs1)
axis(1,c(3,43,83,123,163,203,243,nrow(y),nrow(y)+h),c("1960","1970","1980","1990","2000","2010","2020","",""), col=mcxs1)
axis(2,c(rgdp.range[1],mean(rgdp.range),rgdp.range[2]),c("","rgdp",""), col=mcxs1)
abline(v=246, col=mcxs1)
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[21:1]),
        c(y[246,1],rgdp.interval.f[1,],rgdp.interval.f[2,20:1],y[246,1]),
        col=mcxs1.shade1, border=mcxs1.shade1)

plot(1:(length(y[,2])+h),c(y[,2],p.point.f), type="l", ylim=p.range, axes=FALSE, xlab="", ylab="", lwd=2, col=mcxs2)
axis(1,c(3,43,83,123,163,203,243,nrow(y),nrow(y)+h),c("1960","1970","1980","1990","2000","2010","2020","",""), col=mcxs1)
axis(2,c(p.range[1],mean(p.range),p.range[2]),c("","prices",""), col=mcxs1)
abline(v=246, col=mcxs1)
polygon(c(length(y[,2]):(length(y[,2])+h),(length(y[,2]):(length(y[,2])+h))[21:1]),
        c(y[246,2],p.interval.f[1,],p.interval.f[2,20:1],y[246,2]),
        col=mcxs2.shade1, border=mcxs2.shade1)
dev.off()



