############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 4: Numerical optimization and integration
############################################################

# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"


############################################################
# Empirical illustration for numerical optimization
############################################################

# Sample size T = 100
############################################################
set.seed(123456)

sigma2.true = 2                                     # the true parameter value

T           = 100                                   # set the sample size
X           = as.matrix(rnorm(T, sd=sqrt(10)))      # generate X matrix
Y           = X + rnorm(T, sd=sqrt(sigma2.true))    # generate Y matrix
E           = Y - X                                 # compute residuals for beta=1
EE          = as.numeric(t(E)%*%E)                  # residuals' sum of squares

k.max       = 10                                    # maximum number of iterations in our Newton-Raphson algorithm
e           = 1e-6                                  # required precision of estimation
  
sigma2.k    = 0.5
output      = matrix(NA,k.max,6)
colnames(output) = c("sigma2.k","l","|G|","|G|<e","|diff.sigma2.k|","|diff.sigma2.k|<e")
rownames(output) = 1:k.max

for (k in 1:k.max){
  l.sigma2    = -0.5*T*log(2*pi) - 0.5*T*log(sigma2.k) - 0.5*(1/sigma2.k)*EE
  G.sigma2    = -0.5*T*(1/sigma2.k) + 0.5*(1/sigma2.k^2)*EE
  H.sigma2    = 0.5*T*(1/sigma2.k^2) - (1/sigma2.k^3)*EE
  sigma2.k.old= sigma2.k
  sigma2.k    = sigma2.k - (1/H.sigma2)*G.sigma2
  output[k,]  = c(sigma2.k.old, l.sigma2, abs(G.sigma2), abs(G.sigma2)<e, abs(sigma2.k-sigma2.k.old), abs(sigma2.k-sigma2.k.old)<e)
}
output
sigma2.mle    = EE/T

# Graph
############################################################
log.l.sigma2  = function(sigma2){ -0.5*T*log(2*pi) - 0.5*T*log(sigma2) - 0.5*(1/sigma2)*EE}
from          = 0.45
to            = 3.5
grid          = seq(from=from, to=to, by=0.01)
likel         = log.l.sigma2(grid)
limits        = range(likel)


for (n in 1:9){
  par(col.main=mcxs1,
      col.lab=mcxs2)
  pdf(file=paste("./grphs/newton-raphson-",n,".pdf",sep=""), width=9,height=6)
  plot(grid, likel, type="l", xlim=c(from,to), ylim=c(limits[1],limits[2]+30), axes=FALSE, xlab="", ylab="", main=paste("n =", n), lwd=2,col=mcxs1)
  axis(1,c(output[,1],to),c(rep("",n-1),substitute(sigma^2[(n)], list(n=n)),rep("",k.max-n),expression(sigma^2)), col=mcxs2)
  axis(2,c(output[,2],limits[2]+30),c(rep("",n-1),round(output[n,2],3),rep("",k.max-n),expression(l(sigma^2[(n)]))), col=mcxs2)
  lines(c(from,output[n,1]),rep(output[n,2],2),lwd=2, lty=2, col=mcxs2)
  lines(rep(output[n,1],2),c(limits[1],output[n,2]),lwd=2, lty=2, col=mcxs2)
  abline(a=output[n,2] - output[n,3]*output[n,1], b=output[n,3], lwd=2,col=mcxs1)
  dev.off()
}
round(output[1:9,],3)
round(output[1:9,1],4)


############################################################
# Empirical illustration for numerical integration
# Gibbs sampler for a simple linear regression model
############################################################
rm(list=ls())

# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"

set.seed(123456)

sigma2.true = 2                                     # the true parameter value

T           = 100                                   # set the sample size
X           = as.matrix(rnorm(T, sd=sqrt(10)))      # generate X matrix
Y           = X + rnorm(T, sd=sqrt(sigma2.true))    # generate Y matrix
S1          = 100                                   # detemine the burn-in draws
S2          = 1000                                  # number of draws from the final simulation

prior       = list(
  beta.under        = 0,  
  sigma2.beta.under = 1,
  s.under           = 1,
  nu.under          = 1
)

posterior         = matrix(NA,S1+S2,2)
colnames(posterior) = c("beta","sigma2")
posterior[1,2]    = 0.5                             # starting value for sigma2

set.seed(85561)
for (s in 1:(S1+S2)){
  # draw beta from the normal full conditional posterior
  sigma2.beta.bar   = 1/((1/prior$sigma2.beta.under) + (1/posterior[s,2])*as.numeric(t(X)%*%X))
  beta.bar          = sigma2.beta.bar * ((1/prior$sigma2.beta.under)*prior$beta.under + (1/posterior[s,2])*as.numeric(t(X)%*%Y))
  beta.draw         = rnorm(1, mean=beta.bar, sd=sqrt(sigma2.beta.bar))
  posterior[s,1]    = beta.draw
  
  # draw sigma2 from the inverse gamma 2 full conditional posterior
  if (s!=(S1+S2)){
    s.bar             = prior$s.under + as.numeric(t(Y-posterior[s,1]*X)%*%(Y-posterior[s,1]*X))
    nu.bar            = prior$nu.under + T
    sigma2.draw.tmp   = rchisq(1, df=nu.bar)
    sigma2.draw       = s.bar/sigma2.draw.tmp
    posterior[s+1,2]  = sigma2.draw
  }
}
plot.ts(posterior)
head(posterior)

limits.beta           = range(posterior[,1])
limits.sigma2         = range(posterior[,2])
posterior.plot        = cbind(
  kronecker(posterior[,1],rep(1,2))[1:(2*(S1+S2)-1)],
  c(posterior[1,2],kronecker(posterior[2:(S1+S2),2],rep(1,2)))
  )
beta.mle              = as.numeric(solve(t(X)%*%X)%*%t(X)%*%Y)
sigma2.mle            = as.numeric(t(Y-beta.mle*X)%*%(Y-beta.mle*X))/T


# Gibbs sampler - just the bivariate distribution
############################################################
for (s in 2:241){
  par(col.main=mcxs1)
  pdf(file=paste("./grphs/Gibbs-",s-1,".pdf",sep=""), width=9,height=6)
  plot(posterior.plot[1:s,2], posterior.plot[1:s,1], type="l", xlim=limits.sigma2, ylim=limits.beta, axes=FALSE, xlab="", ylab="", main=paste("s =", floor(s/2)), lwd=2,col=mcxs1)
  axis(1,c(limits.sigma2[1],sigma2.mle,limits.sigma2[2]),c("",expression(sigma^2),""), col=mcxs2)
  axis(2,c(limits.beta[1],beta.mle,limits.beta[2]),c("",expression(beta),""), col=mcxs2)
  dev.off()
}

# # Gibbs sampler - the bivariate distribution and trace plots
# ############################################################
# S   = 500
# for (s in 1:S){
# pdf(file=paste("./grphs/Gibbs3-",s,".pdf",sep=""), width=9,height=6)
#   par(mfrow=c(3,1), mar=rep(3,4),cex.axis=1.5,cex.main=1.5)
#   plot(posterior.plot[1:(2*s),2], posterior.plot[1:(2*s),1], type="l", xlim=limits.sigma2, ylim=limits.beta, axes=FALSE, xlab="", ylab="", main=paste("s =", s), lwd=2,col="black")
#   axis(1,c(limits.sigma2[1],sigma2.mle,limits.sigma2[2]),c("",expression(sigma^2),""))
#   axis(2,c(limits.beta[1],beta.mle,limits.beta[2]),c("",expression(beta),""))
#   plot(1:S,c(posterior[1:s,1],rep(NA,S-s)), type="l", ylim=limits.beta, axes=FALSE, xlab="", ylab="")
#   axis(1,c(1,S/2,S),c("","",""))
#   axis(2,c(limits.beta[1],mean(limits.beta),limits.beta[2]),c("",expression(beta),""))
#   plot(1:S,c(posterior[1:s,2],rep(NA,S-s)), type="l", ylim=limits.sigma2, axes=FALSE, xlab="", ylab="")
#   axis(1,c(1,S/2,S),c("1","s",S))
#   axis(2,c(limits.sigma2[1],mean(limits.sigma2),limits.sigma2[2]),c("",expression(sigma^2),""))
# dev.off()
# }

# trace plots
############################################################
pdf(file=paste("./grphs/bs-trace.pdf",sep=""), width=9,height=6)
par(mfrow=c(2,1), mar=rep(3,4),cex.axis=1.5)
plot(1:S2,posterior[(S1+1):(S1+S2),1], type="l", ylim=limits.beta, axes=FALSE, xlab="", ylab="", col=mcxs1)
axis(1,c(1,S2/2,S2),c("","",""), col=mcxs2, col.lab=mcxs2)
axis(2,c(limits.beta[1],mean(limits.beta),limits.beta[2]),c("",expression(beta),""), col=mcxs2)
plot(1:S2,posterior[(S1+1):(S1+S2),2], type="l", ylim=limits.sigma2, axes=FALSE, xlab="", ylab="", col=mcxs1)
axis(1,c(1,S2/2,S2),c("1","s",S2), col=mcxs2)
axis(2,c(limits.sigma2[1],mean(limits.sigma2),limits.sigma2[2]),c("",expression(sigma^2),""), col=mcxs2)
dev.off()

# marginal posterior densities
############################################################
beta.den      = density(posterior[(S1+1):(S1+S2),1])
sigma2.den    = density(posterior[(S1+1):(S1+S2),2])

limits.beta           = range(posterior[(S1+1):(S1+S2),1])
limits.sigma2         = range(posterior[(S1+1):(S1+S2),2])

library(HDInterval) # a package to compute the highest density intervals
beta.hdi      = hdi(posterior[(S1+1):(S1+S2),1], credMass = 0.90)
sigma2.hdi    = hdi(posterior[(S1+1):(S1+S2),2], credMass = 0.90)

dig2 = function(x,alpha,beta,log=FALSE){ 
  # a function to compute the density of inverse gamma 2
  output = -lgamma(beta/2) + (beta/2)*log(alpha/2) - .5*(beta+2)*log(x) - alpha/(2*x)
  if (log==FALSE){output = exp(output)}
  return(output)
}

beta.prior.d  = dnorm(beta.den$x, mean=prior$beta.under, sd=sqrt(prior$sigma2.beta.under))
sigma2.prior.d= dig2(beta.den$x, alpha=prior$s.under, beta=prior$nu.under)


pdf(file=paste("./grphs/bs-marginal.pdf",sep=""), width=10,height=5)
par(mfrow=c(1,2), mar=c(4,3,1,1))
plot(beta.den, axes=FALSE, xlab="", ylab="",main="", col=mcxs1)
polygon(
  c(beta.hdi[1],beta.hdi[1],beta.den$x[(beta.den$x>beta.hdi[1])&(beta.den$x<beta.hdi[2])],beta.hdi[2],beta.hdi[2],beta.hdi[1]), 
  c(0,mean(beta.den$y[136:137]),beta.den$y[(beta.den$x>beta.hdi[1])&(beta.den$x<beta.hdi[2])],mean(beta.den$y[341:342]),0,0), 
  col=mcxs2,border=mcxs2)
lines(beta.den, col=mcxs1)
lines(beta.den$x,beta.prior.d,lwd=2,lty=2, col=mcxs1)
axis(1,c(limits.beta[1],beta.hdi[1],beta.mle,beta.hdi[2],limits.beta[2]),round(c(NA,beta.hdi[1],beta.mle,beta.hdi[2],NA),2), cex.axis=1, col=mcxs1)
axis(2,c(0,max(beta.den$y)),c("",""), col=mcxs1)
mtext(expression(beta), side = 1, line = 2.5, cex=1.5, col=mcxs1)
mtext("density", side = 2, line = 1, cex=1, col=mcxs1)
plot(sigma2.den, axes=FALSE, xlab="", ylab="",main="",xlim=limits.sigma2, col=mcxs1)
polygon(
  c(sigma2.hdi[1],sigma2.hdi[1],sigma2.den$x[(sigma2.den$x>sigma2.hdi[1])&(sigma2.den$x<sigma2.hdi[2])],sigma2.hdi[2],sigma2.hdi[2],sigma2.hdi[1]), 
  c(0,mean(sigma2.den$y[74:75]),sigma2.den$y[(sigma2.den$x>sigma2.hdi[1])&(sigma2.den$x<sigma2.hdi[2])],mean(sigma2.den$y[239:240]),0,0), 
  col=mcxs2,border=mcxs2)
lines(sigma2.den,lwd=2, col=mcxs1)
lines(sigma2.den$x,sigma2.prior.d,lwd=2,lty=2, col=mcxs1)
axis(1,c(limits.sigma2[1],sigma2.hdi[1],sigma2.mle,sigma2.hdi[2],limits.sigma2[2]),round(c(NA,sigma2.hdi[1],sigma2.mle,sigma2.hdi[2],NA),2), cex.axis=1, col=mcxs1)
axis(2,c(0,max(sigma2.den$y)),c("",""), col=mcxs1)
mtext(expression(sigma^2), side = 1, line = 2.5, cex=1.5, col=mcxs1)
legend(x=2.1, y=1.4, legend=c("posterior density","prior density","highiest density interval"), lwd=c(2,2,8), lty=c(1,2,1),col=c(mcxs1,mcxs1,mcxs2),bty="n")
dev.off()


# joint posterior densities
############################################################
beta.den      = density(posterior[(S1+1):(S1+S2),1])
sigma2.den    = density(posterior[(S1+1):(S1+S2),2])

limits.beta           = range(posterior[(S1+1):(S1+S2),1])
limits.sigma2         = range(posterior[(S1+1):(S1+S2),2])

library(MASS)
bands       = 100
xlime       = limits.sigma2
ylime       = limits.beta
posterior.kernel = kde2d(x=posterior[(S1+1):(S1+S2),2],y=posterior[(S1+1):(S1+S2),1], n =bands, lims=c(xlime,ylime))

pdf(file=paste("./grphs/bs-joint.pdf",sep=""), width=9,height=7)
persp(posterior.kernel, phi = 25, theta =35, xlab="s2", ylab="b", zlab="density", ticktype="simple",ylim=ylime, xlim=xlime, shade = .5, border=NA, col=mcxs2)
dev.off()


# library(FinTS)
# Acf(posterior[(S1+1):(S1+S2),1])
# Acf(posterior[(S1+1):(S1+S2),2])
