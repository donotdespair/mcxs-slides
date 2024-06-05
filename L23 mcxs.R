############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 23: Auxiliary mixture
############################################################

library(mvtnorm)
library(SparseM)
library(Matrix)

# Gibbs sampler for a simple UC model using precision sampler
############################################################
rm(list=ls())
dlogchisq1    = function(x){
  (sqrt(2)/gamma(0.5))*exp(-0.5*x)
}

alpha.st      = c(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000)
sigma.st      = c(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342)
pi.st         = c(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115)

rlogchisq     = log(rnorm(100000)^2)
range.lcs     = range(rlogchisq)
seq.x         = seq(from=range.lcs[1], to=range.lcs[2], length.out=1000)
den.lcs       = density(rlogchisq)
ddnorm        = matrix(NA, 10, length(seq.x))


plot(den.lcs$x,den.lcs$y,type="n")
for (i in 1:10){
  ddnorm[i,]  =  pi.st[i]*dnorm(seq.x, mean=alpha.st[i], sd=sqrt(sigma.st[i]))
  lines(seq.x,ddnorm[i,])
}

pdf(file="aux-mix.pdf", width=9,height=7)
plot(seq.x,apply(ddnorm,2,sum), type="n", axes=FALSE, xlab="", ylab="",cex.lab=1.5,col="gray70",lwd=2)
axis(1,c(-29,0,3),c(-29,0,3))
axis(2,c(0,0.24),c(0,0.24))
abline(h=0)
for (i in 1:10){
  lines(seq.x,ddnorm[i,], col="skyblue1")
}
lines(seq.x,apply(ddnorm,2,sum),col="skyblue4")
lines(den.lcs$x,den.lcs$y,col="maroon4")
legend(-29, 0.23, c(expression(log(chi^2)),expression(Sigma[i]*pi[i]*N(mu[i]*sigma[i]^2)),expression(pi[i]*N(mu[i]*sigma[i]^2))), lwd=c(2,2,2), col=c("maroon4","skyblue4","skyblue1"), box.lwd = "n")
legend(-23, 0.23, c("- exact distribution", "- auxilliary mixture","- components"), box.lwd = "n")
dev.off()


pdf(file="inverse.pdf", width=9,height=7)
plot(0:11, c(0,cumsum(pi.st),1), type="n", axes=FALSE, xlab="m", ylab=expression(Pr(s[t]==m,y)))
axis(1,1:10,c(1,NA,expression(s[t.1]==3),NA,NA,expression(s[t.2]==6),rep(NA,3),10))
axis(2,c(0,0.11,0.67,1),c("0",expression(u[1]==.11),expression(u[2]==.67),"1"))
lines(0:11, c(0,cumsum(pi.st),1), type="s", lwd=3, col="maroon3")
lines(0:3,rep(.11,4),lty=2, col="skyblue4")
lines(rep(3,2),c(0,.11),lty=2, col="skyblue4")
lines(0:6,rep(.67,7),lty=2, col="skyblue4")
lines(rep(6,2),c(0,.67),lty=2, col="skyblue4")
dev.off()
