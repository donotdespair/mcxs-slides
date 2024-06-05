############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 19: Modeling trend inflation
############################################################
# Estimation of the UC-AR(p) model with hierarchical prior
# and time-varying drift parameter
############################################################

# Load data and source code
############################################################
source("L19-codesTVPdrift.R")
load("cpi_au.rda")

# Define colors
############################################################
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
mcxs3.rgb   = col2rgb(mcxs3)
mcxs3.shade1= rgb(mcxs3.rgb[1],mcxs3.rgb[2],mcxs3.rgb[3], alpha=120, maxColorValue=255)

# Format data
############################################################
T           = length(cpi)
Y           = as.matrix(cpi)
series_name = expression(cpi[t])
ts.name     = "cpi-TVPdrift"

estimates   = matrix(NA,0,8)

# A function to check stationarity of an AR model
############################################################
is.ar.stationary = function(alpha.draw){
  max(Mod(eigen(rbind(alpha.draw, cbind(diag(p - 1), rep(0, p - 1))))$values)) < 1
}

# Estimation of the UC-AR(4)
############################################################
p         = 4
kappa1    = 1

H         = diag(T)
sdiag(H, -1) = -1

H.inv     = solve(H)
HH        = t(H)%*%H
HH.inv    = solve(HH)
X.tau     = cbind(rep(1,T),diag(T)[,1])

priors    = list(
  alpha.m = matrix(0,p,1),
  alpha.v = kappa1*diag(p),
  tau0.m  = 0,
  tau0.v  = 10,
  mu0.m   = 0,
  mu0.v   = 1,
  sigma.nu= 3,
  s.s     = 1*var(diff(Y)),
  s.a     = 1,
  H       = H,
  H.inv   = H.inv,
  HH      = HH,
  HH.inv  = HH.inv,
  X.tau   = X.tau
)

epsilon.0 = matrix(rnorm(T,sd=0.01),T,1)
mu.0      = matrix(rnorm(T,sd=0.01),T,1)
starting.values = list(
  Y       = Y,
  tau     = matrix(Y,T,1)-epsilon.0,
  epsilon = epsilon.0,
  mu      = mu.0,
  alpha   = rep(0,p),
  H.alpha = diag(T),
  tau0    = Y[1,],
  mu0     = 0,
  sigma   = rep(var(Y),3),
  sigma.s = var(Y)
)
names(starting.values$sigma) = c("eta","e","m")

t0      = proc.time()
uc.1    = UC.Gibbs.sampler(S=10000, starting.values, priors)
t1      = proc.time()
(t1-t0)[3]
uc.2    = UC.Gibbs.sampler(S=10000, uc.1$last.draw, priors)
t2      = proc.time()
(t2-t1)[3]
save(uc.1,uc.2,priors,cpi,Y,file=paste0("results/",ts.name,"-uc-ar-",p,".RData"))
system("say do not despair")


qqq     = uc.2
tau.m   = apply(qqq$posterior$tau,1,mean)
mu.m    = apply(qqq$posterior$mu,1,mean)
epsilon.m = apply(qqq$posterior$epsilon,1,mean)
plot.ts(tau.m)
plot.ts(mu.m)
plot.ts(epsilon.m)
plot.ts(t(qqq$posterior$sigma))
plot.ts(t(qqq$posterior$alpha))
hist(apply(qqq$posterior$alpha,2,sum),breaks=200)
sum(apply(qqq$posterior$alpha,2,sum)<1)/10000
plot.ts(qqq$posterior$tau0)
plot.ts(qqq$posterior$mu0)

apply(qqq$posterior$sigma,1,mean)
apply(qqq$posterior$alpha,1,mean)
mean(qqq$posterior$tau0)
mean(qqq$posterior$mu0)


load(paste0("results/",ts.name,"-uc-ar-",p,".RData"))
qqq     = uc.2

tau.m   = apply(qqq$posterior$tau,1,mean)
mu.m    = apply(qqq$posterior$mu,1,mean)
epsilon.m = apply(qqq$posterior$epsilon,1,mean)

tau.hdi = apply(qqq$posterior$tau,1,HDInterval::hdi,credMass=0.68)
mu.hdi  = apply(qqq$posterior$mu,1,HDInterval::hdi,credMass=0.68)
epsilon.hdi = apply(qqq$posterior$epsilon,1,HDInterval::hdi,credMass=0.68)

pdf("results/cpi-uc-tvpdrift-tau.pdf", height=5,width=7)
plot.ts(tau.m, ylim=range(tau.hdi), bty="n", ylab="",xlab="", axes=FALSE, main="")
polygon(c(1:T,T:1),c(tau.hdi[1,],tau.hdi[2,T:1]), col=mcxs3.shade1, border=mcxs3.shade1)
lines(tau.m,lwd=2, col=mcxs2)
lines(as.vector(Y), col=mcxs3)
legend("bottomright", c(expression(cpi[t]),expression(tau[t])), lwd=c(1,2),col=c(mcxs2,mcxs1), bty="n")
axis(2, range(tau.hdi),round(range(tau.hdi),2))
axis(1,16+c(1,41,81,121,161),c("1980","1990","2000","2010","2020"))
dev.off()

pdf("results/cpi-uc-tvpdrift-mu.pdf", height=5,width=7)
plot.ts(mu.m, ylim=range(mu.hdi), bty="n", ylab="",xlab="", axes=FALSE, main="")
polygon(c(1:T,T:1),c(mu.hdi[1,],mu.hdi[2,T:1]), col=mcxs3.shade1, border=mcxs3.shade1)
lines(pi/100, col=mcxs1)
lines(c(NA,diff(Y)), col=mcxs1)
lines(mu.m, col=mcxs2, lwd=2)
abline(h=0)
legend("topright", c(expression(Delta*cpi[t]),expression(mu[t])), lwd=c(1,2),col=c(mcxs1,mcxs2), bty="n")
axis(2, range(mu.hdi),round(range(mu.hdi), 2))
axis(1,16+c(1,41,81,121,161),c("1980","1990","2000","2010","2020"))
dev.off()

pdf("results/cpi-uc-tvpdrift-epsilon.pdf", height=5,width=7)
plot.ts(epsilon.m, ylim=range(epsilon.hdi), bty="n", ylab="",xlab="", axes=FALSE, main="")
polygon(c(1:T,T:1),c(epsilon.hdi[1,],epsilon.hdi[2,T:1]), col=mcxs3.shade1, border=mcxs3.shade1)
lines(epsilon.m, col=mcxs2, lwd=2)
lines(as.vector(Y-tau.m), col=mcxs1)
abline(h=0)
legend("topright", c(expression(cpi[t]-tau[t]),expression(epsilon[t])), lwd=c(1,2),col=c(mcxs1,mcxs2), bty="n")
axis(2, range(epsilon.hdi),round(range(epsilon.hdi), 2))
axis(1,16+c(1,41,81,121,161),c("1980","1990","2000","2010","2020"))
dev.off()

ar5.stat  = mean(apply(qqq$posterior$alpha,2,is.ar.stationary))

round(c(
  apply(qqq$posterior$alpha,1,mean),
  mean(ar5.stat),
  mean(qqq$posterior$sigma[1,]/qqq$posterior$sigma[2,]),
  mean(qqq$posterior$tau0),
  mean(qqq$posterior$mu0)
), 3)


round(c(
  apply(qqq$posterior$alpha,1,sd), NA, NA,
  sd(qqq$posterior$tau0),
  sd(qqq$posterior$mu0)
), 3)
