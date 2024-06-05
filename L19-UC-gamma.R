############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 19: Modeling trend inflation
############################################################
# Estimation of the UC-AR(p) model with hierarchical prior,
# uncorrelated shocks, and 
# gamma prior for sigma2_eta
############################################################

# Load data and source code
############################################################
source("L19-codes.R")
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
T           = length(pi)
Y           = as.matrix(pi)
series_name = expression(pi[t])
ts.name     = "pi-gamma"

estimates   = matrix(NA,0,8)

# A function to check stationarity of an AR model
############################################################
is.ar.stationary = function(alpha.draw){
  max(Mod(eigen(rbind(alpha.draw, cbind(diag(p - 1), rep(0, p - 1))))$values)) < 1
}

# Estimation of UC-AR(p) models for p=1,...,5
############################################################
p         = 4
set.seed(1234)
# for (p in 1:5) {
  kappa1    = 1
  kappa2    = 1
  
  H         = diag(T)
  sdiag(H, -1) = -1
  
  
  H.inv     = solve(H)
  HH        = t(H) %*% H
  HH.inv    = solve(HH)
  X.tau     = cbind(rep(1, T), diag(T)[, 1])
  
  priors    = list(
    alpha.m = matrix(0,p,1),
    alpha.v = kappa1 * diag(p),
    beta.m  = matrix(0,2,1),
    beta.v  = kappa2 * diag(2),
    sigma.nu= 3,
    s.s     = 1 * var(Y),
    s.a     = 1,
    H       = H,
    H.inv   = H.inv,
    HH      = HH,
    HH.inv  = HH.inv,
    X.tau   = X.tau
  )
  
  epsilon.0 = matrix(rnorm(T, sd = 0.01), T, 1)
  starting.values = list(
    Y       = Y,
    tau     = matrix(Y,T,1) - epsilon.0,
    epsilon = epsilon.0,
    alpha   = rep(0, p),
    H.alpha = diag(T),
    beta    = c(0, Y[1,]),
    sigma   = rep(var(Y), 2),
    sigma.s = var(Y)
  )
  
  t0      = proc.time()
  uc.1    = UC.AR.Gibbs.sampler.sigma(S=1000, starting.values, priors)
  t1      = proc.time()          
  (t1-t0)/60
  uc.2    = UC.AR.Gibbs.sampler.sigma(S=5000, uc.1$last.draw, priors)
  t2      = proc.time()          
  (t2-t1)/60
  save(uc.1,uc.2,starting.values,priors, file=paste0("results/",ts.name,"-uc-ar-",p,".RData"))
  
  tau.hdi     = apply(uc.2$posterior$tau,1,HDInterval::hdi,probs = 0.68)
  epsilon.hdi = apply(uc.2$posterior$epsilon,1,HDInterval::hdi,probs = 0.68)
  
  pdf(file=paste0("results/",ts.name,"-uc-ar",p,".pdf"), height=7, width=10)
  par(mfrow=c(2,1), mar=c(2,2,2,2))
  plot(1:T,apply(uc.2$posterior$tau,1,mean), ylim=range(tau.hdi,Y), type="n", main="",xlab="",ylab="", axes=FALSE)
  polygon(c(1:T,T:1),c(tau.hdi[1,],tau.hdi[2,T:1]), col=mcxs2.shade1, border=mcxs2.shade1)
  lines(1:T,Y,lwd=3,col=mcxs1)
  lines(1:T,apply(uc.2$posterior$tau,1,quantile,probs=0.5),lwd=2, col=mcxs2)
  axis(1,16+c(1,41,81,121,161),c("1980","1990","2000","2010","2020"))
  axis(2,round(range(tau.hdi,Y),2))
  abline(h=0)
  legend(140,max(tau.hdi,Y),legend=c(series_name, expression(tau[t])),bty="n", col=c(mcxs1,mcxs2), text.col=c(mcxs1,mcxs2),lwd=c(3,2))
  plot(1:T,apply(uc.2$posterior$epsilon,1,mean), ylim=range(epsilon.hdi), type="n", main="",xlab="",ylab="", axes=FALSE)
  polygon(c(1:T,T:1),c(epsilon.hdi[1,],epsilon.hdi[2,T:1]), col=mcxs3.shade1, border=mcxs3.shade1)
  lines(1:T,apply(uc.2$posterior$epsilon,1,quantile,probs=0.5),lwd=2, col=mcxs2)
  abline(h=0)
  axis(1,16+c(1,41,81,121,161),c("1980","1990","2000","2010","2020"))
  axis(2,round(range(epsilon.hdi),2))
  legend(140,max(epsilon.hdi),legend=expression(epsilon[t]),bty="n", col=mcxs2, text.col=mcxs2,lwd=c(2))
  dev.off()
  
  ar5.a.m   = apply(uc.2$posterior$alpha,1,mean)
  ar5.a.sd  = apply(uc.2$posterior$alpha,1,sd)
  ar5.stat  = mean(apply(uc.2$posterior$alpha,2,is.ar.stationary))
  ar5.b.m   = apply(uc.2$posterior$beta,1,mean)
  ar5.b.sd  = apply(uc.2$posterior$beta,1,sd)
  ar5.ntsr  = mean(uc.2$posterior$sigma[1,]/uc.2$posterior$sigma[2,])
  
  estimates     = rbind(estimates,
                        rbind(
                          c(ar5.a.m,rep(NA,5-p),ar5.stat,ar5.b.m[1],ar5.ntsr),
                          c(ar5.a.sd,rep(NA,6-p),ar5.b.sd[1],rep(NA,1))
                        )
  )
# }
system("say This is the end")

estimates
xtable::xtable(estimates, digits=c(1,rep(3,8)))


omega_sign  = rep(1, 5000)
omega_sign[runif(5000)<0.5] = -1
omega = omega_sign * sqrt(uc.2$posterior$sigma[1,])
xx = seq(from = range(omega)[1], to = range(omega)[2], length.out = 200)
yy = rep(NA, length(xx))
for (s in 1:length(xx)) {
    yy[s] = mean(dnorm(xx[s], sd = sqrt(uc.2$posterior$sigma.s)))
}

pdf(file=paste0("results/",ts.name,"-omega.pdf"), height=5, width=7)
hist(omega, breaks=100, col = mcxs1.shade1, border = mcxs1.shade1, xlab = expression(+-sqrt(sigma[eta]^2)), main = "", freq=FALSE)
lines(x = xx, y = yy, lwd = 2, col = mcxs1)
dev.off()

