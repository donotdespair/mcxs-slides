
rm(list=ls())

set.seed(123456)
source("lessThan2-codes.R")
load("lessThan2-data.RData")

yel  = rgb(red=254, green=213, blue=0, maxColorValue = 255)
gre  = rgb(red=78, green=177, blue=173, maxColorValue = 255)
blu  = rgb(red=134, green=91, blue=184, maxColorValue = 255)

data          = list(
  F           = F,
  XF          = XF,
  G           = G,
  XG          = XG,
  tau         = tau,
  Xtau        = Xtau,
  tech        = tech,
  tech.trend  = tech.trend
)

T           = T-1
aux         = list(
  gamma     = as.vector(solve(t(XF)%*%XF)%*%t(XF)%*%F),
  phi       = rep(0.98,N-1),
  betatau   = c(0,.98,rep(0,1+N)),
  sigmaf2   = 0.5,
  sigmag2   = c(0.6, 0.1, 0.6, 2.5, 0.4, 2.1, 6.3, 1, 1.6),
  sigma2    = rep(0.1,N),
  phi.mu    = 0.1,
  phi.sig   = 1,
  delta.mu  = -0.2,
  delta.sig = 3,
  s         = 3,
  s.sig     = 0.1
)

t0          = proc.time()
qqq1        = lessThan2.Gibbs(S=1000, data, starting.values=aux)
t1          = proc.time()
(t1-t0)/60

t0          = proc.time()
qqq2        = lessThan2.Gibbs(S=20000, data, starting.values=qqq1$last.draw)
t1          = proc.time()
(t1-t0)/60
save(qqq2,data,aux,file="lessThan2-MCMC.RData")


load("lessThan2-MCMC.RData")
qqq = qqq2

# # trace plots
# plot.ts(t(qqq$posterior$gamma))
# plot.ts(t(qqq$posterior$phi))
# plot.ts(t(qqq$posterior$betatau[c(1:2,13),]))
# plot.ts(t(qqq$posterior$betatau[3:12,]))
# plot.ts(qqq$posterior$sigmaf2)
# plot.ts(t(qqq$posterior$sigmag2))
# plot.ts(t(qqq$posterior$sigma2))
# plot.ts(qqq$posterior$phi.mu)
# plot.ts(qqq$posterior$phi.sig)
# plot.ts(qqq$posterior$delta.mu)
# plot.ts(qqq$posterior$delta.sig)
# plot.ts(qqq$posterior$s)
# plot.ts(qqq$posterior$s.sig)
# 
# # report estimation results
# # frontier economy
# round(apply(qqq$posterior$gamma,1,mean),3)
# round(apply(qqq$posterior$gamma,1,sd),3)
# round(mean(sqrt(qqq$posterior$sigmaf2)),3)
# round(sd(sqrt(qqq$posterior$sigmaf2)),3)
# # 
# # # other economies
# round(apply(qqq$posterior$phi,1,mean),2)
# round(apply(qqq$posterior$phi,1,sd),3)
# round(apply(sqrt(qqq$posterior$sigmag2),1,mean),2)
# round(apply(sqrt(qqq$posterior$sigmag2),1,sd),3)
# round(mean(qqq$posterior$phi.mu),2)
# round(sd(qqq$posterior$phi.mu),2)
# round(mean(sqrt(qqq$posterior$phi.sig)),2)
# round(sd(sqrt(qqq$posterior$phi.sig)),2)
# round(mean(sqrt(qqq$posterior$s)),3)
# round(sd(sqrt(qqq$posterior$s)),3)
# # 
# # # carbon intensity
# round(apply(qqq$posterior$beta,1,mean),3)
# round(apply(qqq$posterior$beta,1,sd),2)
# round(apply(sqrt(qqq$posterior$sigma2),1,mean),3)
# round(apply(sqrt(qqq$posterior$sigma2),1,sd),3)
# round(mean(qqq$posterior$delta.mu),2)
# round(sd(qqq$posterior$delta.mu),2)
# round(mean(sqrt(qqq$posterior$delta.sig)),2)
# round(sd(sqrt(qqq$posterior$delta.sig)),2)
# round(mean(sqrt(qqq$posterior$s.sig)),3)
# round(sd(sqrt(qqq$posterior$s.sig)),3)


# forecasting
rm(list=ls())
load("lessThan2-data.RData")
load("lessThan2-MCMC.RData")
qqq         = qqq2
S           = dim(qqq$posterior$gamma)[2]
T           = nrow(data$F)
N           = length(aux$sigma2)
h           = 85

# predict F_t
gdp.forecast            = array(NA,c(N,h,S))
for (s in 1:S){
  gdp.forecast[1,,s]    = log(Ys[57,4,1])+  qqq$posterior$gamma[1,s]*(1:h) + sqrt(qqq$posterior$sigmaf2[s])*rnorm(h)
  for (c in 2:N){
    H.phi               = diag(h)
    H.phi[2:h,1:(h-1)]  = H.phi[2:h,1:(h-1)] -qqq$posterior$phi[c-1,s]*diag(h-1)
    H.phi.inv           = solve(H.phi)
    # gdp.forecast[c,,s]  = -((matrix(G,ncol=N-1))[56,c-1]*H.phi.inv[,1] + H.phi.inv%*%(sqrt(qqq$posterior$sigmag2[c-1,s])*rnorm(h))-gdp.forecast[1,,s])
    gdp.forecast[c,,s]  = (matrix(G,ncol=N-1))[56,c-1]*H.phi.inv[,1] + H.phi.inv%*%(sqrt(qqq$posterior$sigmag2[c-1,s])*rnorm(h))+gdp.forecast[1,,s]
  }
}

gdp           = ts(log(Y[,2,]),start=c(1959), frequency=1)
gdp.median    = t(apply(gdp.forecast,1:2,quantile,probs=0.5))
gdp.lower     = t(apply(gdp.forecast,1:2,quantile,probs=0.05))
gdp.lower68   = t(apply(gdp.forecast,1:2,quantile,probs=0.16))
gdp.upper     = t(apply(gdp.forecast,1:2,quantile,probs=0.95))
gdp.upper68   = t(apply(gdp.forecast,1:2,quantile,probs=0.84))
xs            = 1959:2100

for (c in 1:N){
  ylims     = range(gdp[,c],gdp.upper[,c],gdp.lower[,c])
  pdf(paste("gdp-",c,"-",colnames(gdp)[c],".pdf",sep=""),height=5,width=10)
  plot(xs,c(gdp[,c],gdp.median[,c]), type="l", ylim=ylims, axes=FALSE, xlab="time", ylab=paste("GDP per capita for",colnames(gdp)[c],"[$1000] (log-scale)"),col="deepskyblue4")
  polygon(c(2015:2100,2100:2015), c(c(gdp[57,c],gdp.lower[,c]),c(gdp.upper[85:1,c],gdp[57,c])), col="deepskyblue1",border="deepskyblue1")
  polygon(c(2015:2100,2100:2015), c(c(gdp[57,c],gdp.lower68[,c]),c(gdp.upper68[85:1,c],gdp[57,c])), col="deepskyblue3",border="deepskyblue3")
  lines(2015:2100,c(gdp[57,c],gdp.median[,c]),col="deepskyblue4")
  axis(1,seq(from=1960,to=2100,by=10),seq(from=1960,to=2100,by=10))
  axis(2,c(ylims[1],gdp[57,c],ylims[2]),c(round(exp(ylims[1])/1000,0),round(exp(gdp[57,c])/1000,0),round(exp(ylims[2])/1000,0)))
  dev.off()
}

# predict \tau_t 
cie.forecast            = array(NA,c(N,h,S))
for (s in 1:S){
  for (c in 1:N){
    H.beta              = diag(h)
    H.beta[2:h,1:(h-1)] = H.beta[2:h,1:(h-1)] -qqq$posterior$beta[2,s]*diag(h-1)
    H.beta.inv          = solve(H.beta)
    cie.forecast[c,,s]  = H.beta.inv[,1]*(Y[57,3,c]- qqq$posterior$beta[2+c,s]) + H.beta.inv%*%(qqq$posterior$beta[1,s]*(29:113)) + H.beta.inv%*%(rnorm(h) + qqq$posterior$betatau[13,s]*qqq$posterior$sigma2[c,s]* rnorm(h))
  }
}

cie           = ts(Y[,3,],start=c(1959), frequency=1)
cie.median    = t(apply(cie.forecast,1:2,quantile,probs=0.5))
cie.lower     = t(apply(cie.forecast,1:2,quantile,probs=0.05))
cie.lower68   = t(apply(cie.forecast,1:2,quantile,probs=0.16))
cie.upper     = t(apply(cie.forecast,1:2,quantile,probs=0.95))
cie.upper68   = t(apply(cie.forecast,1:2,quantile,probs=0.84))
xs            = 1959:2100

for (c in 1:N){
  ylims     = range(cie[,c],cie.upper[,c],cie.lower[,c])
  pdf(paste("cie-",c,"-",colnames(cie)[c],".pdf",sep=""),height=5,width=10)
  plot(xs,c(cie[,c],cie.median[,c]), type="l", ylim=ylims, axes=FALSE, xlab="time", ylab=paste("carbon intensity for",colnames(cie)[c],"(log-scale)"),col="maroon4")
  polygon(c(2015:2100,2100:2015), c(c(cie[57,c],cie.lower[,c]),c(cie.upper[85:1,c],cie[57,c])), col="maroon1",border="maroon1")
  polygon(c(2015:2100,2100:2015), c(c(cie[57,c],cie.lower68[,c]),c(cie.upper68[85:1,c],cie[57,c])), col="maroon3",border="maroon3")
  lines(2015:2100,c(cie[57,c],cie.median[,c]),col="maroon4")
  axis(1,seq(from=1960,to=2100,by=10),seq(from=1960,to=2100,by=10))
  axis(2,c(ylims[1],cie[57,c],ylims[2]),c(round(exp(ylims[1])/1000,2),round(exp(cie[57,c])/1000,2),round(exp(ylims[2])/1000,2)))
  dev.off()
}
save(cie.forecast,gdp.forecast,file="lessThan2-predictions.RData")


rm(list=ls())
load("lessThan2-data.RData")
load("lessThan2-predictions.RData")
N           = dim(cie.forecast)[1]
h           = dim(cie.forecast)[2]
S           = dim(cie.forecast)[3]

cie.forecast    = cie.forecast[,seq(from=5,to=85,by=5),]
gdp.forecast    = gdp.forecast[,seq(from=5,to=85,by=5),]
logpopulation   = log(population)

CO2.forecast    = cie.forecast + gdp.forecast - log(10000)
  #array(rep(t(logpopulation),S),c(N,17,S))

CO2.median    = t(apply(CO2.forecast,1:2,quantile,probs=0.5))
CO2.lower     = t(apply(CO2.forecast,1:2,quantile,probs=0.05))
CO2.lower68   = t(apply(CO2.forecast,1:2,quantile,probs=0.16))
CO2.upper     = t(apply(CO2.forecast,1:2,quantile,probs=0.95))
CO2.upper68   = t(apply(CO2.forecast,1:2,quantile,probs=0.84))

c=1
ylims         = range(log(Ys[,6,c]),CO2.lower[,c],CO2.upper[,c])
plot(1959:2100,c(log(Ys[,6,c]),rep(NA,h)),type="l",ylim=ylims)
lines(seq(from=2015,to=2100,by=5),c(log(Ys[57,6,c]),CO2.median[,c]),type="l")
lines(seq(from=2015,to=2100,by=5),c(log(Ys[57,6,c]),CO2.lower[,c]),type="l")
lines(seq(from=2015,to=2100,by=5),c(log(Ys[57,6,c]),CO2.lower68[,c]),type="l")
lines(seq(from=2015,to=2100,by=5),c(log(Ys[57,6,c]),CO2.upper[,c]),type="l")
lines(seq(from=2015,to=2100,by=5),c(log(Ys[57,6,c]),CO2.upper68[,c]),type="l")

log(Ys[57,5,1])
log(Ys[57,5,1])





# cie.median    = t(apply(cie.forecast,1:2,quantile,probs=0.5))
# cie.lower     = t(apply(cie.forecast,1:2,quantile,probs=0.05))
# cie.lower68   = t(apply(cie.forecast,1:2,quantile,probs=0.16))
# cie.upper     = t(apply(cie.forecast,1:2,quantile,probs=0.95))
# cie.upper68   = t(apply(cie.forecast,1:2,quantile,probs=0.84))
# 
# gdp.median    = t(apply(gdp.forecast,1:2,quantile,probs=0.5))
# gdp.lower     = t(apply(gdp.forecast,1:2,quantile,probs=0.05))
# gdp.lower68   = t(apply(gdp.forecast,1:2,quantile,probs=0.16))
# gdp.upper     = t(apply(gdp.forecast,1:2,quantile,probs=0.95))
# gdp.upper68   = t(apply(gdp.forecast,1:2,quantile,probs=0.84))

# CO2.median    = gdp.median + cie.median - logpopulation
# CO2.lower     = gdp.lower + cie.lower - logpopulation
# CO2.lower68   = gdp.lower68 + cie.lower68 - logpopulation
# CO2.upper     = gdp.upper + cie.upper - logpopulation
# CO2.upper68   = gdp.upper68 + cie.upper68 - logpopulation

# CO2.median    = gdp.median + cie.median - log(10000)
# CO2.lower     = gdp.lower + cie.lower - log(10000)
# CO2.lower68   = gdp.lower68 + cie.lower68 - log(10000)
# CO2.upper     = gdp.upper + cie.upper - log(10000)
# CO2.upper68   = gdp.upper68 + cie.upper68 - log(10000)
