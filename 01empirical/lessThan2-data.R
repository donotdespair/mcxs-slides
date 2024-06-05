

rm(list=ls())
############################################################
# data preparation
############################################################
load("data_medium.Rda")

# select the countries
countries     = c("USA","AUS","CAN","CHN","IDN","IND","JPN","NZL","PHL","POL")
N             = length(countries)
T             = nrow(data.medium[data.medium$Isocode=="USA",])

# coopy corresponding data
Ys            = array(NA,c(T,6,N), dimnames=list(
  1959:2015,
  names(data.medium)[c(3,5:9)],
  countries
))

for (i in 1:length(countries)){
  Ys[,,i]     = as.matrix(data.medium[data.medium$Isocode==countries[i],c(3,5:9)])
}
Y             = Ys[,c(1,4,6),]

# # plot the data (not logs)
# pdf("grph-GDP.pdf", height=7, width=12)
# plot(ts(Y[,2,],start=c(1959), frequency=1),main="",nc=2, yax.flip=TRUE, col="deepskyblue2",lwd=2)
# dev.off()
# pdf("grph-Intensity.pdf", height=7, width=12)
# plot(ts(Y[,3,],start=c(1959), frequency=1),main="Carbon intensity",nc=2, yax.flip=TRUE,col="maroon2",lwd=2)
# dev.off()

# form model matrices
# log-GDP series
gdp   = ts(log(Y[,2,]),start=c(1959), frequency=1)
F     = as.matrix(diff(gdp[,1]))
XF    = cbind(rep(1,T-1), c(rep(1,14),rep(0,T-1-14)))
G     = matrix((gdp[,2:N] - matrix(rep(gdp[,1],N-1),ncol=N-1))[2:T,], ncol=1)
XG    = matrix(0,nrow(G),N-1)
for (c in 1:(N-1)){
  XG[(1+(c-1)*(T-1)):((T-1)+(c-1)*(T-1)),c]  = (gdp[,2:N] - matrix(rep(gdp[,1],N-1),ncol=N-1))[1:(T-1),c]
}

# identify the Tech peak and select the sub-samples
tau           = matrix(NA,0,1)
Xtau          = matrix(NA,0,2+N)
tech          = vector("list",N)
tech.trend    = vector("list",N)
tech.loess    = matrix(NA,T,N)
which.max.t   = rep(NA,N)

for (c in 1:N){
  gdp.loess           = loess(Tech ~ Year, data=as.data.frame(Y[,,c]), span=0.25)
  tech.loess[,c]      = predict(gdp.loess)
  which.max.t[c]      = which.max(tech.loess[,c])
  tech[[c]]           = log(Y[which.max.t[c]:T,3,c])
  tech.trend[[c]]     = (which.max.t[c]+1):T - mean(1:T)
  tau                 = rbind(tau,as.matrix(tech[[c]][2:length(tech[[c]])]))
  Xdelta.tmp          = t(matrix(rep(diag(N)[c,], length(tech[[c]])-1),ncol=length(tech[[c]])-1))
  Xtau.c              = cbind(tech.trend[[c]],tech[[c]][1:(length(tech[[c]])-1)],Xdelta.tmp)
  Xtau                = rbind(Xtau,Xtau.c)
}

# plot intensity, loess and the cut out point
for (c in 1:N){
  maxx  = ts(rep(NA,T),start=c(1959), frequency=1)
  maxx[which.max.t[c]] = tech.loess[which.max.t[c],c]
  pdf(paste("grph-Intensity-",c,".pdf",sep=""), height=4, width=9)
  plot(ts(Y[,3,c],start=c(1959), frequency=1),main="",ylab=colnames(Y[,3,])[c],xlab="",nc=2, yax.flip=TRUE,col="maroon2",lwd=2)
  lines(ts(tech.loess[,c],start=c(1959), frequency=1),col="maroon4")
  lines(maxx,lwd=3,col="maroon4",type="p")
  dev.off()
}

pop.tmp       = as.data.frame(read.csv("population.csv"))[,2:11]
population    = ts(pop.tmp,start=c(2020),frequency=0.2)

save(population,F,XF,G,XG,tau,Xtau,which.max.t,tech,tech.trend,tech.loess,Y,Ys,N,T,countries,file="lessThan2-data.RData")
