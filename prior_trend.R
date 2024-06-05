
library(MASS)
library(plot3D)


#Construct function for the 
log_cond_var_NC<-function(x,t,rho,sigma.omega.sq){
  (pi*sqrt(sigma.omega.sq*(1-rho^(2*t))/(1-rho^(2))))^(-1)*besselK(abs(x)/sqrt(sigma.omega.sq*(1-rho^(2*t))/(1-rho^(2))),0)
}


#Sampling step. Sample the priors required. 

iter=10000
sigmav.s= 0.5
sigmav.nu= 3
rho=rep(NA,iter)
sigma.omega.sq=rep(NA,iter)
sigma.v.sq=rep(NA,iter)
rho=runif(iter,-1,1)
sigma.omega.sq= unlist(sapply(1:iter, function(ii){runif(1,0,1-(rho[ii]^2))}))
sigma.v.sq=sigmav.s/rchisq(iter,sigmav.nu)

T=10
grid=seq(from=-3, to=3, by=0.1)
N=length(grid)
density=matrix(NA,N,iter)
log_sigma_2_nc=matrix(NA,N,T)
for (t in 1:T){
  for (g in 1:N){
    for (i in 1:iter){
      density[,i]=log_cond_var_NC(grid,t,rho[i],sigma.omega.sq[i])
    }
    log_sigma_2_nc[g,t]=mean(density[g,])
  }
}


log_sigma_2_c=matrix(NA,N,T)
rho_sq_t<-function(rho,time){
  rho^(2*time)
}
for (t in 1:T){
  for (g in 1:N){
    for (i in 1:iter){
      density[,i]=dnorm(grid,mean = 0,sd=sqrt(sigma.v.sq[i]*(1-rho_sq_t(rho[i],t))/(1-rho[i]^(2))))
    }
    log_sigma_2_c[g,t]=mean(density[g,])
  }
}
z= t(log_sigma_2_nc)
zlimabrar=2.5
x= seq(from=1,to=T,by=1)
z[z>zlimabrar]  = zlimabrar

z1              = t(log_sigma_2_c)
z1[z1>zlimabrar]= zlimabrar

phi             = 15.5
theta           = 80
# theta           = 106.5

mcxs1  = "#05386B"
mcxs2  = "#379683" #
mcxs3  = "#5CDB95" #
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"
purple = "#b02442"

pdf(file="prior_trend.pdf", height=7, width=7)
f4    = persp3D(x=x, y=grid, z=z, phi=phi, theta=theta, xlab="time", ylab="log-variances", zlab="density", shade=NA, border=NA, ticktype="detailed", nticks=3,cex.lab=1, col=NA,plot=FALSE, bty="n", zlim=c(0,zlimabrar))
perspbox (x=x, y=grid,z=z, bty="f", lwd.grid = 0.3, lwd.panel = 0.3, col.axis="black", phi=phi, theta=theta, xlab="time", ylab="trend", zlab="density", ticktype="detailed", nticks=3,cex.lab=1, col = NULL, plot = TRUE, zlim=c(0,zlimabrar))
f4.a1 = trans3d(x=x, y=0, z=0, pmat=f4)
lines(f4.a1, lwd=1)

for (i in 1:length(x)){
  f4.l = trans3d(x=x[i], y=grid, z=z[i,], pmat=f4)
  lines(f4.l, lwd=1, col=mcxs3)
}
for (i in 1:length(x)){
  f4.l2 = trans3d(x=x[i], y=grid, z=z1[i,], pmat=f4)
  lines(f4.l2, lwd=1, col=mcxs2)
}
dev.off()

system("say After everything I've done for you! That you didn't ask for")