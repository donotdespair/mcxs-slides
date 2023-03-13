############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 6: Understanding Unit Rooters
############################################################

# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"

############################################################
# Reproduction of the results from Sims, Uhlig (1991)
# Understanding Unit Rooters: A Helicopter Tour
############################################################
rm(list=ls())

library(plot3D)

set.seed(123456)                                # for the reproduction of the results

ar1.ols    = function(x){
  # a function to estimate the autoregressive coefficient of the AR(1) model 
  as.numeric(solve(t(x[1:(length(x)-1)])%*%(x[1:(length(x)-1)]))%*%t(x[1:(length(x)-1)])%*%(x[2:length(x)]))
}

S           = 50000                             # simulation size
T           = 100                               # dimention of the simulated data
sigma2      = 1                                 # fix the variance of the error term in the AR(1) process with zero constant term and autregressive coefficient rho

rho.grid    = seq(from=0.8, to=1.1, by=0.01)    # the grid of values of rho considered by Sims & Uhlig (1991)
R           = length(rho.grid)
wn          = matrix(rnorm(S*T, sd=sigma2),T,S)            # generate the S white noise processes
Y           = array(NA,c(T,S,R))
rho.bins    = c(-Inf,seq(from=0.795, to=1.105, by=0.01),Inf)    # the bounds for the bins of the histogram of rho.ols

for (r in 1:length(rho.grid)){
  H               = diag(T)
  H[2:T,1:(T-1)]  = H[2:T,1:(T-1)] - rho.grid[r]*diag(T-1)
  Y[,,r]          = solve(H)%*%wn               # create the data according to the DGP: y_t = rho_r * y_{t-1} + wn_t
}

rho.ols     = apply(Y,2:3,ar1.ols)               # ols estimation

c.rho       = matrix(NA,R,R)
for (r in 1:R){
  hh        = hist(rho.ols[,r],breaks=rho.bins, plot=FALSE)
  c.rho[r,] = hh$counts[2:32]
}
scaling     = 20
p.rho       = scaling*c.rho/sum(c.rho)


pdf(file="f1.pdf", width=9,height=7)
f1    = persp3D(x=rho.grid[1:21], y=rho.grid[1:21], z=p.rho[1:21,1:21], phi=20, theta=125, zlim=c(0,p.rho[21,21]), xlab="\nrho", ylab="rho.ml", zlab="\n20 x density", col=mcxs4, shade=.01, border="black", ticktype="detailed", nticks=3,cex.lab=1.5, scale=FALSE, bty="f", lwd=.5)
f1.l1 = trans3d(x=rho.grid[1:21], y=rho.grid[1:21], z=diag(p.rho[1:22,1:21]), pmat=f1)
lines(f1.l1, lwd=.5)
f1.l2 = trans3d(x=rho.grid[1:21], y=rho.grid[21], z=p.rho[1:21,21], pmat=f1)
lines(f1.l2, lwd=4, col=mcxs1)
f1.l3 = trans3d(x=rho.grid[21], y=rho.grid[1:21], z=p.rho[21,1:21], pmat=f1)
lines(f1.l3, lwd=4, col=mcxs2)
dev.off()

pdf(file="f2.pdf", width=9,height=7)
f2    = persp3D(x=rho.grid[1:21], y=rho.grid[1:21], z=p.rho[1:21,1:21], phi=15.5, theta=135, zlim=c(0,p.rho[21,21]), xlab="\nrho", ylab="rho.ml", zlab="\n20 x density", col=mcxs4, shade=.01, border="black", ticktype="detailed", nticks=3,cex.lab=1.5, scale=FALSE, bty="f", lwd=.5)
f2.l1 = trans3d(x=rho.grid[1:21], y=rho.grid[1:21], z=diag(p.rho[1:21,1:21]), pmat=f2)
lines(f2.l1, lwd=.5)
f2.l2 = trans3d(x=rho.grid[1:21], y=rho.grid[21], z=p.rho[1:21,21], pmat=f2)
lines(f2.l2, lwd=4, col=mcxs1)
f2.l3 = trans3d(x=rho.grid[21], y=rho.grid[1:21], z=p.rho[21,1:21], pmat=f2)
lines(f2.l3, lwd=4, col=mcxs2)
dev.off()

pdf(file="f3.pdf", width=9,height=7)
f3    = persp3D(x=rho.grid[1:21], y=rho.grid[1:31], z=p.rho[1:21,1:31], phi=12, theta=110, zlim=c(0,p.rho[21,21]), xlab="\nrho", ylab="rho.ml", zlab="\n20 x density", shade=.01, border="black", ticktype="detailed", nticks=3,cex.lab=1.5, col=mcxs4, scale=FALSE, bty="f")
f3.l1 = trans3d(x=rho.grid[1:21], y=rho.grid[1:21], z=diag(p.rho[1:21,1:21]), pmat=f3)
lines(f3.l1, lwd=1)
f3.l2 = trans3d(x=rho.grid[21], y=rho.grid[1:31], z=p.rho[21,1:31], pmat=f3)
lines(f3.l2, lwd=4, col=mcxs2)
f3.l3 = trans3d(x=rho.grid[21], y=rho.grid[21], z=c(0,p.rho[21,21]), pmat=f3)
lines(f3.l3, lwd=2)
f3.l4 = trans3d(x=rho.grid[1:21], y=rho.grid[21], z=rep(p.rho[21,21],2), pmat=f3)
lines(f3.l4, lwd=2)
f3.l5 = trans3d(x=rho.grid[1], y=rho.grid[21], z=c(0,p.rho[21,21]), pmat=f3)
lines(f3.l5, lwd=2)
f3.l6 = trans3d(x=rho.grid[1:21], y=rho.grid[21], z=rep(0,2), pmat=f3)
lines(f3.l6, lwd=2)
dev.off()


# for canvas banner
############################################################
library(plot3D)
jpeg(file="f3.jpg", width=900,height=600)
f3    = persp3D(x=rho.grid[1:21], y=rho.grid[1:31], z=p.rho[1:21,1:31], phi=12, theta=110, zlim=c(0,p.rho[21,21]), xlab="", ylab="", zlab="", shade=.01, border="black", ticktype="simple", nticks=3,cex.lab=1.5, col="white", scale=FALSE)
f3.l1 = trans3d(x=rho.grid[1:21], y=rho.grid[1:21], z=diag(p.rho[1:21,1:21]), pmat=f3)
lines(f3.l1, lwd=1)
f3.l2 = trans3d(x=rho.grid[21], y=rho.grid[1:31], z=p.rho[21,1:31], pmat=f3)
lines(f3.l2, lwd=8, col="maroon1")
f3.l3 = trans3d(x=rho.grid[21], y=rho.grid[21], z=c(0,p.rho[21,21]), pmat=f3)
lines(f3.l3, lwd=2)
f3.l4 = trans3d(x=rho.grid[1:21], y=rho.grid[21], z=rep(p.rho[21,21],2), pmat=f3)
lines(f3.l4, lwd=2)
f3.l5 = trans3d(x=rho.grid[1], y=rho.grid[21], z=c(0,p.rho[21,21]), pmat=f3)
lines(f3.l5, lwd=2)
f3.l6 = trans3d(x=rho.grid[1:21], y=rho.grid[21], z=rep(0,2), pmat=f3)
lines(f3.l6, lwd=2)
f3.l7 = trans3d(x=rho.grid[21], y=rho.grid[31], z=c(0,p.rho[21,21]), pmat=f3)
lines(f3.l7, lwd=2)
f3.l8 = trans3d(x=rho.grid[1:21], y=rho.grid[31], z=rep(p.rho[21,21],2), pmat=f3)
lines(f3.l8, lwd=2)
f3.l9 = trans3d(x=rho.grid[21], y=rho.grid[1:31], z=rep(p.rho[21,21],2), pmat=f3)
lines(f3.l9, lwd=2)
dev.off()






pdf(file="f4.pdf", width=9,height=7)
f4    = persp3D(x=rho.grid[1:31], y=rho.grid[1:21], z=p.rho[1:31,1:21], phi=12, theta=160, zlim=c(0,p.rho[21,21]), xlab="\nrho", ylab="rho.ml", zlab="\n20 x density", shade=.01, border="black", ticktype="detailed", nticks=3,cex.lab=1.5, col=mcxs4, scale=FALSE, bty="f", lwd=.5)
f4.l1 = trans3d(x=rho.grid[1:21], y=rho.grid[1:21], z=diag(p.rho[1:21,1:21]), pmat=f4)
lines(f4.l1, lwd=.5)
f4.l2 = trans3d(x=rho.grid[1:31], y=rho.grid[21], z=p.rho[1:31,21], pmat=f4)
lines(f4.l2, lwd=4, col=mcxs1)
f4.l3 = trans3d(x=rho.grid[21], y=rho.grid[21], z=c(0,p.rho[21,21]), pmat=f4)
lines(f4.l3, lwd=2)
f4.l6 = trans3d(x=rho.grid[21], y=rho.grid[1:21], z=rep(0,2), pmat=f4)
lines(f4.l6, lwd=2)
f4.l7 = trans3d(x=rho.grid[21], y=rho.grid[1], z=c(0,p.rho[21,21]), pmat=f4)
lines(f4.l7, lwd=2)
f4.l8 = trans3d(x=rho.grid[21], y=rho.grid[1:21], z=rep(p.rho[21,21],2), pmat=f4)
lines(f4.l8, lwd=2)
dev.off()







pdf(file="f5.pdf", width=9,height=7)
f5    = persp3D(x=rho.grid[1:31], y=rho.grid[1:16], z=p.rho[1:31,1:16], phi=12, theta=160, zlim=c(0,p.rho[21,21]), xlab="\nrho", ylab="rho.ml", zlab="\n20 x density", shade=.01, border="black", ticktype="detailed", nticks=3,cex.lab=1.5, col=mcxs4, scale=FALSE, bty="f")
f5.l1 = trans3d(x=rho.grid[1:16], y=rho.grid[1:16], z=diag(p.rho[1:16,1:16]), pmat=f5)
lines(f5.l1, lwd=1)
f5.l2 = trans3d(x=rho.grid[1:31], y=rho.grid[16], z=p.rho[1:31,16], pmat=f5)
lines(f5.l2, lwd=4, col=mcxs1)
f5.l3 = trans3d(x=rho.grid[21], y=rho.grid[16], z=c(0,p.rho[21,21]), pmat=f5)
lines(f5.l3, lwd=2)
f5.l6 = trans3d(x=rho.grid[21], y=rho.grid[1:16], z=rep(0,2), pmat=f5)
lines(f5.l6, lwd=2)
f5.l7 = trans3d(x=rho.grid[21], y=rho.grid[1], z=c(0,p.rho[21,21]), pmat=f5)
lines(f5.l7, lwd=2)
f5.l8 = trans3d(x=rho.grid[21], y=rho.grid[1:16], z=rep(p.rho[21,21],2), pmat=f5)
lines(f5.l8, lwd=2)
dev.off()



pdf(file="f6.pdf", width=9,height=7)
plot(x=rho.grid[1:31], y=p.rho[1:31,21]/sum(p.rho[1:31,21]), type="l", main="", xlab=expression(rho),ylab="density", axes=FALSE, lwd=2, col="gray50")
lines(x=rho.grid[1:31], y=p.rho[1:31,21]/sum(p.rho[1:31,21]), lwd=4, col=mcxs1)
lines(x=rho.grid[1:31], y=p.rho[21,1:31]/sum(p.rho[21,1:31]), lwd=4, col=mcxs2)
axis(1, c(rho.grid[1], mean(rho.grid[12:13]), rho.grid[21],rho.grid[31]), c("","","1",""))
axis(2, c(0,p.rho[21,21]/sum(p.rho[1:31,21])), c("",""))
dev.off()



mcxs1.rgb   = col2rgb(mcxs1)
mcxs1.sh1   = rgb(mcxs1.rgb[1],mcxs1.rgb[2],mcxs1.rgb[3], alpha=50, maxColorValue = 255)
mcxs1.sh2   = rgb(mcxs1.rgb[1],mcxs1.rgb[2],mcxs1.rgb[3], alpha=100, maxColorValue = 255)
mcxs1.sh3   = rgb(mcxs1.rgb[1],mcxs1.rgb[2],mcxs1.rgb[3], alpha=200, maxColorValue = 255)
x.super   = 24
m.point   = max(p.rho[1:31,x.super]/sum(p.rho[1:31,x.super]))

pdf(file="f7.pdf", width=9,height=7)
plot(x=rho.grid[1:31], y=p.rho[1:31,21]/sum(p.rho[1:31,21]), type="l", main="", xlab=expression(rho),ylab="density", axes=FALSE, lwd=2, col="gray50", ylim=c(0,m.point))
lines(x=rho.grid[1:31], y=p.rho[1:31,21]/sum(p.rho[1:31,21]), lwd=4, col=mcxs1)
lines(x=rho.grid[1:31], y=p.rho[1:31,x.super]/sum(p.rho[1:31,x.super]), lwd=4, col=mcxs1.sh3)
lines(x=rho.grid[1:31], y=p.rho[1:31,16]/sum(p.rho[1:31,16]), lwd=4, col=mcxs1.sh1)
lines(x=rho.grid[1:31], y=p.rho[1:31,11]/sum(p.rho[1:31,11]), lwd=4, col=mcxs1.sh2)
axis(1, c(rho.grid[1], mean(rho.grid[12:13]), rho.grid[21],rho.grid[31]), c("","","1",""))
axis(2, c(0,m.point), c("",""))
text(x=rho.grid[21], y=m.point*0.5, substitute(p(rho/hat(rho)[ML]==r),list(r=rho.grid[21])), col=mcxs1)
text(x=rho.grid[x.super], y=m.point*0.95, substitute(p(rho/hat(rho)[ML]==r),list(r=rho.grid[x.super])), col=mcxs1.sh3)
text(x=rho.grid[16], y=m.point*0.25, substitute(p(rho/hat(rho)[ML]==r),list(r=rho.grid[16])), col=mcxs1.sh1)
text(x=rho.grid[11], y=m.point*0.2, substitute(p(rho/hat(rho)[ML]==r),list(r=rho.grid[11])), col=mcxs1.sh2)
dev.off()









p.rho       = scaling*c.rho/sum(c.rho)

p.rho = p.rho/1.5

pdf(file="behind.pdf", width=9,height=7)
b1    = persp3D(theta=120, x=rho.grid[1:31], y=rho.grid[1:31], z=p.rho[1:31,1:31], phi=12, zlim=c(0,p.rho[24,24]), xlab="\nrho", ylab="rho.ml", zlab="\n 13.3(3) x density", shade=.01, border="black", ticktype="detailed", nticks=3,cex.lab=1.5, col=mcxs4, scale=FALSE, bty="f")
b1.l1 = trans3d(x=rho.grid[1:31], y=rho.grid[1:31], z=diag(p.rho[1:31,1:31]), pmat=b1)
lines(b1.l1, lwd=1)
b1.l3 = trans3d(x=rho.grid[21], y=rho.grid[21], z=c(0,p.rho[21,21]), pmat=b1)
lines(b1.l3, lwd=2)
b1.l4 = trans3d(x=rho.grid[1:21], y=rho.grid[21], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l4, lwd=2)
b1.l4a = trans3d(x=rho.grid[21:31], y=rho.grid[21], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l4a, lwd=1)
b1.l5 = trans3d(x=rho.grid[1], y=rho.grid[21], z=c(0,p.rho[21,21]), pmat=b1)
lines(b1.l5, lwd=2)
b1.l6 = trans3d(x=rho.grid[1:21], y=rho.grid[21], z=rep(0,2), pmat=b1)
lines(b1.l6, lwd=2)
b1.l6a = trans3d(x=rho.grid[21:31], y=rho.grid[21], z=rep(0,2), pmat=b1)
lines(b1.l6a, lwd=1)
b1.l3a = trans3d(x=rho.grid[21], y=rho.grid[1], z=c(0,p.rho[21,21]), pmat=b1)
lines(b1.l3a, lwd=2)
b1.l4a = trans3d(x=rho.grid[1:21], y=rho.grid[1], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l4a, lwd=2)
b1.l4b = trans3d(x=rho.grid[21:31], y=rho.grid[1], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l4b, lwd=1)
b1.l5a = trans3d(x=rho.grid[1], y=rho.grid[1], z=c(0,p.rho[21,21]), pmat=b1)
lines(b1.l5a, lwd=2)
b1.l6a = trans3d(x=rho.grid[1:21], y=rho.grid[1], z=rep(0,2), pmat=b1)
lines(b1.l6a, lwd=2)
b1.l6b = trans3d(x=rho.grid[21:31], y=rho.grid[1], z=rep(0,2), pmat=b1)
lines(b1.l6b, lwd=1)
b1.l9 = trans3d(x=rho.grid[21], y=rho.grid[1:21], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l9, lwd=2)
b1.l9e = trans3d(x=rho.grid[21], y=rho.grid[1:21], z=rep(0,2), pmat=b1)
lines(b1.l9e, lwd=2)
b1.l9f = trans3d(x=rho.grid[21], y=rho.grid[21:31], z=rep(0,2), pmat=b1)
lines(b1.l9f, lwd=1)
b1.l9b = trans3d(x=rho.grid[1], y=rho.grid[1:21], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l9b, lwd=2)
b1.l9d = trans3d(x=rho.grid[31], y=rho.grid[1:21], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l9d, lwd=1)
b1.l5a = trans3d(x=rho.grid[31], y=rho.grid[21], z=c(0,p.rho[21,21]), pmat=b1)
lines(b1.l5a, lwd=1)
b1.l9a = trans3d(x=rho.grid[21], y=rho.grid[21:31], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l9a, lwd=1)
b1.l9c = trans3d(x=rho.grid[1], y=rho.grid[21:31], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l9c, lwd=1)
b1.l7 = trans3d(x=rho.grid[21], y=rho.grid[31], z=c(0,p.rho[21,21]), pmat=b1)
lines(b1.l7, lwd=1)
b1.l8 = trans3d(x=rho.grid[1:21], y=rho.grid[31], z=rep(p.rho[21,21],2), pmat=b1)
lines(b1.l8, lwd=1)
b1.l2 = trans3d(x=rho.grid[1:31], y=rho.grid[16], z=p.rho[1:31,16], pmat=b1)
lines(b1.l2, lwd=3, col=mcxs1)
b1.l2 = trans3d(x=rho.grid[1:31], y=rho.grid[21], z=p.rho[1:31,21], pmat=b1)
lines(b1.l2, lwd=3, col=mcxs1)
b1.l2 = trans3d(x=rho.grid[21], y=rho.grid[1:31], z=p.rho[21,1:31], pmat=b1)
lines(b1.l2, lwd=3, col=mcxs2)
dev.off()




pdf(file="behind-nopdf.pdf", width=9,height=7)
b1    = persp3D(theta=120, x=rho.grid[1:31], y=rho.grid[1:31], z=p.rho[1:31,1:31], phi=12, zlim=c(0,p.rho[24,24]), xlab="\nrho", ylab="rho.ml", zlab="\n 13.3(3) x density", shade=.01, border="black", ticktype="detailed", nticks=3,cex.lab=1.5, col="white", scale=FALSE)
b1.l1 = trans3d(x=rho.grid[1:31], y=rho.grid[1:31], z=diag(p.rho[1:31,1:31]), pmat=b1)
lines(b1.l1, lwd=1)
b1.l3 = trans3d(x=rho.grid[21], y=rho.grid[21], z=c(0,p.rho[21,21]), pmat=b1)
lines(b1.l3, lwd=2)
dev.off()