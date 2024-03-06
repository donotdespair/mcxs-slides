
# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"

load("honeyhoney-data.RData")
y        = Y[,1]

########################
# OLS
ols      = ar(y, aic=FALSE, order.max=1)
ols.mu   = ols$ar
ols.var  = var(ols$resid[2:192])

# prior
pri      = ar(y, aic=FALSE, order.max=14)
pri.sig  = var(pri$resid[15:192])

pri.mu   = 0.9
pri.var  = (2^2)*pri.sig

# posterior
pos.var  = 1/(1/pri.var + 1/ols.var)
pos.mu   = ((1/pri.var)*pri.mu + (1/ols.var)*ols.mu)*pos.var

(1/pri.var)*pos.var
(1/ols.var)*pos.var

from     = 0.82
to       = 1.03
grid     = seq(from=from, to=to, length.out=300)
prior    = dnorm(grid, mean = pri.mu, sd=sqrt(pri.var))
likel    = dnorm(grid, mean = ols.mu, sd=sqrt(ols.var))
poste    = dnorm(grid, mean = pos.mu, sd=sqrt(pos.var))

limits   = range(c(prior,likel,poste))

pdf(file="learning.pdf", width=9,height=6)
plot(grid, likel, type="l", xlim=c(from,to), ylim=limits, axes=FALSE, xlab="", ylab="", main="", lwd=3,col=mcxs4)
axis(1,c(0.85,pri.mu,pos.mu,ols.mu,1),c("","","","",""), col=mcxs2)
axis(2,c(0,20,40),c("","",""), col=mcxs2)
mtext("parameter", side = 1, line = 2, cex=1.5, col=mcxs2)
mtext("density", side = 2, line = 1, cex=1.5, col=mcxs2)
lines(grid, prior, lwd=3,col=mcxs3)
lines(grid, poste, lwd=3,col=mcxs2)
text(pri.mu,10,"prior", cex=1.5, col=mcxs3)
text(1.01,30,"likelihood", cex=1.5, col=mcxs4)
text(0.94,40,"posterior", cex=1.5, col=mcxs2)
dev.off()







# better learning
########################
# OLS
ols      = ar(y, aic=FALSE, order.max=1)
ols.mu   = ols$ar
ols.var  = var(ols$resid[2:192])

# prior
pri      = ar(y, aic=FALSE, order.max=14)
pri.sig  = var(pri$resid[15:192])

pri.mu   = 1
pri.var  = 0.01

# posterior
pos.var  = 1/(1/pri.var + 1/ols.var)
pos.mu   = ((1/pri.var)*pri.mu + (1/ols.var)*ols.mu)*pos.var

(1/pri.var)*pos.var
(1/ols.var)*pos.var

from     = 0.8
to       = 1.2
grid     = seq(from=from, to=to, length.out=300)
prior    = dnorm(grid, mean = pri.mu, sd=sqrt(pri.var))
likel    = dnorm(grid, mean = ols.mu, sd=sqrt(ols.var))
poste    = dnorm(grid, mean = pos.mu, sd=sqrt(pos.var))

limits   = range(c(prior,likel,poste))

pdf(file="grphs/learning_better.pdf", width=9,height=6)
plot(grid, likel, type="l", xlim=c(from,to), ylim=limits, axes=FALSE, xlab="", ylab="", main="", lwd=3,col=mcxs4)
axis(1,c(0.85,pri.mu,pos.mu,ols.mu,1.15),c("","","","",""), col=mcxs2)
axis(2,c(0,20,40),c("","",""), col=mcxs2)
mtext("parameter", side = 1, line = 2, cex=1.5, col=mcxs2)
mtext("density", side = 2, line = 1, cex=1.5, col=mcxs2)
lines(grid, prior, lwd=3,col=mcxs3)
lines(grid, poste, lwd=3,col=mcxs2)
text(pri.mu,10,"prior", cex=1.5, col=mcxs3)
text(1.03,30,"likelihood", cex=1.5, col=mcxs4)
text(0.94,35,"posterior", cex=1.5, col=mcxs2)
dev.off()





##############################
from     = -1
to       = 1
grid     = seq(from=from, to=to, length.out=300)
prior    = dnorm(grid, mean = 0, sd=0.3)

pdf(file="shrinkage.pdf", width=9,height=6)
plot(grid, prior, type="l", xlim=c(from,to), axes=FALSE, xlab="", ylab="", main="", lwd=3,col=mcxs3)
axis(1,c(-1,0,1),c("","",""),col=mcxs2)
axis(2,c(0,0.6,1.2),c("","",""),col=mcxs2)
mtext("parameter", side = 1, line = 2, cex=1.5,col=mcxs2)
mtext("density", side = 2, line = 1, cex=1.5,col=mcxs2)
text(0,0.2,"shrinkage", cex=1.5,col=mcxs3)
arrows(0.2,0.2,0.58,0.2,lwd=3,col=mcxs3)
arrows(-0.2,0.2,-0.58,0.2,lwd=3,col=mcxs3)
dev.off()




######################
from     = -2
to       = 2
sdev     = 0.3
l.ratio  = 3
mu       = 0.8
grid     = seq(from=from, to=to, length.out=300)
ignorance= rep(1/(to-from),length(grid))
believes = dnorm(grid, mean = 0, sd=sdev)
likelihood= dnorm(grid, mean = mu, sd=l.ratio*sdev)

limits   = range(c(ignorance,believes,likelihood))

pdf(file="priors.pdf", width=9,height=6)
plot(grid, believes, type="l", xlim=c(from,to), ylim=limits, axes=FALSE, xlab="", ylab="", main="", lwd=3,col=mcxs5, cex.lab=1.5)
axis(1,c(from,0,mu,to),c("","",expression(hat(theta)),""), cex.axis=1.5, col=mcxs1)
axis(2,c(0,0.6,1.2),c("","",""), col=mcxs1)
mtext("parameter", side = 1, line = 2, cex=1.5, col=mcxs1)
mtext("density", side = 2, line = 1, cex=1.5, col=mcxs1)
lines(grid, ignorance, lwd=3, col=mcxs3)
lines(grid, likelihood, lwd=3, col=mcxs1)
text(-1.5,0.3,"agnostic view", cex=1.5, col=mcxs3)
text(1.2,0.5,"likelihood", cex=1.5, col=mcxs1)
text(0,0.5,"strong beliefs", cex=1.5, col=mcxs5)
abline(v=mu,lty=2, col=mcxs1)
lines(c(grid[grid<mu],mu), rep(dnorm(mu, mean = 0, sd=sdev),length(grid[grid<mu])+1),lty=2, col=mcxs1)
dev.off()



######################
from     = -2
to       = 2
sdev1    = 0.5
mu1      = 1
sdev2    = 1
mu2      = -1
grid     = seq(from=from, to=to, length.out=300)
modelA   = dnorm(grid, mean = mu1, sd=sdev1)
modelB   = dnorm(grid, mean = mu2, sd=sdev2)
Y        = 0.5
pA       = dnorm(Y, mean = mu1, sd=sdev1)
pB       = dnorm(Y, mean = mu2, sd=sdev2)
limits   = range(c(modelA,modelB))

pdf(file="mdd.pdf", width=9,height=6)
plot(grid, modelA, type="l", xlim=c(from,to), ylim=limits, axes=FALSE, xlab="", ylab="", main="", lwd=3,col="gray50", cex=1.5)
axis(1,c(from,0,Y,to),c("","","Y",""), cex.axis=1.5)
axis(2,c(0,pA,pB,0.8),c("","","",""), cex.axis=1.5)
mtext("data", side = 1, line = 2, cex=1.5)
mtext("density", side = 2, line = 1, cex=1.5)
lines(grid, modelB, lwd=3,col="gray70")
text(mu1,0.6,expression("Model M"[1]), cex=1.5)
text(mu2,0.3,expression("Model M"[2]), cex=1.5)
text(-1.8,pA+0.03,expression("p(Y|M"[1]), cex=1.5)
text(-1.57,pA+0.03,")", cex=1.5)
text(-1.8,pB+0.03,expression("p(Y|M"[2]), cex=1.5)
text(-1.57,pB+0.03,")", cex=1.5)
abline(v=Y,lty=2)
lines(c(grid[grid<Y],0.5), rep(pA,length(grid[grid<Y])+1),lty=2)
lines(c(grid[grid<Y],0.5), rep(pB,length(grid[grid<Y])+1),lty=2)
dev.off()
