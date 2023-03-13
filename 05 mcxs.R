############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Wozniak
# R file for Lecture 5: Unit-root processes
############################################################


# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"


############################################################
# Unit-circle plot
############################################################
library(plotrix)
library(shape)
scale = 1.3
roots.x = c(-.9,-.9,.5) 
roots.y = c(-.9,sqrt(.19),.5)
roots.lab = c("outside of the unit circle","a unit root, on the unit circle","inside of the unit circle") 

pdf(file="unitcircle.pdf",height=7, width=7)
plot(x=c(-scale,scale), y=c(-scale,scale), type="n", axes=FALSE, ylab="", xlab="", asp=1)
draw.circle(0,0,1, nv=2000, lwd=2, border=mcxs2)
arctext(x = "unit circle", center = c(0, 0), radius = 1.05, middle = 7*pi/4, clockwise=FALSE, col=mcxs2)
Arrows(-scale,0,scale,0, arr.type="simple", col=mcxs1)
Arrows(0,-scale,0,scale, arr.type="simple", col=mcxs1)
text(1.05,-.05,"1", col=mcxs1)
text(-.95,-.05,"-1", col=mcxs1)
text(.05,-1.05,"-1", col=mcxs1)
text(.05,.95,"1", col=mcxs1)
text(.95,.05,"real part: a", col=mcxs1)
text(-.05,.85,"imaginary part: b", srt=90, col=mcxs1)
points(roots.x,roots.y, pch=19, col=mcxs3)
text(roots.x,roots.y,roots.lab, pos=4, col=mcxs3)
dev.off()

############################################################
# Empirical illustration for Australian real GDP
############################################################

rgdp.download = read_abs(series_id="A2304402X")

rgdp          = ts(log(rgdp.download[,6]),start=c(1959,3),frequency=4)
drgdp         = 100*diff(rgdp)

library(FinTS)
rgdp.range    = range(rgdp)
rgdp.acf      = Acf(rgdp,type="correlation", plot=FALSE, lag.max=20)
drgdp.range   = range(drgdp)
drgdp.acf     = Acf(drgdp,type="correlation", plot=FALSE, lag.max=20)

pdf(file="rgdp.pdf", width=15,height=10)
par(mfrow=c(2,2), mar=rep(3,4),cex.axis=1.5, col.lab=mcxs1)
plot(1:length(rgdp),rgdp, type="l", ylim=rgdp.range, axes=FALSE, xlab="", ylab="", col=mcxs2)
axis(1,c(3,43,83,123,163,203,243,246),c("1960","1970","1980","1990","2000","2010","2020",""), col=mcxs1)
axis(2,c(rgdp.range[1],mean(rgdp.range),rgdp.range[2]),c("","ln(RGDP)",""), col=mcxs1)
plot(1:20, rgdp.acf$acf[2:21], axes=FALSE, main="", xlab="", ylab="", ylim=c(0,1), type="h", lwd=10, col=mcxs2)
axis(1,c(1,10,20),c("1","10","20"), col=mcxs1)
axis(2,c(0,0.5,1),c("","acf","1"), col=mcxs1)
abline(h=0, col=mcxs1)

plot(1:length(drgdp),drgdp, type="l", ylim=drgdp.range, axes=FALSE, xlab="", ylab="", col=mcxs2)
axis(1,c(2,42,82,122,162,202,242,245),c("1960","1970","1980","1990","2000","2010","2020",""), col=mcxs1)
axis(2,c(-4,0,4),c("-4",expression(100*Delta*ln(RGDP)),"4"), col=mcxs1)
abline(h=0, col=mcxs1)
plot(1:20, drgdp.acf$acf[2:21], axes=FALSE, main="", xlab="", ylab="", type="h", lwd=10, col=mcxs2)
axis(1,c(1,10,20),c("1","10","20"), col=mcxs1)
axis(2,c(-.1,0,.1),c("-0.1","acf","0.1"), col=mcxs1)
abline(h=0, col=mcxs1)
dev.off()


Y         = as.matrix(rgdp[2:length(rgdp)])
X         = as.matrix(rgdp[1:(length(rgdp)-1)])
T         = nrow(Y)
alpha.hat = as.numeric(solve(t(X)%*%X)%*%t(X)%*%Y)
sigma2.hat= as.numeric(t(Y-alpha.hat*X) %*% (Y-alpha.hat*X))/T
alpha.hat.sd = sqrt(alpha.hat*as.numeric(solve(t(X)%*%X)))

set.seed(1234567)
S         = 10000
rw.ols    = function(x){as.numeric(solve(t(x[1:(length(x)-1)])%*%(x[1:(length(x)-1)]))%*%t(x[1:(length(x)-1)])%*%(x[2:length(x)]))}
rw        = matrix(rnorm(T*S,sd=sigma2.hat),T,S)
rw        = apply(rw,2,cumsum)
rw.alpha  = apply(rw,2,rw.ols)
alpha.as  = density(rw.alpha)

alpha.range = range(rw.alpha)
alpha.05    = quantile(rw.alpha, probs=0.05)

likeli.x    = seq(from=alpha.range[1], to=alpha.range[2], by=0.001)
likeli.y    = dnorm(likeli.x, mean=alpha.hat, sd=alpha.hat.sd)
likeli.lb   = qnorm(c(0.05), mean=alpha.hat, sd=alpha.hat.sd)

mcxs1.rgb 	= col2rgb(mcxs1)
mcxs1.trans = rgb(mcxs1.rgb[1], mcxs1.rgb[2], mcxs1.rgb[3], alpha=100, maxColorValue=255)
mcxs2.rgb 	= col2rgb(mcxs2)
mcxs2.trans = rgb(mcxs2.rgb[1], mcxs2.rgb[2], mcxs2.rgb[3], alpha=100, maxColorValue=255)

pdf(file=paste("inference.pdf",sep=""), width=9,height=7)
plot(alpha.as, axes=FALSE, xlab="", ylab="",main="", ylim=c(0,max(likeli.y)),col=mcxs2)
lines(likeli.x,likeli.y,lwd=2, col=mcxs1)
polygon(
  c(likeli.x[likeli.x<likeli.lb],likeli.lb,likeli.lb,likeli.x[1]),
  c(likeli.y[likeli.x<likeli.lb],mean(likeli.y[110:111]),0,0),
  col=mcxs1.trans ,border=mcxs1.trans)
lines(likeli.x,likeli.y,lwd=2, col=mcxs1)
polygon(
  c(alpha.as$x[alpha.as$x<alpha.05],alpha.05,alpha.05,alpha.as$x[1]),
  c(alpha.as$y[alpha.as$x<alpha.05],mean(alpha.as$y[301:302]),0,0),
  col=mcxs2.trans ,border=mcxs2.trans)
lines(alpha.as,lwd=2, col=mcxs2)
axis(1,c(alpha.range[1],alpha.05,likeli.lb,alpha.hat,alpha.range[2]),c("",".96",".99",expression(hat(alpha)[1]),""), cex.axis=1, col=mcxs1)
axis(2,c(0,max(likeli.y)),c("",""), col=mcxs1)
mtext("density", side = 2, line = 1, cex=1, col=mcxs1)
legend(x=.88, y=65, legend=c("asymptotic distribution","likelihood function"), lwd=c(2,2), lty=c(1,1),col=c(mcxs2,mcxs1),bty="n")
dev.off()
