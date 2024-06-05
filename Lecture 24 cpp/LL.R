
rm(list=ls())

library(mgcv)
library(Rcpp)
library(microbenchmark)

source("LL_nothing_special.R")
source("LL_band_precision.R")
source("LL_tridiag_precision.R")
sourceCpp("LL_arma_inv.cpp")
sourceCpp("LL_arma_solve.cpp")
sourceCpp("LL_arma_stochvol.cpp")

load("LL_sim_data.RData")

########################################################
# Compare the speed for sample size 250
########################################################
S = 10

# system.time(LL_tridiag_precision(y250, S=S))
# system.time(LL_nothing_special(y250, S=S))
# system.time(LL_arma_solve(y250, S=S))

t0  = proc.time()
microbenchmark(R.not     = LL_nothing_special(y250, S),
               R.ban     = LL_band_precision(y250, S),
               R.tri     = LL_tridiag_precision(y250, S),
               cpp.inv   = LL_arma_inv(y250, S),
               cpp.sol   = LL_arma_solve(y250, S),
               cpp.sto   = LL_arma_stochvol(y250, S),
               setup  = set.seed(123)
               )
t1  = proc.time()
(t1-t0)[3]/60

c(R.not$last.draw$mu0, 
  R.ban$last.draw$mu0,
  R.tri$last.draw$mu0, 
  cmp.not$last.draw$mu0,
  cmp.ban$last.draw$mu0,
  cmp.tri$last.draw$mu0,
  cpp.inv$last.draw$mu0,
  cpp.sol$last.draw$mu0, 
  cpp.sto$last.draw$mu0)




########################################################
# Compare the speed for sample size 750
########################################################
S=10
system.time(LL_tridiag_precision(y750, S=S))
system.time(LL_band_precision(y750, S=S))
system.time(LL_arma_solve(y750, S=S))

# this should take around 40 minutes
t0  = proc.time()
microbenchmark(tridiag <- LL_tridiag_precision(y750, S=S),
               band    <- LL_band_precision(y750, S=S),
               usual   <- LL_nothing_special(y750, S=S),
               # cpp.inv = rmvnorm_arma_set_seed(y750, S=S),
               cpp.sol <- LL_arma_solve(y750, S=S),
               # check  = "equal",
               setup  = set.seed(123),
               times  = 100
)
t1  = proc.time()
(t1-t0)[3]/60


########################################################
# Compare estimation results
########################################################
S = 11000
set.seed(123)
system.time(qqq.r <- LL_tridiag_precision(y250, S=S))
set.seed(123)
system.time(qqq.b <- LL_band_precision(y250, S=S))
set.seed(123)
system.time(qqq.c <- LL_arma_solve(y250, S=S))

final = 1001:11000
plot.ts(y250, lwd=2)
lines(apply(qqq.r$posterior$mu[,final],1,mean), col="gray")
lines(apply(qqq.b$posterior$mu[,final],1,mean), col="blue")
lines(apply(qqq.c$posterior$mu[,final],1,mean), col="red")
system("say do not despair")

plot.ts(t(qqq.c$posterior$mu[1:10,final]))
plot.ts(cbind(qqq.c$posterior$mu0[final],qqq.c$posterior$sigma2[final],qqq.c$posterior$sigma2m[final]))
mean(qqq.c$posterior$mu0)
mean(qqq.c$posterior$sigma2)
mean(qqq.c$posterior$sigma2m)



mean(qqq.r$posterior$mu0)
mean(qqq.b$posterior$mu0)
mean(qqq.c$posterior$mu0)

mean(qqq.r$posterior$sigma2)
mean(qqq.b$posterior$sigma2)
mean(qqq.c$posterior$sigma2)

mean(qqq.r$posterior$sigma2m)
mean(qqq.b$posterior$sigma2m)
mean(qqq.c$posterior$sigma2m)





S=10
set.seed(123)
qqq.r <- LL_tridiag_precision(y250, S=S)
set.seed(123)
qqq.b <- LL_band_precision(y250, S=S)
set.seed(123)
qqq.c <- LL_arma_solve(y250, S=S)

qqq.r$last.draw$mu0
qqq.b$last.draw$mu0
qqq.c$last.draw$mu0

qqq.r$last.draw$mu[250]
qqq.b$last.draw$mu[250]
qqq.c$last.draw$mu[250]

