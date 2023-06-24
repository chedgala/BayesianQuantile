################################################################################
## Simulated data
############################################################################

library(lqr)
#n<-200
n<-50

p<-0.90
sigmas<-1.3
nu=3
beta<-c(-1,2,3)

x1<-cbind(1,runif(n),rbinom(n,1,0.3))

#y<-x1%*%beta+rSKD(n, mu = 0, sigma = sigmas, p = p, dist = "normal")
y<-x1%*%beta+rSKD(n, mu = 0, sigma = sigmas, p = p, dist = "t", nu =nu)
x_p<- ncol(x1)


modelT = lqr(y~x1-1,dist = "t", p=p,precision = 0.000000001)

library(rstan, quietly = T)
library(shinystan)

#initf1 <- function() {list(beta=modelN$beta, sigma=modelN$sigma)}
initf1 <- function() {list(beta=modelT$beta, sigma=modelT$sigma, u=rep(0.8,n))}
initf1

iters = 5000

fit_stan <- stan(file='SKT.stan', init=  initf1,
                 data = list(y=c(y),w=c(y), x_p=x_p,x=x1,n=n,p=p),
                 chains = 1, iter = iters*1.05, warmup = iters*0.05) #seed = 9955

# library(rstanmulticore)
# 
# initf1 <- function(n) {list(beta=c(-1.051333 ,2.129753 ,2.651769),
#                            sigma=1.384141, u=rep(0.8,200))}
# 
# fit_stan <- pstan(file='SKT.stan', init=initf1,
#                  data = list(y=c(y),w=c(y), x_p=x_p,x=x1,n=n,p=p),
#                  thin = 1000, chains = 4, iter = 10000, warmup = 1000,
#                  seed = 9955,verbose = TRUE)
# 
# library(microbenchmark)
# microbenchmark(stan(file='SKT.stan', init=initf1,
#                                      data = list(y=c(y),w=c(y), x_p=x_p,x=x1,n=n,p=p),
#                                      thin = 1000, chains = 4, iter = 10000, warmup = 1000,
#                                      seed = 9955,verbose = TRUE),
#                                pstan(file='SKT.stan', init=initf1,
#                                      data = list(y=c(y),w=c(y), x_p=x_p,x=x1,n=n,p=p),
#                                      thin = 1000, chains = 4, iter = 10000, warmup = 1000,
#                                      seed = 9955,verbose = TRUE),times = 5)

#qoi <- c("beta","sigma")
qoi <- c("beta","sigma","nu")
#pairs(fit_stan, pars=qoi)
plot(fit_stan,pars = qoi)
print(fit_stan, pars=qoi)

# lqr(y~x1-1,dist = "t", p=p, precision = 0.000000001)

#modelN = lqr(y~x1-1,dist = "slash", p=p, precisio=0.000000001)
##################################################

# length(rstan::extract(fit_stan)$beta[,1])
# length(rstan::extract(fit_stan)$sigma)
# 
# plot.ts(rstan::extract(fit_stan)$beta[,1])
# plot.ts(c(rstan::extract(fit_stan)$sigma))
# plot.ts(c(rstan::extract(fit_stan)$nu))
