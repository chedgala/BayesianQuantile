################################################################################
## Simulated data
############################################################################

library(lqr)
n<-200
p<-0.75
sigmas<-1.3
#nu=3
beta<-c(-1,2,3)

x1<-cbind(1,runif(n),rbinom(n,1,0.3))

y<-x1%*%beta+rSKD(n, mu = 0, sigma = sigmas, p = p, dist = "laplace")
#y<-x1%*%beta+rSKD(n, mu = 0, sigma = sigmas, p = p, dist = "slash", nu =nu)
x_p<- ncol(x1)


modelL = lqr(y~x1-1,dist = "laplace", p=p,precision = 0.000000001)

library(rstan, quietly = T)
library(shinystan)

#initf1 <- function() {list(beta=modelN$beta, sigma=modelN$sigma)}
initf1 <- function() {list(beta=modelN$beta, sigma=modelN$sigma, u=rep(0.8,n))}

fit_stan <- stan(file='SKL.stan',    init=  initf1,
                 data = list(y=c(y),w=c(y), x_p=x_p,x=x1,n=n,p=p),
                 thin = 1000, chains = 2, iter = 20000, warmup = 1000,
                 seed = 9955)

qoi <- c("beta","sigma")
#qoi <- c("beta","sigma","nu")
#pairs(fit_stan, pars=qoi)
plot(fit_stan,pars = qoi)
print(fit_stan, pars=qoi)

modelL = lqr(y~x1-1,dist = "laplace", p=p, precision = 0.000000001)
#modelN = lqr(y~x1-1,dist = "slash", p=p, precisio=0.000000001)
##################################################