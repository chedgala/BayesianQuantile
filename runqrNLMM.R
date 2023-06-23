#################################################################################
## soybean
################################################################################
library("qrNLMM")
data(Soybean)
attach(Soybean)
y = weight #response
x = Time #time
covar = as.numeric(as.factor(Variety))-1 #factor genotype (0=Forrest, 1=Plan Introduction)
cluster = as.numeric(as.factor(Plot))
q1<-3
nind<-max(cluster)
p<-0.75
x_p<-4
n<-length(y)
library(rstan, quietly = T)
library(rstanmulticore)
library(shinystan)

library(nlme)
plot(Soybean)
meanmodel <- nlme(model = weight ~ SSlogis(Time, Asym, xmid, scal),
     fixed = list(Asym ~ Variety, xmid ~ 1, scal ~ 1),
     random = Asym + xmid + scal ~1|Plot,
     start = c(16,3,55,8),
     data=Soybean)

parinit = c(fixed.effects(meanmodel)[c(1,3,4,2)],
            pi,
            as.numeric(var(meanmodel$coefficients$random$Plot)))

fit_stan <- pstan(file='censNLM-Multi-Laplace.stan',
                 data = list(y=c(y), x_p=x_p, x=covar,n=n, p=p, cluster=cluster, nind=nind,  q1=q1, tt=x),
                 thin = 200, chains = 4, iter = 10000, warmup = 1000,
                 seed = 9955,init = parinit,verbose = TRUE)
qoi <- c("beta","sigma","D1")
print(fit_stan, pars=qoi)


exprNL = expression((fixed[1]+(fixed[4]*covar[1])+random[1])/
                      (1 + exp(((fixed[2]+random[2])- x)/(fixed[3]+random[3]))))
#Initial values for fixed effects
initial = c(max(y),0.6*max(y),0.73*max(y),3)
# A quantile regression for the three quartiles
box_reg = QRNLMM(y,x,Plot,initial,exprNL,covar,p=0.5)
