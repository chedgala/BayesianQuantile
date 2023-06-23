library(qrLMM)
data(Orthodont)
attach(Orthodont)
sex = c()
sex[Sex=="Male"]=0
sex[Sex=="Female"]=1
y = distance #response
x = cbind(1,sex,age) #design matrix for fixed effects
z = cbind(1,age) #design matrix for random effects
cluster = sort(as.numeric(as.factor(Subject)))
#A median regression
median_reg = QRLMM(y,x,z,Subject,MaxIter = 500)

## STAN

q1<-ncol(z)
nind<-max(cluster)
p<-0.5
x_p<-ncol(x)
n<-length(y)
library(rstan, quietly = T)
library(shinystan)

fit_stan <- stan(file='censCho-Multi-Laplace.stan',  
                 data = list(y=c(y), x_p=x_p, x=x,n=n, p=p, cluster=cluster, nind=nind, z=z, q1=q1),
                 thin = 50, chains = 1, iter = 5000, warmup = 1000,
                 seed = 9955)
qoi <- c("beta","sigma","D1")
print(fit_stan, pars=qoi)
median_reg = QRLMM(y,x,z,Subject,MaxIter = 500)