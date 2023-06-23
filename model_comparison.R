# Parameters
sigmas<-1.3
beta<-c(-1,2,3)
x_p<- length(beta)

# Set up
M = 100

nvec = c(50,100,200)
pvec = c(50,75,90)/100

i = 3 # n = 200
j = 1 # p = 0.5

nu = 3

set.seed(13)
x1<-cbind(1,runif(n),rbinom(n,1,0.3))


semilla = 9955

# Modelo ------------------------------------------------------------------

modelo = "t"

if(modelo == "normal"){
y <- x1%*%beta + rSKD(nvec[i],mu = 0,sigma = sigmas,
                      p = pvec[j], dist = "t",nu = nu)
}

if(modelo == "t"){
y <- x1%*%beta + rSKD(nvec[i],mu = 0,sigma = sigmas,
                      p = pvec[j], dist = "normal")
}

if(modelo == "laplace"){
y <- x1%*%beta + rSKD(nvec[i],mu = 0,sigma = sigmas,
                      p = pvec[j], dist = "laplace")
}
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------



modelN = lqr(y~x1-1,dist = "normal", p=pvec[j],
             precision = 0.000000001,silent = TRUE)

initf1 <- function() {list(beta=modelN$beta, sigma=modelN$sigma)}

fit_stanN <- stan(file='SKN.stan', init=initf1,
                 data = list(y=c(y),w=c(y), x_p=x_p,x=x1,
                             n=nvec[i],p=pvec[j]),
                 chains = 1, iter = 5000, warmup = 1000,
                 seed = semilla,verbose = FALSE)


modelT = lqr(y~x1-1,dist = "t", p=pvec[j],
             precision = 0.000000001,silent = TRUE)

initf1 <- function() {list(beta=modelT$beta, sigma=modelT$sigma, u=rep(0.8,nvec[i]))}

fit_stanT <- stan(file='SKT.stan', init=initf1,
                  data = list(y=c(y),w=c(y), x_p=x_p,x=x1,n=nvec[i],p=pvec[j]),
                  chains = 1, iter = 5000, warmup = 1000,
                  seed = semilla)


modelL = lqr(y~x1-1,dist = "laplace", p=pvec[j],
             precision = 0.000000001,silent = TRUE)

initf1 <- function() {list(beta=modelL$beta, sigma=modelL$sigma)}

fit_stanL <- stan(file='SKL.stan', init=initf1,
                 data = list(y=c(y),w=c(y), x_p=x_p,x=x1,
                             n=nvec[i],p=pvec[j]),
                 chains = 1, iter = 5000, warmup = 1000,
                 seed = semilla,verbose = FALSE)

criteria = 
cbind(unlist(getloglik(fit_stanN,dist = "normal",p = 0.5)),
unlist(getloglik(fit_stanT,dist = "t",p = 0.5)),
unlist(getloglik(fit_stanL,dist = "laplace",p = 0.5)))
colnames(criteria) = c("normal","t","laplace")

criteria

best.lqr(y~x1-1,p=pvec[j],
         precision = 0.000000001)

# max(rstan::extract(fit_stanT)$lp__)>max(rstan::extract(fit_stanN)$lp__)
# mean(rstan::extract(fit_stanT)$lp__)>mean(rstan::extract(fit_stanN)$lp__)
# quantile(rstan::extract(fit_stanT)$lp__,probs = 0.999)>quantile(rstan::extract(fit_stanN)$lp__,probs = 0.999)

plot(rstan::extract(fit_stanN)$lp__,type="l",
     ylim = range(c(rstan::extract(fit_stanT)$lp__,
                    rstan::extract(fit_stanN)$lp__,
                    rstan::extract(fit_stanL)$lp__)),
     main = modelo)
lines(rstan::extract(fit_stanT)$lp__,type="l",col=2)
lines(rstan::extract(fit_stanL)$lp__,type="l",col=3)
# Agregar leyenda

legend("right", legend = c("Normal","t Student","laplace"),
       col = 1:3, lty = 1,lwd = 2)

# 
# 
# 
# mi_lista = rstan::extract(fit_stan)
#   
# #Utilizamos la función sapply para extraer el quinto objeto de cada vector numérico
# maxx <- which.max(rstan::extract(fit_stan)$lp__)
# maxx_objetos <- sapply(mi_lista, function(x) as.matrix(x)[maxx,])
# 
# # Imprimimos el vector con los quintos objetos
# print(maxx_objetos)
