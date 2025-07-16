# Librerías
library(lqr)
library(rstan)

# Directorio
pathfiles = paste0(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(pathfiles)
source("auxiliar.R")


# Código de ejemplo: Skew-Normal ------------------------------------------
# Valores iniciales
x1  # Matriz de diseño
y   # Variable de respuesta
p   # Cuantil 

modelN = lqr(y~x1-1,dist = "normal", p=p, precision = 0.000000001,silent = TRUE)
initf1 = function() {list(beta=modelN$beta, sigma=modelN$sigma)}
x_p    = length(initf1$beta)
n      = nrow(x1)

# Parameter estimation
fit_stan = stan(file='SKN.stan', init=initf1,
                data=list(x_p=x_p, n=n, x=x1, y=c(y), w=c(y), p=p),
                chains=1, iter=5000, warmup=1000, seed=9955, verbose=FALSE)

# Information criteria
criteria = get_criteria(fit_stan)

# Confidence interval
CIbeta = apply(rstan::extract(fit_stan)$beta, 2, hpd, alpha=0.05)
CIsigm = hpd(rstan::extract(fit_stan)$sigma, alpha=0.05)

# Otros resultados
c(beta = colMeans(rstan::extract(fit_stan)$beta),
  sigma = mean(rstan::extract(fit_stan)$sigma),
  sdbeta = sqrt(diag(var(rstan::extract(fit_stan)$beta))),
  sdsigma = sd(rstan::extract(fit_stan)$sigma),
  criteria = unlist(criteria))


# Código de ejemplo: Skew-t -----------------------------------------------
modelT = lqr(y~x1-1,dist="t", p=p, precision = 0.000000001,silent = TRUE)
initf1 = function() {list(beta=modelT$beta, sigma=modelT$sigma, u=rep(0.8, n))}

# Parameter estimation
fit_stan = try(stan(file='SKT.stan', init=initf1,
                    data = list(x_p=x_p, n=n, x=x1, y=c(y), w=c(y), p=p),
                    chains = 1, iter = 5000, warmup = 1000, seed = 7162), silent=TRUE)


# Código de ejemplo: Skew-Laplace -----------------------------------------
modelL = lqr(y~x1-1,dist="laplace", p=p, precision = 0.000000001,silent = TRUE)
initf1 = function() {list(beta=modelL$beta, sigma=modelL$sigma)}

# Parameter estimation
fit_stan = stan(file='SKL.stan', init=initf1,
                data = list(x_p=x_p, n=n, x=x1, y=c(y), w=c(y), p=p),
                chains = 1, iter = 5000, warmup = 1000, seed = 9966,verbose = FALSE)


