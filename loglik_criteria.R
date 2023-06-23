
getloglik <- function(fit_stan,dist,p){
  
  pbeta = ncol(rstan::extract(fit_stan)$beta)
  
  if(dist == "normal"){
    
    chains = cbind(
      beta = rstan::extract(fit_stan)$beta,
      sigma = rstan::extract(fit_stan)$sigma)
    
    npar = pbeta + 1
    
    logliks = apply(X = chains,MARGIN = 1,FUN = function(x) loglikN(x = y,
                                                                    mu = x1%*%x[1:pbeta],
                                                                    sigma = x[npar],
                                                                    p = p))
  }
  
  if(dist == "t"){
    
    chains = cbind(
      beta = rstan::extract(fit_stan)$beta,
      sigma = rstan::extract(fit_stan)$sigma,
      nu = rstan::extract(fit_stan)$nu
    )
    
    npar = pbeta + 2
    
    logliks = apply(X = chains,MARGIN = 1,FUN = function(x) loglikT(x = y,
                                                                    mu = x1%*%x[1:pbeta],
                                                                    sigma = x[npar - 1],
                                                                    nu = x[npar],
                                                                    p = p))
  }
  
  if(dist == "laplace"){
    
    chains = cbind(
      beta = rstan::extract(fit_stan)$beta,
      sigma = rstan::extract(fit_stan)$sigma)
    
    npar = pbeta + 1
    
    logliks = apply(X = chains,MARGIN = 1,FUN = function(x) loglikL(x = y,
                                                                    mu = x1%*%x[1:pbeta],
                                                                    sigma = x[npar],
                                                                    p = p))
  }
  
  AIC  = -2*logliks +2*npar
  BIC  = -2*logliks +log(n)*npar
  HQ   = -2*logliks +2*log(log(n))*npar
  
  eloglik = mean(logliks)
  eAIC = mean(AIC)
  eBIC = mean(BIC)
  eHQ = mean(HQ)
  
  return(list(eloglik = eloglik,eAIC = eAIC,eBIC = eBIC, eHQ = eHQ))
  
}

# Run ---------------------------------------------------------------------

source("auxiliar.R")
getloglik(fit_stanN,dist = "normal",p = 0.5)

