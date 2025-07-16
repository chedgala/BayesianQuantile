library(lqr)

# Normal Density
densN = function(x,mu=0,sigma=1,p=0.5){
  return(2*ifelse(x-mu<=0,p*dnorm(x,mu,sqrt((sigma^2)/(4*(1-p)^2))),(1-p)*dnorm(x,mu,sqrt((sigma^2)/(4*p^2)))))
}

# Slash Density
densSl = function(x,mu,sigma,nu,p){
  u  = ifelse(x<=mu,yes = rtrunc(n = 1,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(1-p)^2),no = rtrunc(n = 1,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(p)^2))
  gu = densN(x,mu,sigma*(u^(-1/2)),p)
  fu = dbeta(x = u,shape1 = nu,shape2 = 1)
  qu = ifelse(x<=mu,yes = dtrunc(x = u,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(1-p)^2),no = dtrunc(x = u,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(p)^2))
  return(gu*fu/qu)
}

# Log-likelihood : Slash
loglikSl = function(x,mu,sigma,nu,p){
  return( (log(densSl(x,mu,sigma,nu,p))) )
}

# HDP interval
hpd = function(x, alpha) {
  n = length(x)
  m = max(1, ceiling(alpha * n))
  y = sort(x)
  a = y[1:m]
  b = y[(n - m + 1):n]
  i = order(b - a)[1]
  structure(c(a[i], b[i]))
}

# CPO criterion
CPO = function(result){
  ll   = data.frame(rstan::extract(result, pars = "log_invlik"))
  CPOb = sum(log(1/apply(ll,2,mean)))
  return(CPOb)
}

# Information criteria
get_criteria = function(fit){
  eAIC = mean(rstan::extract(fit)$AIC) 
  eBIC = mean(rstan::extract(fit)$BIC) 
  eHQ  = mean(rstan::extract(fit)$HQ)
  eloglik = mean(rstan::extract(fit)$likel)
  
  log_lik_n = loo::extract_log_lik(fit, merge_chains=FALSE)
  loo_n     = loo::loo(log_lik_n)
  lmpl      = loo_n$estimates[1]
  WAIC_n    = loo::waic(log_lik_n)
  critFin   = list(eloglik=eloglik, eAIC=eAIC, eBIC=eBIC, eHQ=eHQ, lmpl=lmpl,
                   loo=loo_n$estimates[[3]], waic=WAIC_n$estimates[[3]], cpo=CPO(fit))
  return(critFin)  
}

# Information criteria : Slash
getloglik = function(fit_stan, dist, p, n, y){
  pbeta = ncol(rstan::extract(fit_stan)$beta)
  
  if(dist == "slash"){
    chains = cbind( beta=rstan::extract(fit_stan)$beta,
                    sigma=rstan::extract(fit_stan)$sigma, 
                    nu=rstan::extract(fit_stan)$nu)
    npar = pbeta + 2
    loglike = t(apply(X=chains, MARGIN=1, FUN=function(x) loglikSl(x=y, mu=x1%*%x[1:pbeta],
                                                                   sigma=x[npar-1], nu=x[npar], p=p)))
  }
  loo_n   = loo::loo(loglike)
  lmpl    = loo_n$estimates[1]
  WAIC_n  = loo::waic(loglike)
  CPOb    = sum(log(1/apply( exp(-loglike),2,mean )))

  logliks = rowSums(loglike)
  eAIC = mean( -2*logliks + 2*npar )
  eBIC = mean( -2*logliks + log(n)*npar )
  eHQ  = mean( -2*logliks + 2*log(log(n))*npar )
  eloglik = mean(logliks)
  
  return(list(eloglik=eloglik, eAIC=eAIC, eBIC=eBIC, eHQ=eHQ, lmpl=lmpl,
              loo=loo_n$estimates[[3]], waic=WAIC_n$estimates[[3]], cpo=CPOb))
}
