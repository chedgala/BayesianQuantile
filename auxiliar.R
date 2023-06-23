########################################################################
#LOG LIKELIHOOD FUNCTIONS
########################################################################

#DENSITIES

densN = function(x,mu=0,sigma=1,p=0.5)
{
  return(2*ifelse(x-mu<=0,p*dnorm(x,mu,sqrt((sigma^2)/(4*(1-p)^2))),(1-p)*dnorm(x,mu,sqrt((sigma^2)/(4*p^2)))))
}

densT = function(x,mu=0,sigma=1,nu=4,p=0.5)
{
  return(
    ifelse(test=x<mu,
           yes=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(1-p)^2 +1)^(-(nu+1)/2),
           no=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(p)^2 +1)^(-(nu+1)/2))
  )
}

densT(0,mu=0,sigma=1,nu=4,p=0.75)



densL = function(x,mu=0,sigma=1,p=0.5)
{
  return(ifelse(test=x<mu,yes=(2*p*(1-p)/sigma)*exp(2*(1-p)*(x-mu)/sigma),no=(2*p*(1-p)/sigma)*exp(-(2*p)*(x-mu)/sigma)))
}

densSl = function(x,mu,sigma,nu,p)
{
  u   = ifelse(x<=mu,yes = rtrunc(n = 1,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(1-p)^2),no = rtrunc(n = 1,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(p)^2))
  gu   = densN(x,mu,sigma*(u^(-1/2)),p)
  fu   = dbeta(x = u,shape1 = nu,shape2 = 1)
  qu   = ifelse(x<=mu,yes = dtrunc(x = u,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(1-p)^2),no = dtrunc(x = u,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(p)^2))
  return(gu*fu/qu)
}

densNC = function(x,mu=0,sigma=1,nu=0.1,gamma=0.1,p=0.5)
{
  return(nu*densN(x,mu,sigma/sqrt(gamma),p) + (1-nu)*densN(x,mu,sigma,p))
}


loglikN = function(x,mu,sigma,p)
{
  return(sum(log(ifelse(x-mu<=0,p*dnorm(x,mu,sqrt((sigma^2)/(4*(1-p)^2))),(1-p)*dnorm(x,mu,sqrt((sigma^2)/(4*p^2)))))))
}

loglikL = function(x,mu,sigma,p)
{
  return(sum(log(densL(x,mu,sigma,p))))
}

loglikT = function(x,mu,sigma,nu,p)
{
  return(sum(log(
    ifelse(test=x<mu,
           yes=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(1-p)^2 +1)^(-(nu+1)/2),
           no=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(p)^2 +1)^(-(nu+1)/2))
  )))
}

loglikSl = function(x,mu,sigma,nu,p)
{
  return(sum(log(
    densSl(x,mu,sigma,nu,p)
  )))
}

loglikNC = function(x,mu,sigma,nu,gamma,p)
{
  return(sum(log(
    densNC(x,mu,sigma,nu,gamma,p)
  )))
}

AUXloglikNC = function(x,mu,sigma,par,p)
{
  return(sum(log(
    densNC(x,mu,sigma,par[1],par[2],p)
  )))
}

