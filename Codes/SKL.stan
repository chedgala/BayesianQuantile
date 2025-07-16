data{
  int x_p;
  int n;
  matrix[n,x_p] x;
  vector[n] y;
  vector[n] w; 
  real<lower=0, upper=1> p;
}

parameters {
  vector[x_p]  beta;
  real<lower=0> sigma;
  vector<lower=0>[n] u; // slash, Student, Laplace
}

transformed parameters{
   vector[n] mu = x*beta;  
}

model{ 
  vector[n] sigmai;
  for (r in 1:x_p){
    beta[r] ~ normal(0, 1000); 
  } 
  sigma ~ cauchy(0, 4);
  
  for (i in 1:n){
    u[i] ~ inv_gamma(1, 0.5);  //laplace
    if (w[i] <= mu[i]){
      sigmai[i] = inv_sqrt(u[i])*sigma/(2*(1-p));  // slash, Student, Laplace
      target += bernoulli_lpmf(1 | p) + normal_lpdf(y[i] | mu[i], sigmai[i]);
    } else {
      sigmai[i] =  inv_sqrt(u[i])*sigma/(2*p); //slash, Student,Laplace
      target += bernoulli_lpmf(0 | p) + normal_lpdf(y[i] | mu[i], sigmai[i]); 
    }                     
  }               
}

generated quantities{
	vector[n] log_lik;
	vector[n] log_invlik;
	for (j in 1:n){
	  if (w[j] < mu[j]){
      log_lik[j] = log(2*p*(1-p)/sigma) + 2*(1-p)*(y[j]-mu[j])/sigma;
    } else {
      log_lik[j] = log(2*p*(1-p)/sigma) - 2*p*(y[j]-mu[j])/sigma; 
    } 
    log_invlik[j]= exp(-log_lik[j]);
	}
  real likel= sum(log_lik);
  real AIC  = -2*likel + 2*(x_p+1);
  real BIC  = -2*likel + log(n)*(x_p+1);
  real HQ   = -2*likel + 2*log(log(n))*(x_p+1);
}

