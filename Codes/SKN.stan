data{
  int x_p;
  int n;
  matrix[n,x_p] x;
  vector[n] y;
  vector[n] w; 
  real<lower=0, upper=1> p;
}

parameters {
  vector[x_p] beta;
  real<lower=0> sigma;
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
    if (w[i] <= mu[i]){
      sigmai[i] = sigma/(2*(1-p)); // normal
      target += bernoulli_lpmf(1 | p) + normal_lpdf(y[i] | mu[i], sigmai[i]);
    } else {
      sigmai[i] = sigma/(2*p); //  normal
      target += bernoulli_lpmf(0 | p) + normal_lpdf(y[i] | mu[i], sigmai[i]); 
    }                     
  }               
}

generated quantities{
	vector[n] log_lik;
	vector[n] log_invlik;
	vector[n] sigmai;
	for (j in 1:n){
	  if (w[j] <= mu[j]){
      sigmai[j]  = sigma/(2*(1-p)); // normal
      log_lik[j] = log(p) + normal_lpdf(y[j] | mu[j], sigmai[j]);
    } else {
      sigmai[j]  = sigma/(2*p); //  normal
      log_lik[j] = log(1-p) + normal_lpdf(y[j] | mu[j], sigmai[j]); 
    } 
    log_invlik[j]= exp(-log_lik[j]);
	}
  real likel= sum(log_lik);
  real AIC  = -2*likel + 2*(x_p+1);
  real BIC  = -2*likel + log(n)*(x_p+1);
  real HQ   = -2*likel + 2*log(log(n))*(x_p+1);
}