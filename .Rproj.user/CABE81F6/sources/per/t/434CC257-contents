data{
   real<lower=0,upper=1> p;
  int<lower=0> n;   //total observations
  int<lower=0> nind;   //total individuos
  int<lower=0> x_p;    // number covariates
  vector[n] y;      // responses
  int<lower=1,upper=nind> cluster[n];
  matrix[n,x_p] x;  //covariates
}


parameters {
  vector[x_p]  beta;
  vector[nind]  b;
  real<lower=0> sigmab;
  real<lower=0> sigma;
  vector<lower =0>[n] w;
}
 
transformed parameters{
    real<lower=0> sigmab2;
    sigmab2 = sigmab^2;
}

model{
    vector[n] mu;
    vector[n] media;
    vector[n] pe;
    
    mu = x*beta;
   	for(i in 1 : nind) {
    b[i]~normal(0,sigmab);
		}
  
  for (r in 1:x_p){
  beta[r] ~ normal(0, 10); //regression parameters
  }
  
  sigmab ~ student_t(4,0,5);
  sigma ~ student_t(4,0,5);         

for(i in 1 : n){
 w[i] ~ exponential(1);
 pe[i] = (2*w[i]*sigma^2)/(p*(1-p));
 media[i] = (1-2*p)*sigma/(p*(1-p))*w[i]+mu[i]+b[cluster[i]]; 
 y[i] ~ normal(media[i],sqrt(pe[i]));   
}
}

