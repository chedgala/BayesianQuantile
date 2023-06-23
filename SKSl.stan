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
  // real<lower=2> nu;      //Student
  real<lower=1> nu;        //Slash
  vector<lower = 0>[n] u; // slash, Student, Laplace
}

transformed parameters{
  vector[n] mu;
  mu = x*beta;  
}

model{ 
  vector[n] sigmai;
  for (r in 1:x_p){
    beta[r] ~ normal(0, 1000); 
  } 
  sigma ~ cauchy(0, 4);
  nu ~ student_t(4,0,5);  // Student and Slash
  for (i in 1:n){
    // u[i] ~ gamma(nu,nu);    //Student
     u[i] ~ beta(nu,1);  //Slash
    //   u[i] ~ inv_gamma(1, 0.5);  //laplace
    if (w[i]<=mu[i]){
      sigmai[i] = inv_sqrt(u[i])*sigma/(2*(1-p));  // slash, Student, Laplace
      //sigmai[i] = sigma/(2*(1-p));              // normal
      target += bernoulli_lpmf(1 | p)+normal_lpdf(y[i] | mu[i], sigmai[i]);
    }
    else {
      sigmai[i] =  inv_sqrt(u[i])*sigma/(2*p); //slash, Student,Laplace
      //sigmai[i] =  sigma/(2*p);  //  normal
      target += bernoulli_lpmf(0 | p)+normal_lpdf(y[i] | mu[i], sigmai[i]); 
    }                     
  }               
}
