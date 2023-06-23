data{
   real<lower=0,upper=1> p;
  int<lower=0> n;   //total observations
  int<lower=0> nind;   //total individuos
  int<lower=0> x_p;    // number covariates
  int<lower=0> q1;    // dimension random effects
  vector[n] y;      // responses
  int<lower=1,upper=nind> cluster[n];
  vector[n] x;  //covariates
  vector[n] tt;//tempo
}


parameters {
  vector[x_p]  beta;
  real<lower=0> sigma;
  vector<lower =0>[n] w;
   vector<lower=0>[q1] ddsqrt; 
    matrix[q1, nind] etavec;
   cholesky_factor_corr[q1] Lcorr;// cholesky factor (L_u matrix for D1R)
}
 
transformed parameters{
  matrix[q1, nind] bvec;
  cov_matrix[q1] D1; // VCV matrix
  D1 = quad_form_diag(multiply_lower_tri_self_transpose(Lcorr), ddsqrt); // quad_form_diag: diag_matrix(ddsqrt) * d1R * diag_matrix(ddsqrt)
  {
    matrix[q1,q1] dL = diag_pre_multiply(ddsqrt, Lcorr);
	for (j in 1:nind){
	  bvec[,j] = dL * etavec[,j];
    }
}
}

model{
    vector[n] mu;
    vector[n] media;
    vector[n] pe;
    ddsqrt ~ student_t(4,0,5);
    Lcorr ~ lkj_corr_cholesky(5.0);
    to_vector(etavec) ~ std_normal();
  for (r in 1:x_p){
  beta[r] ~ normal(0, 10); //regression parameters
  }
  
  sigma ~ student_t(4,0,5);  
       
for(i in 1 : n){
 mu[i]=(beta[1]+beta[4]*x[i]+bvec[1,cluster[i]])/(1+exp(-(tt[i]-(beta[2]+bvec[2,cluster[i]]))/(beta[3]+bvec[3,cluster[i]])));
 w[i] ~ exponential(1);
 pe[i] = (2*w[i]*sigma^2)/(p*(1-p));
 media[i] = (1-2*p)*sigma/(p*(1-p))*w[i]+mu[i]; 
 y[i] ~ normal(media[i],sqrt(pe[i]));   
}
}

