# library(devtools)
# install_github('nathanvan/rstanmulticore')

library(rstan)
## The data to analyze (Yes, it is very little!)
schools_dat <- list(
  J = 8, y = c(28,  8, -3,  7, -1,  1, 18, 12),
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

## The Stan model for the data, stored as a string
schools_code <- 'data {
  int J; // number of schools 
  real y[J]; // estimated treatment effects
  real sigma[J]; // s.e. of effect estimates 
}
parameters {
  real mu; 
  real tau;
  real eta[J];
}
transformed parameters {
  real theta[J];
  for (j in 1:J)
    theta[j] <- mu + tau * eta[j];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}'

## The data to analyze (Yes, it is very little!)
schools_dat <- list(
  J = 8, y = c(28,  8, -3,  7, -1,  1, 18, 12),
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

## The Stan model for the data, stored as a string
schools_code <- 'data {
  int J; // number of schools 
  real y[J]; // estimated treatment effects
  real sigma[J]; // s.e. of effect estimates 
}
parameters {
  real mu; 
  real tau;
  real eta[J];
}
transformed parameters {
  real theta[J];
  for (j in 1:J)
    theta[j] <- mu + tau * eta[j];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}'

## Estimating the model 
fit.serial   <- stan( model_code = schools_code, data = schools_dat, 
                      iter = 1000, chains = 4, seed = 1)


library(rstanmulticore)
## Loading required package: parallel

fit.parallel <- pstan( model_code = schools_code, data = schools_dat, 
                       iter = 1000, chains = 4, seed = 1)
