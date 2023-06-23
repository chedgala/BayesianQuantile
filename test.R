# Instalar R 4.2
# Instalar Rtools

#Sys.getenv("BINPREF") # debe ser igual a ""
#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# Sys.which("make") #probar ruta a Rtools

# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# 
# example(stan_model, package = "rstan", run.dontrun = TRUE)

#install.packages("rstudioapi")


library("rstan") # observe startup messages
library("rstudioapi")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
fit <- stan(file = 'schools.stan', data = schools_dat)

print(fit)
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))

la <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- la$mu 

### return an array of three dimensions: iterations, chains, parameters 
a <- extract(fit, permuted = FALSE)

### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)
