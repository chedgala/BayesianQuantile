# Simulation SKN
library(lqr)
library(rstan)

# Parameters
sigmas<-1.3
beta<-c(-1,2,3)
x_p<- length(beta)

# Set up
M = 100
nvec = c(50,100,200)
pvec = c(50,75,90)/100

total_rows = length(nvec)*length(pvec)*M
total_cols = 2 + 2*(x_p + 1)
out = matrix(NA,nrow = total_rows,ncol = total_cols)

iter = 1


# Testing chains ----------------------------------------------------------

# library(shinystan)
# 
# fit_stan <- stan(file='SKN.stan', init=initf1,
#                  data = list(y=c(y),w=c(y), x_p=x_p,x=x1,
#                              n=nvec[i],p=pvec[j]),
#                  thin = 1000, chains = 1, iter = 10000, warmup = 1000,
#                  seed = 9955,verbose = FALSE)
# 
# launch_shinystan(fit_stan)
# 
# fit_stan <- stan(file='SKN.stan', init=initf1,
#                  data = list(y=c(y),w=c(y), x_p=x_p,x=x1,
#                              n=nvec[i],p=pvec[j]),
#                  chains = 1, iter = 5000, warmup = 1000,
#                  seed = 9955,verbose = FALSE)
# 
# qoi <- c("beta","sigma")
# print(fit_stan, pars=qoi)
# #pairs(fit_stan, pars=qoi)
# 
# acf(rstan::extract(fit_stan)$beta[,1])
# acf(rstan::extract(fit_stan)$beta[,2])
# acf(rstan::extract(fit_stan)$beta[,3])
# plot.ts(rstan::extract(fit_stan)$beta[,3])
# 
# #launch_shinystan(fit_stan1)
# 
# matplot(rstan::extract(fit_stan)$beta,type = "l",lty = c(1,2,3))
# matplot(rstan::extract(fit_stan)$sigma,type = "l")

# Bucle -------------------------------------------------------------------

for(i in seq_along(nvec)){
  
  set.seed(13)
  x1<-cbind(1,runif(nvec[i]),rbinom(nvec[i],1,0.3))
  
  for(j in seq_along(pvec)){
    
    for(k in 1:M){
      
      print(paste0("iteration ", iter, " of ",total_rows))
      
      y <- x1%*%beta + rSKD(nvec[i],mu = 0,sigma = sigmas,
                            p = pvec[j], dist = "normal")
      
      modelN = lqr(y~x1-1,dist = "normal", p=pvec[j],
                   precision = 0.000000001,silent = TRUE)
      
      initf1 <- function() {list(beta=modelN$beta, sigma=modelN$sigma)}
      
      fit_stan <- stan(file='SKN.stan', init=initf1,
                       data = list(y=c(y),w=c(y), x_p=x_p,x=x1,
                                   n=nvec[i],p=pvec[j]),
                       chains = 1, iter = 5000, warmup = 1000,
                       seed = 9955,verbose = FALSE)
      
      out[iter,] = c(n = nvec[i],quantile = pvec[j],
                     beta = colMeans(rstan::extract(fit_stan)$beta),
                     sigma = mean(rstan::extract(fit_stan)$sigma),
                     sdbeta = sqrt(diag(var(rstan::extract(fit_stan)$beta))),
                     sdsigma = sd(rstan::extract(fit_stan)$sigma))
      
      # out[iter,] = c(n = nvec[i],quantile = pvec[j],
      #                beta=modelN$beta,sigma=modelN$sigma,
      #                se = modelN$SE)
      # 
      iter = iter + 1
    }
  }  
}

save(out,file = "simSKN.RData")


outdf <- as.data.frame(out)

library(tidyverse)
colnames(outdf) <- c("n","quantile",
                     "mean_1","mean_2","mean_3","mean_4",
                     "sd_1","sd_2","sd_3","sd_4")

head(outdf)
tail(outdf)

outdf_long <- outdf %>%
  pivot_longer(cols = c("mean_1", "mean_2", "mean_3", "mean_4", "sd_1", "sd_2", "sd_3", "sd_4"),
               names_to = c("type", "parameter"),
               values_to = "value",
               names_pattern = "([a-z]+)_(\\d)") %>%
  mutate(parameter = recode(parameter,
                            "1" = "beta 0",
                            "2" = "beta 1",
                            "3" = "beta 2",
                            "4" = "sigma"))

head(outdf_long)

outdf_long$n <- as.factor(outdf_long$n)
outdf_long$quantile <- as.factor(outdf_long$quantile)

library(ggplot2)

# ggplot(outdf_long, aes(x = parameter, y = value, fill = parameter)) +
#   geom_violin() +
#   labs(title = "Gráfico de violín",
#        x = "parameter",
#        y = "value") +
#   scale_fill_discrete(name = "parameter") +
#     facet_wrap(~type,ncol=2)

ggplot(outdf_long,
       aes(x = parameter, y = value, fill = parameter)) +
  geom_violin() +
  labs(title = "Gráfico de violín",
       x = "parameter",
       y = "value") +
  scale_fill_discrete(name = "parameter") +
  facet_wrap(type~quantile,ncol=length(pvec))

ggplot(outdf_long %>% filter(type == "mean"),
       aes(x = parameter, y = value, fill = parameter)) +
  geom_violin() +
  labs(title = "Gráfico de violín",
       x = "parameter",
       y = "value") +
  scale_fill_discrete(name = "parameter") +
  facet_wrap(~quantile,ncol=length(pvec))


ggplot(outdf_long %>% filter(type == "mean"),
       aes(x = n, y = value, fill = parameter)) +
  geom_violin() +
  labs(title = "Gráfico de violín",
       x = "n",
       y = "value") +
  scale_fill_discrete(name = "parameter") +
  facet_wrap(parameter~quantile,ncol=length(pvec),scales = "free_y")

#
outdf_long_averages <- outdf_long  %>%
  group_by(type,parameter,n,quantile) %>%
  summarise(estimated = mean(value, na.rm = TRUE))


medias <- outdf_long_averages[1:(nrow(outdf_long_averages)/2),]

outdf_long_averages[-(1:(nrow(outdf_long_averages)/2)),"estimated"]^2
medias$par <- rep(c(beta,sigmas),each=length(pvec)*length(nvec))
medias$bias2 <- (medias$estimated - medias$par)^2
vars <- outdf_long_averages[-(1:(nrow(outdf_long_averages)/2)),"estimated"]^2
medias$mse <- unlist(medias$bias2) + unlist(vars)

ggplot(medias,mapping = aes(x = n,y = bias2,fill = parameter)) +
  geom_point() +
  geom_line(aes(color = parameter,group = parameter)) +
  scale_fill_discrete(name = "parameter") +
  facet_wrap(parameter~quantile,ncol=length(pvec),scales = "free_y")

ggplot(medias,mapping = aes(x = n,y = mse,fill = parameter)) +
  geom_point() +
  geom_line(aes(color = parameter,group = parameter)) +
  scale_fill_discrete(name = "parameter") +
  facet_wrap(parameter~quantile,ncol=length(pvec),scales = "free_y")
