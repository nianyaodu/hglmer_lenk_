### Simulate data for Section 5.1 - linear regression case ###
rm(list=ls())
source("hglmer_linear.R")
library(mvtnorm)

n=100;m=10
psi = c(0.2,0.3,0.5)
thetas = list(c(0,0),c(-10,7),c(5,5))
Lambdas = list(matrix(c(1,0,0,1),nrow=2,byrow=T),
               matrix(c(25,9,9,4),nrow=2,byrow=T),
               matrix(c(9,-5,-5,5),nrow=2,byrow=T))
alpha=-1; tau2 = 4

simydata = vector(mode='list', length=n)
simxdata = vector(mode='list', length=n)
simind   = vector(length=n)
simbeta  = matrix(0, nrow=n, ncol=2)

set.seed(876)
for(i in 1:n){
  simind[i] = ind = sample(3,1,prob=psi)
  theta = thetas[[ind]];
  Lambda = Lambdas[[ind]];
  beta = as.vector(rmvnorm(1,mean=theta,sigma=Lambda))
  X = cbind(1,rnorm(m,2,1))
  phi = rnorm(1,alpha,sqrt(tau2))
  y = X%*%beta + rnorm(m,mean=0,sd=sqrt(exp(phi)))
  
  simbeta[i,] = beta
  simxdata[[i]] = X
  simydata[[i]] = y
  
  if(i == 1) print(phi)
}

betas <- list()
for(k in 1:3){
  betas[[k]] <- simbeta[simind==k,]
}
par(mfrow=c(1,1))
plot(betas[[1]][,1],betas[[1]][,2],pch=1,xlim=c(-25,15),ylim=c(-5,15))
points(betas[[2]][,1],betas[[2]][,2],pch=2)
points(betas[[3]][,1],betas[[3]][,2],pch=2)

fit_hglmer_linear <- gibbs_hglmer_linear(ydata=simydata,xdata=simxdata,
                                         K=3,a0=0,d0=sqrt(10),r0=1,s0=1,u0=0,V0=100,w0=1)

#-------------------------------------------------------------------------------
# Results & Visualization

result_psi = do.call(rbind,fit_hglmer_linear$psi)

apply(result_psi,2,mean)

par(mfrow=c(1,3))
plot(result_psi[,1],type='l')
plot(result_psi[,2],type='l')
plot(result_psi[,3],type='l')

result_theta = fit_hglmer_linear$theta
res_theta_1 = do.call(rbind,lapply(result_theta, function(x) x[1,]))
res_theta_2 = do.call(rbind,lapply(result_theta, function(x) x[2,]))
res_theta_3 = do.call(rbind,lapply(result_theta, function(x) x[3,]))

apply(res_theta_1,2,mean)
apply(res_theta_2,2,mean)
apply(res_theta_3,2,mean)

par(mfrow=c(3,2),mar=c(4,4,2,1))
plot(res_theta_1[,1],type='l')
plot(res_theta_1[,2],type='l')

plot(res_theta_2[,1],type='l')
plot(res_theta_2[,2],type='l')

plot(res_theta_3[,1],type='l')
plot(res_theta_3[,2],type='l')

result_phi = fit_hglmer_linear$phi

par(mfrow=c(1,1))
plot(unlist(lapply(result_phi, function(x) x[76])),type='l')



#-------------------------------------------------------------------------------
# Model Selection

Ks = c(1,2,3,4,5,6)
res_gibbs = list()
res_logML = list()

set.seed(3)
for (K in Ks){
  cat("\nFitting mixture model with K =", K, "components\n")
  fit = gibbs_hglmer_linear(ydata=simydata,xdata=simxdata,
                K=K,a0=0,d0=sqrt(10),r0=1,s0=1,u0=0,V0=100,w0=1)
  res_gibbs[[K]] = fit
  
  cat("Calculating log marginal likelihood...\n")
  logML_value = logML_hglmer_linear(post_draws=fit,ydata=simydata,xdata=simxdata,
                K=K,a0=0,d0=sqrt(10),r0=1,s0=1,u0=0,V0=100,w0=1)
  res_logML[[K]] = logML_value
  cat("Log marginal likelihood for K =", K, ":", logML_value, "\n")
}


