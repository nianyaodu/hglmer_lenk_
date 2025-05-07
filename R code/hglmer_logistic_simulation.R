### Simulate data for Section 5.1 - logistic regression case ###
rm(list=ls())
source("hglmer_logistic.R")
library(boot)

n=100;m=40
psi = c(0.4,0.6)
thetas = list(c(0,-1),c(-1,1))
Lambdas = list(matrix(c(0.01,0,0,0.01),nrow=2,byrow=T),
                matrix(c(0.04,-0.028,-0.028,0.04),nrow=2,byrow=T))

simydata = vector(mode='list', length=n)
simxdata = vector(mode='list', length=n)
simind   = vector(length=n)
simbeta  = matrix(0, nrow=n, ncol=2)

set.seed(856)
for(i in 1:n){
  simind[i] = ind = sample(2,1,prob=psi)
  theta = thetas[[ind]];
  Lambda = Lambdas[[ind]];
  beta = as.vector(rmvnorm(1,mean=theta,sigma=Lambda))
  X = cbind(1,rnorm(m))
  y = rbinom(m,1,inv.logit(X%*%beta))
  
  simbeta[i,] = beta
  simxdata[[i]] = X
  simydata[[i]] = y
}

betas <- list()
for(k in 1:2){
  betas[[k]] <- simbeta[simind==k,]
}
par(mfrow=c(1,1))
plot(betas[[1]][,1],betas[[1]][,2],pch=1,xlim=c(-1.5,0.5),ylim=c(-2,2))
points(betas[[2]][,1],betas[[2]][,2],pch=2)

fit_hglmer_logistic <- gibbs_hglmer_logistic(ydata=simydata,xdata=simxdata,
                                             K=2,u0=0,V0=100,w0=1)

#-------------------------------------------------------------------------------
# Results & Visualization

result_psi = do.call(rbind,fit_hglmer_logistic$psi)

apply(result_psi,2,mean)

par(mfrow=c(1,2))
plot(result_psi[,1],type='l')
plot(result_psi[,2],type='l')

result_theta = fit_hglmer_logistic$theta
res_theta_1 = do.call(rbind,lapply(result_theta, function(x) x[1,]))
res_theta_2 = do.call(rbind,lapply(result_theta, function(x) x[2,]))

apply(res_theta_1,2,mean)
apply(res_theta_2,2,mean)

par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(res_theta_1[,1],type='l')
plot(res_theta_1[,2],type='l')

plot(res_theta_2[,1],type='l')
plot(res_theta_2[,2],type='l')









