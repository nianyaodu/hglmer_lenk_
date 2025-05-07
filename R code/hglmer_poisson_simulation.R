setwd("/Users/nydu/Desktop/stat526proj/project/Hglmer_Lenk/R Code")
### Simulate data for Poisson regression case ###
rm(list=ls())
source("hglmer_poisson.R")
library(mvtnorm)

# simulation params
n = 100 # number of subjects
m = 10 # number of observations per subject
psi = c(0.4, 0.6)  # mixing 
thetas = list(c(-0.5, 0.2), c(1.5, 0.8)) # means
Lambdas = list(matrix(c(0.1, 0, 0, 0.1), nrow=2, byrow=T),
               matrix(c(0.2, 0.05, 0.05, 0.2), nrow=2, byrow=T))

simydata = vector(mode='list', length=n)
simxdata = vector(mode='list', length=n)
simind = vector(length=n)
simbeta = matrix(0, nrow=n, ncol=2)

set.seed(123)
for(i in 1:n) {
  simind[i] = ind = sample(2, 1, prob=psi)
  
  # subject-specific parameters from the component
  theta = thetas[[ind]]
  Lambda = Lambdas[[ind]]
  beta = as.vector(rmvnorm(1, mean=theta, sigma=Lambda))
  X = cbind(1, runif(m, -1, 1))
  log_mu = X %*% beta
  mu = exp(log_mu)
  y = rpois(m, lambda=mu)
  simbeta[i,] = beta
  simxdata[[i]] = X
  simydata[[i]] = y
}
betas = list()
for(k in 1:2) {
  betas[[k]] = simbeta[simind==k,]
}

par(mfrow=c(1,1))
plot(betas[[1]][,1], betas[[1]][,2], pch=1, col="red",
     xlim=c(-1.5, 2.5), ylim=c(-0.5, 1.5), 
     xlab="Intercept", ylab="Slope",
     main="True Regression Coefficients")
points(betas[[2]][,1], betas[[2]][,2], pch=2, col="blue")
points(thetas[[1]][1], thetas[[1]][2], pch=16, cex=2, col="red")
points(thetas[[2]][1], thetas[[2]][2], pch=17, cex=2, col="blue")
text(thetas[[1]][1], thetas[[1]][2], paste0("(", round(thetas[[1]][1], 2), 
                                           ",", round(thetas[[1]][2], 2), ")"), 
     pos=3, cex=1.5)
text(thetas[[2]][1], thetas[[2]][2], paste0("(", round(thetas[[2]][1], 2), 
                                           ",", round(thetas[[2]][2], 2), ")"), 
     pos=1, cex=1.5)
legend("topright", c("Component 1", "Component 2"), 
       pch=c(1, 2), col=c("red", "blue"))
fit_hglmer_poisson = gibbs_hglmer_poisson(
  ydata = simydata,
  xdata = simxdata,
  K = 2, # if set to null will automatically determine k 
  u0 = 0,
  V0 = 100,
  w0 = 1,
  mcmc = list(nblow = 1000, smcmc = 1000)
)
result_psi = do.call(rbind, fit_hglmer_poisson$psi)
cat("Estimated mixture proportions:\n")
print(apply(result_psi, 2, mean))
par(mfrow=c(1,2))
plot(result_psi[,1], type='l', main="Component 1 Proportion", 
     xlab="Iteration", ylab="Mixture Probability")
plot(result_psi[,2], type='l', main="Component 2 Proportion",
     xlab="Iteration", ylab="Mixture Probability")

# thetas 
result_theta = fit_hglmer_poisson$theta
res_theta_1 = do.call(rbind, lapply(result_theta, function(x) x[1,]))
res_theta_2 = do.call(rbind, lapply(result_theta, function(x) x[2,]))

cat("\nComponent 1 means (intercept, slope):\n")
print(apply(res_theta_1, 2, mean))
cat("\nComponent 2 means (intercept, slope):\n")
print(apply(res_theta_2, 2, mean))

par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(res_theta_1[,1], type='l', main="Component 1 Intercept",
     xlab="Iteration", ylab="Value")
plot(res_theta_1[,2], type='l', main="Component 1 Slope",
     xlab="Iteration", ylab="Value")

plot(res_theta_2[,1], type='l', main="Component 2 Intercept",
     xlab="Iteration", ylab="Value")
plot(res_theta_2[,2], type='l', main="Component 2 Slope",
     xlab="Iteration", ylab="Value")
true_theta1 = thetas[[1]]
true_theta2 = thetas[[2]]
est_theta1 = apply(res_theta_1, 2, mean)
est_theta2 = apply(res_theta_2, 2, mean)

cat("\nParameter Comparison:\n")
cat("Component 1 - True:", true_theta1, "Estimated:", est_theta1, "\n")
cat("Component 2 - True:", true_theta2, "Estimated:", est_theta2, "\n")


