setwd("/Users/nydu/Desktop/stat526proj/project/Hglmer_Lenk/R Code")
########################
# Sleep Study Analysis #
########################
rm(list=ls())
source("hglmer_linear.R")
library(mvtnorm)
library(lme4)
library(MCMCpack)


data(sleepstudy)

# data prep 
subjects = unique(sleepstudy$Subject)
n = length(subjects)
simydata = vector(mode='list', length=n)
simxdata = vector(mode='list', length=n)

for(i in 1:n) {
  subj_data = sleepstudy[sleepstudy$Subject == subjects[i],]
  X = cbind(1, subj_data$Days)
  y = subj_data$Reaction
  
  simxdata[[i]] = X
  simydata[[i]] = y
}

# approach i:
Kvalues = c(2, 3, 4, 5, 6)
set.seed(123) 
results = list()
logML_values = list()
for (K in Kvalues){
  cat("\nFitting mixture model with K =", K, "components\n")
  fit_sleepstudy = gibbs_hglmer_linear(
    ydata = simydata,
    xdata = simxdata,
    K = K,
    a0 = 0,
    d0 = sqrt(10),
    r0 = 1,
    s0 = 1,
    u0 = 0,
    V0 = 100,
    w0 = 1,
    mcmc = list(nblow = 2000, smcmc = 2000)
    )
  results[[K]] = fit_sleepstudy

  cat("Calculating log marginal likelihood...\n")
  logML_value = logML_hglmer_linear(
      post_draws = fit_sleepstudy,
      ydata = simydata,
      xdata = simxdata,
      K = K,
      a0 = 0,
      d0 = sqrt(10),
      r0 = 1,
      s0 = 1,
      u0 = 0,
      V0 = 100,
      w0 = 1
      )
  logML_values[[K]] = logML_value
  cat("Log marginal likelihood for K =", K, ":", logML_value, "\n")

  # plot psi 
  result_psi = do.call(rbind, fit_sleepstudy$psi)
  psi_means = apply(result_psi, 2, mean)
  cat("Estimated mixture proportions:\n")
  print(psi_means)
  
  png(paste0("./plots/sleep_study/psi_plot_K_", K, ".png"), 
    width=800, height=300*K, res=100)
  par(mfrow=c(K, 1), mar=c(3, 3, 2, 1)) 
  for (i in 1:K) {
    plot(result_psi[,i], type='l')
  }
  dev.off()
  
  # plot theta 
  png(paste0("./plots/sleep_study/theta_plot_K_", K, ".png"))
  result_theta = fit_sleepstudy$theta
  theta_list = list()
  for (i in 1:K) {
    theta_list[[i]] = do.call(rbind, lapply(result_theta, function(x) x[i,]))
  }
  means_by_component = lapply(theta_list, function(x) apply(x, 2, mean))
  for (i in 1:K) {
    cat("Component", i, "means (intercept, days):", means_by_component[[i]], "\n")
  }
  par(mfrow=c(K, 2), mar=c(4,4,2,1))
  for (i in 1:K) {
    plot(theta_list[[i]][,1], type='l') # intercept
    plot(theta_list[[i]][,2], type='l') # slope
  }
  dev.off()

  rm(fit_sleepstudy)
  gc()
}

# approach ii
cat("\nFitting mixture model with 2 stage algorithm\n")
fit_sleepstudy_2stage = gibbs_hglmer_linear(
  ydata = simydata,
  xdata = simxdata,
  K = NULL,
  a0 = 0,
  d0 = sqrt(10),
  r0 = 1,
  s0 = 1,
  u0 = 0,
  V0 = 100,
  w0 = 1,
  mcmc = list(nblow = 2000, smcmc = 2000)
)
K = fit_sleepstudy_2stage$K
# plot psi 
result_psi = do.call(rbind, fit_sleepstudy_2stage$psi)
psi_means = apply(result_psi, 2, mean)
cat("Estimated mixture proportions:\n")
print(psi_means)

png(paste0("./plots/sleep_study/psi_plot_approach2.png"), 
  width=800, height=300*K, res=100)
par(mfrow=c(K, 1), mar=c(3, 3, 2, 1)) 
for (i in 1:K) {
  plot(result_psi[,i], type='l')
}
dev.off()

# plot theta 
png(paste0("./plots/sleep_study/theta_plot_approach2.png"))
result_theta = fit_sleepstudy_2stage$theta
theta_list = list()
for (i in 1:K) {
  theta_list[[i]] = do.call(rbind, lapply(result_theta, function(x) x[i,]))
}
means_by_component = lapply(theta_list, function(x) apply(x, 2, mean))
for (i in 1:K) {
  cat("Component", i, "means (intercept, days):", means_by_component[[i]], "\n")
}
par(mfrow=c(K, 2), mar=c(4,4,2,1))
for (i in 1:K) {
  plot(theta_list[[i]][,1], type='l') # intercept
  plot(theta_list[[i]][,2], type='l') # slope
}
dev.off()

# TODO: random effect model 
