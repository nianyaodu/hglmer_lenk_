setwd("/Users/nydu/Desktop/stat526proj/project/Hglmer_Lenk/R Code")
#########################
# InstEval Data Analysis #
#########################
rm(list=ls())
source("hglmer_logistic.R")
library(mvtnorm)
library(lme4)
library(MCMCpack)

data(InstEval)

# sample student with sufficient data 
min_ratings = 5  # minimum number of ratings per student
students = names(which(table(InstEval$s) >= min_ratings))
# select a sample if more than 100 
set.seed(123)
if(length(students) > 100) {
  students = sample(students, 100)
}

# data prep 
simydata = vector(mode='list', length=0)
simxdata = vector(mode='list', length=0)
student_ids = c()

for(i in 1:length(students)) {
  student_id = students[i]
  student_data = InstEval[InstEval$s == student_id,]
  
  # skip students with too few observations or only one service type
  if(nrow(student_data) < min_ratings) next
  if(length(unique(student_data$service)) < 2) next
  
  # convert to binary outcome (ratings >= 4 means positive)
  y = as.numeric(student_data$y >= 4)
  
  X = cbind(1, as.numeric(student_data$service == "1"))
  simydata = c(simydata, list(y))
  simxdata = c(simxdata, list(X))
  student_ids = c(student_ids, student_id)
}
n = length(simydata)
cat("Analysis will use", n, "students\n")


u0 = 0 
V0 = 100 
w0 = 1 

Kvalues = c(1, 2, 3, 4) 
set.seed(123)
results = list()
logML_values = list()

for (K in Kvalues) {
  cat("\nFitting mixture model with K =", K, "components\n")
  fit_insteval = gibbs_hglmer_logistic(
    ydata = simydata,
    xdata = simxdata,
    K = K,
    u0 = u0,
    V0 = V0,
    w0 = w0,
    mcmc = list(nblow = 2000, smcmc = 2000)
  )
  results[[K]] = fit_insteval

  
  # plot psi 
  result_psi = do.call(rbind, fit_insteval$psi)
  psi_means = apply(result_psi, 2, mean)
  cat("Estimated mixture proportions:\n")
  print(psi_means)
  
  png(paste0("./plots/insteval/psi_plot_K_", K, ".png"),
    width=800, height=300*K, res=100)
  par(mfrow=c(K, 1), mar=c(3, 3, 2, 1)) 
  for (i in 1:K) {
    plot(result_psi[,i], type='l', main=paste("Component", i),
         xlab="Iteration", ylab="Mixture Probability")
  }
  dev.off()
  
  # plot thetas 
  png(paste0("./plots/insteval/theta_plot_K_", K, ".png"))
  
  result_theta = fit_insteval$theta
  theta_list = list()
  for (i in 1:K) {
    theta_list[[i]] = do.call(rbind, lapply(result_theta, function(x) x[i,]))
  }
  means_by_component = lapply(theta_list, function(x) apply(x, 2, mean))
  for (i in 1:K) {
    cat("Component", i, "means (intercept, service):", means_by_component[[i]], "\n")
  }
  
  par(mfrow=c(K, 2), mar=c(4,4,2,1))
  for (i in 1:K) {
    plot(theta_list[[i]][,1], type='l', main=paste("Component", i, "Intercept"),
         xlab="Iteration", ylab="Value")  # intercept
    plot(theta_list[[i]][,2], type='l', main=paste("Component", i, "Service Effect"),
         xlab="Iteration", ylab="Value")  # effect
  }
  dev.off()
  
  rm(fit_insteval)
  gc()
}

# approach ii
cat("\nFitting mixture model with 2 stage algorithm\n")
fit_insteval_2stage = gibbs_hglmer_logistic(
    ydata = simydata,
    xdata = simxdata,
    K = NULL,
    u0 = u0,
    V0 = V0,
    w0 = w0,
    mcmc = list(nblow = 2000, smcmc = 2000)
  )

K = fit_insteval_2stage$K
# plot psi 
result_psi = do.call(rbind, fit_insteval_2stage$psi)
psi_means = apply(result_psi, 2, mean)
cat("Estimated mixture proportions:\n")
print(psi_means)

png(paste0("./plots/insteval/psi_plot_approach2.png"),
  width=800, height=300*K, res=100)
par(mfrow=c(K, 1), mar=c(3, 3, 2, 1)) 
for (i in 1:K) {
  plot(result_psi[,i], type='l', main=paste("Component", i),
        xlab="Iteration", ylab="Mixture Probability")
}
dev.off()

# plot thetas 
png(paste0("./plots/insteval/theta_plot_approach2.png"))

result_theta = fit_insteval_2stage$theta
theta_list = list()
for (i in 1:K) {
  theta_list[[i]] = do.call(rbind, lapply(result_theta, function(x) x[i,]))
}
means_by_component = lapply(theta_list, function(x) apply(x, 2, mean))
for (i in 1:K) {
  cat("Component", i, "means (intercept, service):", means_by_component[[i]], "\n")
}

par(mfrow=c(K, 2), mar=c(4,4,2,1))
for (i in 1:K) {
  plot(theta_list[[i]][,1], type='l', main=paste("Component", i, "Intercept"),
        xlab="Iteration", ylab="Value")  # intercept
  plot(theta_list[[i]][,2], type='l', main=paste("Component", i, "Service Effect"),
        xlab="Iteration", ylab="Value")  # effect
}
dev.off()

# comparison with baseline
insteval_subset = InstEval[InstEval$s %in% student_ids,]
insteval_subset$binary_y = as.numeric(insteval_subset$y >= 4)

standard_model = glmer(binary_y ~ service + (service|s), 
                       data=insteval_subset, family=binomial)
print(summary(standard_model))

# TODO: need to test for k = 1
# TODO: random effect model
## glmm
## brms - stan

# TODO: fix poisson model and update bike sharing example