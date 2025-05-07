#########################
# Bike Sharing Analysis #
#########################
rm(list=ls())
source("hglmer_poisson.R")
library(MASS)
library(MCMCpack)
library(mclust)
library(lme4)

calculate_pred_logscore = function(fit, ydata, xdata, K) {
  cat("Calculating predictive log score for K =", K, "components...\n")
  
  thetas = fit$theta
  Lambdas = fit$Lambda
  psis = fit$psi
  n_samples = length(thetas)
  
  n = length(ydata)
  log_scores = numeric(n)
  
  for (i in 1:n) {
    y_test = ydata[[i]]
    x_test = xdata[[i]]
    
    sample_loglik = numeric(n_samples)
    
    for (s in 1:n_samples) {
      theta_s = thetas[[s]]
      Lambda_s = Lambdas[[s]]
      psi_s = psis[[s]]
      
      mix_pred_dens = 0
      
      for (k in 1:K) {
        beta_pred = mvrnorm(1, mu=theta_s[k,], Sigma=Lambda_s[,,k])
        
        eta = x_test %*% beta_pred
        log_lik = sum(y_test * eta - exp(eta) - lfactorial(y_test))

        mix_pred_dens = mix_pred_dens + psi_s[k] * exp(log_lik)
      }
      
      sample_loglik[s] = log(mix_pred_dens + 1e-10) 
    }
    
    log_scores[i] = mean(sample_loglik)
  }
  
  mean_logscore = mean(log_scores)
  cat("Mean predictive log score:", mean_logscore, "\n")
  return(mean_logscore)
}




bike_data = read.csv(".././data/bike+sharing+dataset/hour.csv")

cat("Dataset structure:\n")
str(bike_data)

days = unique(bike_data$dteday)

orig_n = length(days)
cat("Number of unique days:", orig_n, "\n")

simydata = vector(mode='list', length=0)
simxdata = vector(mode='list', length=0)
days_kept = c()

cat("Starting data filtering and preparation...\n")
for(i in 1:length(days)) {
  day_data = bike_data[bike_data$dteday == days[i],]
  
  if(nrow(day_data) < 8) {
    if(i %% 50 == 0) cat("Skipping day", i, "due to insufficient hours (<8):", days[i], "\n")
    next
  }
  
  y = day_data$cnt
  
  if(var(y) < 0.1) {
    if(i %% 50 == 0) cat("Skipping day", i, "due to zero variance in counts:", days[i], "\n")
    next
  }
  
  X = cbind(1, 
            as.numeric(scale(day_data$temp, center=TRUE, scale=TRUE)), 
            as.numeric(scale(day_data$hum, center=TRUE, scale=TRUE)), 
            as.numeric(scale(day_data$windspeed, center=TRUE, scale=TRUE)))
  
  if(any(is.na(X)) || any(is.infinite(X)) || any(is.nan(X))) {
    if(i %% 50 == 0) cat("Skipping day", i, "due to NA/Inf values in predictors:", days[i], "\n")
    next
  }
  
  max_abs_val = max(abs(X[,-1]))
  if(max_abs_val > 5) {
    if(i %% 50 == 0) cat("Skipping day", i, "due to extreme predictor values (", max_abs_val, "):", days[i], "\n")
    next
  }
  
  simxdata = c(simxdata, list(X))
  simydata = c(simydata, list(y))
  days_kept = c(days_kept, days[i])
  
  if(i %% 50 == 0) cat("Processed", i, "days, kept", length(simydata), "so far\n")
}

n = length(simydata)
cat("Number of usable days after filtering:", n, "\n")

cat("Standardizing observation lengths...\n")
max_obs = max(sapply(simydata, length))
min_obs = min(sapply(simydata, length))
cat("Current observation lengths - min:", min_obs, "max:", max_obs, "\n")

standard_length = 24 

std_simydata = vector(mode='list', length=n)
std_simxdata = vector(mode='list', length=n)

for(i in 1:n) {
  curr_length = length(simydata[[i]])
  
  if(curr_length > standard_length) {
    std_simydata[[i]] = simydata[[i]][(curr_length-standard_length+1):curr_length]
    std_simxdata[[i]] = simxdata[[i]][(curr_length-standard_length+1):curr_length,]
  } else if(curr_length < standard_length) {
    cat("Removing day", i, "with insufficient observations:", curr_length, "\n")
    std_simydata[[i]] = NULL
    std_simxdata[[i]] = NULL
  } else {
    std_simydata[[i]] = simydata[[i]]
    std_simxdata[[i]] = simxdata[[i]]
  }
}

valid_indices = which(!sapply(std_simydata, is.null))
std_simydata = std_simydata[valid_indices]
std_simxdata = std_simxdata[valid_indices]
days_kept = days_kept[valid_indices]

n = length(std_simydata)
cat("Dataset size after standardization:", n, "\n")

obs_lengths = sapply(std_simydata, length)
cat("All observations now have length", unique(obs_lengths), "\n")

simydata = std_simydata
simxdata = std_simxdata

cat("Running verification check on beta values...\n")
problem_days = c()
for(i in 1:n) {
  y_i = simydata[[i]]
  x_i = simxdata[[i]]
  
  if(min(y_i) == 0 && mean(y_i) < 5) {
    problem_days = c(problem_days, i)
    cat("Potential problem with day", i, ": low counts with zeros\n")
  }
  
  var_values = apply(x_i[, -1, drop=FALSE], 2, var)
  if(any(var_values < 0.01)) {
    problem_days = c(problem_days, i)
    cat("Potential problem with day", i, ": low variance in predictors\n")
  }
}

if(length(problem_days) > 0) {
  cat("Found", length(problem_days), "potentially problematic days. Consider removing if errors persist.\n")
  cat("Removing problematic days from dataset...\n")
  
  for(i in sort(problem_days, decreasing=TRUE)) {
    simydata = simydata[-i]
    simxdata = simxdata[-i]
    days_kept = days_kept[-i]
  }
  
  n = length(simydata)
  cat("Dataset size after removing problematic days:", n, "\n")
}

cat("Performing additional data quality checks...\n")
to_remove = c()
for(i in 1:n) {
  y_i = simydata[[i]]
  x_i = simxdata[[i]]
  
  skew = mean((y_i - mean(y_i))^3) / (sd(y_i)^3)
  if(abs(skew) > 5) {
    cat("Day", i, "has extreme skewness in counts:", skew, "\n")
    to_remove = c(to_remove, i)
    next
  }
  for(j in 2:ncol(x_i)) {
    corr = cor(x_i[,j], y_i)
    if(abs(corr) > 0.95) {
      cat("Day", i, "has near-perfect correlation between predictor", j, "and response:", corr, "\n")
      to_remove = c(to_remove, i)
      break
    }
  }
  
  xcov = cov(x_i[,-1])
  eigen_vals = eigen(xcov)$values
  if(any(eigen_vals < 0.001) || (max(eigen_vals)/min(eigen_vals) > 1000)) {
    cat("Day", i, "has potentially unstable covariance matrix\n")
    to_remove = c(to_remove, i)
  }
}

if(length(to_remove) > 0) {
  cat("Found", length(to_remove), "additional problematic days to remove\n")
  for(i in sort(unique(to_remove), decreasing=TRUE)) {
    simydata = simydata[-i]
    simxdata = simxdata[-i]
    days_kept = days_kept[-i]
  }
  n = length(simydata)
  cat("Dataset size after additional filtering:", n, "\n")
}

count_means = sapply(simydata, mean)
count_vars = sapply(simydata, var)
cat("Response variable summary:\n")
cat("Min mean count:", min(count_means), "Max mean count:", max(count_means), "\n")
cat("Min variance:", min(count_vars), "Max variance:", max(count_vars), "\n")
cat("Variance/Mean ratio range:", min(count_vars/count_means), "-", max(count_vars/count_means), "\n")

u0 = 0  
V0 = 100 
w0 = 1 

p = ncol(simxdata[[1]])
cat("Number of predictors (p):", p, "\n")

beta_init = matrix(0, nrow=n, ncol=p)
regression_failures = 0
cat("Starting beta initialization...\n")

for (i in 1:n) {
  y_i = simydata[[i]]
  x_i = simxdata[[i]]
  
  if (max(y_i) > 1000 || min(y_i) < 0) {
    cat("Warning: Extreme count values for day", i, ": min =", min(y_i), ", max =", max(y_i), "\n")
  }
  
  log_y = log(y_i + 0.1)
  
  tryCatch({
    fit = lm(log_y ~ x_i[,-1])
    if (!is.null(fit$coefficients) && length(fit$coefficients) == p && 
        sum(is.na(fit$coefficients)) == 0) {
      beta_init[i,] = as.numeric(fit$coefficients)
    } else {
      regression_failures = regression_failures + 1
      cat("Warning: Regression produced invalid coefficients for day", i, "\n")
      beta_init[i,1] = log(mean(y_i) + 0.1)
      beta_init[i,2:p] = 0
    }
  }, error = function(e) {
    regression_failures = regression_failures + 1
    cat("Error in regression for day", i, ":", conditionMessage(e), "\n")
    beta_init[i,1] = log(mean(y_i) + 0.1)
    beta_init[i,2:p] = 0
  })
  
  if (i %% 20 == 0) {
    cat("Processed", i, "beta initializations\n")
  }
}

cat("Total regression failures:", regression_failures, "out of", n, "days\n")

na_count_before = sum(is.na(beta_init))
inf_count_before = sum(is.infinite(beta_init))
extreme_count_before = sum(abs(beta_init) > 10)

if (na_count_before > 0) cat("Found", na_count_before, "NA values in beta_init\n")
if (inf_count_before > 0) cat("Found", inf_count_before, "Inf values in beta_init\n")
if (extreme_count_before > 0) cat("Found", extreme_count_before, "extreme values in beta_init\n")

beta_init[is.na(beta_init)] = 0
beta_init[is.infinite(beta_init)] = 0
beta_init[abs(beta_init) > 10] = sign(beta_init[abs(beta_init) > 10]) * 10

cat("Beta initialization summary after cleaning:\n")
print(summary(as.vector(beta_init)))

cat("Checking numerical stability of beta covariance...\n")
beta_cov = cov(beta_init)
eigen_vals = eigen(beta_cov)$values
cat("Eigenvalues of beta covariance matrix:", eigen_vals, "\n")
if (any(eigen_vals < 1e-10)) {
  cat("Warning: Near-zero eigenvalues detected in beta covariance matrix\n")
}

Kvalues = c(2, 3, 4, 5) 
results = list()
logML_values = list()

patch_hglmer_poisson = function() {
  if (!exists("original_gibbs_hglmer_poisson")) {
    original_gibbs_hglmer_poisson <<- gibbs_hglmer_poisson
    
    patched_gibbs_hglmer_poisson <- function(...) {
      args <- list(...)
      
      ydata <- args$ydata
      xdata <- args$xdata
      K <- args$K
      u0 <- args$u0
      V0 <- args$V0
      w0 <- args$w0
      mcmc <- args$mcmc
      show_pbar <- if("show_pbar" %in% names(args)) args$show_pbar else TRUE
      init <- args$init

      n = length(ydata)
      p = ncol(xdata[[1]])
      
      print(paste("Number of subjects:", n))
      print(paste("Number of predictors:", p))
      
      mcvals = list(nblow=1000, smcmc=1000)
      mcvals[names(mcmc)] = mcmc
      nblow = mcvals$nblow 
      smcmc = mcvals$smcmc  
      nmcmc = nblow + smcmc 
      
      f0 = p + 1
      G0 = p*diag(rep(1,p))
      
      if (!is.null(init)) {
        beta = init$beta
        z = init$z
        theta = init$theta
        Lambda = init$Lambda
        psi = init$psi
        beta_hat = beta 
        
        cat("Using user-provided initialization\n")
      } else {
        beta = matrix(0, nrow=n, ncol=p)
        beta_hat = matrix(0, nrow=n, ncol=p)
        
        for (i in 1:n) {
          y_i = ydata[[i]]
          x_i = xdata[[i]]
          
          df = data.frame(y = y_i)
          for(j in 2:p) {
            df[[paste0("x", j-1)]] = x_i[,j]
          }
          
          if(p > 2) {
            formula_str = paste("y ~", paste(paste0("x", 1:(p-1)), collapse=" + "))
            formula_obj = as.formula(formula_str)
          } else {
            formula_obj = as.formula("y ~ x1")
          }
          
          fit_result = tryCatch({
            fit_poisson = glm(formula_obj, data=df, family=poisson())
            coef(fit_poisson)
          }, error = function(e) {
            c(log(mean(y_i) + 0.1), rep(0, p-1))
          })
          
          beta[i,] = beta_hat[i,] = fit_result
        }
        
        init_clusters = kmeans(beta, centers=K, nstart=10)
        z = init_clusters$cluster
        
        theta = matrix(0, nrow=K, ncol=p)
        Lambda = array(0, dim=c(p, p, K))
        psi = numeric(K)
        
        for (k in 1:K) {
          idx = which(z == k)
          if (length(idx) > 0) {
            theta[k,] = colMeans(beta[idx, , drop=FALSE])
            if (length(idx) > 1) {
              Lambda[,,k] = cov(beta[idx, , drop=FALSE])
              eigen_vals = eigen(Lambda[,,k])$values
              if (any(eigen_vals <= 0)) {
                Lambda[,,k] = Lambda[,,k] + diag(0.1 + abs(min(eigen_vals)), p)
              }
            } else {
              Lambda[,,k] = diag(0.5, p)
            }
            psi[k] = length(idx) / n
          } else {
            theta[k,] = colMeans(beta)
            Lambda[,,k] = diag(0.5, p)
            psi[k] = 1/K
          }
        }
        
        ord = order(psi)
        psi = psi[ord]
        theta = theta[ord,]
        Lambda_temp = Lambda
        for (k in 1:K) {
          Lambda[,,k] = Lambda_temp[,,ord[k]]
        }
        z = match(z, ord)
      }
      
      beta.save = list()
      beta.hat.save = list()
      z.save = list()
      theta.save = list()
      Lambda.save = list()
      psi.save = list()
      
      cat(sprintf("Generating hglmer posterior samples \n"))
      if(show_pbar) { pb <- txtProgressBar(min = 1, max = nmcmc, style=3) }
      
      for (iter in 1:nmcmc) {
        for (i in 1:n) {
          lpsi_i = numeric(K)
          for (k in 1:K) {
            if (psi[k] <= 0 || is.na(psi[k])) {
              psi[k] = 1e-10
            }
            
            tryCatch({
              lpsi_i[k] = log(psi[k]) + dmvnorm(beta[i,], mean=theta[k,], sigma=Lambda[,,k], log=TRUE)
            }, error = function(e) {
              lpsi_i[k] = -1000
            })
          }
          
          lpsi_i[is.na(lpsi_i) | is.infinite(lpsi_i)] = -1000
          
          lpsi_max = max(lpsi_i)
          lpsi_shifted = lpsi_i - lpsi_max
          psi_i = exp(lpsi_shifted)
          
          if (sum(psi_i) == 0) {
            psi_i = rep(1/K, K)
          } else {
            psi_i = psi_i / sum(psi_i)
          }
          
          if (any(is.na(psi_i))) {
            cat("WARNING: NA values in probabilities for subject", i, "- using uniform distribution\n")
            psi_i = rep(1/K, K)
          }
          
          z[i] = sample(K, 1, prob=psi_i)
        }
        
        for (i in 1:n) {
          y_i = ydata[[i]]
          x_i = xdata[[i]]
          beta_hat_i = beta_hat[i,]
          
          if (z[i] < 1 || z[i] > K || is.na(z[i])) {
            cat("WARNING: Invalid z value for subject", i, ":", z[i], "- resetting to 1\n")
            z[i] = 1
          }
          
          Lamb_inv = tryCatch({
            solve(Lambda[,,z[i]])
          }, error = function(e) {
            cat("WARNING: Error solving Lambda for subject", i, "component", z[i], "\n")
            diag(1, p) 
          })
          
          log_mu = numeric(length(y_i))
          mu = numeric(length(y_i))
          
          for (j in 1:length(y_i)) {
            log_mu[j] = sum(x_i[j,] * beta_hat_i)
            if (log_mu[j] > 20) log_mu[j] = 20 
            mu[j] = exp(log_mu[j])
          }
          

          V_i_sum = matrix(0, nrow=p, ncol=p)
          for (j in 1:length(mu)) {
            w_j = mu[j]
            x_j = x_i[j,]
            V_i_sum = V_i_sum + w_j * (x_j %*% t(x_j))
          }
          

          V_i_sum = V_i_sum + diag(1e-6, p)
          
          V_i = tryCatch({
            solve(V_i_sum + Lamb_inv)
          }, error = function(e) {
            cat("WARNING: Error solving V_i for subject", i, "\n")
            diag(0.1, p) 
          })
          
          resids = y_i - mu
          grad_llik_i = numeric(p)
          for (j in 1:length(resids)) {
            grad_llik_i = grad_llik_i + resids[j] * x_i[j,]
          }
          
          grad_prior_i = -Lamb_inv %*% (beta_hat_i - theta[z[i],])
          grad_i = grad_llik_i + grad_prior_i
          
          beta_hat[i,] = as.vector(beta_hat_i + V_i %*% grad_i)
          
          beta_i_prop = as.vector(mvrnorm(n=1, mu=beta_hat[i,], Sigma=V_i))
          
          log_mu_prop = numeric(length(y_i))
          log_post_prop = 0
          
          for (j in 1:length(y_i)) {
            log_mu_prop[j] = sum(x_i[j,] * beta_i_prop)
            if (log_mu_prop[j] > 20) log_mu_prop[j] = 20
            log_post_prop = log_post_prop + (y_i[j] * log_mu_prop[j] - exp(log_mu_prop[j]))
          }
          
          prior_term_prop = -0.5 * t(beta_i_prop - theta[z[i],]) %*% Lamb_inv %*% (beta_i_prop - theta[z[i],])
          proposal_term_prop = 0.5 * t(beta_i_prop - beta_hat[i,]) %*% solve(V_i) %*% (beta_i_prop - beta_hat[i,])
          log_post_prop = log_post_prop + prior_term_prop + proposal_term_prop
          

          log_mu_curr = numeric(length(y_i))
          log_post_curr = 0
          
          for (j in 1:length(y_i)) {
            log_mu_curr[j] = sum(x_i[j,] * beta[i,])
            if (log_mu_curr[j] > 20) log_mu_curr[j] = 20
            log_post_curr = log_post_curr + (y_i[j] * log_mu_curr[j] - exp(log_mu_curr[j]))
          }
          
          prior_term_curr = -0.5 * t(beta[i,] - theta[z[i],]) %*% Lamb_inv %*% (beta[i,] - theta[z[i],])
          proposal_term_curr = 0.5 * t(beta[i,] - beta_hat[i,]) %*% solve(V_i) %*% (beta[i,] - beta_hat[i,])
          log_post_curr = log_post_curr + prior_term_curr + proposal_term_curr
          
          if (is.na(log_post_prop) || is.infinite(log_post_prop)) log_post_prop = -1e10
          if (is.na(log_post_curr) || is.infinite(log_post_curr)) log_post_curr = -1e10
          
          accprob = min(0, log_post_prop - log_post_curr)
          
          if (is.na(accprob) || is.infinite(accprob)) {
            cat("WARNING: Invalid acceptance probability for subject", i, "\n")
          } else {
            if(log(runif(1)) < accprob) {
              beta[i,] = beta_i_prop
            }
          }
        }
        
        for (k in 1:K) {
          n_k = sum(z==k)
          if (n_k > 0) { 
            Vn_k = tryCatch({
              solve(n_k*solve(Lambda[,,k]) + solve(V0*diag(rep(1,p))))
            }, error = function(e) {
              cat("WARNING: Error solving Vn_k for component", k, "\n")
              diag(1, p) 
            })
            
            subjects_in_k = which(z==k)
            beta_sum = rep(0, p)
            for (i in subjects_in_k) {
              beta_sum = beta_sum + beta[i,]
            }
            
            un_k = Vn_k %*% (solve(Lambda[,,k]) %*% beta_sum + 
                            solve(V0*diag(rep(1,p))) %*% rep(u0,p))
            
            theta[k,] = as.vector(mvrnorm(n=1, mu=un_k, Sigma=Vn_k))
          } else {
            theta[k,] = rep(u0, p) + rnorm(p, 0, sqrt(V0))
          }
        }
        
        for (k in 1:K) {
          n_k = sum(z==k)
          if (n_k > 0) {
            fn_k = f0 + n_k
            M = sweep(beta[z==k,,drop=FALSE], 2, theta[k,], `-`)
            Gn_k = G0 + t(M) %*% M
            
            min_eig = min(eigen(Gn_k)$values)
            if (min_eig <= 0) {
              Gn_k = Gn_k + diag(0.1 + abs(min_eig), p)
            }
            
            Lambda[,,k] = riwish(fn_k, Gn_k)
          } else {
            Lambda[,,k] = riwish(f0, G0)
          }
          
          eigen_vals = eigen(Lambda[,,k])$values
          if (any(eigen_vals <= 0)) {
            cat("WARNING: Lambda matrix for component", k, "not positive definite after update\n")
            Lambda[,,k] = Lambda[,,k] + diag(0.1 + abs(min(eigen_vals)), p)
          }
        }
        
        n_K = numeric(K)
        for (k in 1:K) { n_K[k] = sum(z==k) }
        shape = rep(w0, K) + n_K
        
        psi_new = tryCatch({
        as.vector(rng_orderdiri(1, shape))
        }, error = function(e) {
        cat("WARNING: Error in rng_orderdiri:", conditionMessage(e), "\n")
        temp = MCMCpack::rdirichlet(1, shape)
        sort(temp)
        })

        if (any(is.na(psi_new)) || any(psi_new <= 0) || abs(sum(psi_new) - 1) > 0.01) {
        cat("WARNING: Invalid psi values generated:", psi_new, "\n")
        } else {
        psi = psi_new / sum(psi_new)
        }
        
        if (iter > nblow) {
          beta.save[[iter-nblow]] = beta
          z.save[[iter-nblow]] = z
          theta.save[[iter-nblow]] = theta
          Lambda.save[[iter-nblow]] = Lambda
          psi.save[[iter-nblow]] = psi
          beta.hat.save[[iter-nblow]] = beta_hat
        }
        
        if (show_pbar) { setTxtProgressBar(pb, iter) }
      }
      
      cat(sprintf("\nMCMC is done! \n"))
      
      out = NULL
      out$beta = beta.save
      out$beta_hat = beta.hat.save
      out$theta = theta.save
      out$Lambda = Lambda.save
      out$psi = psi.save
      out$z = z.save
      
      return(out)
    }
    
    gibbs_hglmer_poisson <<- patched_gibbs_hglmer_poisson
    cat("Patched the hglmer_poisson function to improve stability\n")
  }
}


patch_hglmer_poisson()


sample_size = 200 
if (n > sample_size) {
  cat("Using a smaller subset of", sample_size, "days for analysis\n")
  random_indices = sample(1:n, sample_size)
  simydata = simydata[random_indices]
  simxdata = simxdata[random_indices]
  beta_init = beta_init[random_indices,]
  n = sample_size
}



#################################################
# FIT BASELINE MODELS
#################################################
fit_baseline_models = function(simydata, simxdata) {
  cat("\n===== FITTING BASELINE MODELS =====\n")
  
  n = length(simydata)
  p = ncol(simxdata[[1]])
  hourly_obs = 24 
  
  combined_data = data.frame(
    response = numeric(n * hourly_obs),
    temp = numeric(n * hourly_obs),
    humidity = numeric(n * hourly_obs),
    windspeed = numeric(n * hourly_obs),
    hour = numeric(n * hourly_obs),
    day_id = numeric(n * hourly_obs)
  )
  
  idx = 1
  for (i in 1:n) {
    y_i = simydata[[i]]
    x_i = simxdata[[i]]
    
    for (j in 1:length(y_i)) {
      combined_data$response[idx] = y_i[j]
      combined_data$temp[idx] = x_i[j, 2]  
      combined_data$humidity[idx] = x_i[j, 3] 
      combined_data$windspeed[idx] = x_i[j, 4] 
      combined_data$hour[idx] = j
      combined_data$day_id[idx] = i
      idx = idx + 1
    }
  }
  
  combined_data$hour_factor = as.factor(combined_data$hour)
  
  k_folds = 5
  set.seed(123)
  fold_indices = sample(1:k_folds, n, replace=TRUE)
  
  models = list(
    "Poisson" = function(train_data) {
      glm(response ~ temp + humidity + windspeed, 
          family = poisson(), data = train_data)
    },
    "PoissonWithHour" = function(train_data) {
      glm(response ~ temp + humidity + windspeed + hour_factor, 
          family = poisson(), data = train_data)
    },
    "PoissonInteraction" = function(train_data) {
      glm(response ~ temp * humidity * windspeed + hour_factor, 
          family = poisson(), data = train_data)
    },
    "NegativeBinomial" = function(train_data) {
      suppressWarnings(
        MASS::glm.nb(response ~ temp + humidity + windspeed + hour_factor, 
                     data = train_data)
      )
    },
    "RandomEffects" = function(train_data) {
      suppressWarnings(
        lme4::glmer(response ~ temp + humidity + windspeed + hour_factor + (1|day_id), 
                   family = poisson(), data = train_data)
      )
    }
  )
  
  model_results = list()
  
  for (model_name in names(models)) {
    cat("Fitting", model_name, "model with cross-validation...\n")
    
    fold_log_scores = numeric(k_folds)
    
    for (fold in 1:k_folds) {
      train_days = which(fold_indices != fold)
      test_days = which(fold_indices == fold)
      
      train_idx = which(combined_data$day_id %in% train_days)
      test_idx = which(combined_data$day_id %in% test_days)
      
      train_data = combined_data[train_idx, ]
      test_data = combined_data[test_idx, ]
      
      model_fit = tryCatch({
        models[[model_name]](train_data)
      }, error = function(e) {
        cat("Error fitting", model_name, ":", conditionMessage(e), "\n")
        return(NULL)
      })
      
      if (!is.null(model_fit)) {
        if (model_name == "RandomEffects") {
          pred_means = rep(NA, nrow(test_data))
          
          unique_test_days = unique(test_data$day_id)
          
          for (day in unique_test_days) {
            day_idx = which(test_data$day_id == day)
            day_data = test_data[day_idx, ]
            
            mm <- model.matrix(~ temp + humidity + windspeed + hour_factor, data = day_data)
            beta_fixed <- lme4::fixef(model_fit)

            mm_cols <- colnames(mm)
            coef_names <- names(beta_fixed)
            matching_cols <- mm_cols %in% coef_names
            
            if(!all(matching_cols)) {
              mm <- mm[, matching_cols, drop=FALSE]
            }
            
            matching_coefs <- coef_names %in% colnames(mm)
            if(!all(matching_coefs)) {
              beta_fixed <- beta_fixed[matching_coefs]
            }

            linear_pred <- mm %*% beta_fixed[match(colnames(mm), names(beta_fixed))]
            pred_means[day_idx] <- exp(linear_pred)
          }
          
          log_lik = sum(dpois(test_data$response, lambda = pred_means, log = TRUE))
          
          fold_log_scores[fold] = log_lik / length(test_idx)
        } else {
          pred_means = predict(model_fit, newdata = test_data, type = "response")
          log_lik = sum(dpois(test_data$response, lambda = pred_means, log = TRUE))
          fold_log_scores[fold] = log_lik / length(test_idx)
        }
      } else {
        fold_log_scores[fold] = NA
      }
    }
    mean_log_score = mean(fold_log_scores, na.rm = TRUE)
    cat("Mean predictive log score for", model_name, ":", mean_log_score, "\n")
    
    model_results[[model_name]] = list(
      log_score = mean_log_score
    )
  }
  
  cat("Fitting K-means + Poisson model...\n")
  
  day_features = matrix(0, nrow = n, ncol = p)
  for (i in 1:n) {
    day_features[i, ] = colMeans(simxdata[[i]])
  }
  
  kmeans_log_scores = numeric(k_folds)
  
  for (fold in 1:k_folds) {
    train_days = which(fold_indices != fold)
    test_days = which(fold_indices == fold)
    
    kmeans_k = 2  
    kmeans_result = kmeans(day_features[train_days, ], centers = kmeans_k)
    
    test_features = day_features[test_days, ]
    test_clusters = numeric(length(test_days))
    
    for (i in 1:length(test_days)) {
      distances = apply(kmeans_result$centers, 1, function(center) {
        sum((test_features[i, ] - center)^2)
      })
      test_clusters[i] = which.min(distances)
    }
    
    cluster_models = list()
    for (j in 1:kmeans_k) {
      train_idx = which(combined_data$day_id %in% train_days & 
                         rep(kmeans_result$cluster == j, each = 24))
      
      if (length(train_idx) > 0) {
        cluster_models[[j]] = glm(
          response ~ temp + humidity + windspeed + hour_factor, 
          family = poisson(),
          data = combined_data[train_idx, ]
        )
      } else {
        cluster_models[[j]] = NULL
      }
    }
    
    test_log_lik = 0
    test_count = 0
    
    for (i in 1:length(test_days)) {
      day_id = test_days[i]
      cluster = test_clusters[i]
      
      if (!is.null(cluster_models[[cluster]])) {
        test_idx = which(combined_data$day_id == day_id)
        
        pred_means = predict(
          cluster_models[[cluster]], 
          newdata = combined_data[test_idx, ], 
          type = "response"
        )
        
        log_lik = sum(dpois(combined_data$response[test_idx], lambda = pred_means, log = TRUE))
        test_log_lik = test_log_lik + log_lik
        test_count = test_count + length(test_idx)
      }
    }
    
    if (test_count > 0) {
      kmeans_log_scores[fold] = test_log_lik / test_count
    } else {
      kmeans_log_scores[fold] = NA
    }
  }
  
  mean_kmeans_log_score = mean(kmeans_log_scores, na.rm = TRUE)
  cat("Mean predictive log score for K-means + Poisson:", mean_kmeans_log_score, "\n")
  
  model_results[["K-means + Poisson"]] = list(
    log_score = mean_kmeans_log_score
  )
  
  return(model_results)
}


#################################################
# FIT HGLMER 
#################################################

for (K in Kvalues) {
  cat("\nFitting mixture model with K =", K, "components\n")
  
  cat("Starting manual clustering...\n")
  set.seed(123)
  
  cat("Beta init dimensions:", dim(beta_init), "\n")
  cat("Checking for NA/Inf values in beta_init before clustering:\n")
  cat("NA values:", sum(is.na(beta_init)), "\n")
  cat("Inf values:", sum(is.infinite(beta_init)), "\n")
  
  init_clusters = kmeans(beta_init, centers=K, nstart=10)
  z_init = init_clusters$cluster
  
  cat("Cluster sizes:", table(z_init), "\n")
  
  cat("Creating initial parameter values...\n")
  theta_init = matrix(0, nrow=K, ncol=p)
  Lambda_init = array(0, dim=c(p, p, K))
  psi_init = numeric(K) 
  
  for (k in 1:K) {
    cat("Processing component", k, "...\n")
    cluster_size = sum(z_init == k)
    cat("Cluster", k, "size:", cluster_size, "\n")
    
    if (cluster_size > 0) {
      theta_init[k,] = colMeans(beta_init[z_init == k, , drop=FALSE])
      cat("Component", k, "theta:", theta_init[k,], "\n")
      
      if (cluster_size > 1) {
        cov_matrix = tryCatch({
          cov_result = cov(beta_init[z_init == k, , drop=FALSE])
          if (any(is.na(cov_result)) || any(is.infinite(cov_result))) {
            cat("Warning: Invalid values in covariance calculation for component", k, "\n")
            diag(0.5, p)
          } else {
            cov_result
          }
        }, error = function(e) {
          cat("Error in covariance calculation for component", k, ":", conditionMessage(e), "\n")
          diag(0.5, p)
        })
        
        Lambda_init[,,k] = cov_matrix
      } else {
        cat("Only one observation in component", k, ", using diagonal matrix\n")
        Lambda_init[,,k] = diag(0.5, p)
      }
      
      psi_init[k] = cluster_size / n
    } else {
      cat("Empty cluster", k, ", using defaults\n")
      theta_init[k,] = colMeans(beta_init)
      Lambda_init[,,k] = diag(0.5, p)
      psi_init[k] = 1/K
    }
  }
  
  psi_init = psi_init / sum(psi_init)
  
  cat("Checking positive definiteness of Lambda matrices...\n")
  for (k in 1:K) {
    eigen_vals = eigen(Lambda_init[,,k])$values
    cat("Component", k, "eigenvalues:", eigen_vals, "\n")
    
    if (any(eigen_vals <= 0)) {
      cat("Fixing non-positive definite matrix for component", k, "\n")
      Lambda_init[,,k] = Lambda_init[,,k] + diag(0.1 + abs(min(eigen_vals)), p)
      # Verify fix
      eigen_vals_fixed = eigen(Lambda_init[,,k])$values
      cat("Fixed eigenvalues:", eigen_vals_fixed, "\n")
    }
  }
  
  if (K > 1) {
    cat("Ordering components by size...\n")
    cat("Initial psi values:", psi_init, "\n")
    
    sort_idx = order(psi_init)
    psi_init = psi_init[sort_idx]
    theta_init = theta_init[sort_idx,]
    
    Lambda_init_temp = Lambda_init
    for (k in 1:K) {
      Lambda_init[,,k] = Lambda_init_temp[,,sort_idx[k]]
    }
    
    z_init = match(z_init, sort_idx)
    cat("Ordered psi values:", psi_init, "\n")
  }
  
  cat("\nFinal initialization summary:\n")
  cat("Component proportions (psi):", psi_init, "\n")
  for (k in 1:K) {
    cat("Component", k, "mean (theta):", theta_init[k,], "\n")
    cat("Component", k, "cluster size:", sum(z_init == k), "\n")
  }
  
  cat("\nStarting MCMC sampling...\n")
  fit_bike = tryCatch({
    for(i in 1:n) {
      if(length(simydata[[i]]) != standard_length) {
        cat("ERROR: Observation", i, "has incorrect length:", length(simydata[[i]]), "\n")
        stop("Inconsistent observation length detected")
      }
      if(nrow(simxdata[[i]]) != standard_length) {
        cat("ERROR: Design matrix", i, "has incorrect rows:", nrow(simxdata[[i]]), "\n")
        stop("Inconsistent design matrix detected")
      }
    }
    
    gibbs_hglmer_poisson(
      ydata = simydata,
      xdata = simxdata,
      K = K,
      u0 = u0,
      V0 = V0,
      w0 = w0,
      mcmc = list(nblow = 1000, smcmc = 1000),
      init = list(
        beta = beta_init,
        z = z_init,
        theta = theta_init,
        Lambda = Lambda_init,
        psi = psi_init
      )
    )
  }, error = function(e) {
    cat("ERROR in MCMC sampling:", conditionMessage(e), "\n")
    cat("Dimensions of key variables:\n")
    cat("- simydata length:", length(simydata), "\n")
    cat("- simxdata length:", length(simxdata), "\n") 
    cat("ERROR in MCMC sampling:", conditionMessage(e), "\n")
    cat("Dimensions of key variables:\n")
    cat("- simydata length:", length(simydata), "\n")
    cat("- simxdata length:", length(simxdata), "\n") 
    cat("- beta_init dimensions:", dim(beta_init)[1], "x", dim(beta_init)[2], "\n")
    cat("- z_init length:", length(z_init), "\n")
    return(NULL)
  })
  
  if (is.null(fit_bike)) {
    cat("MCMC failed, skipping remaining analysis\n")
    next
  }
  
  results[[K]] = fit_bike

if (!is.null(fit_bike)) {
  pred_logscore = calculate_pred_logscore(fit_bike, simydata, simxdata, K)
  results[[K]]$pred_logscore = pred_logscore
  }
  
  cat("Calculating log marginal likelihood...\n")
  logML_value = tryCatch({
    logML_hglmer_poisson(
      post_draws = fit_bike,
      ydata = simydata,
      xdata = simxdata,
      K = K,
      u0 = u0,
      V0 = V0,
      w0 = w0
    )
  }, error = function(e) {
    cat("Error calculating log marginal likelihood:", conditionMessage(e), "\n")
    return(NA)
  })
  
  if (!is.na(logML_value)) {
    logML_values[[K]] = logML_value
    cat("Log marginal likelihood for K =", K, ":", logML_value, "\n")
  }
  
  result_psi = do.call(rbind, fit_bike$psi)
  psi_means = apply(result_psi, 2, mean)
  cat("Estimated mixture proportions:\n")
  print(psi_means)
  
  dir.create("./plots", showWarnings = FALSE) 
  
  png(paste0("./plots/bike_psi_plot_K_", K, ".png"))
  par(mfrow=c(1, K))
  for (i in 1:K) {
    plot(result_psi[,i], type='l', main=paste("Component", i),
         xlab="Iteration", ylab="Mixture Probability")
  }
  dev.off()
  
  png(paste0("./plots/bike_theta_plot_K_", K, ".png"))
  
  result_theta = fit_bike$theta
  theta_list = list()
  for (i in 1:K) {
    theta_list[[i]] = do.call(rbind, lapply(result_theta, function(x) x[i,]))
  }
  
  means_by_component = lapply(theta_list, function(x) apply(x, 2, mean))
  for (i in 1:K) {
    cat("Component", i, "means (intercept, temp, humidity, windspeed):", means_by_component[[i]], "\n")
  }
  
  par(mfrow=c(K, 4), mar=c(4,4,2,1))
  param_names = c("Intercept", "Temperature", "Humidity", "Wind Speed")
  
  for (i in 1:K) {
    for (j in 1:4) {
      plot(theta_list[[i]][,j], type='l', 
           main=paste("Component", i, param_names[j]),
           xlab="Iteration", ylab="Value")
    }
  }
  dev.off()
  
  rm(fit_bike)
  gc()
}


baseline_results = fit_baseline_models(simydata, simxdata)
cat("Predictive Log Scores:\n")

all_log_scores = c()

for (K in Kvalues) {
  if (!is.null(results[[K]]) && !is.null(results[[K]]$pred_logscore)) {
    all_log_scores = rbind(all_log_scores, 
                          data.frame(Model = paste("hglmer-K", K), 
                                     LogScore = results[[K]]$pred_logscore))
  }
}

for (model_name in names(baseline_results)) {
  all_log_scores = rbind(all_log_scores, 
                        data.frame(Model = model_name, 
                                   LogScore = baseline_results[[model_name]]$log_score))
}
all_log_scores = all_log_scores[order(-all_log_scores$LogScore), ]
print(all_log_scores)
