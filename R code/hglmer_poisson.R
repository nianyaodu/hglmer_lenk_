#-------------------------------------------------------------------------------
#' Finite mixtures of generalized linear models with random effects
#' : hierarchical Poisson regression model

#' @param ydata List; vector of count dependent variables
#' @param xdata List; vector of independent variables
#' @param K Integer; Number of mixture components
#' @param u0 Prior mean for theta
#' @param V0 Prior variance for theta
#' @param w0 Prior concentration parameter for Dirichlet
#' @param mcmc List; settings for MCMC
#' @param show_pbar Logical; whether to show the progress bar with ... Default to be TRUE.
#' @param init List; optional initialization values (beta, z, theta, Lambda, psi)
#'
#' @return Posterior samples of beta, theta, Lambda, psi, z

source("ordered_dirichlet.R")
library(MASS)
library(MCMCpack)

gibbs_hglmer_poisson <- function(ydata, xdata, K=NULL, u0, V0, w0, mcmc=list(), show_pbar=TRUE, init=NULL) {
  
  # Preliminary settings
  n = length(ydata)
  p = ncol(xdata[[1]])
  
  print(paste("Number of subjects:", n))
  print(paste("Number of predictors:", p))
  
  mcvals = list(nblow=1000, smcmc=1000)
  mcvals[names(mcmc)] = mcmc
  nblow = mcvals$nblow    # number of initial MCMC parameters
  smcmc = mcvals$smcmc    # number of MCMC parameters to save
  nmcmc = nblow + smcmc   # total number of MCMC parameters
  
  # Initialize sampler
  f0 = p + 1
  G0 = p*diag(rep(1,p))
  
  # Check if we have user-provided initialization
  if (!is.null(init)) {
    # Use provided initialization values
    beta = init$beta
    z = init$z
    theta = init$theta
    Lambda = init$Lambda
    psi = init$psi
    beta_hat = beta  # Initialize beta_hat
    
    cat("Using user-provided initialization\n")
  } else {
    # First stage: Initialize beta using glm with Poisson family
    beta = matrix(0, nrow=n, ncol=p)
    beta_hat = matrix(0, nrow=n, ncol=p)

    for (i in 1:n) {
      y_i = ydata[[i]]
      x_i = xdata[[i]]
      
      # Create proper data frame with all variables
      df = data.frame(y = y_i)
      for(j in 2:p) {
        df[[paste0("x", j-1)]] = x_i[,j]
      }
      
      # Construct formula dynamically based on number of predictors
      if(p > 2) {
        formula_str = paste("y ~", paste(paste0("x", 1:(p-1)), collapse=" + "))
        formula_obj = as.formula(formula_str)
      } else {
        formula_obj = as.formula("y ~ x1")
      }
      
      # Safely fit Poisson GLM with tryCatch
      fit_result = tryCatch({
        fit_poisson = glm(formula_obj, data=df, family=poisson())
        coef(fit_poisson)
      }, error = function(e) {
        # Fallback to reasonable values if GLM fails
        c(log(mean(y_i) + 0.1), rep(0, p-1))
      })
      
      beta[i,] = beta_hat[i,] = fit_result
    }
    
    # Second stage: Apply EM algorithm using Mclust
    Mclust = mclust:::Mclust
    mclustBIC = mclust:::mclustBIC
    
    suppressMessages({
      if (is.null(K)) {
        # If K is NULL, try K-means instead of Mclust
        kmeans_results = lapply(1:5, function(k) {
          kmeans(beta, centers=k, nstart=10)
        })
        wss = sapply(kmeans_results, function(km) sum(km$withinss))
        K = min(which(diff(wss)/wss[-length(wss)] < 0.2)) + 1
        if (is.infinite(K) || K < 2) K = 2
        
        init_clusters = kmeans_results[[K-1]]
        z = init_clusters$cluster
        
        # Create initial values
        theta = matrix(0, nrow=K, ncol=p)
        Lambda = array(0, dim=c(p, p, K))
        psi = numeric(K)
        
        for (k in 1:K) {
          idx = which(z == k)
          if (length(idx) > 0) {
            theta[k,] = colMeans(beta[idx, , drop=FALSE])
            if (length(idx) > 1) {
              Lambda[,,k] = cov(beta[idx, , drop=FALSE])
              # Ensure positive definiteness
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
        
        # Ensure psi values are ordered
        ord = order(psi)
        psi = psi[ord]
        theta = theta[ord,]
        Lambda_temp = Lambda
        for (k in 1:K) {
          Lambda[,,k] = Lambda_temp[,,ord[k]]
        }
        z = match(z, ord)
      } else {
        # Try Mclust but fall back to K-means if it fails
        tryCatch({
          fit_EM_beta = Mclust(beta, G=K, verbose=FALSE)
          
          # Resolve label switching by ordering components
          ord = order(fit_EM_beta$parameters$pro)
          theta = t(fit_EM_beta$parameters$mean)[ord,]
          Lambda = fit_EM_beta$parameters$variance$sigma[,,ord]
          psi = fit_EM_beta$parameters$pro[ord]
          z = match(fit_EM_beta$classification, ord)
        }, error = function(e) {
          cat("Mclust failed, falling back to K-means clustering\n")
          init_clusters = kmeans(beta, centers=K, nstart=10)
          z = init_clusters$cluster
          
          # Create initial values
          theta = matrix(0, nrow=K, ncol=p)
          Lambda = array(0, dim=c(p, p, K))
          psi = numeric(K)
          
          for (k in 1:K) {
            idx = which(z == k)
            if (length(idx) > 0) {
              theta[k,] = colMeans(beta[idx, , drop=FALSE])
              if (length(idx) > 1) {
                Lambda[,,k] = cov(beta[idx, , drop=FALSE])
                # Ensure positive definiteness
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
          
          # Ensure psi values are ordered
          ord = order(psi)
          psi = psi[ord]
          theta = theta[ord,]
          Lambda_temp = Lambda
          for (k in 1:K) {
            Lambda[,,k] = Lambda_temp[,,ord[k]]
          }
          z = match(z, ord)
        })
      }
    })
  }
  
  # Check if initialization is valid
  if (any(is.na(beta)) || any(is.infinite(beta))) {
    cat("Warning: Invalid values in beta initialization\n")
    beta[is.na(beta)] = 0
    beta[is.infinite(beta)] = 0
  }
  
  for (k in 1:K) {
    if (any(is.na(Lambda[,,k])) || any(is.infinite(Lambda[,,k])) || 
        any(eigen(Lambda[,,k])$values <= 0)) {
      cat("Warning: Invalid Lambda matrix for component", k, "\n")
      Lambda[,,k] = diag(0.5, p)
    }
  }
  
  # Create storage for MCMC samples
  beta.save = list()
  beta.hat.save = list()
  z.save = list()
  theta.save = list()
  Lambda.save = list()
  psi.save = list()
  
  # Gibbs Updates
  cat(sprintf("Generating hglmer posterior samples \n"))
  if(show_pbar) { pb <- txtProgressBar(min = 1, max = nmcmc, style=3) }
  
  for (iter in 1:nmcmc) {
    # Update z (component assignments)
    for (i in 1:n) {
      lpsi_i = log(psi) + 
        sapply(1:K, function(k) dmvnorm(beta[i,], mean=theta[k,], sigma=Lambda[,,k], log=TRUE))
      lpsi_i1 = lpsi_i - max(lpsi_i)
      psi_i1 = exp(lpsi_i1)
      psi_i = psi_i1/sum(psi_i1)
      z[i] = sample(K, 1, prob=psi_i)
    }
    
    # Update beta (regression coefficients)
    for (i in 1:n) {
      y_i = ydata[[i]]
      x_i = xdata[[i]]
      beta_hat_i = beta_hat[i,]
      Lamb_inv = solve(Lambda[,,z[i]])
      
      # Calculate gradient and Hessian for Poisson regression
      # Ensure dimensions match by explicitly calculating one-by-one
      log_mu = numeric(length(y_i))
      mu = numeric(length(y_i))
      
      for (j in 1:length(y_i)) {
        log_mu[j] = sum(x_i[j,] * beta_hat_i)
        mu[j] = exp(log_mu[j])
      }
      
      # Hessian - build it using a more robust approach
      V_i_sum = matrix(0, nrow=p, ncol=p)
      for (j in 1:length(mu)) {
        w_j = mu[j]
        x_j = x_i[j,]
        V_i_sum = V_i_sum + w_j * (x_j %*% t(x_j))
      }
      V_i = solve(V_i_sum + Lamb_inv)
      
      # Gradient - calculate each component individually
      resids = y_i - mu
      grad_llik_i = numeric(p)
      for (j in 1:length(resids)) {
        grad_llik_i = grad_llik_i + resids[j] * x_i[j,]
      }
      grad_prior_i = -Lamb_inv %*% (beta_hat_i - theta[z[i],])
      grad_i = grad_llik_i + grad_prior_i
      
      # Metropolis-Hastings step
      beta_hat[i,] = as.vector(beta_hat_i + V_i %*% grad_i)
      
      beta_i_prop = as.vector(mvrnorm(n=1, mu=beta_hat[i,], Sigma=V_i))
      
      # Calculate log-posterior for proposed value
      log_mu_prop = numeric(length(y_i))
      log_post_prop = 0
      for (j in 1:length(y_i)) {
        log_mu_prop[j] = sum(x_i[j,] * beta_i_prop)
        log_post_prop = log_post_prop + (y_i[j] * log_mu_prop[j] - exp(log_mu_prop[j]))
      }
      log_post_prop = log_post_prop - 0.5 * t(beta_i_prop - theta[z[i],]) %*% Lamb_inv %*% (beta_i_prop - theta[z[i],]) + 
        0.5 * t(beta_i_prop - beta_hat[i,]) %*% solve(V_i) %*% (beta_i_prop - beta_hat[i,])
      
      # Calculate log-posterior for current value
      log_mu_curr = numeric(length(y_i))
      log_post_curr = 0
      for (j in 1:length(y_i)) {
        log_mu_curr[j] = sum(x_i[j,] * beta[i,]) 
        log_post_curr = log_post_curr + (y_i[j] * log_mu_curr[j] - exp(log_mu_curr[j]))
      }
      log_post_curr = log_post_curr - 0.5 * t(beta[i,] - theta[z[i],]) %*% Lamb_inv %*% (beta[i,] - theta[z[i],]) + 
        0.5 * t(beta[i,] - beta_hat[i,]) %*% solve(V_i) %*% (beta[i,] - beta_hat[i,])
      
      # Accept/reject
      accprob = min(0, log_post_prop - log_post_curr)
      
      if(log(runif(1)) < accprob) {
        beta[i,] = beta_i_prop
      }
    }
    
    # Update theta (component means)
    for (k in 1:K) {
      n_k = sum(z==k)
      if (n_k > 0) { # Ensure there's at least one subject in this component
        Vn_k = solve(n_k*solve(Lambda[,,k]) + solve(V0*diag(rep(1,p))))
        un_k = Vn_k %*% (solve(Lambda[,,k]) %*% apply(beta[z==k,,drop=FALSE], 2, sum) + 
                        solve(V0*diag(rep(1,p))) %*% rep(u0,p))
        theta[k,] = as.vector(mvrnorm(n=1, mu=un_k, Sigma=Vn_k))
      } else {
        # Handle empty components
        theta[k,] = rep(u0, p) + rnorm(p, 0, sqrt(V0))
      }
    }
    
    # Update Lambda (component covariances)
    for (k in 1:K) {
      n_k = sum(z==k)
      if (n_k > 0) { # Ensure there's at least one subject in this component
        fn_k = f0 + n_k
        M = sweep(beta[z==k,,drop=FALSE], 2, theta[k,], `-`)
        Gn_k = G0 + t(M) %*% M
        Lambda[,,k] = riwish(fn_k, Gn_k)
      } else {
        # Handle empty components
        Lambda[,,k] = riwish(f0, G0)
      }
    }
    
    # Update psi (mixture probabilities)
    n_K = numeric(K)
    for (k in 1:K) { n_K[k] = sum(z==k) }
    shape = rep(w0, K) + n_K
    psi = as.vector(rng_orderdiri(1, shape))
    
    # Record Posterior Samples
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

# Log marginal likelihood calculation for model selection
logML_hglmer_poisson <- function(post_draws, ydata, xdata, K, u0, V0, w0) {
  ## Preliminary settings
  n = length(ydata)
  p = ncol(xdata[[1]])
  niter = length(post_draws$beta)
  f0 = p + 1
  G0 = p*diag(rep(1,p))
  
  ## Calculate posterior estimates
  
  # 1. beta
  betas = lapply(1:n, function(i) {
    do.call(rbind, lapply(post_draws$beta, function(M) M[i,,drop=FALSE]))})
  MU_beta = lapply(betas, colMeans)
  Sigma_beta = lapply(betas, cov)
  
  # 2. z
  pi_z = colMeans(do.call(rbind, post_draws$psi))
  
  # 3. theta
  thetas = lapply(1:K, function(k) {
    do.call(rbind, lapply(post_draws$theta, function(M) M[k,,drop=FALSE]))})
  MU_theta = lapply(thetas, colMeans)
  Sigma_theta = lapply(thetas, cov)
  
  # 4. lambda
  MU_psi = lapply(1:K, function(k) {
    mean(sapply(post_draws$psi, function(x) x[k]))})
  MU_Lambdainv = lapply(1:K, function(k) { 
    apply(simplify2array(lapply(post_draws$Lambda, function(M) solve(M[,,k]))), c(1,2), mean)})
  h_Lambda = lapply(1:K, function(k) f0 + n * MU_psi[[k]])
  H_Lambda = lapply(1:K, function(k) h_Lambda[[k]]*solve(MU_Lambdainv[[k]]))
  
  # 5. psi
  MU_psi = lapply(1:K, function(k) {
    mean(sapply(post_draws$psi, function(x) x[k]))})
  Sigma_psi = lapply(1:K, function(k) {
    var(sapply(post_draws$psi, function(x) x[k]))})
  V_est = (1-sum(unlist(MU_psi)^2))/sum(unlist(Sigma_psi))-1
  V_psi = lapply(1:K, function(k) V_est*MU_psi[[k]])
  
  ## marginal likelihood
  g_vec = numeric(niter)
  
  for (u in 1:niter) {
    beta = post_draws$beta[[u]]
    for (i in 1:n) {
      g_vec[u] = g_vec[u] + dmvnorm(beta[i,], MU_beta[[i]], Sigma_beta[[i]], log=TRUE)
    }
    
    z_mat <- t(sapply(post_draws$z[[u]], function(k) {
      out <- numeric(K)
      out[k] <- 1
      out
    }))
    
    for (i in 1:n) {
      g_vec[u] = g_vec[u] + dmultinom(z_mat[i,], prob=pi_z, log=TRUE)
    }
    
    theta = post_draws$theta[[u]]
    for (k in 1:K) {
      g_vec[u] = g_vec[u] + dmvnorm(theta[k,], MU_theta[[k]], Sigma_theta[[k]], log=TRUE)
    }
    
    Lambda = post_draws$Lambda[[u]]
    for (k in 1:K) {
      g_vec[u] = g_vec[u] + log(diwish(Lambda[,,k], h_Lambda[[k]], H_Lambda[[k]]))
    }
    
    psi = post_draws$psi[[u]]
    g_vec[u] = g_vec[u] + log(ddirichlet(psi, unlist(V_psi)))
  }
  
  llik_vec = numeric(niter)
  
  for (u in 1:niter) {
    beta = post_draws$beta[[u]]
    z = post_draws$z[[u]]
    theta = post_draws$theta[[u]]
    Lambda = post_draws$Lambda[[u]]
    psi = post_draws$psi[[u]]
    
    # likihood for Poisson observations
    for (i in 1:n) {
      y_i = ydata[[i]]
      x_i = xdata[[i]]
      log_mu = x_i %*% beta[i,]
      llik_vec[u] = llik_vec[u] + sum(dpois(y_i, exp(log_mu), log=TRUE))
    }
    
    # likelihood for subject-level parameters
    for (i in 1:n) {
      llik_vec[u] = llik_vec[u] + dmvnorm(beta[i,], theta[z[i],], Lambda[,,z[i]], log=TRUE) +
        log(psi[z[i]])
    }
    
    # prior for component parameters
    for (k in 1:K) {
      llik_vec[u] = llik_vec[u] + dmvnorm(theta[k,], rep(u0,p), V0*diag(1,p), log=TRUE)
    }
    
    for (k in 1:K) {
      llik_vec[u] = llik_vec[u] + log(diwish(Lambda[,,k], f0, G0))
    }
    
    llik_vec[u] = llik_vec[u] + log(ddirichlet(psi, rep(w0,K))) + lfactorial(K)
  }
  
  logint = g_vec - llik_vec
  logint_m = max(logint)
  logML = -(logint_m + log(mean(exp(logint-logint_m))))
  
  return(logML)
}