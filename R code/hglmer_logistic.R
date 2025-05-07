#-------------------------------------------------------------------------------
#' Finite mixtures of generalized linear models with random effects
#' : hierarchical logistic regression model

#' @param ydata List; vector of dependent variable
#' @param xdata List; vector of independent variables
#' @param K Integer; Number of mixture components
#' @u0 
#' @V0 
#' @w0 
#' @param n_iter Integer; number of iterations
#' @param show_pbar Logical; whether to show the progress bar with ... Default to be TRUE.
#' 
#' @return Posterior samples of beta, theta, Lambda, psi, z

source("ordered_dirichlet.R")
library(MASS)
library(MCMCpack)
library(emulator)

"gibbs_hglmer_logistic" <- function(ydata, xdata, K=NULL, u0, V0, w0, mcmc=list(), show_pbar=TRUE){
  
  # Preliminary settings
  n = length(ydata)
  p = ncol(xdata[[1]])
  
  mcvals = list(nblow=1000,smcmc=1000)
  mcvals[names(mcmc)]=mcmc
  nblow     = mcvals$nblow    # number of initial MCMC parameters
  smcmc     = mcvals$smcmc    # number of MCMC parameters to save
  nmcmc     = nblow + smcmc   # total number of MCMC parameters
  
  # Initialize sampler
  f0       = p + 1
  G0       = p*diag(rep(1,p))
  beta     = matrix(0, nrow=n, ncol=p)
  beta_hat = matrix(0, nrow=n, ncol=p)
  for (i in 1:n) {
    simdata   = data.frame(y=ydata[[i]], x=xdata[[i]][,2])
    fit_logit = glm(y ~ x, data=simdata, family=binomial())
    beta[i,]  = beta_hat[i,] = coef(fit_logit)
  }
  Mclust    = mclust:::Mclust
  mclustBIC = mclust:::mclustBIC
  
  suppressMessages({
    if (is.null(K)) {
      fit_EM_beta = Mclust(beta, verbose=FALSE)
      K = fit_EM_beta$G
    } else {
      fit_EM_beta = Mclust(beta, G=K, verbose=FALSE)
    }
  })
  
  ord    = order(fit_EM_beta$parameters$pro)
  theta  = t(fit_EM_beta$parameters$mean)[ord,,drop=FALSE]
  Lambda = fit_EM_beta$parameters$variance$sigma[,,ord,drop=FALSE]
  psi    = fit_EM_beta$parameters$pro[ord]
  z      = match(fit_EM_beta$classification, ord)
  
  # Making objects
  beta.save     = list()
  beta.hat.save = list()
  z.save        = list()
  theta.save    = list()
  Lambda.save   = list()
  psi.save      = list()
  
  # Gibbs Updates
  cat(sprintf("Generating hglmer posterior samples \n"))
  if(show_pbar) { pb <- txtProgressBar(min = 1, max = nmcmc, style=3) }
  
  for (iter in 1:nmcmc) {
    # update z
    for (i in 1:n) {
      lpsi_i  = log(psi) + 
        sapply(1:K, function(k) dmvnorm(beta[i,],mean=theta[k,],sigma=Lambda[,,k],log=TRUE))
      lpsi_i1 = lpsi_i - max(lpsi_i)
      psi_i1  = exp(lpsi_i1)
      psi_i   = psi_i1/sum(psi_i1)
      z[i]    = sample(K,1,prob=psi_i)
    }
    
    # update beta
    for (i in 1:n) {
      y_i        = ydata[[i]]
      x_i        = xdata[[i]]
      beta_hat_i = beta_hat[i,]
      Lamb_inv   = solve(Lambda[,,z[i]])
      
      # hessian
      pi  = 1/(1+exp(-x_i%*%beta_hat_i))
      w   = as.vector(pi*(1-pi))
      Hs  = lapply(seq_along(w), function(k) { w[k]*x_i[k,]%*%t(x_i[k,]) })
      V_i = solve(Reduce("+", Hs) + Lamb_inv)
      
      # gradient
      resids       = as.vector(y_i-pi)
      grad_llik_i  = t(x_i)%*%resids
      grad_prior_i = -Lamb_inv%*%(beta_hat_i-theta[z[i],])
      grad_i       = grad_llik_i + grad_prior_i
      
      # Metropolis-Hastings step
      beta_hat[i,] = as.vector(beta_hat_i + V_i %*% grad_i)
      beta_i_prop  = mvrnorm(n=1, mu = beta_hat[i,], Sigma = V_i)
      
      lpost_prop = sum(y_i*x_i%*%beta_i_prop-log1p(exp(x_i%*%beta_i_prop))) - 
        0.5*t(beta_i_prop-theta[z[i],])%*%Lamb_inv%*%(beta_i_prop-theta[z[i],]) + 
        0.5*t(beta_i_prop-beta_hat[i,])%*%solve(V_i)%*%(beta_i_prop-beta_hat[i,])
      
      lpost_curr = sum(y_i*x_i%*%beta[i,]-log1p(exp(x_i%*%beta[i,]))) - 
        0.5*t(beta[i,]-theta[z[i],])%*%Lamb_inv%*%(beta[i,]-theta[z[i],]) + 
        0.5*t(beta[i,]-beta_hat[i,])%*%solve(V_i)%*%(beta[i,]-beta_hat[i,])
      
      accprob = min(0, lpost_prop-lpost_curr)
      
      if(log(runif(1)) < accprob){
        beta[i,] = beta_i_prop
      }
    }
    
    # update theta
    for (k in 1:K) {
      n_k       = sum(z==k)
      Vn_k      = solve(n_k*solve(Lambda[,,k])+solve(V0*diag(rep(1,p))))
      un_k      = Vn_k%*%(solve(Lambda[,,k])%*%apply(beta[z==k,,drop=FALSE],2,sum)+solve(V0*diag(rep(1,p)))%*%rep(u0,p))
      theta[k,] = rmvnorm(1,mean=un_k,sigma=Vn_k)
    }
    
    # update Lambda
    for (k in 1:K) {
      n_k         = sum(z==k)
      fn_k        = f0 + n_k
      M           = sweep(beta[z==k,,drop=FALSE], 2, theta[k,], `-`)
      Gn_k        = G0 + t(M)%*%M
      Lambda[,,k] = riwish(fn_k,Gn_k)
    }
    
    # update psi
    if (K > 1) {
      n_K = vector(length=K)
      for (k in 1:K) { n_K[k] = sum(z==k) }
      shape = rep(w0,K) + n_K
      psi   = as.vector(rng_orderdiri(1,shape))
    }
    
    # Record Posterior Samples
    if (iter > nblow) {
      beta.save[[iter-nblow]]     = beta
      z.save[[iter-nblow]]        = z
      theta.save[[iter-nblow]]    = theta
      Lambda.save[[iter-nblow]]   = Lambda
      psi.save[[iter-nblow]]      = psi
      beta.hat.save[[iter-nblow]] = beta_hat
    }
    
    if (show_pbar) { setTxtProgressBar(pb,iter) }
  }
  
  cat(sprintf("\nMCMC is done! \n"))
  
  out = NULL
  out$beta     = beta.save
  out$beta_hat = beta.hat.save
  out$theta    = theta.save
  out$Lambda   = Lambda.save
  out$psi      = psi.save
  out$z        = z.save
  out$K        = K
  
  return(out)
}




