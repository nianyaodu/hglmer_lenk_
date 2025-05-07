#-------------------------------------------------------------------------------
#' Finite mixtures of generalized linear models with random effects
#' : hierarchical linear regression model

#' @param ydata List; vector of dependent variables
#' @param xdata List; vector of independent variables
#' @param K Integer; Number of mixture components
#' @ a0
#' @ d0
#' @ r0
#' @ s0
#' @ u0
#' @ V0
#' @ w0
#' @param n_iter Integer; number of iterations
#' @param show_pbar Logical; whether to show the progress bar with ... Default to be TRUE.
#'
#' @return Posterior samples of beta, theta, Lambda, psi, z

source("ordered_dirichlet.R")
library(MASS)
library(MCMCpack)

gibbs_hglmer_linear <- function(ydata, xdata, K=NULL, a0, d0, r0, s0, u0, V0, w0, mcmc=list(), show_pbar=TRUE) {
  
  # Preliminary settings
  n = length(ydata)
  p = ncol(xdata[[1]])
  
  mcvals = list(nblow=1000,smcmc=1000)
  mcvals[names(mcmc)]=mcmc
  nblow     = mcvals$nblow    # number of initial MCMC parameters
  smcmc     = mcvals$smcmc    # number of MCMC parameters to save
  nmcmc     = nblow + smcmc   # total number of MCMC parameters
  
  # Initialize sampler
  f0      = p + 1
  G0      = p*diag(rep(1,p))
  alpha   = rnorm(1,a0,d0)
  tau2    = rinvgamma(1,shape=r0/2,scale=s0/2)
  
  beta    = matrix(0, nrow=n, ncol=p)
  phi     = numeric(n)
  phi_hat = numeric(n)
  for (i in 1:n) {
    simdata  = data.frame(y=ydata[[i]], x=xdata[[i]][,2])
    fit      = lm(y ~ x, data=simdata)
    beta[i,] = fit$coefficients
    phi[i]   = phi_hat[i] = log(summary(fit)$sigma^2)
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
  phi.save      = list()
  phi.hat.save  = list() 
  z.save        = list()
  theta.save    = list()
  Lambda.save   = list()
  psi.save      = list()
  alpha.save    = list() 
  tau2.save     = list()
  
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
      y_i      = ydata[[i]]
      x_i      = xdata[[i]]
      sig2_i   = exp(phi[i])
      D_i      = solve(t(x_i)%*%x_i/sig2_i+solve(Lambda[,,z[i]]))
      b_i      = D_i%*%(t(x_i)%*%y_i/sig2_i + solve(Lambda[,,z[i]])%*%theta[z[i],])
      beta[i,] = rmvnorm(1, mean=b_i, sigma=D_i)
    }
    
    # update theta
    for (k in 1:K) {
      n_k       = sum(z==k)
      Vn_k      = solve(n_k*solve(Lambda[,,k])+solve(V0*diag(rep(1,p))))
      un_k      = Vn_k%*%(solve(Lambda[,,k])%*%apply(beta[z==k,,drop=FALSE],2,sum)+solve(V0*diag(rep(1,p)))%*%rep(u0,p))
      theta[k,] = rmvnorm(1, mean=un_k, sigma=Vn_k)
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
    
    # update phi
    for (i in 1:n) {
      y_i        = ydata[[i]]
      x_i        = xdata[[i]]
      beta_i     = beta[i,]
      phi_hat_i  = phi_hat[i]
      SE_i       = sum(y_i*as.vector(x_i%*%beta_i)-as.vector(x_i%*%beta_i)^2/2)
      v_i        = 1/((sum(y_i^2)/2-SE_i)*exp(-phi_hat_i)+1/tau2)
      phi_hat[i] = phi_hat_i + v_i*((sum(y_i^2)/2-SE_i)*exp(-phi_hat_i)
                                    -length(y_i)/2-(phi_hat_i-alpha)/tau2)
      phi_i_prop = rnorm(1, phi_hat[i], sqrt(v_i))
      
      lpost_prop = (SE_i-sum(y_i^2)/2)*exp(-phi_i_prop) - 
        length(y_i)*phi_i_prop/2 - (phi_i_prop-alpha)^2/(2*tau2) + 
        (phi_i_prop-phi_hat[i])^2/(2*v_i)
      
      lpost_curr = (SE_i-sum(y_i^2)/2)*exp(-phi[i]) - 
        length(y_i)*phi[i]/2 - (phi[i]-alpha)^2/(2*tau2) + 
        (phi[i]-phi_hat[i])^2/(2*v_i)
      
      accprob = min(0, lpost_prop-lpost_curr)
      
      if(log(runif(1)) < accprob) {
        phi[i] = phi_i_prop
      }
    }
    
    # update alpha
    d_n2  = 1/(n/tau2+1/d0^2)
    a_n   = d_n2*(sum(phi)/tau2+a0/d0^2)
    alpha = rnorm(1,a_n,sqrt(d_n2))
    
    # update tau2
    r_n  = r0 + n
    s_n  = s0 + sum((phi-alpha)^2)
    tau2 = rinvgamma(1, shape=r_n/2, scale=s_n/2)
    
    # Record Posterior Samples
    if (iter > nblow) {
      beta.save[[iter-nblow]]    = beta
      phi.save[[iter-nblow]]     = phi
      phi.hat.save[[iter-nblow]] = phi_hat
      z.save[[iter-nblow]]       = z
      theta.save[[iter-nblow]]   = theta
      Lambda.save[[iter-nblow]]  = Lambda
      psi.save[[iter-nblow]]     = psi
      alpha.save[[iter-nblow]]   = alpha
      tau2.save[[iter-nblow]]    = tau2
    }
    
    if (show_pbar) { setTxtProgressBar(pb,iter) }
  }
  
  cat(sprintf("\nMCMC is done! \n"))
  
  out = NULL
  out$beta    = beta.save
  out$phi     = phi.save
  out$phi_hat = phi.hat.save
  out$theta   = theta.save
  out$Lambda  = Lambda.save
  out$psi     = psi.save
  out$alpha   = alpha.save
  out$tau2    = tau2.save
  out$z       = z.save
  out$K       = K
  
  return(out)
}


logML_hglmer_linear <- function(post_draws, ydata, xdata, K, a0, d0, r0, s0, u0, V0, w0) {
  ## Preliminary settings
  n      = length(ydata)
  p      = ncol(xdata[[1]])
  niter  = length(post_draws$beta)
  f0     = p + 1
  G0     = p*diag(rep(1,p))
  
  ## Calculate posterior estimates
  
  # 1. beta
  betas      = lapply(1:n, function(i) {
    do.call(rbind, lapply(post_draws$beta, function(M) M[i,,drop=FALSE]))})
  MU_beta    = lapply(betas, colMeans)
  Sigma_beta = lapply(betas, cov)
  
  # 2. z
  pi_z = colMeans(do.call(rbind, post_draws$psi))
  
  # 3. theta
  thetas      = lapply(1:K, function(k) {
    do.call(rbind, lapply(post_draws$theta, function(M) M[k,,drop=FALSE]))})
  MU_theta    = lapply(thetas, colMeans)
  Sigma_theta = lapply(thetas, cov)
  
  # 4. phi
  phis      = lapply(1:n, function(i) {
    sapply(post_draws$phi, function(x) x[i])})
  MU_phi    = lapply(phis, mean)
  Sigma_phi = lapply(phis, var)
  
  # 5. alpha and log-tau2
  alpha_ltau2s      = cbind(unlist(post_draws$alpha), log(unlist(post_draws$tau2)))
  MU_alpha_ltau2    = colMeans(alpha_ltau2s)
  Sigma_alpha_ltau2 = cov(alpha_ltau2s)
  
  # 6. Lambda
  MU_psi       = lapply(1:K, function(k) {
    mean(sapply(post_draws$psi, function(x) x[k]))})
  MU_Lambdainv = lapply(1:K, function(k) { apply(simplify2array(lapply(post_draws$Lambda, 
                                                                       function(M) solve(M[,,k]))), c(1,2), mean)})
  h_Lambda     = lapply(1:K, function(k) f0 + n * MU_psi[[k]])
  H_Lambda     = lapply(1:K, function(k) h_Lambda[[k]]*solve(MU_Lambdainv[[k]]))
  
  # 7. psi
  MU_psi    = lapply(1:K, function(k) {
    mean(sapply(post_draws$psi, function(x) x[k]))})
  Sigma_psi = lapply(1:K, function(k) {
    var(sapply(post_draws$psi, function(x) x[k]))})
  V_est     = (1-sum(unlist(MU_psi)^2))/sum(unlist(Sigma_psi))-1
  V_psi     = lapply(1:K, function(k) V_est*MU_psi[[k]])
  
  
  ## Calculate marginal likelihood
  g_vec = numeric(niter)
  
  for (u in 1:niter) {
    beta = post_draws$beta[[u]]
    for (i in 1:n) {
      g_vec[u]=g_vec[u]+dmvnorm(beta[i,],MU_beta[[i]],Sigma_beta[[i]],log=TRUE)
    }
    
    if (K > 1) {
      z_mat <- t(sapply(post_draws$z[[u]], function(k) {
        out <- numeric(K)
        out[k] <- 1
        out
      }))
      
      for (i in 1:n) {
        g_vec[u]=g_vec[u]+dmultinom(z_mat[i,],prob=pi_z,log=TRUE)
      }
    }
    
    theta = post_draws$theta[[u]]
    for (k in 1:K) {
      g_vec[u]=g_vec[u]+dmvnorm(theta[k,],MU_theta[[k]],Sigma_theta[[k]],log=TRUE)
    }
    
    phi = post_draws$phi[[u]]
    for (i in 1:n) {
      g_vec[u]=g_vec[u]+dnorm(phi[i],MU_phi[[i]],sqrt(Sigma_phi[[i]]),log=TRUE)
    }
    
    alpha_ltau2 = c(post_draws$alpha[[u]], log(post_draws$tau2[[u]]))
    g_vec[u]=g_vec[u]+dmvnorm(alpha_ltau2,MU_alpha_ltau2,Sigma_alpha_ltau2,log=TRUE)
    
    Lambda = post_draws$Lambda[[u]]
    for (i in 1:K) {
      g_vec[u]=g_vec[u]+log(diwish(Lambda[,,k],h_Lambda[[k]],H_Lambda[[k]]))
    }
    
    if (K > 1) {
      psi = post_draws$psi[[u]]
      g_vec[u]=g_vec[u]+log(ddirichlet(psi,unlist(V_psi)))
    }
  }
  
  llik_vec = numeric(niter)
  
  for (u in 1:niter) {
    beta   = post_draws$beta[[u]]
    z      = post_draws$z[[u]]
    phi    = post_draws$phi[[u]]
    psi    = post_draws$psi[[u]]
    theta  = post_draws$theta[[u]]
    Lambda = post_draws$Lambda[[u]]
    psi    = post_draws$psi[[u]]
    alpha  = post_draws$alpha[[u]]
    tau2   = post_draws$tau2[[u]]
    
    for (i in 1:n) {
      llik_vec[u]=llik_vec[u]+sum(dnorm(as.vector(ydata[[i]]),as.vector(xdata[[i]]%*%beta[i,]),
                                        sqrt(exp(phi[i])),log=TRUE))
    }
    
    for (i in 1:n) {
      llik_vec[u]=llik_vec[u]+dmvnorm(beta[i,],theta[z[i],],Lambda[,,z[i]],log=TRUE) +
        log(psi[z[i]])
    }
    
    for (k in 1:K) {
      llik_vec[u]=llik_vec[u]+dmvnorm(theta[k,],rep(u0,p),V0*diag(1,p),log=TRUE)
    }
    
    for (k in 1:K) {
      llik_vec[u]=llik_vec[u]+log(diwish(Lambda[,,k],f0,G0))
    }
    
    if (K > 1) {
      llik_vec[u]=llik_vec[u]+log(ddirichlet(psi,rep(w0,K)))+lfactorial(K)
    }
    
    for (i in 1:n) {
      llik_vec[u]=llik_vec[u]+dnorm(phi[i],alpha,sqrt(tau2),log=TRUE)
    }
    
    llik_vec[u]=llik_vec[u]+dnorm(alpha,a0,d0,log=TRUE)
    llik_vec[u]=llik_vec[u]+log(dinvgamma(tau2,shape=r0/2,scale=s0/2))
  }
  
  logint   = g_vec-llik_vec
  logint_m = max(logint)
  logML    = -(logint_m+log(mean(exp(logint-logint_m))))
  
  return(logML)
}












