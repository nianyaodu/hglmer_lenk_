setwd("/Users/nydu/Desktop/stat526proj/project/Hglmer_Lenk/R Code")

########################
# Sleep Study Analysis 
########################
rm(list=ls())
source("hglmer_linear.R")
library(mvtnorm)
library(lme4)
library(MCMCpack)
library(boot) 

data(sleepstudy)

set.seed(456)
train_prop <- 0.8 
subjects <- unique(sleepstudy$Subject)
n_subjects <- length(subjects)
n_train <- round(n_subjects * train_prop)
train_subjects <- sample(subjects, n_train)
test_subjects <- setdiff(subjects, train_subjects)

train_data <- sleepstudy[sleepstudy$Subject %in% train_subjects, ]
test_data <- sleepstudy[sleepstudy$Subject %in% test_subjects, ]
prepare_data <- function(data) {
  subjects <- unique(data$Subject)
  n <- length(subjects)
  ydata <- vector(mode='list', length=n)
  xdata <- vector(mode='list', length=n)
  
  for(i in 1:n) {
    subj_data <- data[data$Subject == subjects[i],]
    X <- cbind(1, subj_data$Days)
    y <- subj_data$Reaction
    
    xdata[[i]] <- X
    ydata[[i]] <- y
  }
  
  return(list(ydata=ydata, xdata=xdata, subjects=subjects))
}
train <- prepare_data(train_data)
test <- prepare_data(test_data)
cat("Data split summary:\n")
cat("Training set:", length(train$subjects), "subjects,", nrow(train_data), "observations\n")
cat("Testing set:", length(test$subjects), "subjects,", nrow(test_data), "observations\n\n")

cat("\nFitting Pooled OLS model...\n")
model_pooled <- lm(Reaction ~ Days, data = train_data)
summary_pooled <- summary(model_pooled)
cat("Pooled OLS coefficients:\n")
print(coef(model_pooled))
cat("Residual standard error:", summary_pooled$sigma, "\n")
cat("R-squared:", summary_pooled$r.squared, "\n")
cat("\nFitting Random Effects model...\n")
model_random <- lmer(Reaction ~ Days + (Days | Subject), data = train_data)
summary_random <- summary(model_random)
cat("Random Effects fixed coefficients:\n")
print(fixef(model_random))
cat("Random Effects variance components:\n")
print(VarCorr(model_random))
cat("Residual standard error:", attr(VarCorr(model_random), "sc"), "\n")

n_bootstrap <- 1000

cat("\nGenerating bootstrap samples for OLS model...\n")
boot_ols <- boot::boot(data = train_data, statistic = function(data, indices) {
  model_boot <- lm(Reaction ~ Days, data = data[indices,])
  c(coef(model_boot), summary(model_boot)$sigma)
}, R = n_bootstrap)
ols_coef_boot <- boot_ols$t[,1:2] 
ols_sigma_boot <- boot_ols$t[,3] 

cat("\nGenerating bootstrap samples for Random Effects model...\n")
re_coef_boot <- matrix(0, nrow = n_bootstrap, ncol = 2)
re_var_boot <- array(0, dim = c(n_bootstrap, 2, 2))
re_sigma_boot <- numeric(n_bootstrap)

for (b in 1:n_bootstrap) {
  if (b %% 100 == 0) cat("  Bootstrap sample", b, "of", n_bootstrap, "\n")
  sim_data <- simulate(model_random, newdata = train_data, re.form = NULL)
  names(sim_data) <- "Reaction"
  boot_data <- cbind(train_data[,c("Subject", "Days")], sim_data)
  tryCatch({
    boot_model <- lmer(Reaction ~ Days + (Days | Subject), data = boot_data)
    re_coef_boot[b,] <- fixef(boot_model)
    vcov_comps <- VarCorr(boot_model)
    re_var_boot[b,1,1] <- as.numeric(vcov_comps$Subject[1,1]) 
    re_var_boot[b,2,2] <- as.numeric(vcov_comps$Subject[2,2]) 
    re_var_boot[b,1,2] <- re_var_boot[b,2,1] <- as.numeric(vcov_comps$Subject[1,2]) 
    re_sigma_boot[b] <- attr(vcov_comps, "sc")
  }, error = function(e) {
    re_coef_boot[b,] <- fixef(model_random)
    
    vcov_comps <- VarCorr(model_random)
    re_var_boot[b,1,1] <- as.numeric(vcov_comps$Subject[1,1])
    re_var_boot[b,2,2] <- as.numeric(vcov_comps$Subject[2,2])
    re_var_boot[b,1,2] <- re_var_boot[b,2,1] <- as.numeric(vcov_comps$Subject[1,2])
    
    re_sigma_boot[b] <- attr(vcov_comps, "sc")
  })
}

calculate_bootstrap_logscore <- function(model_type, test_data, n_bootstrap) {
  test_df <- data.frame()
  for (i in 1:length(test_data$subjects)) {
    subj_data <- data.frame(
      Subject = rep(test_data$subjects[i], nrow(test_data$xdata[[i]])),
      Days = test_data$xdata[[i]][,2],
      Reaction = test_data$ydata[[i]]
    )
    test_df <- rbind(test_df, subj_data)
  }
  
  X_test <- cbind(1, test_df$Days)
  y_test <- test_df$Reaction
  n_test <- nrow(test_df)
  log_scores <- numeric(n_test)
  
  if (model_type == "pooled") {
    for (i in 1:n_test) {
      x_i <- X_test[i,]
      y_i <- y_test[i]
      
      densities <- numeric(n_bootstrap)
      
      for (b in 1:n_bootstrap) {
        beta_b <- ols_coef_boot[b,]
        sigma_b <- ols_sigma_boot[b]
        mu_b <- sum(x_i * beta_b)

        densities[b] <- dnorm(y_i, mean = mu_b, sd = sigma_b)
      }
      
      log_scores[i] <- log(max(mean(densities), .Machine$double.eps))
    }
    
  } else if (model_type == "random") {
    for (i in 1:n_test) {
      x_i <- X_test[i,]
      y_i <- y_test[i]
      densities <- numeric(n_bootstrap)
      
      for (b in 1:n_bootstrap) {
        beta_b <- re_coef_boot[b,]
        mu_b <- sum(x_i * beta_b)
        var_b <- re_sigma_boot[b]^2 +  
                re_var_boot[b,1,1] +   
                re_var_boot[b,2,2] * x_i[2]^2 +  
                2 * re_var_boot[b,1,2] * x_i[2]  
        densities[b] <- dnorm(y_i, mean = mu_b, sd = sqrt(var_b))
      }
      log_scores[i] <- log(max(mean(densities), .Machine$double.eps))
    }
  }
  if (model_type == "pooled") {
    y_pred <- X_test %*% colMeans(ols_coef_boot)
  } else {
    y_pred <- X_test %*% colMeans(re_coef_boot)
  }
  rmse <- sqrt(mean((y_test - y_pred)^2))
  
  return(list(
    log_scores = log_scores,
    mean_log_score = mean(log_scores),
    rmse = rmse
  ))
}
cat("\nCalculating predictive metrics for Pooled OLS...\n")
pooled_metrics <- calculate_bootstrap_logscore("pooled", test, n_bootstrap)
cat("Pooled OLS - RMSE:", pooled_metrics$rmse, "\n")
cat("Pooled OLS - Mean Log Score:", pooled_metrics$mean_log_score, "\n")

cat("\nCalculating predictive metrics for Random Effects...\n")
random_metrics <- calculate_bootstrap_logscore("random", test, n_bootstrap)
cat("Random Effects - RMSE:", random_metrics$rmse, "\n")
cat("Random Effects - Mean Log Score:", random_metrics$mean_log_score, "\n")

calculate_hglmer_logscore <- function(model, test_data) {
  K <- model$K
  n_iter <- length(model$theta)
  n_subjects <- length(test_data$subjects)

  all_preds <- list()
  all_log_scores <- list()

  for (i in 1:n_subjects) {
    x_i <- test_data$xdata[[i]]
    y_i <- test_data$ydata[[i]]
    n_obs <- length(y_i)
    
    preds_i <- matrix(0, nrow=n_obs, ncol=1)
    log_scores_i <- numeric(n_obs)
    pred_densities <- matrix(0, nrow=n_obs, ncol=n_iter)
    for (j in 1:n_obs) {
      x_ij <- x_i[j,]
      y_ij <- y_i[j]
      for (iter in 1:n_iter) {
        theta_samples <- model$theta[[iter]] 
        Lambda_samples <- model$Lambda[[iter]]
        psi_samples <- model$psi[[iter]] 
        phi_samples <- model$phi[[iter]][i] 
        iter_mean <- 0
        iter_density <- 0
        
        for (k in 1:K) {
          theta_k <- theta_samples[k,]
          Lambda_k <- Lambda_samples[,,k]
          psi_k <- psi_samples[k]
          mu_ijk <- sum(x_ij * theta_k)
          iter_mean <- iter_mean + psi_k * mu_ijk
          var_ijk <- t(x_ij) %*% Lambda_k %*% x_ij + exp(phi_samples)
          comp_density <- dnorm(y_ij, mean=mu_ijk, sd=sqrt(var_ijk))
          iter_density <- iter_density + psi_k * comp_density
        }
        preds_i[j] <- preds_i[j] + iter_mean/n_iter
        pred_densities[j, iter] <- iter_density
      }
      max_log_density <- max(log(pred_densities[j,]))
      log_scores_i[j] <- max_log_density + 
        log(mean(exp(log(pred_densities[j,]) - max_log_density)))
    }
    
    all_preds[[i]] <- preds_i
    all_log_scores[[i]] <- log_scores_i
  }
  y_true <- unlist(test_data$ydata)
  y_pred <- unlist(all_preds)
  log_scores <- unlist(all_log_scores)
  rmse <- sqrt(mean((y_true - y_pred)^2))
  mean_log_score <- mean(log_scores)
  
  return(list(
    predictions = y_pred,
    true_values = y_true,
    log_scores = log_scores,
    rmse = rmse,
    mean_log_score = mean_log_score
  ))
}

Kvalues <- list(1, 2, 3, 4, NULL)
set.seed(123)
results <- list()
logML_values <- list()
pred_metrics <- list()
model_labels <- c() 

for (i in 1:length(Kvalues)) {
  K <- Kvalues[[i]]
  if (is.null(K)) {
    model_label <- "HGLMER (Auto)"
    cat("\nFitting mixture model with automatic K selection\n")
  } else {
    model_label <- paste0("HGLMER (K=", K, ")")
    cat("\nFitting mixture model with K =", K, "components\n")
  }
  model_labels <- c(model_labels, model_label)
  fit_sleepstudy <- gibbs_hglmer_linear(
    ydata = train$ydata,
    xdata = train$xdata,
    K = K,
    a0 = 0,
    d0 = sqrt(10),
    r0 = 1,
    s0 = 1,
    u0 = 0,
    V0 = 100,
    w0 = 1,
    mcmc = list(nblow = 2000, smcmc = 1000)
  )
  actual_K <- fit_sleepstudy$K
  if (is.null(K)) {
    cat("Automatic selection chose K =", actual_K, "\n")
    model_labels[i] <- paste0("HGLMER (Auto, K=", actual_K, ")")
    model_tag <- "Auto"
  } else {
    model_tag <- paste0("K", K)
  }
  results[[i]] <- fit_sleepstudy
  cat("Calculating log marginal likelihood...\n")
  logML_value <- logML_hglmer_linear(
    post_draws = fit_sleepstudy,
    ydata = train$ydata,
    xdata = train$xdata,
    K = actual_K, 
    a0 = 0,
    d0 = sqrt(10),
    r0 = 1,
    s0 = 1,
    u0 = 0,
    V0 = 100,
    w0 = 1
  )
  logML_values[[i]] <- logML_value
  cat("Log marginal likelihood for", model_labels[i], ":", logML_value, "\n")
  cat("Calculating predictive metrics...\n")
  pred_metrics[[i]] <- calculate_hglmer_logscore(fit_sleepstudy, test)
  cat(model_labels[i], "- RMSE:", pred_metrics[[i]]$rmse, "\n")
  cat(model_labels[i], "- Mean Log Score:", pred_metrics[[i]]$mean_log_score, "\n")
  dir.create("./plots/sleep_study", recursive = TRUE, showWarnings = FALSE)
  result_psi <- do.call(rbind, fit_sleepstudy$psi)
  psi_means <- apply(result_psi, 2, mean)
  cat("Estimated mixture proportions:\n")
  print(psi_means)
  
  png(paste0("./plots/sleep_study/psi_plot_", model_tag, ".png"), 
    width=800, height=300*actual_K, res=100)
  par(mfrow=c(actual_K, 1), mar=c(3, 3, 2, 1)) 
  for (j in 1:actual_K) {
    plot(result_psi[,j], type='l', 
         main=paste0("Component ", j, " Weight (Mean = ", round(psi_means[j], 3), ")"),
         xlab="MCMC Iteration", ylab="Mixture Weight")
  }
  dev.off()
  result_theta <- fit_sleepstudy$theta
  theta_list <- list()
  for (j in 1:actual_K) {
    theta_list[[j]] <- do.call(rbind, lapply(result_theta, function(x) x[j,]))
  }
  means_by_component <- lapply(theta_list, function(x) apply(x, 2, mean))
  for (j in 1:actual_K) {
    cat("Component", j, "means (intercept, days):", means_by_component[[j]], "\n")
  }
  
  png(paste0("./plots/sleep_study/theta_plot_", model_tag, ".png"), 
      width=900, height=250*actual_K, res=100)
  par(mfrow=c(actual_K, 2), mar=c(4,4,2,1))
  for (j in 1:actual_K) {
    int_mean <- mean(theta_list[[j]][,1])
    plot(theta_list[[j]][,1], type='l', 
         main=paste0("Component ", j, " Intercept (Mean = ", round(int_mean, 2), ")"),
         xlab="MCMC Iteration", ylab="Intercept")
    abline(h=int_mean, col="red", lty=2)
    slope_mean <- mean(theta_list[[j]][,2])
    plot(theta_list[[j]][,2], type='l', 
         main=paste0("Component ", j, " Slope (Mean = ", round(slope_mean, 2), ")"),
         xlab="MCMC Iteration", ylab="Slope")
    abline(h=slope_mean, col="red", lty=2)
  }
  dev.off()
  rm(fit_sleepstudy)
  gc()
}

model_names <- c("Pooled OLS", "Random Effects", model_labels)

rmse_values <- c(
  pooled_metrics$rmse,
  random_metrics$rmse,
  sapply(1:length(Kvalues), function(i) pred_metrics[[i]]$rmse)
)
log_score_values <- c(
  pooled_metrics$mean_log_score,
  random_metrics$mean_log_score,
  sapply(1:length(Kvalues), function(i) pred_metrics[[i]]$mean_log_score)
)

comparison_table <- data.frame(
  Model = model_names,
  RMSE = rmse_values,
  Mean_Log_Score = log_score_values
)

cat("\n=========== Model Comparison ==========\n")
print(comparison_table)
write.csv(comparison_table, "model_comparison_results.csv", row.names = FALSE)