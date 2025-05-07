setwd("/Users/nydu/Desktop/stat526proj/project/Hglmer_Lenk/R Code")
#########################
# InstEval Data Analysis #
#########################
rm(list=ls())
source("hglmer_logistic.R")
library(mvtnorm)
library(lme4)
library(MCMCpack)
library(boot)

data(InstEval)

min_ratings = 5 
students = names(which(table(InstEval$s) >= min_ratings))
set.seed(123)
# if(length(students) > 100) {
#   students = sample(students, 100)
# }

set.seed(456)
train_prop <- 0.8
n_students <- length(students)
n_train <- round(n_students * train_prop)
train_students <- sample(students, n_train)
test_students <- setdiff(students, train_students)

prepare_data <- function(student_ids, min_ratings=5) {
  ydata = vector(mode='list', length=0)
  xdata = vector(mode='list', length=0)
  ids = c()
  
  for(i in 1:length(student_ids)) {
    student_id = student_ids[i]
    student_data = InstEval[InstEval$s == student_id,]
    
    if(nrow(student_data) < min_ratings) next
    if(length(unique(student_data$service)) < 2) next
    
    y = as.numeric(student_data$y >= 4)
    
    X = cbind(1, as.numeric(student_data$service == "1"))
    ydata = c(ydata, list(y))
    xdata = c(xdata, list(X))
    ids = c(ids, student_id)
  }
  
  return(list(ydata=ydata, xdata=xdata, ids=ids))
}

train <- prepare_data(train_students)
test <- prepare_data(test_students)

insteval_train <- InstEval[InstEval$s %in% train$ids,]
insteval_train$binary_y <- as.numeric(insteval_train$y >= 4)

insteval_test <- InstEval[InstEval$s %in% test$ids,]
insteval_test$binary_y <- as.numeric(insteval_test$y >= 4)

cat("Training set:", length(train$ids), "students,", length(unlist(train$ydata)), "observations\n")
cat("Testing set:", length(test$ids), "students,", length(unlist(test$ydata)), "observations\n\n")

#################################################
# FIT BASELINE MODELS
#################################################

cat("Fitting baseline models...\n")
cat("Fitting Fixed Effects model...\n")
fixed_model <- glm(binary_y ~ service, data=insteval_train, family=binomial)
cat("Fixed Effects - AIC:", AIC(fixed_model), "BIC:", BIC(fixed_model), "\n")
cat("Fitting Random Intercept model...\n")
random_int_model <- glmer(binary_y ~ service + (1|s), 
                          data=insteval_train, family=binomial,
                          control=glmerControl(optimizer="bobyqa"))
cat("Random Intercept - AIC:", AIC(random_int_model), "BIC:", BIC(random_int_model), "\n")
cat("Fitting Random Slope model...\n")
tryCatch({
  random_slope_model <- glmer(binary_y ~ service + (service|s), 
                             data=insteval_train, family=binomial,
                             control=glmerControl(optimizer="bobyqa"))
  cat("Random Slope - AIC:", AIC(random_slope_model), "BIC:", BIC(random_slope_model), "\n")
}, error = function(e) {
  cat("Could not fit Random Slope model. Using simpler structure...\n")
  random_slope_model <<- glmer(binary_y ~ service + (1|s) + (0+service|s), 
                             data=insteval_train, family=binomial,
                             control=glmerControl(optimizer="bobyqa"))
  cat("Random Slope (uncorrelated) - AIC:", AIC(random_slope_model), "BIC:", BIC(random_slope_model), "\n")
})

calculate_baseline_logscore <- function(model_type, test_data) {
  if (model_type == "fixed") {
    pred_probs <- predict(fixed_model, newdata=test_data, type="response")
    pred_class <- as.numeric(pred_probs > 0.5)
  } else if (model_type == "random_int") {
    pred_probs <- predict(random_int_model, newdata=test_data, type="response", allow.new.levels=TRUE)
    pred_class <- as.numeric(pred_probs > 0.5)
  } else if (model_type == "random_slope") {
    pred_probs <- predict(random_slope_model, newdata=test_data, type="response", allow.new.levels=TRUE)
    pred_class <- as.numeric(pred_probs > 0.5)
  }
  
  actual <- test_data$binary_y
  rmse <- sqrt(mean((actual - pred_probs)^2))
  
  log_scores <- log(pred_probs*actual + (1-pred_probs)*(1-actual))
  mean_log_score <- mean(log_scores)
  accuracy <- mean(pred_class == actual)
  
  return(list(
    predictions = pred_probs,
    log_scores = log_scores,
    rmse = rmse,
    mean_log_score = mean_log_score,
    accuracy = accuracy
  ))
}
calculate_hglmer_logscore <- function(model, test_data) {
  K <- model$K
  n_iter <- length(model$theta)
  n_subjects <- length(test_data$ydata)
  all_preds <- list()
  all_log_scores <- list()
  
  for (i in 1:n_subjects) {
    x_i <- test_data$xdata[[i]]
    y_i <- test_data$ydata[[i]]
    n_obs <- length(y_i)
    
    preds_i <- numeric(n_obs)
    log_scores_i <- numeric(n_obs)

    for (j in 1:n_obs) {
      x_ij <- x_i[j,]
      y_ij <- y_i[j]

      pred_prob_sum <- 0
      
      for (iter in 1:n_iter) {
        theta_samples <- model$theta[[iter]] 
        Lambda_samples <- model$Lambda[[iter]]  
        psi_samples <- model$psi[[iter]] 
        
        iter_prob <- 0
        
        for (k in 1:K) {
          theta_k <- theta_samples[k,]
          Lambda_k <- Lambda_samples[,,k]
          psi_k <- psi_samples[k]
          
          beta_samples <- rmvnorm(10, mean = theta_k, sigma = Lambda_k)

          comp_probs <- numeric(10)
          for (s in 1:10) {
            beta_s <- beta_samples[s,]
            linear_pred <- sum(x_ij * beta_s)
            prob_s <- 1 / (1 + exp(-linear_pred))
            comp_probs[s] <- prob_s
          }
          comp_prob <- mean(comp_probs)
        
          iter_prob <- iter_prob + psi_k * comp_prob
        }
        
        pred_prob_sum <- pred_prob_sum + iter_prob
      }
      pred_prob <- pred_prob_sum / n_iter
      preds_i[j] <- pred_prob
      log_scores_i[j] <- log(pred_prob*y_ij + (1-pred_prob)*(1-y_ij))
    }
    
    all_preds[[i]] <- preds_i
    all_log_scores[[i]] <- log_scores_i
  }
  
  y_true <- unlist(test_data$ydata)
  y_pred <- unlist(all_preds)
  log_scores <- unlist(all_log_scores)
  
  rmse <- sqrt(mean((y_true - y_pred)^2))
  mean_log_score <- mean(log_scores)
  accuracy <- mean((y_pred > 0.5) == y_true)
  
  return(list(
    predictions = y_pred,
    true_values = y_true,
    log_scores = log_scores,
    rmse = rmse,
    mean_log_score = mean_log_score,
    accuracy = accuracy
  ))
}

cat("\nCalculating predictive metrics for baseline models...\n")
fixed_metrics <- calculate_baseline_logscore("fixed", insteval_test)
cat("Fixed Effects - RMSE:", fixed_metrics$rmse, "\n")
cat("Fixed Effects - Mean Log Score:", fixed_metrics$mean_log_score, "\n")
cat("Fixed Effects - Accuracy:", fixed_metrics$accuracy, "\n")

random_int_metrics <- calculate_baseline_logscore("random_int", insteval_test)
cat("Random Intercept - RMSE:", random_int_metrics$rmse, "\n")
cat("Random Intercept - Mean Log Score:", random_int_metrics$mean_log_score, "\n")
cat("Random Intercept - Accuracy:", random_int_metrics$accuracy, "\n")

random_slope_metrics <- calculate_baseline_logscore("random_slope", insteval_test)
cat("Random Slope - RMSE:", random_slope_metrics$rmse, "\n")
cat("Random Slope - Mean Log Score:", random_slope_metrics$mean_log_score, "\n")
cat("Random Slope - Accuracy:", random_slope_metrics$accuracy, "\n")


u0 = 0 
V0 = 100 
w0 = 1 
Kvalues <- list(1, 2, 3, 4, NULL)
set.seed(123)
results <- list()
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
  fit_insteval <- gibbs_hglmer_logistic(
    ydata = train$ydata,
    xdata = train$xdata,
    K = K,
    u0 = u0,
    V0 = V0,
    w0 = w0,
    mcmc = list(nblow = 2000, smcmc = 2000)
  )
  
  actual_K <- fit_insteval$K
  if (is.null(K)) {
    cat("Automatic selection chose K =", actual_K, "\n")
    model_labels[i] <- paste0("HGLMER (Auto, K=", actual_K, ")")
    model_tag <- "Auto"
  } else {
    model_tag <- paste0("K", K)
  }
  
  results[[i]] <- fit_insteval
  
  cat("Calculating predictive metrics...\n")
  pred_metrics[[i]] <- calculate_hglmer_logscore(fit_insteval, test)
  cat(model_labels[i], "- RMSE:", pred_metrics[[i]]$rmse, "\n")
  cat(model_labels[i], "- Mean Log Score:", pred_metrics[[i]]$mean_log_score, "\n")
  cat(model_labels[i], "- Accuracy:", pred_metrics[[i]]$accuracy, "\n")
  dir.create("./plots/insteval", recursive = TRUE, showWarnings = FALSE)
  
  result_psi <- do.call(rbind, fit_insteval$psi)
  psi_means <- apply(result_psi, 2, mean)
  cat("Estimated mixture proportions:\n")
  print(psi_means)
  
  png(paste0("./plots/insteval/psi_plot_", model_tag, ".png"), 
    width=800, height=300*actual_K, res=100)
  par(mfrow=c(actual_K, 1), mar=c(3, 3, 2, 1)) 
  for (j in 1:actual_K) {
    plot(result_psi[,j], type='l', 
         main=paste0("Component ", j, " Weight (Mean = ", round(psi_means[j], 3), ")"),
         xlab="MCMC Iteration", ylab="Mixture Weight")
  }
  dev.off()
  result_theta <- fit_insteval$theta
  theta_list <- list()
  for (j in 1:actual_K) {
    theta_list[[j]] <- do.call(rbind, lapply(result_theta, function(x) x[j,]))
  }
  means_by_component <- lapply(theta_list, function(x) apply(x, 2, mean))
  for (j in 1:actual_K) {
    cat("Component", j, "means (intercept, service):", means_by_component[[j]], "\n")
  }
  
  png(paste0("./plots/insteval/theta_plot_", model_tag, ".png"), 
      width=900, height=250*actual_K, res=100)
  par(mfrow=c(actual_K, 2), mar=c(4,4,2,1))
  for (j in 1:actual_K) {
    int_mean <- mean(theta_list[[j]][,1])
    plot(theta_list[[j]][,1], type='l', 
         main=paste0("Component ", j, " Intercept (Mean = ", round(int_mean, 2), ")"),
         xlab="MCMC Iteration", ylab="Intercept")
    abline(h=int_mean, col="red", lty=2)
    
    effect_mean <- mean(theta_list[[j]][,2])
    plot(theta_list[[j]][,2], type='l', 
         main=paste0("Component ", j, " Service Effect (Mean = ", round(effect_mean, 2), ")"),
         xlab="MCMC Iteration", ylab="Service Effect")
    abline(h=effect_mean, col="red", lty=2)
  }
  dev.off()
  rm(fit_insteval)
  gc()
}

model_names <- c("Fixed Effects", "Random Intercept", "Random Slope", model_labels)

rmse_values <- c(
  fixed_metrics$rmse,
  random_int_metrics$rmse,
  random_slope_metrics$rmse,
  sapply(1:length(Kvalues), function(i) pred_metrics[[i]]$rmse)
)

log_score_values <- c(
  fixed_metrics$mean_log_score,
  random_int_metrics$mean_log_score,
  random_slope_metrics$mean_log_score,
  sapply(1:length(Kvalues), function(i) pred_metrics[[i]]$mean_log_score)
)

accuracy_values <- c(
  fixed_metrics$accuracy,
  random_int_metrics$accuracy,
  random_slope_metrics$accuracy,
  sapply(1:length(Kvalues), function(i) pred_metrics[[i]]$accuracy)
)

comparison_table <- data.frame(
  Model = model_names,
  RMSE = rmse_values,
  Mean_Log_Score = log_score_values,
  Accuracy = accuracy_values
)

cat("\n=========== Model Comparison ==========\n")
print(comparison_table)
write.csv(comparison_table, "insteval_model_comparison.csv", row.names = FALSE)
