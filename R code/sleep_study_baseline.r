# Baseline models for sleep study comparison
library(lme4) 
library(flexmix) 
library(ggplot2) 
library(MASS) 

data(sleepstudy)

cat("Number of subjects:", length(unique(sleepstudy$Subject)), "\n")
cat("Number of observations:", nrow(sleepstudy), "\n")

######################################
# Baseline Model 1: Pooled OLS
######################################
model_pooled <- lm(Reaction ~ Days, data = sleepstudy)
cat("\n==========================================\n")
cat("Baseline 1: Pooled OLS (ignoring subject effects)\n")
cat("==========================================\n")
print(summary(model_pooled))
cat("AIC:", AIC(model_pooled), "\n")
cat("BIC:", BIC(model_pooled), "\n")

######################################
# Baseline Model 2: Standard Random Effects
######################################
model_random <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
cat("\n==========================================\n")
cat("Baseline 2: Standard Random Effects Model\n")
cat("==========================================\n")
print(summary(model_random))
cat("AIC:", AIC(model_random), "\n")
cat("BIC:", BIC(model_random), "\n")

var_comps <- VarCorr(model_random)
cat("\nRandom Effects Variance Components:\n")
print(var_comps)

cat("\nFixed Effects (Population Average):\n")
print(fixef(model_random))

ranef_vals <- ranef(model_random)$Subject
cat("\nRandom Effects Ranges:\n")
print(apply(ranef_vals, 2, range))

subj_coefs <- coef(model_random)$Subject
colnames(subj_coefs) <- c("Intercept", "Slope")
cat("\nIndividual Subject Coefficients (first 5):\n")
print(head(subj_coefs, 5))

######################################
# Baseline Model 3: Classical Finite Mixture
######################################
cat("\n==========================================\n")
cat("Baseline 3: Classical Finite Mixture Models\n")
cat("==========================================\n")

set.seed(123) 
flexmix_models <- list()
flexmix_BICs <- numeric(5)

for (k in 1:5) {
  if (k == 1) {
    # For k=1, just use standard linear model
    flexmix_BICs[k] <- BIC(model_pooled)
    cat("FlexMix K=1 (equivalent to pooled OLS) - BIC:", flexmix_BICs[k], "\n")
  } else {
    flexmix_models[[k]] <- flexmix(Reaction ~ Days, data = sleepstudy, k = k)
    flexmix_BICs[k] <- BIC(flexmix_models[[k]])
    cat("FlexMix K=", k, "- BIC:", flexmix_BICs[k], "\n")
    
    cat("Component parameters:\n")
    print(parameters(flexmix_models[[k]]))
    
    cluster_sizes <- table(clusters(flexmix_models[[k]]))
    cat("Component sizes:", cluster_sizes, "\n\n")
  }
}

best_k <- which.min(flexmix_BICs)
cat("Best FlexMix model has", best_k, "components based on BIC\n")
