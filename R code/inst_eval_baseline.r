# Basic Baseline Models for InstEval Data
library(lme4)
library(glmmTMB)
library(GLMMadaptive)
library(mclogit)
library(brms)
library(rstanarm)

data(InstEval)

min_ratings = 5
students = names(which(table(InstEval$s) >= min_ratings))
set.seed(123)
if(length(students) > 100) {
  students = sample(students, 100)
}

# Filter data for selected students and create binary outcome
insteval_subset = InstEval[InstEval$s %in% students,]
insteval_subset$binary_y = as.numeric(insteval_subset$y >= 4)

#-------------------------------------------------------
# Baseline 1: Fixed Effects Logistic Regression
#-------------------------------------------------------
fixed_model = glm(binary_y ~ service, 
                  data=insteval_subset, family=binomial)
summary(fixed_model)
AIC(fixed_model)
BIC(fixed_model)

#-------------------------------------------------------
# Baseline 2: Random Intercept Model
#-------------------------------------------------------
random_int_model = glmer(binary_y ~ service + (1|s), 
                         data=insteval_subset, family=binomial)
summary(random_int_model)
AIC(random_int_model)
BIC(random_int_model)

#-------------------------------------------------------
# Baseline 3: Random Intercept and Slope Model
#-------------------------------------------------------
random_slope_model = glmer(binary_y ~ service + (service|s), 
                           data=insteval_subset, family=binomial)
summary(random_slope_model)
AIC(random_slope_model)
BIC(random_slope_model)

#-------------------------------------------------------
# Baseline 4: Single-Component Model (K=1)
#-------------------------------------------------------
# This is essentially your hglmer model with K=1
# Run your existing code with K forced to 1
# Compare this with your mixture models with K>1

#-------------------------------------------------------
# Optional: Bayesian Logistic Mixed Model with brms
#-------------------------------------------------------
# Note: This can take time to run
bayesian_model = brm(
  binary_y ~ service + (service|s),
  data = insteval_subset,
  family = bernoulli(),
  prior = c(
    prior(normal(0, 3), class = "b"),
    prior(normal(0, 3), class = "sd"),
    prior(lkj(2), class = "cor")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4
)
summary(bayesian_model)

