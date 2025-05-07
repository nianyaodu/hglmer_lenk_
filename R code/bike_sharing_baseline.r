#########################
# Baseline Models for Bike Sharing Dataset
#########################
library(MASS)
library(glmnet) 
library(mgcv) 
library(randomForest)
library(xgboost)

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

combined_data$is_morning = ifelse(combined_data$hour >= 6 & combined_data$hour <= 10, 1, 0)
combined_data$is_evening = ifelse(combined_data$hour >= 16 & combined_data$hour <= 20, 1, 0)
combined_data$is_night = ifelse(combined_data$hour >= 22 | combined_data$hour <= 5, 1, 0)

set.seed(123)
train_days = sample(1:n, 0.7*n)
test_days = setdiff(1:n, train_days)

train_idx = which(combined_data$day_id %in% train_days)
test_idx = which(combined_data$day_id %in% test_days)

train_data = combined_data[train_idx, ]
test_data = combined_data[test_idx, ]

cat("Fitting Poisson Regression model...\n")
poisson_model = glm(response ~ temp + humidity + windspeed, 
                   family = poisson(),
                   data = train_data)

cat("Fitting Poisson Regression with hour factor...\n")
poisson_hour_model = glm(response ~ temp + humidity + windspeed + hour_factor, 
                        family = poisson(),
                        data = train_data)

cat("Fitting Poisson Regression with interactions...\n")
poisson_interact_model = glm(response ~ temp * humidity * windspeed + hour_factor, 
                           family = poisson(),
                           data = train_data)

cat("Fitting Negative Binomial model...\n")
nb_model = glm.nb(response ~ temp + humidity + windspeed + hour_factor, 
                 data = train_data)

cat("Fitting GAM model...\n")
gam_model = gam(response ~ s(temp) + s(humidity) + s(windspeed) + hour_factor, 
               family = poisson(),
               data = train_data)


cat("Fitting Random Forest model...\n")
train_data$log_response = log(train_data$response + 1)
test_data$log_response = log(test_data$response + 1)

rf_model = randomForest(log_response ~ temp + humidity + windspeed + hour_factor + 
                          is_morning + is_evening + is_night,
                       data = train_data,
                       ntree = 100)
cat("Fitting XGBoost model...\n")
features = c("temp", "humidity", "windspeed", "hour", "is_morning", "is_evening", "is_night")
dtrain = xgb.DMatrix(data = as.matrix(train_data[, features]), 
                    label = train_data$response)
dtest = xgb.DMatrix(data = as.matrix(test_data[, features]), 
                   label = test_data$response)

xgb_params = list(
  objective = "count:poisson",
  eta = 0.1,
  max_depth = 6,
  gamma = 1,
  subsample = 0.8,
  colsample_bytree = 0.8
)

xgb_model = xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = 100,
  watchlist = list(train = dtrain, test = dtest),
  early_stopping_rounds = 10,
  verbose = 0
)

evaluate_model = function(model, test_data, model_type) {
  if (model_type == "xgboost") {
    predictions = predict(model, dtest)
  } else if (model_type == "rf") {
    log_pred = predict(model, test_data)
    predictions = exp(log_pred) - 1
  } else {
    predictions = predict(model, test_data, type = "response")
  }
  
  rmse = sqrt(mean((predictions - test_data$response)^2))
  mae = mean(abs(predictions - test_data$response))
  null_dev = sum((test_data$response - mean(test_data$response))^2 / mean(test_data$response))
  res_dev = sum((test_data$response - predictions)^2 / predictions)
  explained_deviance = 1 - res_dev/null_dev
  
  return(data.frame(
    Model = model_type,
    RMSE = rmse,
    MAE = mae,
    ExplainedDeviance = explained_deviance
  ))
}
results = rbind(
  evaluate_model(poisson_model, test_data, "Poisson"),
  evaluate_model(poisson_hour_model, test_data, "PoissonWithHour"),
  evaluate_model(poisson_interact_model, test_data, "PoissonInteract"),
  evaluate_model(nb_model, test_data, "NegativeBinomial"),
  evaluate_model(gam_model, test_data, "GAM"),
  evaluate_model(rf_model, test_data, "RandomForest") #, 
  # evaluate_model(xgb_model, test_data, "XGBoost")
)
print(results)

# best model based on RMSE
best_model_name = results$Model[which.min(results$RMSE)]
cat("Best model based on RMSE:", best_model_name, "\n")

pred = predict(poisson_interact_model, test_data, type = "response")
plot(test_data$response[1:100], type = "l", col = "blue", 
     ylim = range(c(test_data$response[1:100], pred[1:100])),
     xlab = "Observation", ylab = "Bike Count",
     main = paste("Actual vs Predicted:", best_model_name))
lines(pred[1:100], col = "red")
legend("topright", legend = c("Actual", "Predicted"), col = c("blue", "red"), lty = 1)
