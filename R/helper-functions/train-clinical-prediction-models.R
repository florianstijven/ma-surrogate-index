# Function to emulate a prediction that just retuns the surrogate. 
train_clinical_prediction_model_surrogate = function(simulated_data) {
  return(function(newdata) {
    return(newdata$surrogate)
  })
}

# Function to train the clinical prediction model based on linear regression for
# the proof-of-concept scenario.
train_clinical_prediction_model_lm <- function(simulated_data) {
  # Train a correctly specified linear regression model for predicting the
  # clinical endpoint given the surrogate and baseline covariate
  clinical_model <- glm(
    formula = clinical ~ surrogate*X1 + surrogate ^ 2,
    data = simulated_data,
    x = FALSE,
    y = FALSE,
    family = gaussian()
  )
  
  return(function(newdata) {
    return(predict(clinical_model, newdata = newdata))
  })
}

# Function to train the clinical prediction mode based on logistic regression
# for the vaccine scenario.
train_clinical_prediction_model_logistic = function(simulated_data) {
  # Train a correctly specified logistic regression model for predicting the
  # clinical endpoint given the surrogate and baseline covariate
  clinical_model <- glm(
    formula = clinical ~ surrogate*X1 + surrogate*X2 + surrogate*X3,
    data = simulated_data,
    x = FALSE,
    y = FALSE,
    family = binomial()
  )
  
  return(function(newdata) {
    return(predict(clinical_model, newdata = newdata, type = "response"))
  })
}

# Function to train the clinical prediction mode based on GAM logistic regression
# for the vaccine scenario.
train_clinical_prediction_model_gam = function(simulated_data) {
  # Fit GAM logistic model.
  clinical_model <- mgcv::gam(
    formula = clinical ~ s(
      surrogate,
      by = X3,
      bs = "cr",
      k = 4
    ) +
      s(
        surrogate,
        by = 1 - X3,
        bs = "cr",
        k = 4
      ) +
      s(X1, bs = "cr", k = 4) + s(X2, bs = "cr", k = 4) + X3,
    data = simulated_data,
    family = binomial()
  )
  
  return(function(newdata) {
    return(mgcv::predict.gam(clinical_model, newdata = newdata, type = "response"))
  })
}

# Function to train the clinical prediction mode based on HAL logistic regression
# for the vaccine scenario.
train_clinical_prediction_model_hal = function(simulated_data) {
  # Fit a GLM to get the correct model matrix.
  glm_model <- glm(
    formula = clinical ~ surrogate + X1 + X2 + X3,
    data = simulated_data,
    family = binomial()
  )
  # Set tuning parameters for HAL. We follow the recommended settings for a fast
  # runtime (see documentation).
  smoothness_orders = 0
  max_degree = 3
  num_knots = c(100, 50, 25)
  
  # Train a correctly specified logistic regression model for predicting the
  # clinical endpoint given the surrogate and baseline covariate
  clinical_model <- hal9001::fit_hal(
    Y = simulated_data$clinical,
    X = model.matrix(glm_model)[, -1],
    family = "binomial",
    smoothness_orders = smoothness_orders,
    max_degree = max_degree,
    num_knots = num_knots,
    return_lasso = FALSE,
    return_x_basis = FALSE
  )
  
  print(clinical_model)
  
  return(function(newdata) {
    return(predict(clinical_model, new_data = newdata, type = "response"))
  })
}

# Function to train the clinical prediction mode based on random forests.
train_clinical_prediction_model_rf = function(simulated_data) {
  # # Fit random forest
  # clinical_model <- randomForest::randomForest(
  #   formula = clinical ~ surrogate + X1 + X2 + X3,
  #   data = simulated_data,
  #   ntree = 500,
  #   mtry = 2,
  #   nodesize = 3
  # )
  
  # Fit random forest with function from ranger package.
  clinical_model <- ranger::ranger(
    formula = clinical ~ surrogate + X1 + X2 + X3,
    data = simulated_data,
    num.trees = 100,
    m.try = 1,
    max.depth = 6,
    replace = FALSE,
    save.memory = TRUE,
    num.threads = 1,
    probability = TRUE
  )
  
  return(function(newdata) {
    return(predict(clinical_model, num.threads = 1, data = newdata, type = "response")$predictions[, 2])
  })
}


# Function that is just a wrapper to select one of the above functions.
train_clinical_prediction_model = function(simulated_data, method) {
  # Select the correct method-specific function.
  switch(
    method,
    "surrogate" = train_clinical_prediction_model_surrogate(simulated_data),
    "lm" = train_clinical_prediction_model_lm(simulated_data),
    "logistic" = train_clinical_prediction_model_logistic(simulated_data),
    "gam" = train_clinical_prediction_model_gam(simulated_data),
    "hal" = train_clinical_prediction_model_hal(simulated_data),
    "rf" = train_clinical_prediction_model_rf(simulated_data)
  )
}




