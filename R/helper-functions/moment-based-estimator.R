# Function to estimate the residual variance under the assumption of an identity
# relation between treatment effects on the clinical endpoint and treatment
# effects on the surrogate endpoint.
residual_variance_identity_line = function(alpha_hat, # Vector of treatment effect estimates on the surrogate.
                                           beta_hat, # Vector of treatment effect estimates on the clinical endpoint.
                                           vcov_list, # List of variance-covariance matrices for the treatment effect estimates.
                                           weights = rep(1, length(alpha_hat))) {
  # Total number of independent units.
  N = length(alpha_hat)
  
  # Make sure that the weights sum to 1.
  weights = weights / sum(weights)
  
  # Compute the variance of the estimated treatment effect differences. 
  total_variance = sum(weights * (beta_hat - alpha_hat) ^ 2)
  
  # Compute mean of the estimated variances due to sampling variability.
  var_alpha = purrr::map_dbl(vcov_list, ~.x[1, 1])
  var_beta = purrr::map_dbl(vcov_list, ~.x[2, 2])
  var_alpha_beta = purrr::map_dbl(vcov_list,~.x[1, 2])
  mean_est_var = sum(weights * (var_alpha + var_beta - 2 * var_alpha_beta))
  
  # Compute corrected residual variance.
  residual_variance = total_variance - mean_est_var
  
  return(residual_variance)
}


# The function below implements the so-called moment based estimator for the
# mean and covariance parameters of the trial-level treatment effects in a
# meta-analytic sample. 
moment_estimator = function(
    alpha_hat, # Vector of treatment effect estimates on the surrogate.
    beta_hat, # Vector of treatment effect estimates on the clinical endpoint.
    vcov_list, # List of variance-covariance matrices for the treatment effect estimates.
    estimator_adjustment = "N - 1", # Finite-sample adjustment to the estimator
    sandwich_adjustment = "N -1", # Finite-sample adjustment to the sandwich estimator
    weights = rep(1, length(alpha_hat)),
    SE = TRUE,
    nearest_PD = FALSE) 
  {
  # Total number of independent units.
  N = length(alpha_hat)
  
  # Make sure that the weights sum to 1.
  weights = weights / sum(weights)
  
  # Estimate the overall means.
  mu_alpha_hat = sum(weights * alpha_hat)
  mu_beta_hat = sum(weights * beta_hat)
  
  # Estimate the unadjusted covariance matrix.
  total_residual = cbind(alpha_hat - mu_alpha_hat, beta_hat - mu_beta_hat) %>%
    as.matrix()
  S =  t(weights * total_residual) %*% total_residual
  
  # If required, do finite-sample adjustment to the estimated covariance matrix.
  if (estimator_adjustment == "N - 1") {
    # S = (N / (N - 1)) * S
    S = (1 / (1 - sum(weights ** 2))) * S
  }

  # Estimate the mean of the within-trial covariance matrices.
  # Compute the element-wise mean of a list of matrices.
  weighted_vcov_list = purrr::map2(
    .x = vcov_list,
    .y = weights,
    .f = function(vcov_matrix, weight) weight * vcov_matrix
  )
  mean_est_error_vcov = Reduce("+", weighted_vcov_list)
  
  # Estimate the adjusted covariance matrix.
  S_adj = S - mean_est_error_vcov
  
  
  # Compute sandwich estimator if required.
  if (SE) {
    # Compute bread matrix for sandwich estimator. This is a diagonal matrix with
    # the first two diagonal elements ones, and the remaining three minus ones.
    bread = matrix(0, nrow = 5, ncol = 5)
    diag(bread) = c(1, 1, -1, -1, -1)
    
    # Compute the ham matrix for the sandwich estimator. This is the outer product
    # of the estimation equations evaluated at the estimated parameters.
    U = est_function(alpha_hat = alpha_hat,
                     beta_hat = beta_hat,
                     vcov_list = vcov_list,
                     params = c(mu_alpha_hat, mu_beta_hat, S_adj[1, 1], S_adj[2, 2], S_adj[1, 2]),
                     estimator_adjustment = estimator_adjustment,
                     weights = weights)
    
    # Compute the outer product of estimating equations.
    ham = t(weights * U) %*% U
    
    # Compute the sandwich estimator.
    sandwich = bread %*% ham %*% bread
    
    # Correct the sandwich estimate if required.
    if (sandwich_adjustment == "N - 1") {
      sandwich = (1 / (1 - sum(weights ** 2))) * sandwich
    }
  } else {
    sandwich = NA
  }
  
  residual_var = residual_variance_identity_line(
    alpha_hat = alpha_hat,
    beta_hat = beta_hat,
    vcov_list = vcov_list,
    weights = weights
  )
  
  # Find the nearest positive definite matrix, if this is asked by the user.
  if (nearest_PD) {
    S_adj = tryCatch({
      Matrix::nearPD(S_adj)$mat
    }, error = function(e) {
      S_adj = matrix(NA, nrow = 2, ncol = 2)
    })
    
  }
  

  # Return the estimated mean and covariance parameters and the corresponding 
  # sandwich estimator.
  return(list(
    coefs = c(mu_alpha_hat, mu_beta_hat, S_adj[1, 1], S_adj[2, 2], S_adj[1, 2]),
    vcov = sandwich / N,
    residual_var = residual_var
  ))
}


# Function that evaluates the estimating function in set of data given the parameters.
est_function = function(alpha_hat, beta_hat, vcov_list, params, estimator_adjustment, weights) {
  # Total number of independent units.
  N = length(alpha_hat)
  
  # Make sure that the weights sum to 1.
  weights = weights / sum(weights)
  
  # Extract the parameters from params.
  mu_alpha = params[1]
  mu_beta = params[2]
  d_alpha = params[3]
  d_beta = params[4]
  d_alpha_beta = params[5]
  
  # Extract within-trial covariance estimates.
  sigma_alpha = purrr::map_dbl(vcov_list, ~.x[1, 1])
  sigma_beta = purrr::map_dbl(vcov_list, ~.x[2, 2])
  sigma_alpha_beta = purrr::map_dbl(vcov_list, ~.x[1, 2])
  
  # Estimating function evaluated in the data and parameters for the mean
  # parameters.
  ee_mu_alpha = alpha_hat - mu_alpha
  ee_mu_beta = beta_hat - mu_beta
  
  # Estimating function evaluated in the data and parameters for the covariance
  # parameters.
  if (estimator_adjustment == "N - 1") {
    # Finite sample adjustment taking weights into account.
    finite_sample_adj = (1 / (1 - sum(weights ** 2))) 
    
    ee_d_alpha = finite_sample_adj * (alpha_hat - mu_alpha) ^ 2 - sigma_alpha - d_alpha
    ee_d_beta = finite_sample_adj * (beta_hat - mu_beta) ^ 2 - sigma_beta - d_beta
    ee_d_alpha_beta = finite_sample_adj * (alpha_hat - mu_alpha) * (beta_hat - mu_beta) - sigma_alpha_beta - d_alpha_beta
  } else if (estimator_adjustment == "none") {
    ee_d_alpha = (alpha_hat - mu_alpha) ^ 2 - sigma_alpha - d_alpha
    ee_d_beta = (beta_hat - mu_beta) ^ 2 - sigma_beta - d_beta
    ee_d_alpha_beta = (alpha_hat - mu_alpha) * (beta_hat - mu_beta) - sigma_alpha_beta - d_alpha_beta
  }

  
  # Return the evaluated estimating functions as a matrix with each unit's
  # values as a single row.
  U_unweighted = matrix(
    c(ee_mu_alpha, ee_mu_beta, ee_d_alpha, ee_d_beta, ee_d_alpha_beta),
    nrow = length(alpha_hat),
    ncol = 5,
    byrow = FALSE
  )
  U = U_unweighted
  return(U)
}
