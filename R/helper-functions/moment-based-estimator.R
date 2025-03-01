# The function below implements the so-called moment based estimator for the
# mean and covariance parameters of the trial-level treatment effects in a
# meta-analytic sample. 
moment_estimator = function(
    alpha_hat, # Vector of treatment effect estimates on the surrogate.
    beta_hat, # Vector of treatment effect estimates on the clinical endpoint.
    vcov_list # List of variance-covariance matrices for the treatment effect estimates.
                            ) 
  {
  # Total number of independent units.
  N = length(alpha_hat)
  
  # Estimate the overall means.
  mu_alpha_hat = mean(alpha_hat)
  mu_beta_hat = mean(beta_hat)
  
  # Estimate the unadjusted covariance matrix.
  S = cov(cbind(alpha_hat, beta_hat))
  
  # Estimate the mean of the within-trial covariance matrices.
  # Compute the element-wise mean of a list of matrices.
  mean_est_error_vcov = Reduce("+", vcov_list) / N
  
  # Estimate the adjusted covariance matrix.
  S_adj = S - mean_est_error_vcov
  
  # Compute bread matrix for sandwich estimator. This is a diagonal matrix with
  # the first two diagonal elements ones, and the remaining three minus ones.
  bread = matrix(0, nrow = 5, ncol = 5)
  diag(bread) = c(1, 1, -1, -1, -1)
  
  # Compute the ham matrix for the sandwich estimator. This is the outer product
  # of the estimation equations evaluated at the estimated parameters.
  U = est_function(alpha_hat,
                   beta_hat,
                   vcov_list,
                   c(mu_alpha_hat, mu_beta_hat, S_adj[1, 1], S_adj[2, 2], S_adj[1, 2]))
  
  # Compute the outer product of estimating equations.
  ham = (t(U) %*% U) * (1 / N)
  
  # Compute the sandwich estimator.
  sandwich = bread %*% ham %*% bread
  
  # Return the estimated mean and covariance parameters and the corresponding 
  # sandwich estimator.
  return(list(
    coefs = c(mu_alpha_hat, mu_beta_hat, S_adj[1, 1], S_adj[2, 2], S_adj[1, 2]),
    vcov = sandwich / N
  ))
}


# Function that evaluates the estimating function in set of data given the parameters.
est_function = function(alpha_hat, beta_hat, vcov_list, params) {
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
  ee_d_alpha = (alpha_hat - mu_alpha) ^ 2 - sigma_alpha - d_alpha
  ee_d_beta = (beta_hat - mu_beta) ^ 2 - sigma_beta - d_beta
  ee_d_alpha_beta = (alpha_hat - mu_alpha) * (beta_hat - mu_beta) - sigma_alpha_beta - d_alpha_beta
  
  # Return the evaluated estimating functions as a matric with each unit's
  # values as a single row.
  return(
    matrix(
      c(ee_mu_alpha, ee_mu_beta, ee_d_alpha, ee_d_beta, ee_d_alpha_beta),
      nrow = length(alpha_hat),
      ncol = 5,
      byrow = FALSE
    )
  )
}
