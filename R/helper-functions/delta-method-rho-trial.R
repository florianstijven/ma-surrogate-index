# Function that computes the a confidence interval for the Pearson correlation
# parameter based on the estimated variances and covariance, the corresponding
# variance matrix. The confidence interval is based on the delta method and
# Fisher's Z transformation.
rho_delta_method = function(coefs, vcov, alpha = 0.05, method = "t-adjustment", N = NULL){
  # Extract the estimated parameters.
  mu_alpha = coefs[1]
  mu_beta = coefs[2]
  d_alpha = coefs[3]
  d_beta = coefs[4]
  d_alpha_beta = coefs[5]
  
  # Estimate the correlation
  rho = d_alpha_beta / sqrt(d_alpha * d_beta)
  
  # Compute the standard error for rho_hat (without any transformation) based on
  # the elements of the variance matrix and the delta method.

  # Vector of partial derivative of rho wrt d_alpha, d_beta, and d_alpha_beta.
  partial_rho = c(
    (-1 / 2) * d_alpha_beta / (sqrt(d_alpha * d_beta) * d_alpha),
    (-1 / 2) * d_alpha_beta / (sqrt(d_alpha * d_beta) * d_beta),
    1 / sqrt(d_alpha * d_beta)
  )
  # This vector of partial derivatives should be a column matrix.
  partial_rho = matrix(partial_rho, nrow = 3, ncol = 1)

  # Compute the SE of rho.
  se_rho_hat = sqrt(t(partial_rho) %*% vcov[3:5, 3:5] %*% partial_rho) %>%
    as.vector()
  
  # Compute the SE for the fisher's z transformation of rho.
  se_fisher_z = se_rho_hat / ((1 + rho) * (1 - rho))
  
  if (method == "no adjustment") {
    # Compute the Z-score for the confidence interval.
    z = qnorm(1 - alpha / 2)
  } else if (method == "t-adjustment") {
    if (is.null(N)) stop("Number of independent trials N required for t-adjustment.")
    z = qt(1 - alpha / 2, df = N - 1)
  }

  
  # Compute the confidence interval for the Pearson correlation parameter on the
  # Fisher's Z scale.
  ci_fisher_z = atanh(rho) + c(-1, 1) * z * se_fisher_z
  
  # Transform the limits back to the rho scale.
  ci_rho_hat = tanh(ci_fisher_z)
  
  return(list(rho = rho, se = se_rho_hat, ci = ci_rho_hat))
}
