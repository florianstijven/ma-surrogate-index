# Function that generates bootstrap replications based on weights sampled from
# unit exponential distributions.
multiplier_bootstrap_sampling <- function(data, statistic, B) {
  # Number of rows in the data set that is to be bootstrapped
  n <- nrow(data)
  statistic(data, weights = rep(1, n))
  
  # Function that samples weights and evaluates statistic in reweighted data.
  bootstrap_replicate_f = function(x) {
    # Generate random weights from the unit exponential distribution
    weights <- rexp(n = n, rate = 1)
    
    # Compute statistic on the weighted data set
    estimate = statistic(data, weights)
    return(estimate)
  }
  # Compute estimate on reweighted data. A list of estimates is returned.
  bootstrap_est_list = purrr::map(.x = 1:B, .f = bootstrap_replicate_f)
  
  # Convert the list of estimates to a matrix with one row per bootstrap
  # replicate.
  bootstrap_matrix = as.matrix(do.call(rbind, bootstrap_est_list))
  return(bootstrap_matrix)
}

# Compute the percentile confidence intervals
multiplier_bootstrap_ci_help = function(samples_matrix, alpha = 0.05) {
  if (any(is.na(samples_matrix))) {
    warning("Some bootstrap replicates have missing values. These are removed for computing percentile confidence intervals.")
  }
  ci_lower = apply(samples_matrix, 2, quantile, probs = alpha / 2, na.rm = TRUE)
  ci_upper = apply(samples_matrix, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  
  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}

multiplier_bootstrap_ci = function(data, statistic, B, alpha = 0.05 ) {
  samples_matrix = multiplier_bootstrap_sampling(data, statistic, B)
  multiplier_bootstrap_ci_help(samples_matrix, alpha)
}


