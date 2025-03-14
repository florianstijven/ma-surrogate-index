# Function that generates bootstrap replications based on weights sampled from
# unit exponential distributions.
multiplier_bootstrap_sampling <- function(data, statistic, B) {
  # Number of rows in the data set that is to be bootstrapped
  n <- nrow(data)
  
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
  
  # Extract estimates from the list of estimates.
  bootstrap_estimates = purrr::map_dbl(.x = bootstrap_est_list, 1)
  
  # Extraxt SEs from the list of estimates.
  bootstrap_ses = purrr::map_dbl(.x = bootstrap_est_list, 2)
  
  # Return list with bootstrap replicates and standard errors.
  return(list(
    bootstrap_estimates = bootstrap_estimates,
    bootstrap_ses = bootstrap_ses
  ))
}



multiplier_bootstrap_ci = function(data, statistic, B, alpha = 0.05, type = "BCa") {
  bootstrap_replications_list = multiplier_bootstrap_sampling(data, statistic, B)
  # Compute the required type of bootstrap CI.
  if (type == "BCa") {
    # Compute BCa interval
    return(
      BCa_CI(
        boot_replicates = bootstrap_replications_list$bootstrap_estimates,
        estimate = statistic(data, weights = rep(1, nrow(data)))$estimate,
        statistic = statistic,
        data = data,
        alpha = alpha
      )
    )
  } else if (type == "percentile") {
    # Compute percentile interval
    return(percentile_CI(bootstrap_replications_list$bootstrap_estimates, alpha))
  } else {
    stop("Invalid type. Must be 'BCa' or 'percentile'.")
  }
}

# Compute the percentile confidence intervals
percentile_CI = function(boot_replicates, alpha = 0.05) {
  # Check for missing values or NaN in the bootstrap replicates. A warning is
  # raised if there are any missing values. The warning also details the number
  # of problematic values such that the user can decide whether to ignore this.
  if (any(is.na(boot_replicates)) | any(is.nan(boot_replicates))) {
    warning_message = paste(
      "Some bootstrap replicates have missing values. These are removed for computing BCa confidence intervals.",
      "\nNumber of missing values in bootstrap replicates: ",
      sum(is.na(boot_replicates) | is.nan(boot_replicates)))
    # Remove missing values from the bootstrap replicates and standard errors.
    boot_replicates = boot_replicates[!is.na(boot_replicates) & !is.nan(boot_replicates)] 
  }
  ci_lower = quantile(x = boot_replicates, probs = alpha / 2, na.rm = TRUE)
  ci_upper = quantile(x = boot_replicates, probs = 1 - (alpha / 2), na.rm = TRUE)
  
  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}



# Function that implements the BCa bootstrap.
BCa_CI <- function(boot_replicates,
                   estimate,
                   data,
                   statistic,
                   alpha = 0.05) {
  # Check for missing values or NaN in the bootstrap replicates. A warning is
  # raised if there are any missing values. The warning also details the number
  # of problematic values such that the user can decide whether to ignore this.
  if (any(is.na(boot_replicates)) | any(is.nan(boot_replicates))) {
    warning_message = paste(
      "Some bootstrap replicates have missing values. These are removed for computing BCa confidence intervals.",
      "\nNumber of missing values in bootstrap replicates: ",
      sum(is.na(boot_replicates) | is.nan(boot_replicates)))
    # Remove missing values from the bootstrap replicates and standard errors.
    boot_replicates = boot_replicates[!is.na(boot_replicates) & !is.nan(boot_replicates)] 
  }
  
  
  # Calculate the number of bootstrap replicates
  B <- length(boot_replicates)
  # Number of independent replicates in the data
  N <- nrow(data)
  
  
  # Calculate the z0 (bias correction) term
  median_bias = sum(boot_replicates < estimate) / B
  if (median_bias == 0 | median_bias == 1) {
    warning("Extreme median bias. Interpret results with care and increase number of bootstrap replications.")
    if (median_bias == 0) {
      median_bias = 1 / B
    } else {
      median_bias = (B - 1) / B
    }
  }
  z0 <- qnorm(sum(boot_replicates < estimate) / B)
  
  # Calculate the acceleration term (a)
  jackknife_est <- sapply(1:N, function(i) {
    statistic(data[-i, ], weights = rep(1, N - 1))$estimate
  })
  if (any(is.na(jackknife_est)) | any(is.nan(jackknife_est))) {
    warning("Some jacknknife estimates could not be computed. These are ignored in the computation of the BCa interval. The percentile CI may be advised.")
    jackknife_est = jackknife_est[!is.na(jackknife_est)]
    }
  jackknife_mean <- mean(jackknife_est)
  a <- sum((jackknife_mean - jackknife_est)^3) / (6 * sum((jackknife_mean - jackknife_est)^2)^1.5)
  
  # Calculate the adjusted alpha levels
  alpha1 <- pnorm(z0 + (z0 + qnorm(alpha / 2)) / (1 - a * (
    z0 + qnorm(alpha / 2)
  )))
  alpha2 <- pnorm(z0 + (z0 + qnorm(1 - alpha / 2)) / (1 - a * (
    z0 + qnorm(1 - alpha / 2)
  )))
  
  # Calculate the BCa intervals
  ci_lower <- quantile(boot_replicates, alpha1)
  ci_upper <- quantile(boot_replicates, alpha2)
  
  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}


