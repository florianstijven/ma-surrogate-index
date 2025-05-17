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
  if (type != "double") {
    estimate = statistic(data, weights = rep(1, nrow(data)))
    original_se = estimate$se
    estimate = estimate$estimate
    
    bootstrap_replications_list = multiplier_bootstrap_sampling(data, statistic, B)
    bootstrap_estimates = bootstrap_replications_list$bootstrap_estimates
    # If there are any NA or NaNs in the boostrap estimates, we assume that they
    # are +Inf or -Inf. 
    bootstrap_estimates[is.na(bootstrap_estimates) |
                          is.nan(bootstrap_estimates)] = ifelse(estimate > 0, +Inf, -Inf)
  }
  # Compute the required type of bootstrap CI.
  if (type == "BCa") {
    # Compute BCa interval
    return(
      BCa_CI(
        boot_replicates = bootstrap_estimates,
        estimate = statistic(data, weights = rep(1, nrow(data)))$estimate,
        statistic = statistic,
        data = data,
        alpha = alpha
      )
    )
  } else if (type == "percentile") {
    # Compute percentile interval
    return(percentile_CI(bootstrap_estimates, alpha))
  } else if (type == "BC percentile") {
    estimate = statistic(data, weights = rep(1, nrow(data)))$estimate
    return(BC_percentile_CI(estimate, bootstrap_estimates, alpha))
  } else if (type == "studentized") {
    bootstrap_ses = bootstrap_replications_list$bootstrap_ses
    bootstrap_ses[is.infinite(bootstrap_estimates)] = 1
    return(
      studentized_CI(
        bootstrap_estimates,
        bootstrap_ses,
        estimate,
        original_se,
        alpha
      )
    )
  } else if (type == "double") {
    double_bootstrap_CI(data = data, statistic = statistic, alpha = alpha, B = B)
  } else if (type == "basic") {
    estimate = statistic(data, weights = rep(1, nrow(data)))$estimate
    basic_CI(estimate = estimate, boot_replicates = bootstrap_estimates, alpha = alpha)
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
      "Some bootstrap replicates have missing values. These are removed for computing percentile confidence intervals.",
      "\nNumber of missing values in bootstrap replicates: ",
      sum(is.na(boot_replicates) | is.nan(boot_replicates)))
    warning(warning_message)
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

basic_CI = function(estimate, boot_replicates, alpha = 0.05) {
  # Check for missing values or NaN in the bootstrap replicates. A warning is
  # raised if there are any missing values. The warning also details the number
  # of problematic values such that the user can decide whether to ignore this.
  if (any(is.na(boot_replicates)) | any(is.nan(boot_replicates))) {
    warning_message = paste(
      "Some bootstrap replicates have missing values. These are removed for computing percentile confidence intervals.",
      "\nNumber of missing values in bootstrap replicates: ",
      sum(is.na(boot_replicates) | is.nan(boot_replicates)))
    warning(warning_message)
    # Remove missing values from the bootstrap replicates and standard errors.
    boot_replicates = boot_replicates[!is.na(boot_replicates) & !is.nan(boot_replicates)] 
  }
  ci_lower = quantile(x = boot_replicates - estimate, probs = alpha / 2, na.rm = TRUE)
  ci_upper = quantile(x = boot_replicates - estimate, probs = 1 - (alpha / 2), na.rm = TRUE)
  
  return(list(
    ci_lower = estimate - ci_upper,
    ci_upper = estimate - ci_lower
  ))
}

# Compute the Bias-Corrected percentile confidence interval
BC_percentile_CI = function(estimate, boot_replicates, alpha = 0.05) {
  # Check for missing values or NaN in the bootstrap replicates. A warning is
  # raised if there are any missing values. The warning also details the number
  # of problematic values such that the user can decide whether to ignore this.
  if (any(is.na(boot_replicates)) | any(is.nan(boot_replicates))) {
    warning_message = paste(
      "Some bootstrap replicates have missing values. These are removed for computing percentile confidence intervals.",
      "\nNumber of missing values in bootstrap replicates: ",
      sum(is.na(boot_replicates) | is.nan(boot_replicates)))
    warning(warning_message)
    # Remove missing values from the bootstrap replicates and standard errors.
    boot_replicates = boot_replicates[!is.na(boot_replicates) & !is.nan(boot_replicates)] 
  }
  # Compute bias-correction value.
  p0 = mean(boot_replicates < estimate) + 0.5 * mean(boot_replicates == estimate)
  z0 = qnorm(p0)
  
  alpha_lower = pnorm(2 * z0 + qnorm(alpha / 2))
  alpha_upper = pnorm(2 * z0 + qnorm(1 - alpha / 2))
  
  ci_lower = quantile(x = boot_replicates, probs = alpha_lower, na.rm = TRUE)
  ci_upper = quantile(x = boot_replicates, probs = alpha_upper, na.rm = TRUE)
  
  return(list(ci_lower = ci_lower, ci_upper = ci_upper))
}



# Function that implements the BCa bootstrap.
BCa_CI <- function(boot_replicates,
                   estimate,
                   data,
                   statistic,
                   alpha = 0.05) {
  # Check whether estimate is valid. If it is not valid, we return NA is
  # confidence intervals.
  if (is.na(estimate) | is.nan(estimate)) {
    warning("`estimate` is NA or NaN. BCa interval cannot be computed.")
    return(list(
      ci_lower = NA,
      ci_upper = NA
    ))
  }
  
  # Check for missing values or NaN in the bootstrap replicates. A warning is
  # raised if there are any missing values. The warning also details the number
  # of problematic values such that the user can decide whether to ignore this.
  if (any(is.na(boot_replicates)) | any(is.nan(boot_replicates))) {
    warning_message = paste(
      "Some bootstrap replicates have missing values. These are removed for computing BCa confidence intervals.",
      "\nNumber of missing values in bootstrap replicates: ",
      sum(is.na(boot_replicates) | is.nan(boot_replicates)))
    warning(warning_message)
    # Remove missing values from the bootstrap replicates and standard errors.
    boot_replicates = boot_replicates[!is.na(boot_replicates) & !is.nan(boot_replicates)] 
  }
  
  
  # Calculate the number of bootstrap replicates
  B <- length(boot_replicates)
  # Check length of bootstrap replicates.
  if (B <= 2) {
    stop("BCa could not be computed because there are too few valid bootstrap replicates.")
  }
  # Number of independent replicates in the data
  N <- nrow(data)
  
  
  # Calculate the z0 (bias correction) term
  p0 = mean(boot_replicates < estimate) + 0.5 * mean(boot_replicates == estimate)
  
  if (p0 == 0 | p0 == 1) {
    warning("Extreme median bias. Interpret results with care and increase number of bootstrap replications.")
    if (p0 == 0) {
      p0 = 1 / B
    } else {
      p0 = (B - 1) / B
    }
  }
  z0 <- qnorm(p0)
  
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

studentized_CI = function(boot_replicates,
                          se,
                          estimate,
                          original_se,
                          alpha = 0.05) {
  
  # Check whether estimate is valid. If it is not valid, we return NA is
  # confidence intervals.
  if (is.na(estimate) | is.nan(estimate)) {
    warning("`estimate` is NA or NaN. BCa interval cannot be computed.")
    return(list(
      ci_lower = NA,
      ci_upper = NA
    ))
  }
  
  if (any(is.na(boot_replicates)) | any(is.nan(boot_replicates)) | any(is.na(se)) | any(is.nan(se))) {
    warning_message = paste(
      "Some bootstrap replicates have missing values. These are removed for computing studentized confidence intervals.",
      "\nNumber of missing values in bootstrap replicates: ",
      sum(is.na(boot_replicates) | is.nan(boot_replicates) | is.na(se) | is.nan(se)))
    warning(warning_message)
    # Remove missing values from the bootstrap replicates and standard errors.
    good_indices = !is.na(boot_replicates) & !is.nan(boot_replicates) & !is.na(se) & !is.nan(se)
    boot_replicates = boot_replicates[good_indices] 
    se = se[good_indices] 
  }
  
  # Check length of bootstrap replicates.
  if (B <= 2) {
    stop("BCa could not be computed because there are too few valid bootstrap replicates.")
  }
  
  # Compute studentized statistics.
  t = (boot_replicates - estimate) / se
  
  # Determine percentiles of the studentized distribution
  alpha <- 0.05
  lower_percentile <- quantile(t, probs = alpha / 2)
  upper_percentile <- quantile(t, probs = 1 - alpha / 2)
  
  # Step 5: Transform back to original scale
  ci_lower <- estimate - upper_percentile * original_se
  ci_upper <- estimate - lower_percentile * original_se
  
  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}


double_bootstrap_CI <- function(data, statistic, alpha = 0.05, B = 1000, M = 500) {
  n <- nrow(data)
  theta_hat <- statistic(data, weights = rep(1, n))$estimate
  
  # Step 1: Outer bootstrap
  theta_star_star = matrix(0, nrow = B, ncol = M)
  theta_star = rep(0, B)
  for (b in 1:B) {  
    # First-level bootstrap sample
    weights_star <- rexp(n = n, rate = 1)
    theta_star[b] <- statistic(data, weights_star)$estimate
    
    # Step 2: Inner bootstrap
    theta_star_star[b, ] <- replicate(M, {
      weights_star_star <- weights_star * rexp(n = n, rate = 1)
      statistic(data, weights_star_star)$estimate
    })
    
  }
  
  # Step 4: Adjust alpha by searching for alpha' such that coverage ≈ 1 - alpha
  adjust_alpha <- function(alpha_try, alpha_goal) {
    # Compute the confidence intervals based on the level 2 bootstrap.
    ci_limit = apply(
      X = theta_star_star,
      FUN = quantile,
      probs = alpha_try,
      MARGIN = 1, 
      na.rm = TRUE
    )
    
    coverage =  as.numeric(theta_hat <= ci_limit)
    return(mean(coverage - alpha_goal, na.rm = TRUE))
  }
  
  
  # Use bisection method to find alpha' where adjusted coverage matches alpha /2
  # to obtain the adjusted upper limit alpha.
  if (adjust_alpha(0.35, 1 - alpha / 2) * adjust_alpha(1, 1 - alpha / 2) > 0) {
    alpha_upper_adjusted = 1 - alpha / 2
  } else {
    alpha_upper_adjusted <- uniroot(
      adjust_alpha,
      lower = 0.35,
      upper = 1,
      alpha_goal = 1 - alpha / 2,
      tol = 1 / B
    )$root
  }

  # Repeat the above for the lower limit with alpha / 2.
  if (adjust_alpha(0, alpha / 2) * adjust_alpha(0.65, alpha / 2) > 0) {
    alpha_lower_adjusted = alpha / 2
  } else {
    alpha_lower_adjusted <- uniroot(
      adjust_alpha,
      lower = 0,
      upper = 0.65,
      alpha_goal = alpha / 2,
      tol = 1 / B
    )$root
  }
  
  ci_lower <- quantile(theta_star, alpha_lower_adjusted, na.rm = TRUE)
  ci_upper <- quantile(theta_star, alpha_upper_adjusted, na.rm = TRUE)
  
  c(lower = ci_lower, upper = ci_upper)
}


