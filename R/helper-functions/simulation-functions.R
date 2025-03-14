# Function to generate trial-level random parameters. In the current setup, we
# only have three trial-level random parameters: p1, beta_surr_treatment, and
# beta_clin_treatment. p1 is the proportion for the baseline binary covariate.
# beta_surr_treatment is the mean treatment effect on the surrogate endpoint.
# beta_clin_treatment is the "additional" treatment effect on the clinical
# endpoint in a model conditional on the baseline covariate and surrogate.
generate_random_coefficients_proof_of_concept <- function(sd_beta_clin_treatment = 0.10, sd_beta_clin_surrogate = 0.10) {
  # Sample trial-level proportion parameter for baseline covariate from uniform
  # distribution. 
  p1 <- runif(1, min = 0.3, max = 0.7)
  
  # Sample random treatment effect on surrogate. 
  beta_surr_treatment <- rnorm(1, mean = 0.5, sd = 0.3)
  
  # Sample random treatment effect on clinical endpoint.
  beta_clin_treatment <- rnorm(1, mean = 0, sd = sd_beta_clin_treatment)
  beta_clin_surrogate <- rnorm(1, mean = 0, sd = sd_beta_clin_surrogate)
  
  # Return the list of sampled parameters.
  return(
    list(
      p1 = p1,
      beta_surr_treatment = beta_surr_treatment,
      beta_clin_treatment = beta_clin_treatment,
      beta_clin_surrogate = beta_clin_surrogate
    )
  )
}


generate_random_coefficients_vaccine <- function(sd_beta_clin_treatment = 0.10, sd_beta_clin_surrogate = 0.10) {
  # Sample trial-level proportion parameter for X3 from uniform
  # distribution. 
  p1 <- runif(1, min = 0.25, max = 0.75)
  # Trial-level mean parameters for X1 and X2 are also sampled from uniform
  # distributions.
  mu1 = runif(1, min = -1, max = 1)
  mu2 = runif(1, min = -1, max = 1)
  
  # Sample random treatment effect on surrogate. 
  beta_surr_treatment <- runif(1, min = 0.25, max = 3)
  
  # Sample random treatment effect on clinical endpoint.
  beta_clin_treatment <- rnorm(1, mean = 0, sd = sd_beta_clin_treatment)
  beta_clin_surrogate <- rnorm(1, mean = 0, sd = sd_beta_clin_surrogate)
  
  # Return the list of sampled parameters.
  return(
    list(
      mu1 = mu1,
      mu2 = mu2,
      p1 = p1,
      beta_surr_treatment = beta_surr_treatment,
      beta_clin_treatment = beta_clin_treatment,
      beta_clin_surrogate = beta_clin_surrogate
    )
  )
}


# Function to simulate a single trial with passed trial-level parameters.
simulate_trial_proof_of_concept <- function(n, coefficients) {
  # Randomization: Treatment assignment (1 or 0). Simple randomization is
  # assumed.
  treatment <- rbinom(n, 1, 0.5)  # 1 = treatment, 0 = control
  
  # Baseline covariates: One covariate from a Bernoulli distribution with
  # probability p1.
  X1 <- rbinom(n, 1, coefficients$p1)
  
  # Surrogate endpoint: Linear relationship with treatment and random noise. The
  # distribution of the surrogate does not depend on the baseline covariate.
  surrogate <- coefficients$beta_surr_treatment * treatment + rnorm(n)
  
  # Simulate clinical endpoint with no random coefficients. We're assuming that
  # the regression of the clinical endpoint on the surrogate and covariate is
  # the same across trials modulus some small random treatment effect.
  clinical <- -1 * surrogate * X1 + surrogate * (1 - X1) + 
    coefficients$beta_clin_treatment * treatment +
    coefficients$beta_clin_surrogate * surrogate ^ 2 +
    rnorm(n)
  
  # Data frame for a single trial
  trial_data <- data.frame(treatment, X1, surrogate, clinical)
  
  return(trial_data)
}

# Function to simulate a single vaccine trial. 
simulate_trial_vaccine = function(n, coefficients) {
  # Randomization: Treatment assignment (1 or 0). Simple randomization is
  # assumed.
  treatment <- rbinom(n, 1, 0.5)  # 1 = treatment, 0 = control
  
  # Baseline covariates: One covariate from a Bernoulli distribution with
  # probability p1.
  X1 = rnorm(n, mean = coefficients$mu1)
  X2 = rnorm(n, mean = coefficients$mu2)
  X3 <- rbinom(n, 1, coefficients$p1)
  
  # Surrogate endpoint: Linear relationship with treatment and random noise. The
  # distribution of the surrogate does not depend on the baseline covariates.
  surrogate <- coefficients$beta_surr_treatment * treatment + rnorm(n)
  
  # Simulate clinical endpoint with no random coefficients. 
  eta = coefficients$beta_clin_treatment * treatment - 
    (1 + coefficients$beta_clin_surrogate + 1 * X3) * surrogate +
    X1 - 2 * X3
  prob_clin = (1 / (1 + exp(-1 * (eta)))) * 0.10 * exp(sin(1.5 * X2) + 1)
  clinical <- rbinom(n = n, size = 1, prob = prob_clin)
  
  # Data frame for a single trial
  trial_data <- data.frame(treatment, X1, X2, X3, surrogate, clinical)
  
  return(trial_data)
}

# Function to simulate data for multiple trials with a common within-trial
# sample size. The random trial-level parameters are automatically sampled
# within this function.
simulate_trials_with_random_coefficients <- function(N, n, sd_beta, scenario) {
  
  if (scenario == "proof-of-concept") {
    simulate_trial = simulate_trial_proof_of_concept
    generate_random_coefficients = generate_random_coefficients_proof_of_concept
  } else if (scenario == "vaccine") {
    simulate_trial = simulate_trial_vaccine
    generate_random_coefficients = generate_random_coefficients_vaccine
  } else {
    stop("`scenario` is not valid.")
  }
  
  
  trials <- lapply(1:N, function(i,
                                 sd_beta_clin_treatment,
                                 sd_beta_clin_surrogate) {
    coefficients <- generate_random_coefficients(sd_beta_clin_treatment, sd_beta_clin_surrogate)  # generate new random coefficients for each trial
    trial_data <- simulate_trial(n = n, coefficients = coefficients)
    trial_data$trial <- i  # Add a trial-level variable
    return(trial_data)
  }, sd_beta_clin_treatment = sd_beta[1], sd_beta_clin_surrogate = sd_beta[2])
  
  # Combine all trial data into a single tibble
  all_trials_data <- bind_rows(trials) %>%
    as_tibble()  # Ensure the result is a tibble
  
  return(all_trials_data)
}

# Function to compute treatment effects with robust standard errors and covariance matrix
compute_treatment_effects <- function(trial_data, method = "adjusted", formula = NULL, measure = rep("mean difference", "mean difference")) {
  # If adjusted treatment effects are requires, we fit two univariate linear
  # regression models and estimate the covariance of the treatment effect
  # estimators through the sandwich estimators.
  if (method == "adjusted") {
    # Adjusted method cannot be combined with the VE measure.
    if (any(measure != c("mean difference", "mean difference"))) {
      stop("The adjusted method can only be combined with the mean difference measure.")
    }
    
    # Combine the coefficients into a single model to compute the covariance matrix for both treatment effects
    combined_model <- lm(formula, data = trial_data)
    
    # Compute the covariance matrix using the sandwich estimator
    cov_matrix <- sandwich::vcovHC(combined_model, type = "HC0")
    
    # Extract the covariance matrix between the treatment effects on surrogate and clinical
    cov_surr_clin <- cov_matrix[c("surrogate:treatment", "clinical:treatment"), c("surrogate:treatment", "clinical:treatment")]  # Covariance between treatment effects
    
    # Return a tibble with the treatment effects, standard errors, and covariance matrix
    treatment_effects <- tibble(
      treatment_effect_surr = coef(combined_model)["treatment", "surrogate"],  # Treatment effect on surrogate
      se_surr = sqrt(cov_surr_clin["surrogate:treatment", "surrogate:treatment"]),  # Robust standard error for surrogate treatment effect
      treatment_effect_clin = coef(combined_model)["treatment", "clinical"],  # Treatment effect on clinical endpoint
      se_clin = sqrt(cov_surr_clin["clinical:treatment", "clinical:treatment"]),  # Robust standard error for clinical treatment effect
      covariance_matrix = list(cov_surr_clin)  # Store covariance matrix as a list
    )
  } else if (method == "mean") {
    # Compute treatment and outcome specific means.
    surr_mean_0 = mean(trial_data$surrogate[trial_data$treatment == 0])
    surr_mean_1 = mean(trial_data$surrogate[trial_data$treatment == 1])
    clin_mean_0 = mean(trial_data$clinical[trial_data$treatment == 0])
    clin_mean_1 = mean(trial_data$clinical[trial_data$treatment == 1])
    
    # Sample covariance matrices in each treatment group.
    vcov_0 = var(trial_data[trial_data$treatment == 0, c("surrogate", "clinical")])
    vcov_1 = var(trial_data[trial_data$treatment == 1, c("surrogate", "clinical")])
    
    # Treatment group sample sizes
    n0 = sum(trial_data$treatment == 0)
    n1 = sum(trial_data$treatment == 1)
    
    # Compute covariance matrix of (surr_mean_0, clin_mean_0, surr_mean_1, clin_mean_1)'.
    covariance_matrix_means = Matrix::bdiag(list(vcov_0 / n0, vcov_1 / n1)) %>%
      as.matrix()
    
    # Estimate treatment effect on the surrogate and compute corresponding
    # gradient (for subsequent use in the delta method).
    if (measure[1] == "mean difference") {
      treatment_effect_surr = surr_mean_1 - surr_mean_0
      
      # Vector of partial derivatives of surr_mean_1 - surr_mean_0 with
      # respect to (surr_mean_0, clin_mean_0, surr_mean_1, clin_mean_1)'.
      partial_s = c(1, 0, 1, 0)
    } else if (measure[1] == "VE") {
      treatment_effect_surr = 1 - (surr_mean_1 / surr_mean_0)
      
      # Vector of partial derivatives of 1 - (surr_mean_1 / surr_mean_0) with
      # respect to (surr_mean_0, clin_mean_0, surr_mean_1, clin_mean_1)'.
      partial_s = c(surr_mean_1 * surr_mean_0 ^ (-2), 0, -1 / surr_mean_0, 0)
    } else if (measure[1] == "log RR") {
      treatment_effect_surr = log(surr_mean_1 / surr_mean_0)
      
      # Vector of partial derivatives of log(surr_mean_1 / surr_mean_0) with
      # respect to (surr_mean_0, clin_mean_0, surr_mean_1, clin_mean_1)'.
      partial_s = c(-1 / surr_mean_0, 0, 1 / surr_mean_1, 0)
    } else if (measure[1] == "exp mean difference") {
      treatment_effect_surr = exp(surr_mean_1) - exp(surr_mean_0)
      
      # Vector of partial derivatives of exp(surr_mean_1) - exp(surr_mean_0) with
      # respect to (surr_mean_0, clin_mean_0, surr_mean_1, clin_mean_1)'.
      partial_s = c(- exp(surr_mean_0), 0, exp(surr_mean_1), 0)
    } else {
      stop("`measure` is not valid.")
    }
    
    # Estimate treatment effect on the clinical endpoint and compute corresponding
    # gradient (for subsequent use in the delta method).
    if (measure[2] == "mean difference") {
      treatment_effect_clin = clin_mean_1 - clin_mean_0
      
      # Vector of partial derivatives of clin_mean_1 - clin_mean_0 with
      # respect to (surr_mean_0, clin_mean_0, surr_mean_1, clin_mean_1)'.
      partial_c = c(0, 1, 0, 1)
    } else if (measure[2] == "VE") {
      treatment_effect_clin = 1 - (clin_mean_1 / clin_mean_0)
      
      # Vector of partial derivatives of 1 - (clin_mean_1 / clin_mean_0) with
      # respect to (surr_mean_0, clin_mean_0, surr_mean_1, clin_mean_1)'.
      partial_c = c(0, clin_mean_1 * clin_mean_0 ^ (-2), 0, -1 / clin_mean_0)
    } else if (measure[2] == "log RR") {
      treatment_effect_clin = log(clin_mean_1 / clin_mean_0)
      
      # Vector of partial derivatives of log(clin_mean_1 / clin_mean_0) with
      # respect to (surr_mean_0, clin_mean_0, surr_mean_1, clin_mean_1)'.
      partial_c = c(0, -1 / clin_mean_0, 0, 1 / clin_mean_1)
    } else if (measure[2] == "exp mean difference") {
      treatment_effect_clin = exp(clin_mean_1) - exp(clin_mean_0)
      
      # Vector of partial derivatives of exp(surr_mean_1) - exp(surr_mean_0) with
      # respect to (surr_mean_0, clin_mean_0, surr_mean_1, clin_mean_1)'.
      partial_c = c(0, - exp(clin_mean_0), 0, exp(clin_mean_1))
    } else {
      stop("`measure` is not valid.")
    }

    # Jacobian matrix with the gradients on the columns.
    J = matrix(c(partial_s, partial_c), ncol = 2, byrow = FALSE)
    
    # Compute covariance matrix of the treatment effects based on the delta method.
    covariance_matrix = t(J) %*% covariance_matrix_means %*% J
    
    # Return a tibble with the treatment effects, standard errors, and covariance matrix
    treatment_effects <- tibble(
      treatment_effect_surr = treatment_effect_surr,  # Treatment effect on surrogate
      se_surr = sqrt(covariance_matrix[1, 1]),  # Standard error for surrogate treatment effect
      treatment_effect_clin = treatment_effect_clin,  # Treatment effect on clinical endpoint
      se_clin = sqrt(covariance_matrix[2, 2]),  # Standard error for clinical treatment effect
      covariance_matrix = list(covariance_matrix)  # Store covariance matrix as a list
    )
  }
  
  return(treatment_effects)
}

# Function that generates a meta-analytic data set, adds the estimated surrogate
# index as a variable, and computes the treatment effects on (i) the surrogate
# and clinical endpoints and (ii) the surrogate index and clinical endpoint.
# This function uses the previously defined functions.
generate_meta_analytic_data <- function(N, n, surrogate_index_estimators, sd_beta, scenario) {
  # Simulate data for N trials with random coefficients
  trials_data <- simulate_trials_with_random_coefficients(N = N, n = n, sd_beta = sd_beta, scenario = scenario) 
  
  # List for MA data set.
  MA_data_list = as.list(rep(0L, length(surrogate_index_estimators)))
  # List for prediction functions.
  surrogate_index_f_list = as.list(rep(0L, length(surrogate_index_estimators)))
  
  # For every surrogate index estimator, train the clinical prediction model and
  # estimate the treatment effects on the corresponding surrogate index.
  for (i in seq_along(1:length(surrogate_index_estimators))) {
    if (scenario == "proof-of-concept") {
      # For the proof of concept scenario, we use adjusted treatment effect
      # estimators for the mean difference.
      method = "adjusted"
      formula = formula(cbind(surrogate, clinical) ~ treatment + X1)
      measure = c("mean difference", "mean difference")
    } else if (scenario == "vaccine") {
      # For the vaccine scenario, we use the mean treatment effect estimators for
      # the log RR. Adjusted estimators would be difficult to implement here. 
      method = "mean"
      formula = NULL
      measure = c("log RR", "log RR")
      # IF the surrogate index is just the surrogate, we still have to use the
      # mean difference for the surrogate index.
      if (surrogate_index_estimators[i] == "surrogate") {
        measure[1] = "mean difference"
      }
      
    } else {
      stop("`scenario` is not valid.")
    }
    
    # Train the clinical prediction model and add the surrogate index to the data
    surrogate_index_f = train_clinical_prediction_model(trials_data, method = surrogate_index_estimators[i])
    trials_data$surr_index <- surrogate_index_f(trials_data)
    # Compute treatment effects for each trial
    MA_data_list[[i]] <- trials_data %>%
      group_by(trial) %>%
      reframe(
        pick(everything()) %>% 
          select(-surrogate) %>% 
          rename(surrogate = surr_index) %>%
          compute_treatment_effects(formula = formula, measure = measure, method = method)
      ) %>%
      ungroup()  # Ensure the result is a tibble without grouping
    surrogate_index_f_list[[i]] = surrogate_index_f
  }
  
  # Return a list with the computed treatment effects (and associated relevant
  # information).
  return(
    list(
      MA_data_list,
      surrogate_index_f_list
    )
  )
}

# Function that approximates the true trial-level correlation given a surrogate
# index function.
rho_MC_approximation = function(f = NULL,
                                N_approximation_MC,
                                n_approximation_MC,
                                sd_beta,
                                scenario) {
  # We follow a different strategy here for simulating trial-level treatment
  # effects. Essentially, we generate data for a single trial and estimate the
  # corresponding treatment effects, all trial-by-trial. This avoids excessive
  # memory usage; at any time, we only have IPD data from a single trial in
  # memory (instead of IPD data from N_approximation_MC trials).
  
  if (scenario == "proof-of-concept") {
    # If f is NULL, the surrogate index is the surrogate.
    if (is.null(f)) {
      f = function(x) x$surrogate
    }
    # For the proof of concept scenario, we again use mean differences as effect
    # measure.
    measure = c("mean difference", "mean difference")
  } else if (scenario == "vaccine") {
    # For the vaccine scenario, we use the log RR as effect measure, except when
    # f is just the surrogate (which is true when f is NULL). In that case, we use the mean difference.
    if (is.null(f)) {
      measure = c("mean difference", "log RR")
      f = function(x) x$surrogate
    } else {
      measure = c("log RR", "log RR")
    }
    # The covariance estimator for log RR is not unbiased (but it is
    # consistent). This means that `n_approximation_MC` should be large enough.
    if (n_approximation_MC < 2000) {
      stop(
        "The number of MC samples should be at least 2000 for the vaccine scenario to ensure a good approximation."
      )
    }
  } else {
    stop("`scenario` is not valid.")
  }
  

  
  # Initialize a list to put the simulated estimated treatment effects etc in.
  # By defining this list before the for loop, we can speed up the for loop.
  treatment_effect_surr = rep(0L, N_approximation_MC)
  treatment_effect_clin = rep(0L, N_approximation_MC)
  vcov_list = as.list(rep(0L, N_approximation_MC))
  
  # Generate estimated treatment effects etc and add them to the already defined
  # list. By using a for loop and the simulate_and_compute_trt_effects function,
  # we avoid the IPD data to accumulate in memory. Indeed, for each iteration in
  # the for loop, an IPD data set is generated and implicitly deleted from
  # memory.
  for (i in seq_along(1:N_approximation_MC)) {
    # Simulate data for one trial with random coefficients
    treatment_effects_index_row = simulate_trials_with_random_coefficients(
      N = 1,
      n = n_approximation_MC,
      sd_beta = sd_beta,
      scenario = scenario
    ) %>%
      # Compute the surrogate index given the function f.
      mutate(surr_index = f(pick(everything()))) %>%
      # Compute treatment effects for each trial (only considering the surrogate
      # index and clinical endpoint).
      select(-surrogate) %>%
      rename(surrogate = surr_index) %>%
      compute_treatment_effects(method = "mean", measure = measure)
    
    treatment_effect_surr[i] = treatment_effects_index_row$treatment_effect_surr
    treatment_effect_clin[i] = treatment_effects_index_row$treatment_effect_clin
    vcov_list[[i]] = treatment_effects_index_row$covariance_matrix[[1]]
    
    # Garbage clean every 100 iterations.
    if (i %% 100 == 0) {
      gc()
    }
  }
  
  # Check for NAs and NaNs. These are removed while returning a warning.
  treatment_effect_surr_nas = is.na(treatment_effect_surr) | is.nan(treatment_effect_surr)
  treatment_effect_clin_nas = is.na(treatment_effect_clin) | is.nan(treatment_effect_clin)
  vcov_nas = sapply(
    X = vcov_list,
    FUN = function(x) {
      any(is.na(x)) | any(is.nan(x))
    }
  )
  total_nas = treatment_effect_surr_nas | treatment_effect_clin_nas | vcov_nas
  if (any(total_nas)) {
    warning("NAs or NaNs were generated while approximating the true correlation. These are removed.")
    treatment_effect_surr = treatment_effect_surr[!total_nas]
    treatment_effect_clin = treatment_effect_clin[!total_nas]
    vcov_list = vcov_list[!total_nas]
  }

  # Estimate the meta-analytic parameters using the moment-based estimator.
  temp = moment_estimator(
    treatment_effect_surr,
    treatment_effect_clin,
    vcov_list, 
    SE = FALSE
  )
  
  
  # Estimate the trial-level Pearson correlation using the delta method. We only
  # need the point estimate in this MC approximation.
  return(temp$coefs[5] / sqrt(temp$coefs[3] * temp$coefs[4]))
}
