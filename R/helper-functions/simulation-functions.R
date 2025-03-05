# Function to generate trial-level random parameters. In the current setup, we
# only have three trial-level random parameters: p1, beta_surr_treatment, and
# beta_clin_treatment. p1 is the proportion for the baseline binary covariate.
# beta_surr_treatment is the mean treatment effect on the surrogate endpoint.
# beta_clin_treatment is the "additional" treatment effect on the clinical
# endpoint in a model conditional on the baseline covariate and surrogate.
generate_random_coefficients <- function(sd_beta_clin_treatment = 0.10, sd_beta_clin_surrogate_sq = 0.10) {
  # Sample trial-level proportion parameter for baseline covariate from uniform
  # distribution. 
  p1 <- runif(1, min = 0.3, max = 0.7)
  
  # Sample random treatment effect on surrogate. 
  beta_surr_treatment <- rnorm(1, mean = 0.5, sd = 0.3)
  
  # Sample random treatment effect on clinical endpoint.
  beta_clin_treatment <- rnorm(1, mean = 0, sd = sd_beta_clin_treatment)
  beta_clin_surrogate_sq <- rnorm(1, mean = 0, sd = sd_beta_clin_surrogate_sq)
  
  # Return the list of sampled parameters.
  return(
    list(
      p1 = p1,
      beta_surr_treatment = beta_surr_treatment,
      beta_clin_treatment = beta_clin_treatment,
      beta_clin_surrogate_sq = beta_clin_surrogate_sq
    )
  )
}

# Function to simulate a single trial with passed trial-level parameters.
simulate_trial <- function(n, coefficients) {
  # Randomization: Treatment assignment (1 or 0). Simple randomization is
  # assumed.
  treatment <- rbinom(n, 1, 0.5)  # 1 = treatment, 0 = control
  
  # Baseline covariates: One covariate from a Bernoulli distribution with
  # probability p1.
  covariate <- rbinom(n, 1, coefficients$p1)
  
  # Surrogate endpoint: Linear relationship with treatment and random noise. The
  # distribution of the surrogate does not depend on the baseline covariate.
  surrogate <- coefficients$beta_surr_treatment * treatment + rnorm(n)
  
  # Simulate clinical endpoint with no random coefficients. We're assuming that
  # the regression of the clinical endpoint on the surrogate and covariate is
  # the same across trials modulus some small random treatment effect.
  clinical <- -1 * surrogate * covariate + surrogate * (1 - covariate) + 
    coefficients$beta_clin_treatment * treatment +
    coefficients$beta_clin_surrogate_sq * surrogate ^ 2 +
    rnorm(n)
  
  # Data frame for a single trial
  trial_data <- data.frame(treatment, covariate, surrogate, clinical)
  
  return(trial_data)
}

# Function to simulate data for multiple trials with a common within-trial
# sample size. The random trial-level parameters are automatically sampled
# within this function.
simulate_trials_with_random_coefficients <- function(N, n, sd_beta_clin_treatment, sd_beta_clin_surrogate_sq) {
  
  trials <- lapply(1:N, function(i,
                                 sd_beta_clin_treatment,
                                 sd_beta_clin_surrogate_sq) {
    coefficients <- generate_random_coefficients(sd_beta_clin_treatment, sd_beta_clin_surrogate_sq)  # generate new random coefficients for each trial
    trial_data <- simulate_trial(n = n, coefficients = coefficients)
    trial_data$trial <- i  # Add a trial-level variable
    return(trial_data)
  }, sd_beta_clin_treatment = sd_beta_clin_treatment, sd_beta_clin_surrogate_sq = sd_beta_clin_surrogate_sq)
  
  # Combine all trial data into a single tibble
  all_trials_data <- bind_rows(trials) %>%
    as_tibble()  # Ensure the result is a tibble
  
  return(all_trials_data)
}

# Function to compute treatment effects with robust standard errors and covariance matrix
compute_treatment_effects <- function(trial_data, method = "adjusted") {
  # If adjusted treatment effects are requires, we fit two univariate linear
  # regression models and estimate the covariance of the treatment effect
  # estimators through the sandwich estimators.
  if (method == "adjusted") {
    # Combine the coefficients into a single model to compute the covariance matrix for both treatment effects
    combined_model <- lm(cbind(surrogate, clinical) ~ treatment + covariate, data = trial_data)
    
    # Compute the covariance matrix using the sandwich estimator
    cov_matrix <- sandwich::vcovHC(combined_model, type = "HC3")
    
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
    # Unadjusted treatment effect estimates are computed through differences of
    # sample means.
    treatment_effect_surr = mean(trial_data$surrogate[trial_data$treatment == 1]) - mean(trial_data$surrogate[trial_data$treatment == 0])
    treatment_effect_clin = mean(trial_data$clinical[trial_data$treatment == 1]) - mean(trial_data$clinical[trial_data$treatment == 0])
    
    # Sample covariance matrices in each treatment group.
    vcov_0 = var(trial_data[trial_data$treatment == 0, c("surrogate", "clinical")])
    vcov_1 = var(trial_data[trial_data$treatment == 1, c("surrogate", "clinical")])
    
    # Treatment group sample sizes
    n0 = sum(trial_data$treatment == 0)
    n1 = sum(trial_data$treatment == 1)
    
    # Estimated covariance matrix for the treatment effect estimators.
    covariance_matrix = (1 / n1) * vcov_1 + (1 / n0) * vcov_0
    
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
generate_meta_analytic_data <- function(N, n, sd_beta_clin_treatment, sd_beta_clin_surrogate_sq, train_clinical_prediction_model) {
  # Simulate data for N trials with random coefficients
  trials_data <- simulate_trials_with_random_coefficients(N = N, n = n, sd_beta_clin_treatment, sd_beta_clin_surrogate_sq) 
  
  # Train the clinical prediction model and add the surrogate index to the data
  surrogate_index_f = train_clinical_prediction_model(trials_data)
  trials_data$surr_index <- surrogate_index_f(trials_data)
  
  # Compute treatment effects for each trial
  treatment_effects <- trials_data %>%
    group_by(trial) %>%
    do(compute_treatment_effects(.)) %>%
    ungroup()  # Ensure the result is a tibble without grouping
  
  # Compute treatment effects for each trial
  treatment_effects_index <- trials_data %>%
    group_by(trial) %>%
    reframe(
      pick(everything()) %>% 
        select(-surrogate) %>% 
        rename(surrogate = surr_index) %>%
        compute_treatment_effects()
    ) %>%
    ungroup()  # Ensure the result is a tibble without grouping
  
  # Return a list with the computed treatment effects (and associated relevant
  # information).
  return(
    list(
      treatment_effects,
      treatment_effects_index,
      surrogate_index_f
    )
  )
}

# Function that approximates the true trial-level correlation given a surrogate
# index function.
rho_MC_approximation = function(f,
                                N_approximation_MC,
                                n_approximation_MC,
                                sd_beta_clin_treatment,
                                sd_beta_clin_surrogate_sq) {
  # We follow a different strategy here for simulating trial-level treatment
  # effects. Essentially, we generate data for a single trial and estimate the
  # corresponding treatment effects, all trial-by-trial. This avoids excessive
  # memory usage; at any time, we only have IPD data from a single trial in
  # memory (instead of IPD data from N_approximation_MC trials).
  
  # Function that simulates a single trial and computes the treatment effects
  # and covariance matrix based on sample means.
  simulate_and_compute_trt_effects = function() {
    # Simulate data for one trial with random coefficients
    treatment_effects_index_row = simulate_trials_with_random_coefficients(
      N = 1,
      n = n_approximation_MC,
      sd_beta_clin_treatment = sd_beta_clin_treatment,
      sd_beta_clin_surrogate_sq = sd_beta_clin_surrogate_sq
    ) %>%
      # Compute the surrogate index given the function f.
      mutate(surr_index = f(pick(everything()))) %>%
      # Compute treatment effects for each trial (only considering the surrogate
      # index and clinical endpoint).
      select(-surrogate) %>%
      rename(surrogate = surr_index) %>%
      compute_treatment_effects(method = "mean")
    
    return(treatment_effects_index_row)
  }
  
  # Initialize a list to put the simulated estimated treatment effects etc in.
  # By defining this list before the for loop, we can speed up the for loop.
  treatment_effects_index_list = as.list(
    rep(0L, N_approximation_MC)
  )
  
  # Generate estimated treatment effects etc and add them to the already defined
  # list. By using a for loop and the simulate_and_compute_trt_effects function,
  # we avoid the IPD data to accumulate in memory. Indeed, for each iteration in
  # the for loop, an IPD data set is generated and implicitly deleted from
  # memory.
  for (i in seq_along(1:N_approximation_MC)) {
    treatment_effects_index_list[[i]] = simulate_and_compute_trt_effects()
  }
  
  treatment_effects_index = bind_rows(treatment_effects_index_list)

  # Estimate the meta-analytic parameters using the moment-based estimator.
  temp = moment_estimator(
    treatment_effects_index$treatment_effect_surr,
    treatment_effects_index$treatment_effect_clin,
    treatment_effects_index$covariance_matrix
  )
  
  
  
  # Estimate the trial-level Pearson correlation using the delta method. We only
  # need the point estimate in this MC approximation.
  rho_delta_method(temp$coefs, temp$vcov)$rho
}
