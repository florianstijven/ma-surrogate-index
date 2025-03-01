# Function to generate trial-level random parameters. In the current setup, we
# only have three trial-level random parameters: p1, beta_surr_treatment, and
# beta_clin_treatment. p1 is the proportion for the baseline binary covariate.
# beta_surr_treatment is the mean treatment effect on the surrogate endpoint.
# beta_clin_treatment is the "additional" treatment effect on the clinical
# endpoint in a model conditional on the baseline covariate and surrogate.
generate_random_coefficients <- function() {
  # Sample trial-level proportion parameter for baseline covariate from uniform
  # distribution. 
  p1 <- runif(1, min = 0.3, max = 0.7)
  
  # Sample random treatment effect on surrogate. 
  beta_surr_treatment <- rnorm(1, mean = 0.5, sd = 0.3)
  
  # Sample random treatment effect on clinical endpoint.
  beta_clin_treatment <- rnorm(1, mean = 0, sd = 0.05)
  
  # Return the list of sampled parameters.
  return(
    list(
      p1 = p1,
      beta_surr_treatment = beta_surr_treatment,
      beta_clin_treatment = beta_clin_treatment
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
    rnorm(n)
  
  # Data frame for a single trial
  trial_data <- data.frame(treatment, covariate, surrogate, clinical)
  
  return(trial_data)
}

# Function to simulate data for multiple trials with a common within-trial
# sample size. The random trial-level parameters are automatically sampled
# within this function.
simulate_trials_with_random_coefficients <- function(N = 10, n = 100) {
  
  trials <- lapply(1:N, function(i) {
    coefficients <- generate_random_coefficients()  # generate new random coefficients for each trial
    trial_data <- simulate_trial(n = n, coefficients = coefficients)
    trial_data$trial <- i  # Add a trial-level variable
    return(trial_data)
  })
  
  # Combine all trial data into a single tibble
  all_trials_data <- bind_rows(trials) %>%
    as_tibble()  # Ensure the result is a tibble
  
  return(all_trials_data)
}

# Function to compute treatment effects with robust standard errors and covariance matrix
compute_treatment_effects <- function(trial_data) {
  # Fit a linear model for the surrogate endpoint (surrogate ~ treatment + covariates)
  surrogate_model <- lm(surrogate ~ treatment + covariate, data = trial_data)
  
  # Compute robust standard errors for the surrogate using the sandwich estimator
  se_surr <- sqrt(diag(vcovHC(surrogate_model, type = "HC3")))["treatment"]
  treatment_effect_surr <- coef(surrogate_model)["treatment"]
  
  # Fit a linear model for the clinical endpoint (clinical ~ treatment + covariates)
  clinical_model <- lm(clinical ~ treatment + covariate, data = trial_data)
  
  # Compute robust standard errors for the clinical using the sandwich estimator
  se_clin <- sqrt(diag(vcovHC(clinical_model, type = "HC3")))["treatment"]
  treatment_effect_clin <- coef(clinical_model)["treatment"]
  
  # Compute the covariance matrix between the treatment effects on the surrogate and clinical endpoints
  # For this, we need to extract the covariance matrix that jointly captures the uncertainty of the treatment effects
  # across both models.
  # Here we compute the covariance matrix between the two treatment effects.
  
  # Combine the coefficients into a single model to compute the covariance matrix for both treatment effects
  combined_model <- lm(cbind(surrogate, clinical) ~ treatment + covariate, data = trial_data)
  
  # Compute the covariance matrix using the sandwich estimator
  cov_matrix <- sandwich::vcovHC(combined_model, type = "HC3")
  
  # Extract the covariance matrix between the treatment effects on surrogate and clinical
  cov_surr_clin <- cov_matrix[c("clinical:treatment", "surrogate:treatment"), c("clinical:treatment", "surrogate:treatment")]  # Covariance between treatment effects
  
  # Return a tibble with the treatment effects, standard errors, and covariance matrix
  treatment_effects <- tibble(
    treatment_effect_surr = treatment_effect_surr,  # Treatment effect on surrogate
    se_surr = se_surr,  # Robust standard error for surrogate treatment effect
    treatment_effect_clin = treatment_effect_clin,  # Treatment effect on clinical endpoint
    se_clin = se_clin,  # Robust standard error for clinical treatment effect
    covariance_matrix = list(cov_surr_clin)  # Store covariance matrix as a list
  )
  
  return(treatment_effects)
}

# Function that generates a meta-analytic data set, adds the estimated surrogate
# index as a variable, and computes the treatment effects on (i) the surrogate
# and clinical endpoints and (ii) the surrogate index and clinical endpoint.
# This function uses the previously defined functions.
generate_meta_analytic_data <- function(N = 10, n = 100, train_clinical_prediction_model) {
  # Simulate data for N trials with random coefficients
  trials_data <- simulate_trials_with_random_coefficients(N = N, n = n) 
  
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