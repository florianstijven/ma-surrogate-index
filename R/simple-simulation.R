library(tidyverse)

set.seed(123)  # for reproducibility

# Function to generate random coefficients for each trial
generate_random_coefficients <- function() {
  # Ra
  mu_covar1 <- rnorm(1, mean = 1, sd = 1)
  mu_covar2 <- rnorm(1, mean = 1, sd = 1)
  
  # Random coefficients for surrogate model (from normal distribution)
  beta_surr_covar1 <- rnorm(1, mean = 1, sd = 1)
  beta_surr_covar2 <- rnorm(1, mean = 0.5, sd = 1)
  beta_surr_treatment <- rnorm(1, mean = 0, sd = 1)
  beta_surr_interaction <- rnorm(1, mean = 0.5, sd = 0.2)
  beta_surr_treatment_interaction <- rnorm(1, mean = 0.3, sd = 0.05)
  
  # Random coefficients for clinical model (from normal distribution)
  beta_clin_covar1 <- rnorm(1, mean = 2, sd = sd_multiplier * 1)
  beta_clin_covar2 <- rnorm(1, mean = 1, sd =  sd_multiplier * 0.5)
  beta_clin_surrogate <- rnorm(1, mean = 0.8, sd = sd_multiplier * 0.01)
  beta_clin_treatment <- rnorm(1, mean = true_effect, sd = sd_multiplier * 0.01)
  beta_clin_interaction <- rnorm(1, mean = 0.4, sd = sd_multiplier * 0.3)
  beta_clin_treatment_interaction <- rnorm(1, mean = 0.2, sd = sd_multiplier * 0.1)
  beta_clin_surrogate_interaction <- rnorm(1, mean = 0.6, sd = sd_multiplier * 0.3)
  
  # Return the list of coefficients
  return(list(
    mu_covar1 = mu_covar1,
    mu_covar2 = mu_covar2,
    
    beta_surr_covar1 = beta_surr_covar1, 
    beta_surr_covar2 = beta_surr_covar2,
    beta_surr_treatment = beta_surr_treatment, 
    beta_surr_interaction = beta_surr_interaction, 
    beta_surr_treatment_interaction = beta_surr_treatment_interaction,
    
    beta_clin_covar1 = beta_clin_covar1, 
    beta_clin_covar2 = beta_clin_covar2,
    beta_clin_surrogate = beta_clin_surrogate,
    beta_clin_treatment = beta_clin_treatment,
    beta_clin_interaction = beta_clin_interaction,
    beta_clin_treatment_interaction = beta_clin_treatment_interaction,
    beta_clin_surrogate_interaction = beta_clin_surrogate_interaction
  ))
}

# Function to simulate a single trial with passed coefficients and interaction terms
simulate_trial <- function(n = 100, p1 = 0.5, coefficients) {
  # Randomization: Treatment assignment (1 or 0) based on probability p1
  treatment <- rbinom(n, 1, p1)  # 1 = treatment, 0 = control
  
  # Baseline covariates: Two covariates from normal distribution
  covar1 <- rnorm(n, mean = coefficients$mu_covar1)
  covar2 <- rnorm(n, mean = coefficients$mu_covar2)
  
  # Surrogate endpoint with passed random coefficients
  surrogate <- coefficients$beta_surr_covar1 * covar1 + coefficients$beta_surr_covar2 * covar2 + 
    coefficients$beta_surr_treatment * treatment + 
    (covar1 * covar2) * coefficients$beta_surr_interaction + 
    (treatment * covar1) * coefficients$beta_surr_treatment_interaction + rnorm(n)
  
  # Clinical endpoint with passed random coefficients
  clinical <- coefficients$beta_clin_covar1 * covar1 + coefficients$beta_clin_covar2 * covar2 + 
    coefficients$beta_clin_surrogate * surrogate + coefficients$beta_clin_treatment * treatment + 
    (covar1 * covar2) * coefficients$beta_clin_interaction + 
    (treatment * covar1) * coefficients$beta_clin_treatment_interaction +
    (surrogate * covar1) * coefficients$beta_clin_surrogate_interaction + rnorm(n)
  
  # Data frame for a single trial
  trial_data <- data.frame(treatment, covar1, covar2, surrogate, clinical)
  
  return(trial_data)
}

# Function to simulate data for multiple trials with random coefficients passed as arguments
simulate_trials_with_random_coefficients <- function(N = 10, n = 100, p1 = 0.5, true_effect = 0, sd_multiplier = 1) {
  trials <- lapply(1:N, function(i) {
    coefficients <- generate_random_coefficients(true_effect = true_effect, sd_multiplier = sd_multiplier)  # generate new random coefficients for each trial
    trial_data <- simulate_trial(n = n, p1 = p1, coefficients = coefficients)
    trial_data$trial <- i  # Add a trial-level variable
    return(trial_data)
  })

  # Combine all trial data into a single tibble
  all_trials_data <- bind_rows(trials) %>%
    as_tibble()  # Ensure the result is a tibble
  
  return(all_trials_data)
}

library(sandwich)
library(lmtest)

# Function to compute treatment effects with robust standard errors and covariance matrix
compute_treatment_effects <- function(trial_data) {
  # Fit a linear model for the surrogate endpoint (surrogate ~ treatment + covariates)
  surrogate_model <- lm(surrogate ~ treatment + covar1 + covar2, data = trial_data)
  
  # Compute robust standard errors for the surrogate using the sandwich estimator
  se_surr <- sqrt(diag(vcovHC(surrogate_model, type = "HC3")))["treatment"]
  treatment_effect_surr <- coef(surrogate_model)["treatment"]
  
  # Fit a linear model for the clinical endpoint (clinical ~ treatment + covariates)
  clinical_model <- lm(clinical ~ treatment + covar1 + covar2, data = trial_data)
  
  # Compute robust standard errors for the clinical using the sandwich estimator
  se_clin <- sqrt(diag(vcovHC(clinical_model, type = "HC3")))["treatment"]
  treatment_effect_clin <- coef(clinical_model)["treatment"]
  
  # Compute the covariance matrix between the treatment effects on the surrogate and clinical endpoints
  # For this, we need to extract the covariance matrix that jointly captures the uncertainty of the treatment effects
  # across both models.
  # Here we compute the covariance matrix between the two treatment effects.
  
  # Combine the coefficients into a single model to compute the covariance matrix for both treatment effects
  combined_model <- lm(cbind(surrogate, clinical) ~ treatment + covar1 + covar2, data = trial_data)
  
  # Compute the covariance matrix using the sandwich estimator
  cov_matrix <- sandwich::vcovHC(combined_model, type = "HC3")
  
  # Extract the covariance between the treatment effects on surrogate and clinical
  cov_surr_clin <- cov_matrix["clinical:treatment", "surrogate:treatment"]  # Covariance between treatment effects
  
  # Return a tibble with the treatment effects, standard errors, and covariance matrix
  treatment_effects <- tibble(
    treatment_effect_surr = treatment_effect_surr,  # Treatment effect on surrogate
    se_surr = se_surr,  # Robust standard error for surrogate treatment effect
    treatment_effect_clin = treatment_effect_clin,  # Treatment effect on clinical endpoint
    se_clin = se_clin,  # Robust standard error for clinical treatment effect
    cov_surr_clin = cov_surr_clin,  # Covariance between surrogate and clinical treatment effects
    covariance_matrix = list(cov_matrix)  # Store full covariance matrix as a list
  )
  
  return(treatment_effects)
}


# Function to train the clinical prediction model and add predictions to the dataset
train_clinical_prediction_model <- function(simulated_data) {
  # Train a linear regression model for clinical endpoint prediction using surrogate and covariates
  clinical_model <- mgcv::gam(clinical ~ s(surrogate) + s(covar1) + s(covar2), data = simulated_data)
  
  # Generate predictions from the trained model for all observations
  simulated_data$surr_index <- predict(clinical_model, newdata = simulated_data)
  
  return(simulated_data)
}


# Simulate data for 5 trials with random coefficients
trials_data <- simulate_trials_with_random_coefficients(N = 50, n = 1000, sd_multiplier = 0) 

trials_data = train_clinical_prediction_model(trials_data)

plot(trials_data$surr_index, trials_data$clinical)


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


# Create a scatter plot with treatment effects on surrogate vs. clinical, including 95% CI error bars
ggplot(treatment_effects, aes(x = treatment_effect_surr, y = treatment_effect_clin)) +
  # Scatter plot of treatment effects
  geom_point(color = "purple", size = 3) +
  
  # Horizontal error bars for surrogate treatment effect
  geom_errorbarh(aes(xmin = treatment_effect_surr - 1.96 * se_surr, 
                     xmax = treatment_effect_surr + 1.96 * se_surr), 
                 height = 0.2, color = "purple") +
  
  # Vertical error bars for clinical treatment effect
  geom_errorbar(aes(ymin = treatment_effect_clin - 1.96 * se_clin, 
                    ymax = treatment_effect_clin + 1.96 * se_clin), 
                width = 0.2, color = "purple") +
  
  # Labels and title
  labs(
    x = "Treatment Effect on Surrogate Endpoint",
    y = "Treatment Effect on Clinical Endpoint",
    title = "Treatment Effects on Surrogate vs. Clinical Endpoint with 95% Confidence Intervals"
  ) +
  
  # Customize the theme for better visualization
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1)
  )

# Create a scatter plot with treatment effects on surrogate vs. clinical, including 95% CI error bars
ggplot(treatment_effects_index, aes(x = treatment_effect_surr, y = treatment_effect_clin)) +
  # Scatter plot of treatment effects
  geom_point(color = "purple", size = 3) +
  
  # Horizontal error bars for surrogate treatment effect
  geom_errorbarh(aes(xmin = treatment_effect_surr - 1.96 * se_surr, 
                     xmax = treatment_effect_surr + 1.96 * se_surr), 
                 height = 0.2, color = "purple") +
  
  # Vertical error bars for clinical treatment effect
  geom_errorbar(aes(ymin = treatment_effect_clin - 1.96 * se_clin, 
                    ymax = treatment_effect_clin + 1.96 * se_clin), 
                width = 0.2, color = "purple") +
  
  # Labels and title
  labs(
    x = "Treatment Effect on Surrogate Endpoint",
    y = "Treatment Effect on Clinical Endpoint",
    title = "Treatment Effects on Surrogate vs. Clinical Endpoint with 95% Confidence Intervals"
  ) +
  
  # Customize the theme for better visualization
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1)
  ) +
  geom_abline(slope = 1)



  

