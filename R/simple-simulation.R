# Preliminaries -----------------------------------------------------------

# Load all R packages
library(tidyverse)
library(sandwich)
library(lmtest)

# This script simulates data. We set a seed to ensure reproducibility.
set.seed(123) 

# Helper Functions --------------------------------------------------------

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


# Function to train the clinical prediction model and add predictions to the dataset
train_clinical_prediction_model <- function(simulated_data) {
  # Train a correctly specified linear regression model for predicting the
  # clinical endpoint given the surrogate and baseline covariate
  clinical_model <- lm(clinical ~ surrogate*covariate, data = simulated_data)

  # Generate predictions from the trained model for all observations
  simulated_data$surr_index <- predict(clinical_model, newdata = simulated_data)
  
  return(simulated_data)
}


# Data Simulation ---------------------------------------------------------



# Simulate data for 5 trials with random coefficients
trials_data <- simulate_trials_with_random_coefficients(N = 10, n = 3e3) 

trials_data = train_clinical_prediction_model(trials_data)

plot(trials_data$surr_index, trials_data$clinical)
plot(trials_data$surrogate, trials_data$clinical)

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

# Estimate the meta-analytic parameters using the moment-based estimator.
temp = moment_estimator(treatment_effects$treatment_effect_surr, 
                 treatment_effects$treatment_effect_clin, 
                 treatment_effects$covariance_matrix)
rho_delta_method(temp$coefs, temp$vcov)


temp = moment_estimator(treatment_effects_index$treatment_effect_surr, 
                        treatment_effects_index$treatment_effect_clin, 
                        treatment_effects_index$covariance_matrix)
rho_delta_method(temp$coefs, temp$vcov)

  

