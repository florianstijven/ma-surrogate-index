# Preliminaries -----------------------------------------------------------

# Load all R packages
library(tidyverse)
library(sandwich)
library(lmtest)
library(future)
library(furrr)

# This script simulates data. We set a seed to ensure reproducibility.
set.seed(123) 

# Set of parameters controlling the simulations.
N_MC = 100  # Number of Monte Carlo simulations
n_approximation_MC = 1000  # Number of Monte Carlo simulations for the approximation of the true trial-level correlation

# Source helper functions.
source("R/helper-functions/simulation-functions.R")
source("R/helper-functions/moment-based-estimator.R")
source("R/helper-functions/delta-method-rho-trial.R")

# Set up parallel computing
plan(cluster)  # Adjust the number of workers as needed
# Helper Functions --------------------------------------------------------

# Function to train the clinical prediction model and return prediction function.
train_clinical_prediction_model <- function(simulated_data) {
  # Train a correctly specified linear regression model for predicting the
  # clinical endpoint given the surrogate and baseline covariate
  clinical_model <- lm(clinical ~ surrogate*covariate, data = simulated_data)

  return(
    function(newdata) {
      return(predict(clinical_model, newdata = newdata))
    }
  )
}




# Simulation Study ---------------------------------------------------------

# Generate N meta-analytic data sets and compute the trial-level Pearson
# correlation for (i) he treatment effects on the surrogate and clinical
# endpoints and (ii) the treatment effects on the surrogate index and clinical
# endpoint.
data_set_indicator = 1:N_MC
meta_analytic_data_list = future_map(data_set_indicator, function(i) {
  generate_meta_analytic_data(N = 10, n = 3e3, train_clinical_prediction_model)
})
# Extract generated meta-analytic data sets with treatment effects on the
# original surrogate.
treatment_effects_list = lapply(meta_analytic_data_list, function(x) x[[1]])
# Extract generated meta-analytic data sets with treatment effects on the
# surrogate index.
treatment_effects_index_list = lapply(meta_analytic_data_list, function(x) x[[2]])
# Extract surrogate index prediction functions.
surrogate_index_f_list = lapply(meta_analytic_data_list, function(x) x[[3]])
# For each surrogate index prediction function, approximate the true trial-level
# correlation through MC approximation.
rho_true = future_map_dbl(surrogate_index_f_list, function(f) {
  # Simulate data for N_ trials with random coefficients
  trials_data <- simulate_trials_with_random_coefficients(N = n_approximation_MC, n = 1e3) 
  
  # Compute the surrogate index given the function f.
  trials_data$surr_index <- f(trials_data)
  
  # Compute treatment effects for each trial (only considering the surrogate
  # index and clinical endpoint).
  treatment_effects_index <- trials_data %>%
    group_by(trial) %>%
    reframe(
      pick(everything()) %>% 
        select(-surrogate) %>% 
        rename(surrogate = surr_index) %>%
        compute_treatment_effects()
    ) %>%
    ungroup()  # Ensure the result is a tibble without grouping
 
  # Estimate the meta-analytic parameters using the moment-based estimator.
  temp = moment_estimator(treatment_effects_index$treatment_effect_surr, 
                          treatment_effects_index$treatment_effect_clin, 
                          treatment_effects_index$covariance_matrix)
  
  # Estimate the trial-level Pearson correlation using the delta method. We only
  # need the point estimate in this MC approximation.
  rho_delta_method(temp$coefs, temp$vcov)$rho
})
true_rho_surrogate_index_tbl = 
  tibble(
    data_set = data_set_indicator,
    rho_true = rho_true
  )

# Combine the generated meta-analytic data sets with treatment effects in a
# tibble in tidy format.
meta_analytic_data_simulated =
  bind_rows(
    tibble(
      data_set = 1:N_MC,
      treatment_effects = treatment_effects_list
    ) %>% mutate(analysis_type = "Original surrogate"),
    tibble(
      data_set = 1:N_MC,
      treatment_effects = treatment_effects_index_list
    ) %>% mutate(analysis_type = "Surrogate index")
  )

# Add approximated true trial-level rhos.
meta_analytic_data_simulated = meta_analytic_data_simulated %>%
  left_join(true_rho_surrogate_index_tbl, by = "data_set")

# Estimate trial-level Pearson correlation for the generated meta-analytic data
# sets of treatment effects and add everything to a tibble.
meta_analytic_data_simulated = meta_analytic_data_simulated %>%
  rowwise() %>%
  mutate(
    # Estimate the meta-analytic parameters using the moment-based estimator.
    moment_estimates = list(
      moment_estimator(treatment_effects$treatment_effect_surr, 
                       treatment_effects$treatment_effect_clin, 
                       treatment_effects$covariance_matrix)
    ),
    # Estimate the trial-level Pearson correlation using the delta method.
    rho_delta_method = list(
      rho_delta_method(
        moment_estimates$coefs,
        moment_estimates$vcov
      )
    ),
    rho_est = rho_delta_method$rho,
    rho_se = rho_delta_method$se,
    rho_ci_lower = rho_delta_method$ci[1],
    rho_ci_upper = rho_delta_method$ci[2]
  )

# Compute the coverage of the estimated confidence intervals.
meta_analytic_data_simulated %>%
  group_by(analysis_type) %>%
  summarize(
    coverage = mean(rho_true >= rho_ci_lower & rho_true <= rho_ci_upper, na.rm = TRUE)
  )


# Simulate data for 5 trials with random coefficients
trials_data <- meta_analytic_data_simulated$treatment_effects[[1]]

# Create a scatter plot with treatment effects on surrogate vs. clinical, including 95% CI error bars
ggplot(trials_data, aes(x = treatment_effect_surr, y = treatment_effect_clin)) +
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
ggplot(meta_analytic_data_simulated$treatment_effects[[101]], aes(x = treatment_effect_surr, y = treatment_effect_clin)) +
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

  

