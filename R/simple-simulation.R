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
N_MC = 10  # Number of Monte Carlo simulations

N_approximation_MC = 1e4  # Number of Monte Carlo trial replications for the approximation of the true trial-level correlation
n_approximation_MC = 1e3  # Number of patients in each trial for the approximation of the true trial-level correlation

N = 10  # Number of trials in each meta-analytic data set
n = 2e3  # Number of patients in each trial

sd_beta_clin_treatment = 0.10
sd_beta_clin_surrogate_sq = 0.10

# Set up parallel computing
plan(cluster)  # Adjust the number of workers as needed
# Helper Functions --------------------------------------------------------

# Source helper functions.
source("R/helper-functions/simulation-functions.R")
source("R/helper-functions/moment-based-estimator.R")
source("R/helper-functions/delta-method-rho-trial.R")

# Function to train the clinical prediction model and return prediction function.
train_clinical_prediction_model <- function(simulated_data) {
  # Train a correctly specified linear regression model for predicting the
  # clinical endpoint given the surrogate and baseline covariate
  clinical_model <- lm(clinical ~ surrogate*covariate + surrogate ^ 2, data = simulated_data)

  return(
    function(newdata) {
      return(predict(clinical_model, newdata = newdata))
    }
  )
}




# Simulation Study ---------------------------------------------------------

## Simulate IPD data and compute trial-level estimates  --------------------

# Generate N meta-analytic data sets and compute the trial-level Pearson
# correlation for (i) he treatment effects on the surrogate and clinical
# endpoints and (ii) the treatment effects on the surrogate index and clinical
# endpoint.
data_set_indicator = 1:N_MC
meta_analytic_data_list = future_map(data_set_indicator, function(i) {
  generate_meta_analytic_data(N = 10,
                              n = 3e3,
                              train_clinical_prediction_model = train_clinical_prediction_model,
                              sd_beta_clin_treatment = sd_beta_clin_treatment,
                              sd_beta_clin_surrogate_sq = sd_beta_clin_surrogate_sq)
}, .options = furrr_options(seed = TRUE))
# Extract generated meta-analytic data sets with treatment effects on the
# original surrogate.
treatment_effects_list = lapply(meta_analytic_data_list, function(x) x[[1]])
# Extract generated meta-analytic data sets with treatment effects on the
# surrogate index.
treatment_effects_index_list = lapply(meta_analytic_data_list, function(x) x[[2]])
# Extract surrogate index prediction functions.
surrogate_index_f_list = lapply(meta_analytic_data_list, function(x) x[[3]])

## Approximate the true trial-level correlation parameter ------------------

# For each surrogate index prediction function, approximate the true trial-level
# correlation through MC approximation.
rho_true = future_map_dbl(surrogate_index_f_list, function(f) {
  # Simulate data for N_ trials with random coefficients
  trials_data <- simulate_trials_with_random_coefficients(N = n_approximation_MC, n = n_approximation_MC)
  
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
}, .options = furrr_options(seed = TRUE))

true_rho_surrogate_index_tbl = 
  tibble(
    data_set = data_set_indicator,
    rho_true = rho_true
  )

# Approximate the true trial-level correlation using the surrogate endpoint.
trials_data_MC <- simulate_trials_with_random_coefficients(N = N_approximation_MC, n = n_approximation_MC)
treatment_effects_MC <- trials_data_MC %>%
  group_by(trial) %>%
  reframe(compute_treatment_effects(pick(everything()))) %>%
  ungroup()
moment_estimates_MC = moment_estimator(
  treatment_effects_MC$treatment_effect_surr,
  treatment_effects_MC$treatment_effect_clin,
  treatment_effects_MC$covariance_matrix
)
rho_true_surrogate = rho_delta_method(moment_estimates_MC$coefs, moment_estimates_MC$vcov)$rho


# Combine the generated meta-analytic data sets with treatment effects in a
# tibble in tidy format.
meta_analytic_data_simulated =
  bind_rows(
    tibble(data_set = 1:N_MC, treatment_effects = treatment_effects_list) %>% mutate(analysis_type = "Original surrogate", rho_true = rho_true_surrogate),
    tibble(data_set = 1:N_MC, treatment_effects = treatment_effects_index_list) %>% mutate(analysis_type = "Surrogate index") %>%
      left_join(true_rho_surrogate_index_tbl, by = "data_set")
  )

## Moment-based estimator --------------------------------------------------

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

# Save Results -------------------------------------------------------------

# Save the simulated MA data sets with the corresponding rho estimates and
# confidence intervals. 
saveRDS(meta_analytic_data_simulated, "results/raw-results/simple-simulation/meta_analytic_data_simulated.rds")



  

