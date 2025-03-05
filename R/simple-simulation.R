# Preliminaries -----------------------------------------------------------

a = Sys.time()
# Load all R packages
library(tidyverse)
library(sandwich)
library(lmtest)
library(future)
library(furrr)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}

# This script simulates data. We set a seed to ensure reproducibility.
set.seed(123) 

# Set of parameters controlling the simulations.
N_MC = 50  # Number of Monte Carlo simulations

N_approximation_MC = 5e4  # Number of Monte Carlo trial replications for the approximation of the true trial-level correlation
n_approximation_MC = 1e2  # Number of patients in each trial for the approximation of the true trial-level correlation

N = 10  # Number of trials in each meta-analytic data set
n = 2e3  # Number of patients in each trial

sd_beta_clin_treatment = c(0.03, 0.10)
sd_beta_clin_surrogate_sq = c(0.03, 0.10)

B = 2e3

# Tibble with simulation parameters. Each row corresponds to a data-generating
# mechanism.
dgm_param_tbl = tibble(sd_beta_clin_treatment, sd_beta_clin_surrogate_sq)


# Helper Functions --------------------------------------------------------

# Source helper functions.
source("R/helper-functions/simulation-functions.R")
source("R/helper-functions/moment-based-estimator.R")
source("R/helper-functions/delta-method-rho-trial.R")
source("R/helper-functions/multiplier-bootstrap.R")

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

meta_analytic_data = expand_grid(data_set_indicator,
                                 dgm_param_tbl) %>%
  mutate(
    list_of_ma_data_objects = future_map2(
      sd_beta_clin_treatment, 
      sd_beta_clin_surrogate_sq, 
      function(sd_beta_clin_treatment,
               sd_beta_clin_surrogate_sq) {
        generate_meta_analytic_data(
          N = N,
          n = n,
          train_clinical_prediction_model = train_clinical_prediction_model,
          sd_beta_clin_treatment = sd_beta_clin_treatment,
          sd_beta_clin_surrogate_sq = sd_beta_clin_surrogate_sq
        )
      }, .options = furrr_options(seed = TRUE)),
    # Extract generated meta-analytic data sets with treatment effects on the
    # original surrogate.
    treatment_effects = lapply(list_of_ma_data_objects, function(x)
      x[[1]]),
    # Extract generated meta-analytic data sets with treatment effects on the
    # surrogate index.
    treatment_effects_index = lapply(list_of_ma_data_objects, function(x)
      x[[2]]),
    # Extract surrogate index prediction functions.
    surrogate_index_f_list = lapply(list_of_ma_data_objects, function(x)
      x[[3]])
  ) %>%
  select(-list_of_ma_data_objects)



## Approximate the true trial-level correlation parameter ------------------

# For each surrogate index prediction function, approximate the true trial-level
# correlation through MC approximation.
meta_analytic_data = meta_analytic_data %>%
  mutate(
    rho_true =  future_pmap_dbl(
      .l = list(
        f = surrogate_index_f_list,
        sd_beta_clin_treatment = sd_beta_clin_treatment,
        sd_beta_clin_surrogate_sq = sd_beta_clin_surrogate_sq
      ),
      .f = rho_MC_approximation,
      N_approximation_MC = N_approximation_MC,
      n_approximation_MC = n_approximation_MC,
      .options = furrr_options(seed = TRUE)
    )
  )

# Approximate the true trial-level correlation using the surrogate endpoint.
rho_true_surrogate_tbl = dgm_param_tbl %>%
  rowwise(everything()) %>%
  summarise(
    rho_true_surrogate = rho_MC_approximation(
      f = function(x)
        x$surrogate,
      N_approximation_MC = N_approximation_MC,
      n_approximation_MC = n_approximation_MC,
      sd_beta_clin_treatment = sd_beta_clin_treatment,
      sd_beta_clin_surrogate_sq = sd_beta_clin_surrogate_sq
    )
  )

# Add the true trial-level rho parameters to the tibble with the simulation data
# sets.
meta_analytic_data = meta_analytic_data %>%
  left_join(rho_true_surrogate_tbl)



# Convert the simulated data into long format with one row per simulated data
# set of trial-level treatment effects. We also add a variable indicating
# whether these treatment effects correspond to treatment effects on the
# surrogate or the surrogate index.
meta_analytic_data_simulated = meta_analytic_data %>%
  pivot_longer(
    cols = c(treatment_effects, treatment_effects_index),
    values_to = "treatment_effects_tbl",
    names_to = "analysis_type"
  ) %>%
  # Add informative variable indicating the type of surrogate.
  mutate(
    analysis_type = ifelse(
      analysis_type == "treatment_effects",
      "surrogate",
      "surrogate index"
    )
  ) %>%
  # Only keep the true trial-level correlation parameter that corresponds to the analysis type.
  mutate(rho_true = ifelse(analysis_type == "surrogate", rho_true_surrogate, rho_true))

# Remove redundant tibble from memory.
rm("meta_analytic_data")

# Surrogate index prediction function for the setting where we simply use the
# surrogate endpoint.
surrogate_f = function(x)
  x$surrogate

# If the analysis type is the surrogate, we change the surrogate index function to a function that simpy return the surrogate.
meta_analytic_data_simulated = bind_rows(
  meta_analytic_data_simulated %>%
    filter(analysis_type == "surrogate index"),
  meta_analytic_data_simulated %>%
    filter(analysis_type == "surrogate") %>%
    mutate(surrogate_index_f_list = list(surrogate_f))
) %>% # Drop columns that have become superfluous.
  select(-rho_true_surrogate)


## Moment-based estimator --------------------------------------------------

# Estimate trial-level Pearson correlation for the generated meta-analytic data
# sets of treatment effects and add everything to a tibble.
meta_analytic_data_simulated = meta_analytic_data_simulated %>%
  cross_join(expand_grid(
    estimator_adjustment = c("none", "N - 1"),
    sandwich_adjustment = c("none", "N - 1")
  )) %>%
  rowwise(everything()) %>%
  summarise(
    # Estimate the meta-analytic parameters using the moment-based estimator.
    moment_estimates = list(
      moment_estimator(
        treatment_effects_tbl$treatment_effect_surr,
        treatment_effects_tbl$treatment_effect_clin,
        treatment_effects_tbl$covariance_matrix,
        estimator_adjustment = estimator_adjustment,
        sandwich_adjustment = sandwich_adjustment
      )
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
  ) %>% 
  ungroup() %>%
  # Drop superfluous variables
  select(-moment_estimates, -rho_delta_method)

statistic_function_factory = function(estimator_adjustment) {
  statistic_f = function(data, weights) {
    moment_estimate = moment_estimator(
      alpha_hat = data$treatment_effect_surr,
      beta_hat = data$treatment_effect_clin,
      vcov_list = data$covariance_matrix,
      estimator_adjustment = estimator_adjustment,
      weights = weights
    )
    rho_est = rho_delta_method(moment_estimate$coefs, moment_estimate$vcov)
    
    return(rho_est$rho)
  }
  return(statistic_f)
}



# Compute multiplier bootstrap percentile confidence intervals.
meta_analytic_data_simulated =  bind_rows(
  meta_analytic_data_simulated %>%
    mutate(CI_type = "sandwich"),
  meta_analytic_data_simulated %>%
    filter(sandwich_adjustment == "none") %>%
    mutate(
      # Compute CIs for rho based on the multiplier bootstrap.
      rho_ci_bs = furrr::future_map2(
        .x = treatment_effects_tbl,
        .y = estimator_adjustment,
        .f = function(treatment_effects_tbl, estimator_adjustment) {
          multiplier_bootstrap_ci(
            data = treatment_effects_tbl,
            statistic = statistic_function_factory(estimator_adjustment),
            B = B,
            alpha = 0.05
          )
        }
      ),
      rho_ci_lower = purrr::map_dbl(rho_ci_bs, function(x)
        x[[1]]),
      rho_ci_upper = purrr::map_dbl(rho_ci_bs, function(x)
        x[[2]])
    ) %>%
    mutate(CI_type = "multiplier bootstrap") %>%
    # Drop redundant variables.
    select(-rho_ci_bs)
)


# Save Results -------------------------------------------------------------

# Save the simulated MA data sets with the corresponding rho estimates and
# confidence intervals. 
saveRDS(meta_analytic_data_simulated, "results/raw-results/simple-simulation/meta_analytic_data_simulated.rds")
 
print(Sys.time() - a)

  

