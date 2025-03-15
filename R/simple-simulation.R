# Preliminaries -----------------------------------------------------------

a = Sys.time()
# Load all R packages
library(tidyverse)
library(sandwich)
library(lmtest)
library(future)
library(furrr)
library(mgcv)
library(ranger)
library(splines)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}
# Extract arguments for analysis. 
args = commandArgs(trailingOnly=TRUE)

# Under which scenario should the simulations be conducted? 
scenario = args[1]
# Are the simulations small or large sample? Large sample is just to illustrate
# problems under some asymptotic regimes.
regime = args[2]
# Number of Monte Carlo simulations.
N_MC = as.numeric(args[3])

# This script simulates data. We set a seed to ensure reproducibility.
set.seed(123) 

# Number of bootstrap replications for the multiplier bootstrap
B = 5e2

 
# Set within-trial sample size depending on the scenario
if (scenario == "proof-of-concept") {
  n = 2e3
  N = c(6, 12, 24)  # Number of trials in each meta-analytic data set
  
  # Number of Monte Carlo trial replications for the approximation of the true
  # trial-level correlation.
  N_approximation_MC = 2e4 
  # Number of patients in each trial for the approximation of the true
  # trial-level correlation. This can be small for the proof-of-concept scenario
  # because that is based on differences in sample means, for the standard
  # covariance estimator is unbiased.
  n_approximation_MC = 2e2  
  
} else if (scenario == "vaccine") {
  n = 4e3
  N = c(6, 12)  # Number of trials in each meta-analytic data set
  
  
  # Number of Monte Carlo trial replications for the approximation of the true
  # trial-level correlation.
  N_approximation_MC = 5e3 
  # Number of patients in each trial for the approximation of the true
  # trial-level correlation. This should be large in the vaccine scenario
  # because the standard covariance estimator (based on the delta method) for
  # log RR estimators is consistent but biased. 
  n_approximation_MC = 3e3  
}


if (regime == "small") {
  # sd_beta represents the variance of important trial-level parameters that
  # control the surrogacy and comparability assumptions. A non-zero variance
  # implies some violation of these assumptions.
  if (scenario == "vaccine") {
    sd_beta = list(c(0.10, 0.10), c(0.25, 0.25))
  } else if (scenario == "proof-of-concept") {
    sd_beta = list(c(0.03, 0.03), c(0.1, 0.1))
  }
  SI_violation = c("slight", "moderate")
} else if (regime == "large") {
  if (scenario == "vaccine") {
    stop("The large sample regime is not suitable for the vaccine scenario.")
  }
  # This regime corresponds to the large N small n asymptotic regime. This
  # scenario is artifical and just meant to show when methods can break. 
  N = c(5e2, 1e3, 2e3, 5e3)  # Number of trials in each meta-analytic data set
  # The within-trial sample size is set to something small. 
  n = 40
  
  # The approximation accuracy for the true rho is increased for the large N
  # setting. In this setting, the SD of the estimators will be much smaller; so,
  # the MC error in the true rho is relatively more important. 
  N_approximation_MC = 5e4
  n_approximation_MC = 5e2
  
  sd_beta = list(c(0.1, 0.1))
  SI_violation = c("moderate")
}




# Helper Functions --------------------------------------------------------

# Source helper functions.
source("R/helper-functions/simulation-functions.R")
source("R/helper-functions/moment-based-estimator.R")
source("R/helper-functions/delta-method-rho-trial.R")
source("R/helper-functions/multiplier-bootstrap.R")
source("R/helper-functions/train-clinical-prediction-models.R")

# Formula for simple surrogate index estimator based on linear model.
if (scenario == "proof-of-concept") {
  # Set surrogate index estimators
  surrogate_index_estimator = c("surrogate", "lm")
  # We only look at the original surrogate for the large regime. This saves some
  # computational time.
  if (regime == "large") {
    surrogate_index_estimator = c("surrogate")
  }
} else if (scenario == "vaccine") {
  # Set surrogate index estimators
  surrogate_index_estimator = c("surrogate", "logistic", "gam", "rf")
}


# Simulation Study ---------------------------------------------------------

# Tibble with simulation parameters. Each row corresponds to a data-generating
# mechanism.
dgm_param_tbl = expand_grid(tibble(sd_beta, SI_violation), N, n) %>%
  mutate(surrogate_index_estimators = list(surrogate_index_estimator))

## Simulate IPD data and compute trial-level estimates  --------------------

# Generate N meta-analytic data sets and compute the trial-level Pearson
# correlation for (i) he treatment effects on the surrogate and clinical
# endpoints and (ii) the treatment effects on the surrogate index and clinical
# endpoint.
data_set_indicator = 1:N_MC

meta_analytic_data = expand_grid(data_set_indicator,
                                 dgm_param_tbl) %>%
  mutate(
    list_of_ma_data_objects = future_pmap(
      .l = list(
        sd_beta = sd_beta,
        n = n,
        N = N,
        surrogate_index_estimators = surrogate_index_estimators
      ),
      .f = generate_meta_analytic_data,
      scenario = scenario,
      regime = regime,
      .options = furrr_options(seed = TRUE)),
   ) %>%
  rowwise(everything()) %>%
  reframe(
    tibble(
      treatment_effects = list_of_ma_data_objects[[1]],
      surrogate_index_f_list = list_of_ma_data_objects[[2]],
      surrogate_index_estimator = surrogate_index_estimators
    )
  ) %>%
  ungroup() %>%
  select(-list_of_ma_data_objects, -surrogate_index_estimators)

print(Sys.time() - a)

## Approximate the true trial-level correlation parameter ------------------


# For each surrogate index prediction function, approximate the true trial-level
# correlation through MC approximation.
meta_analytic_data$rho_true = future_pmap_dbl(
  .l = list(
    f = meta_analytic_data$surrogate_index_f_list,
    sd_beta = meta_analytic_data$sd_beta,
    surrogate_index_estimator = meta_analytic_data$surrogate_index_estimator
  ),
  .f = function(f,
                sd_beta,
                surrogate_index_estimator,
                N_approximation_MC,
                n_approximation_MC,
                scenario,
                regime) {
    if (surrogate_index_estimator == "surrogate") {
      # We will not compute the true rho for the surrogate because this is the
      # same true rho across simulations. This value is computed later on.
      return(NA)
    } else {
      return(
        # Compute the true rho given the surrogate index f.
        rho_MC_approximation(
          sd_beta = sd_beta,
          f = f,
          N_approximation_MC = N_approximation_MC,
          n_approximation_MC = n_approximation_MC,
          scenario = scenario,
          regime = regime
        )
      )
    }
  },
  N_approximation_MC = N_approximation_MC,
  n_approximation_MC = n_approximation_MC,
  scenario = scenario,
  regime = regime,
  .options = furrr_options(seed = TRUE, packages = c("ranger", "mgcv"))
)

# Drop surrogate index function as this function is no longer needed.
meta_analytic_data = meta_analytic_data %>%
  select(-surrogate_index_f_list)

print(Sys.time() - a)

# Approximate the true trial-level correlation using the surrogate endpoint.
rho_true_surrogate_tbl = tibble(
  sd_beta, 
  SI_violation
) %>%
  mutate(surrogate_index_estimator = "surrogate") %>%
  mutate(
    rho_true = future_map_dbl(
      .x = sd_beta,
      .f = rho_MC_approximation,
      N_approximation_MC = N_approximation_MC,
      n_approximation_MC = n_approximation_MC,
      f = NULL,
      scenario = scenario,
      regime = regime,
      .options = furrr_options(seed = TRUE)
    )
  )



print(Sys.time() - a)
# Add the true trial-level rho parameters to the tibble with the simulation data
# sets.
meta_analytic_data = meta_analytic_data %>%
  left_join(
    rho_true_surrogate_tbl %>%
      select(-sd_beta),
    by = c("SI_violation", "surrogate_index_estimator"), 
    relationship = "many-to-one"
  ) %>%
  # rho_true is present in both tibbles before joining. Hence, we get rho_true.x
  # and .y in the joined data set. We keep the non-missing value.
  mutate(rho_true = coalesce(rho_true.x, rho_true.y)) %>%
  select(-rho_true.x, -rho_true.y)


print(Sys.time() - a)

## Moment-based estimator --------------------------------------------------

# Estimate trial-level Pearson correlation for the generated meta-analytic data
# sets of treatment effects and add everything to a tibble.
meta_analytic_data_simulated = meta_analytic_data %>%
  cross_join(expand_grid(
    estimator_adjustment = c("N - 1"),
    sandwich_adjustment = c("N - 1")
  )) %>%
  rowwise(everything()) %>%
  summarise(
    # Estimate the meta-analytic parameters using the moment-based estimator.
    moment_estimates = list(
      moment_estimator(
        treatment_effects$treatment_effect_surr,
        treatment_effects$treatment_effect_clin,
        treatment_effects$covariance_matrix,
        estimator_adjustment = estimator_adjustment,
        sandwich_adjustment = sandwich_adjustment
      )
    ),
    # Estimate the trial-level Pearson correlation using the delta method.
    rho_delta_method = list(
      rho_delta_method(
        moment_estimates$coefs,
        moment_estimates$vcov,
        method = "t-adjustment",
        N = N
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

print(Sys.time() - a)

statistic_function_factory = function(estimator_adjustment) {
  statistic_f = function(data, weights) {
    moment_estimate = moment_estimator(
      alpha_hat = data$treatment_effect_surr,
      beta_hat = data$treatment_effect_clin,
      vcov_list = data$covariance_matrix,
      estimator_adjustment = estimator_adjustment,
      weights = weights
    )
    rho = rho_delta_method(
      coefs = moment_estimate$coefs,
      vcov = moment_estimate$vcov,
      method = "t-adjustment",
      # N is only used for the t-adjustment, it doesn't matter for the estimate
      # or SE.
      N = 5
    )
    
    return(list(estimate = rho$rho, se = rho$se))
  }
  return(statistic_f)
}



# Compute multiplier bootstrap percentile confidence intervals. We don't compute
# bootstrap intervals for the large N regime.
if (regime == "small") {
  meta_analytic_data_simulated =  bind_rows(
    meta_analytic_data_simulated %>%
      mutate(CI_type = "sandwich"),
    meta_analytic_data_simulated %>%
      mutate(
        # Compute CIs for rho based on the multiplier bootstrap.
        rho_ci_bs = furrr::future_map2(
          .x = treatment_effects,
          .y = estimator_adjustment,
          .f = function(treatment_effects, estimator_adjustment) {
            multiplier_bootstrap_ci(
              data = treatment_effects,
              statistic = statistic_function_factory(estimator_adjustment),
              B = B,
              alpha = 0.05, 
              type = "BCa"
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
}



# Drop columns that are not needed.
meta_analytic_data_simulated = meta_analytic_data_simulated %>%
  select(-treatment_effects)

print(Sys.time() - a)

# Save Results -------------------------------------------------------------

# Save the simulated MA data sets with the corresponding rho estimates and
# confidence intervals. 
saveRDS(meta_analytic_data_simulated, paste0("results/raw-results/simple-simulation/ma-sim-results-", 
        scenario, "-", regime, ".rds"))
 
print(Sys.time() - a)

  

