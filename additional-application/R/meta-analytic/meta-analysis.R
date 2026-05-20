# Setup ------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(scales)
library(splines)
library(future)
library(furrr)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}

# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

ma_trt_effects_tbl_location = "additional-application/results/raw-results/ma_trt_effects_tbl.rds"

# Specify options for saving the plots to files
figures_dir = "additional-application/results/figures"
tables_dir = "additional-application/results/tables"


## Analysis Parameters -------------------------------------------------- 

set.seed(1)
# Number of bootstrap replications for the multiplier bootstrap for the
# meta-analytic parameters.
B_multiplier = 1e5


## Intermediate Results  --------------------------------------------------

# Load data with trial-level treatment effects. 
ma_trt_effects_tbl = readRDS(ma_trt_effects_tbl_location)

ma_trt_effects_tbl = ma_trt_effects_tbl %>%
  rename(treatment_effect_surr = "trt_effect_surrogate_index_est",
         treatment_effect_clin = "trt_effect_clinical_est",
         covariance_matrix = "vcov")

# Duplicate original data where the vcov matrices are divided by 2, to get
# results that better illustrate the methods.
ma_trt_effects_tbl = bind_rows(
  ma_trt_effects_tbl %>%
    mutate(vcov_multiplier = 1),
  ma_trt_effects_tbl %>%
    mutate(
      covariance_matrix = purrr::map(covariance_matrix, ~ .x / 2),
      vcov_multiplier = 0.5
    )
)
  

# Formal Meta-Analysis -----------------------------------------------------

source("R/helper-functions/moment-based-estimator.R")
source("R/helper-functions/multiplier-bootstrap.R")
source("R/helper-functions/delta-method-rho-trial.R")

# Helper function that returns the trial-level correlation estimate given a
# weighted data set.
statistic_f_rho = function(data, weights) {
  moment_estimate = moment_estimator(
    alpha_hat = data$treatment_effect_surr,
    beta_hat = data$treatment_effect_clin,
    vcov_list = data$covariance_matrix,
    estimator_adjustment = "N - 1",
    weights = weights,
    nearest_PD = FALSE
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

statistic_f_residual_var = function(data, weights) {
  moment_estimate = moment_estimator(
    alpha_hat = data$treatment_effect_surr,
    beta_hat = data$treatment_effect_clin,
    vcov_list = data$covariance_matrix,
    estimator_adjustment = "N - 1",
    weights = weights,
    nearest_PD = FALSE,
    SE = FALSE
  )
  # Residual variance
  residual_var = moment_estimate$residual_var
  
  
  return(list(estimate = residual_var, se = NA))
}

statistic_f_residual_var_prop = function(data, weights) {
  moment_estimate = moment_estimator(
    alpha_hat = data$treatment_effect_surr,
    beta_hat = data$treatment_effect_clin,
    vcov_list = data$covariance_matrix,
    estimator_adjustment = "N - 1",
    weights = weights,
    nearest_PD = TRUE,
    SE = FALSE
  )
  # Residual variance
  residual_var = max(moment_estimate$residual_var, 1e-5)
  var_beta = max(moment_estimate$coefs[4], 1e-5)
  
  # Proportion of variance in beta explained by the identity line
  prop_explained = 1 - (residual_var / var_beta)
  
  return(list(estimate = prop_explained, se = NA))
}

# Estimate the surrogacy parameters on each data set of trial-level treatment
# effect estimates.
surrogate_results_tbl = ma_trt_effects_tbl %>%
  group_by(landmark_time, endpoint, method, vcov_multiplier) %>%
  summarise(data_tbl = list(pick(everything())), N = nrow(data_tbl[[1]])) %>%
  ungroup() %>%
  mutate(
    moment_estimate = purrr::map(
      .x = data_tbl,
      .f = function(data_tbl) {
        moment_estimator(
          alpha_hat = data_tbl$treatment_effect_surr,
          beta_hat = data_tbl$treatment_effect_clin,
          vcov_list = data_tbl$covariance_matrix,
          estimator_adjustment = "N - 1",
          sandwich_adjustment = "N - 1",
          nearest_PD = TRUE
        )
      }
    ),
    bootstrap_ci = future_map(
      .x = data_tbl,
      .f = function(data_tbl) {
        multiplier_bootstrap_ci(
          data = data_tbl,
          statistic = statistic_f_rho,
          B = B_multiplier,
          type = "BCa",
          alpha = 0.05
        )
      }
    ),
    bootstrap_ci_residual_var = future_map(
      .x = data_tbl,
      .f = function(data_tbl) {
        multiplier_bootstrap_ci(
          data = data_tbl,
          statistic = statistic_f_residual_var,
          B = B_multiplier,
          type = "BCa",
          alpha = 0.05
        )
      }
    )
  ) %>%
  mutate(rho_sandwich_inference = purrr::map2(
    .x = moment_estimate,
    .y = N,
    .f = function(moment_estimate, N) {
      rho_delta_method(
        coefs = moment_estimate$coefs,
        vcov = moment_estimate$vcov,
        method = "t-adjustment",
        N = N
      )
    }
  ))


surrogate_results_tbl = surrogate_results_tbl %>%
  mutate(
    d_alpha = map_dbl(moment_estimate, function(x)
      x$coefs[3]),
    d_beta = map_dbl(moment_estimate, function(x)
      x$coefs[4]),
    d_alphabeta = map_dbl(moment_estimate, function(x)
      x$coefs[5]),
    rho_trial = d_alphabeta / sqrt(d_alpha * d_beta),
    residual_var = map_dbl(moment_estimate, function(x)
      x$residual_var),
    CI_lower_bs = purrr::map_dbl(bootstrap_ci, "ci_lower"),
    CI_upper_bs = purrr::map_dbl(bootstrap_ci, "ci_upper"),
    CI_lower_bs_residual_var = purrr::map_dbl(bootstrap_ci_residual_var, "ci_lower"),
    CI_upper_bs_residual_var = purrr::map_dbl(bootstrap_ci_residual_var, "ci_upper"),
    # CI_lower_bs_residual_var_prop = purrr::map_dbl(bootstrap_ci_residual_var_prop, "ci_lower"),
    # CI_upper_bs_residual_var_prop = purrr::map_dbl(bootstrap_ci_residual_var_prop, "ci_upper"),
    CI_lower_sandwich = map_dbl(rho_sandwich_inference, function(x)
      x$ci[[1]]),
    CI_upper_sandwich = map_dbl(rho_sandwich_inference, function(x)
      x$ci[[2]]),
    rho_se = map_dbl(rho_sandwich_inference, function(x)
      x$se)
  )

# Summarize inferences in a table.
surrogate_results_tbl = surrogate_results_tbl %>%
  select(-moment_estimate, -bootstrap_ci, -bootstrap_ci_residual_var, -rho_sandwich_inference, -data_tbl) 

surrogate_results_tbl %>%
  write.csv(file = paste0(tables_dir, "/surrogacy-inferences.csv"))
