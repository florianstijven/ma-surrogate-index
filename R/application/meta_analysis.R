# Setup ------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(scales)
library(splines)
library(WeightedROC)
library(mgcv)
library(future)
library(furrr)
library(survival)
library(RColorBrewer)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}

# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

# The first argument indicates whether the analysis should be conducted on the
# original data or on the synthetic data.
data_set = args[1]
if (data_set == "real") {
  ma_trt_effects_tbl_location = "results/raw-results/application/ma_trt_effects_tbl.rds"
  
  # Specify options for saving the plots to files
  figures_dir = "results/figures/application/meta-analysis"
  tables_dir = "results/tables/application/meta-analysis"
} else if (data_set == "synthetic") {
  ma_trt_effects_tbl_location = "results/raw-results/application-synthetic/ma_trt_effects_tbl.rds"
  
  # Specify options for saving the plots to files
  figures_dir = "results/figures/application-synthetic/meta-analysis"
  tables_dir = "results/tables/application-synthetic/meta-analysis"
}

## Analysis Parameters -------------------------------------------------- 

set.seed(1)
# Number of bootstrap replications for the multiplier bootstrap for the
# meta-analytic parameters.
B_multiplier = 1e5

time_cumulative_incidence = 80

# Set of trials corresponding to each analysis set.

trials_first_four = c("Moderna", "AstraZeneca", "Janssen", "Novavax")
trials_first_four_fct = trials_first_four
trials_naive_only = c("Moderna",
                      "AstraZeneca",
                      "Janssen",
                      "Novavax",
                      "Sanofi-1 naive",
                      "Sanofi-2 naive")
trials_naive_only_fct = c("Moderna",
                          "AstraZeneca",
                          "Janssen",
                          "Novavax",
                          "Sanofi 1 (naive)",
                          "Sanofi 2 (naive)")
trials_mixed = c(
  "Moderna",
  "AstraZeneca",
  "Janssen",
  "Novavax",
  "Sanofi-1 naive",
  "Sanofi-1 non-naive",
  "Sanofi-2 naive",
  "Sanofi-2 non-naive"
)
trials_mixed_fct = c(
  "Moderna",
  "AstraZeneca",
  "Janssen",
  "Novavax",
  "Sanofi 1 (naive)",
  "Sanofi 1 (non-naive)",
  "Sanofi 2 (naive)",
  "Sanofi 2 (non-naive)"
)


## Intermediate Results  --------------------------------------------------

# Load data with trial-level treatment effects. 
ma_trt_effects_tbl = readRDS("results/raw-results/application/ma_trt_effects_tbl.rds")

# Add proper name of the surrogates.
ma_trt_effects_tbl = ma_trt_effects_tbl %>%
  mutate(surrogate_name = case_when(
    surrogate == "bindSpike" ~ "IgG Spike",
    surrogate == "pseudoneutid50" ~ "nAb ID50",
    surrogate == "pseudoneutid50_adjusted" ~ "adjusted nAb ID50"
  ))

# Add indicator for whether only naive individuals are in a given trial. 
ma_trt_effects_tbl = ma_trt_effects_tbl %>%
  mutate(
    only_naive = !(trial %in% c("Sanofi-1 non-naive", "Sanofi-2 non-naive"))
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

# Construct data set for the treatment effect estimates on the untransformed
# surrogates. We duplicate rows such that we have a set of rows corresponding to
# each analysis set.
ma_trt_effects_tbl_untransformed = bind_rows(
  ma_trt_effects_tbl %>%
    filter(method == "untransformed surrogate") %>%
    filter(trial %in% trials_naive_only) %>%
    mutate(analysis_set = "naive_only")
  # ma_trt_effects_tbl %>%
  #   filter(method == "untransformed surrogate") %>%
  #   filter(trial %in% trials_mixed) %>%
  #   mutate(analysis_set = "mixed")
)

# Construct data set with a set of rows corresponding to each analysis set.
ma_trt_effects_tbl_modified = ma_trt_effects_tbl %>%
  filter(method != "untransformed surrogate") %>%
  bind_rows(ma_trt_effects_tbl_untransformed) %>%
  # Make sure that only the correct trials are included for each analysis set.
  filter(ifelse(
    analysis_set == "first_four",
    trial %in% trials_first_four,
    ifelse(analysis_set == "naive_only", trial %in% trials_naive_only, TRUE)
  ))

# Add pseudo-real data by (i) cloning each trial four times and (ii) dividing
# the within-trial variance matrices by four.
ma_trt_effects_tbl_modified = bind_rows(
  ma_trt_effects_tbl_modified %>%
    mutate(scenario = "real data"),
  bind_rows(
    ma_trt_effects_tbl_modified,
    ma_trt_effects_tbl_modified,
    ma_trt_effects_tbl_modified,
    ma_trt_effects_tbl_modified
  ) %>%
    mutate(scenario = "four clones", analysis_set == "naive_only"),
  ma_trt_effects_tbl_modified %>%
    filter(analysis_set == "naive_only") %>%
    mutate(
      vcov = purrr::map(
        .x = vcov,
        .f = function(x) {
          x / 4
        }
      ),
      scenario = "precise trials"
    )
)


# Estimate the surrogacy parameters on each data set of trial-level treatment
# effect estimates.
surrogate_results_tbl = ma_trt_effects_tbl_modified %>%
  group_by(surrogate, method, weighting, analysis_set, include_risk_score, scenario) %>%
  summarise(data_tbl = list(pick(everything())), N = nrow(data_tbl[[1]])) %>%
  ungroup() %>%
  mutate(
    moment_estimate = purrr::map(
      .x = data_tbl,
      .f = function(data_tbl) {
        moment_estimator(
          alpha_hat = data_tbl$trt_effect_surrogate_index_est,
          beta_hat = data_tbl$log_RR_est,
          vcov_list = data_tbl$vcov,
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
          data = data_tbl %>%
            rename(
              treatment_effect_surr = "trt_effect_surrogate_index_est",
              treatment_effect_clin = 'log_RR_est',
              covariance_matrix = "vcov"
            ),
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
          data = data_tbl %>%
            rename(
              treatment_effect_surr = "trt_effect_surrogate_index_est",
              treatment_effect_clin = 'log_RR_est',
              covariance_matrix = "vcov"
            ),
          statistic = statistic_f_residual_var,
          B = B_multiplier,
          type = "BCa",
          alpha = 0.05
        )
      }
    ),
    # bootstrap_ci_residual_var_prop = future_map(
    #   .x = data_tbl,
    #   .f = function(data_tbl) {
    #     multiplier_bootstrap_ci(
    #       data = data_tbl %>%
    #         rename(
    #           treatment_effect_surr = "trt_effect_surrogate_index_est",
    #           treatment_effect_clin = 'log_RR_est',
    #           covariance_matrix = "vcov"
    #         ),
    #       statistic = statistic_f_residual_var_prop,
    #       B = B_multiplier,
    #       type = "BCa",
    #       alpha = 0.05
    #     )
    #   }
    # )
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

