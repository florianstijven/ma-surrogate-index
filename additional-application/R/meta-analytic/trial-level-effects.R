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

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}

# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

surr_indices_tbl_location = "additional-application/R/meta-analytic/intermediate-objects/landmark_predictions_tbl.rds"
demographic_tbl_location = "additional-application/R/data/intermediate-objects/full_data_tbl_wide.rds"
out_file = "additional-application/results/raw-results/ma_trt_effects_tbl.rds"

# Specify options for saving the plots to files
figures_dir = "additional-application/results/figures"
tables_dir = "additional-application/results/tables"

# clustering variable used in the analysis. 
clustering_variable_chr = "COU1A"


## Analysis Parameters -------------------------------------------------- 

# Number of bootstrap replications for computing within-trial covariance
# matrices.
B_within_trial = 5e2

time_cumulative_incidence = 2.5 * 365.25

## Intermediate Results ----------------------------------------------------

# Load data set with IPD and estimated surrogate index for every subject.
ipd_surr_indices_tbl = readRDS(file = surr_indices_tbl_location)
# Load data with demographic information, including the clustering variable.
full_data_new_endpoints_landmark_tbl = readRDS(file = demographic_tbl_location)

# Add the clustering information to ``ipd_surr_indices_tbl``.
ipd_surr_indices_tbl = ipd_surr_indices_tbl %>%
  left_join(full_data_new_endpoints_landmark_tbl %>%
              select(SID1A, any_of(clustering_variable_chr)),
            by = "SID1A")

# Change `ipd_surr_indices_tbl` to long format with one estimated surrogate
# index per row.
ipd_surr_indices_tbl = ipd_surr_indices_tbl %>%
  pivot_longer(
    cols = c("predicted_prob_cox", "predicted_prob_sl"),
    names_to = "method",
    values_to = "surrogate_index"
  )

## IPCW -------------------------------------------------------------------

# To estimate the treatment effect on the estimated surrogate index, we need to
# account for missing surrogate indices. Under independent censoring, we can use
# a IPCW mean, where the weights are w_i = 1 / P(C > min(y_i, t)) where C is the censoring
# time, y_i is the observed event time, and t is the landmark time. We can
# estimate the censoring distribution using the Kaplan-Meier estimator, where we
# treat events as censored and censoring times as events.

ipcw_estimator  = function(time_to_event, event, landmark_time) {
  survfit_object = survfit(Surv(time_to_event, 1 - event) ~ 1)
  
  time = unique(pmin(time_to_event, landmark_time))
  surv_probs_tbl = summary(survfit_object, times =  time)[c("surv", "time")] %>%
    as_tibble()
  
  censoring_probs = tibble(time = pmin(time_to_event, landmark_time))  %>%
    left_join(surv_probs_tbl) %>%
    pull(surv)
  
  return(1 / censoring_probs)
}

# Add IPCW weights to the data set.
ipd_surr_indices_tbl = ipd_surr_indices_tbl %>%
  group_by(landmark_time, endpoint, TRTREG1C, COU1A, method) %>%
  mutate(ipcw = ipcw_estimator(time, 1 - censored, landmark_time[1])) %>%
  ungroup()

# Trial-Level Treatment Effects --------------------------------------------

## Helper Functions --------------------------------------------------------


estimate_treatment_effect_surrogate_index = function(data) {
  # Estimate treatment effect on surrogate index using IPCW mean.
  trt_effect_surrogate_index_est = data %>%
    group_by(TRTREG1C) %>%
    summarise(weighted_mean = weighted.mean(surrogate_index, w = ipcw, na.rm = TRUE))
  
  return(trt_effect_surrogate_index_est %>%
           pivot_wider(names_from = TRTREG1C, values_from = weighted_mean) %>%
           mutate(trt_effect_surrogate_index_est = B - A) %>%
           pull(trt_effect_surrogate_index_est))
}



estimate_treatment_effect_clinical_km = function(data) {
  # Fit KM curve stratified by treatment.
  km_A = survfit(Surv(time, 1 - censored) ~ 1,
                       data = data %>%
                         filter(TRTREG1C == "A"))
  km_B = survfit(Surv(time, 1 - censored) ~ 1, data = data %>%
                         filter(TRTREG1C == "B"))
  
  surv_prob_A = summary(km_A, times = time_cumulative_incidence, extend = TRUE)$surv
  surv_prob_B = summary(km_B, times = time_cumulative_incidence, extend = TRUE)$surv
  
  return(surv_prob_B - surv_prob_A)
}



# Function to estimate the treatment effects in a given trial together with a
# covariance matrix based on the non-parametric bootstrap.
estimate_treatment_effects = function(data, B = 2e2) {
  # Estimate treatment effect on surrogate index.
  trt_effect_surrogate_index_est = estimate_treatment_effect_surrogate_index(data)
  
  # Estimate treatment effect on clinical endpoint.
  trt_effect_clinical_est = estimate_treatment_effect_clinical_km(data)
  
  # Perform non-parametric bootstrap if B > 0.
  vcov_est = NA
  if (B > 0) {
    estimates_list = lapply(
      X = 1:B,
      FUN = function(x) {
        # Resample data with replacement.
        data_bs = data %>%
          slice_sample(n = nrow(data), replace = TRUE)
        
        # Estimate treatment effects on surrogate index and clinical endpoint.
        trt_effect_surrogate_index_est = estimate_treatment_effect_surrogate_index(data_bs)
        # Estimate treatment effect on clinical endpoint.
        trt_effect_clinical_est = estimate_treatment_effect_clinical_km(data_bs)
        
        return(c(
          trt_effect_surrogate_index_est,
          trt_effect_clinical_est
        ))
      }
    )
    vcov_est = var(as.matrix(do.call(rbind, estimates_list)), na.rm = TRUE)
  }
  # Return treatment effect estimates.
  return(list(
    estimates = c(
      trt_effect_surrogate_index_est = trt_effect_surrogate_index_est,
      trt_effect_clinical_est = trt_effect_clinical_est
    ),
    vcov = vcov_est
  ))
}

## Estimate Trial-Level Treatment Effects --------------------------------

# Estimate treatment effects on the estimated surrogate indices/surrogates and
# clinical endpoint.
ma_trt_effects_tbl =
  ipd_surr_indices_tbl %>%
  group_by(landmark_time, endpoint, COU1A, method) %>%
  summarize(data = list(pick(everything()))) %>%
  ungroup() %>%
  mutate(
    est_list = future_pmap(
      .l = list(
        data = data
      ),
      .f = estimate_treatment_effects,
      B = B_within_trial,
      .options = furrr_options(seed = TRUE)
    )
  )

ma_trt_effects_tbl = ma_trt_effects_tbl %>%
  rowwise(everything()) %>%
  summarise(
    trt_effect_surrogate_index_est = est_list$estimates["trt_effect_surrogate_index_est"],
    trt_effect_clinical_est = est_list$estimates["trt_effect_clinical_est"],
    vcov = list(est_list$vcov),
    trt_effect_surrogate_index_est_se = sqrt(vcov[1, 1]),
    trt_effect_clinical_est_se = sqrt(vcov[2, 2])
  ) %>%
  ungroup() %>%
  select(-data)

# Saving Results ----------------------------------------------------------

# Save data with trial-level treatment effects to file.
saveRDS(ma_trt_effects_tbl, file = out_file)
