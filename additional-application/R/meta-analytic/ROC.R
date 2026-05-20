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

# Load data set with IPD and estimated surrogate index for every subject.
ipd_surr_indices_tbl = readRDS(file = surr_indices_tbl_location)
# Load data with demographic information, including the clustering variable.
full_data_new_endpoints_landmark_tbl = readRDS(file = demographic_tbl_location)

# Add the clustering information to ``ipd_surr_indices_tbl``.
ipd_surr_indices_tbl = ipd_surr_indices_tbl %>%
  left_join(full_data_new_endpoints_landmark_tbl %>%
              select(SID1A, any_of(clustering_variable_chr), TRTREG1C),
            by = "SID1A")

time_cumulative_incidence = 2.5 * 365.25

## IPCW -------------------------------------------------------------------

# Function to estimate inverse probability of censoring weights. It takes a
# vector of event times and event indicators (1 if event, 0 if censored). It
# returns the IPCWs for each subject.
ipcw_estimator  = function(time_to_event, event) {
  survfit_object = survfit(Surv(time_to_event, event) ~ 1)
  
  time = unique(pmin(
    time_to_event,
    time_cumulative_incidence
  ))
  surv_probs_tbl = summary(survfit_object, times =  time)[c("surv", "time")] %>%
    as_tibble()
  
  censoring_probs = tibble(time = pmin(
    time_to_event,
    time_cumulative_incidence
  ))  %>%
    left_join(surv_probs_tbl) %>%
    pull(surv)
  
  weights = 1 / censoring_probs
  
  if (any(weights == Inf |
          weights == 0, na.rm = TRUE))
    stop("invalid weights computed")
  
  return(weights)
}

# Compute the IPCW for each patient. These weights are based on the KM estimate
# by country and treatment arm. We immediately also do this for each of the new
# composite endpoints separately to take into account the possibility that
# the censoring time is not the same across all endpoints.
ipd_surr_indices_tbl = ipd_surr_indices_tbl %>%
  group_by(COU1A, TRTREG1C, endpoint) %>%
  mutate(
    ipcw_roc = ipcw_estimator(time, censored),
    ipcw_roc = ifelse(time < time_cumulative_incidence &
                    censored == 1, 0, ipcw_roc)
  ) %>%
  ungroup()


## Prediction Model Performance -------------------------------------------
roc_tbl = ipd_surr_indices_tbl %>%
  filter(!is.na(predicted_prob) & !(is.na(censored))) %>%
  filter(ipcw_roc != 0) %>%
  group_by(model, model_type, COU1A, endpoint, landmark_time) %>%
  reframe(
    WeightedROC(
      guess = predicted_prob,
      label = censored,
      weight = ipcw_roc
    )
  )

# Make and save the plots for all combinations of `surrogate` and `weighting`.
roc_tbl %>%       
  ggplot(aes(color = model, linetype = model_type)) +
  geom_path(aes(FPR, TPR)) +
  facet_grid(COU1A ~ landmark_time)

roc_tbl %>%       
  ggplot(aes(color = model, linetype = model_type)) +
  geom_path(aes(FPR, TPR)) +
  facet_wrap(~ landmark_time)
  

roc_tbl %>%
  group_by(COU1A, landmark_time, endpoint, model, model_type) %>%
  summarise(AUC = WeightedAUC(pick(c(
    TPR, FPR, FP, FN, threshold
  )))) %>%
  ggplot(aes(color = model, x = model_type, y = AUC)) +
  geom_point() +
  facet_grid(COU1A ~ landmark_time)

roc_tbl %>%
  group_by(COU1A, landmark_time, endpoint, model, model_type) %>%
  summarise(AUC = WeightedAUC(pick(c(
    TPR, FPR, FP, FN, threshold
  )))) %>%
  ggplot(aes(color = model, shape = model_type, y = AUC, x = landmark_time)) +
  geom_point() +
  facet_wrap(~COU1A)

roc_tbl %>%
  group_by(COU1A, landmark_time, endpoint, model, model_type) %>%
  summarise(AUC = WeightedAUC(pick(c(
    TPR, FPR, FP, FN, threshold
  )))) %>%
  ggplot(aes(color = model, x = model_type, y = AUC, shape = COU1A)) +
  geom_point() +
  facet_wrap(~ landmark_time)
