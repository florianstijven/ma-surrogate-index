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

## Analysis Parameters -------------------------------------------------- 

# Number of bootstrap replications for computing within-trial covariance
# matrices.
B_within_trial = 2e2

time_cumulative_incidence = 80

## Intermediate Results ----------------------------------------------------

# Load data set with IPD and estimated surrogate index for every subject.
ipd_surr_indices_tbl = readRDS(file = "R/application/ipd_surr_indices_tbl.rds")

# Trial-Level Treatment Effects --------------------------------------------

## Helper Functions --------------------------------------------------------

estimate_treatment_effect_surrogate_index = function(data, VE_surr) {
  # Fit linear regression model adjusted for baseline covariates and weighted by
  # the case-cohort-sampling weights.
  lm_fit = lm(
    formula = surrogate_index ~ treatment + treatment * risk_score,
    data = data,
    weights = sample_weight
  )
  # Predict outcome probabilities given the observed baseline covariates
  # combining both treatment groups and setting the treatment to placebo and
  # vaccine. We quantify the treatment effect on the same scale as VE.
  trt_effect_tbl = data %>%
    mutate(
      predicted_prob_placebo = predict(lm_fit, newdata = pick(everything()) %>%
                                         mutate(treatment = 0)),
      predicted_prob_vax = predict(lm_fit, newdata = pick(everything()) %>%
                                     mutate(treatment = 1))
    ) %>%
    summarise(
      mean_prob_placebo = mean(predicted_prob_placebo),
      mean_prob_vax = mean(predicted_prob_vax)
    )
  if (VE_surr) {
    trt_effect_surrogate_index_est = trt_effect_tbl %>%
      mutate(VE = 1 - (mean_prob_vax / mean_prob_placebo)) %>%
      pull(VE)
  } else {
    trt_effect_surrogate_index_est = trt_effect_tbl %>%
      mutate(trt_effect = mean_prob_vax - mean_prob_placebo) %>%
      pull(trt_effect)
  }
  
  return(trt_effect_surrogate_index_est)
}

estimate_treatment_effect_clinical_km = function(data) {
  # Fit KM curve stratified by treatment.
  km_placebo = survfit(Surv(time_to_event, event) ~ 1,
                       data = data %>%
                         filter(treatment == 0))
  km_vaccine = survfit(Surv(time_to_event, event) ~ 1, data = data %>%
                         filter(treatment == 1))
  
  risk_placebo = 1 - summary(km_placebo, times = time_cumulative_incidence, extend = TRUE)$surv
  risk_vaccine = 1 - summary(km_vaccine, times = time_cumulative_incidence, extend = TRUE)$surv
  
  VE_est = 1 - risk_vaccine / risk_placebo
  
  return(VE_est)
}

estimate_treatment_effect_clinical = function(data) {
  # Fit logistic regression model adjusted for baseline covariates.
  glm_fit = glm(
    formula = infection_120d ~ treatment + treatment*risk_score,
    data = data,
    family = binomial()
  )
  
  # Predict outcome probabilities given the observed baseline covariates
  # combining both treatment groups and setting the treatment to placebo and
  # vaccine.
  VE_est = data %>%
    mutate(
      predicted_prob_placebo = predict(
        glm_fit,
        newdata = pick(everything()) %>%
          mutate(treatment = 0),
        type = "response"
      ),
      predicted_prob_vax = predict(
        glm_fit,
        newdata = pick(everything()) %>%
          mutate(treatment = 1),
        type = "response"
      )
    ) %>%
    summarise(
      mean_prob_placebo = mean(predicted_prob_placebo),
      mean_prob_vax = mean(predicted_prob_vax)
    ) %>%
    mutate(VE = 1 - (mean_prob_vax / mean_prob_placebo)) %>%
    pull(VE)
  
  return(VE_est)
}

# Function to estimate the treatment effects in a given trial together with a
# covariance matrix based on the non-parametric bootstrap.
estimate_treatment_effects = function(data,
                                      B = 2e2,
                                      VE_surr = TRUE,
                                      log_RR_alpha = TRUE,
                                      log_RR_beta = TRUE) {
  # Estimate treatment effect on surrogate index.
  trt_effect_surrogate_index_est = estimate_treatment_effect_surrogate_index(data, VE_surr)
  
  # Estimate treatment effect on clinical endpoint.
  VE_est = estimate_treatment_effect_clinical_km(data)
  
  if (log_RR_alpha) {
    trt_effect_surrogate_index_est = log(1 - trt_effect_surrogate_index_est)
  }
  if (log_RR_beta) {
    VE_est = log(1 - VE_est)
  }
  
  # Perform non-parametric bootstrap if B > 0.
  vcov_est = NA
  if (B > 0) {
    estimates_list = lapply(
      X = 1:B,
      FUN = function(x) {
        # Resample data with replacement.
        data_bs = data %>%
          slice_sample(n = nrow(data), replace = TRUE)
        # # Re-estimate the case-cohort sampling weights.
        # data_bs = data_bs %>%
        #   select(-weight) %>%
        #   left_join(
        #     data_bs %>%
        #       group_by(CC_stratum) %>%
        #       summarize(weight = 1 / mean(Delta)) %>%
        #       ungroup(),
        #     by = "CC_stratum"
        #   )
        
        # # The weights should be set to one for placebo patients. Their titers were set
        # # to the lower limits of the assays.
        # data_bs = data_bs %>%
        #   mutate(weight = ifelse(treatment == 0, 1, weight))
        # 
        # # Set weights to zero for patients in the vaccine groups for whom the titers
        # # were not measured.
        # data_bs = data_bs %>%
        #   mutate(weight = ifelse((treatment == 1) &
        #                            (Delta == 0), 0, weight))
        
        # Estimate treatment effects on surrogate index and clinical endpoint.
        trt_effect_surrogate_index_est = estimate_treatment_effect_surrogate_index(data_bs, VE_surr)
        VE_est = estimate_treatment_effect_clinical_km(data_bs)
        
        if (log_RR_alpha) {
          trt_effect_surrogate_index_est = log(1 - trt_effect_surrogate_index_est)
        }
        if (log_RR_beta) {
          VE_est = log(1 - VE_est)
          # If any estimates are not finite, we just set them to NA and ignore NA in
          # the computation of the variance.
          VE_est = ifelse(is.finite(VE_est), VE_est, NA)
        }
        
        return(c(trt_effect_surrogate_index_est, VE_est))
      }
    )
    vcov_est = var(as.matrix(do.call(rbind, estimates_list)), na.rm = TRUE)
  }
  # Return treatment effect estimates.
  return(list(
    estimates = c(trt_effect_surrogate_index_est = trt_effect_surrogate_index_est, VE_est = VE_est),
    vcov = vcov_est
  ))
}

## Estimate Trial-Level Treatment Effects --------------------------------

# Estimate treatment effects on the estimated surrogate indices/surrogates and
# clinical endpoint.
ma_trt_effects_tbl =
  ipd_surr_indices_tbl %>%
  filter(surrogate != "none") %>%
  group_by(method, surrogate, trial, weighting, analysis_set) %>%
  summarize(data = list(pick(everything()))) %>%
  ungroup() %>%
  # Add variables that determine which treatment effects are estimated. For the original surrogates,
  # the treatment effects are mean differences. 
  mutate(log_RR_alpha = ifelse(method == "none", FALSE, TRUE),
         log_RR_beta = TRUE,
         VE_surr = ifelse(method == "none", FALSE, TRUE)) %>%
  mutate(
    est_list = future_pmap(
      .l = list(
        data = data,
        VE_surr = VE_surr,
        log_RR_alpha = log_RR_alpha,
        log_RR_beta = log_RR_beta
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
    log_RR_est = est_list$estimates["VE_est"],
    vcov = list(est_list$vcov),
    trt_effect_surrogate_index_est_se = sqrt(vcov[1, 1]),
    log_RR_est_se = sqrt(vcov[2, 2])
  ) %>%
  ungroup() %>%
  select(-data)

# Saving Results ----------------------------------------------------------

# Save data with trial-level treatment effects to file.
saveRDS(ma_trt_effects_tbl, file = "R/application/ma_trt_effects_tbl.rds")
