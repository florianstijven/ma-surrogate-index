# Load libraries
library(tidyverse)
library(splines)
library(WeightedROC)
library(mgcv)
library(hal9001)
library(future)
library(furrr)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}

## Analysis Parameters -------------------------------------------------- 

# Number of knots for HAL. 
num_knots = c(25, 10)
# Number of bootstrap replications for computing within-trial covariance
# matrices.
B_within_trial = 5e2
# Number of bootstrap replications for the multiplier bootstrap for the
# meta-analytic parameters.
B_multiplier = 5e3

## Data Preparation ----------------------------------------------------

# Load data. We immediately start working with tibbles instead of classic data
# frames.
ipd_tbl = read.csv("data/processed_data.csv") %>%
  tibble() %>%
  # The clinical outcome Y should be an integer.
  mutate(Y = as.integer(Y)) %>%
  # Code the BMI dummy variables into a single factor variable.
  mutate(
    CC_stratum = as.factor(CC_stratum),
    BMI_stratum = case_when(
      BMI_underweight == 1 ~ "Underweight",
      BMI_normal == 1 ~ "Normal",
      BMI_overweight == 1 ~ "Overweight",
      BMI_obese == 1 ~ "Obese",
      .default = NA
    ),
    BMI_stratum = as.factor(BMI_stratum)
  )

# Recode all J&J subunits into a single trial.
ipd_tbl = ipd_tbl %>%
  mutate(
    trial.lbl = ifelse(str_detect(trial.lbl, "J&J"), "J&J", trial.lbl)
  )

# Estimate case-cohort sampling weights and add these to the data set.
ipd_tbl = ipd_tbl %>%
  left_join(ipd_tbl %>%
              group_by(CC_stratum) %>%
              summarize(weight = 1 / mean(Delta)) %>%
              ungroup())

# The weights should be set to one for placebo patients. Their titers were set
# to the lower limits of the assays.
ipd_tbl = ipd_tbl %>%
  mutate(weight = ifelse(vax == 0, 1, weight))
 
# Set weights to zero for patients in the vaccine groups for whom the titers
# were not measured.
ipd_tbl = ipd_tbl %>%
  mutate(weight = ifelse((vax == 1) & (Delta == 0), 0, weight))

# Make sure that each trial-treatment group combination gets the same total
# weight. We don't want one trial to have an undue influence on the predictions.
ipd_tbl = ipd_tbl %>%
  left_join(ipd_tbl %>%
              group_by(trial.lbl, vax) %>%
              summarize(total_weight = sum(weight)) %>%
              ungroup()) %>%
  mutate(weight = weight / total_weight)

# Ensure that the weights sum to the total number of observations. 
ipd_tbl = ipd_tbl %>%
  mutate(weight = weight / mean(weight))

# Every trial now receives the same total weight.
ipd_tbl %>%
  group_by(trial.lbl, vax) %>%
  summarise(sum(weight))

# Add variable that indicates whether the titers were at or below the limit of
# detection.
lod_bindSpike = ipd_tbl$bindSpike[ipd_tbl$vax == 0][1]
lod_pseudoneutid50 = ipd_tbl$pseudoneutid50[ipd_tbl$vax == 0][1]
ipd_tbl = ipd_tbl %>%
  mutate(
    detected_bindSpike = ifelse(bindSpike == lod_bindSpike, "Not Detected", "Detected"),
    detected_pseudoneut50 = ifelse(
      pseudoneutid50 == lod_pseudoneutid50,
      "Not Detected",
      "Detected"
    )
  )


# Prediction Models -------------------------------------------------------

## Parametric Prediction Model --------------------------------------------

# Estimate logistic regression model for the probability of infection given the
# baseline covariates and the surrogate. Note that we missing predictors: the
# antibody titers are missing according to the case-cohort sampling mechanism. 
glm_fit_baseline = glm(
  formula = Y ~ bs(risk_score) + BMI_stratum + bs(Age),
  data = ipd_tbl, 
  weights = weight,
  subset = (ipd_tbl$Delta == 1) | (ipd_tbl$vax == 0),
  family = quasibinomial()
)
glm_fit_spike = glm(
  formula = Y ~ detected_bindSpike + bs(bindSpike) +
    bs(risk_score) + BMI_stratum + bs(Age),
  data = ipd_tbl, 
  weights = weight,
  subset = (ipd_tbl$Delta == 1) | (ipd_tbl$vax == 0),
  family = quasibinomial()
)
glm_fit_neut = glm(
  formula = Y ~ detected_pseudoneut50 + bs(pseudoneutid50) +
    bs(risk_score) + BMI_stratum + bs(Age),
  data = ipd_tbl, 
  weights = weight,
  subset = (ipd_tbl$Delta == 1) | (ipd_tbl$vax == 0),
  family = quasibinomial()
)
glm_fit_both = glm(
  formula = Y ~ detected_bindSpike + bs(bindSpike) +
    detected_pseudoneut50 + bs(pseudoneutid50) +
    bs(risk_score) + BMI_stratum + bs(Age),
  data = ipd_tbl, 
  weights = weight,
  subset = (ipd_tbl$Delta == 1) | (ipd_tbl$vax == 0),
  family = quasibinomial()
)

surrogate_index_models_tbl = tibble(
  fitted_model = list(glm_fit_spike, glm_fit_neut, glm_fit_both, glm_fit_baseline),
  method = rep("glm", 4),
  surrogate = c("bindSpike", "pseudoneutid50", "bindSpike + pseudoneutid50", "none")
)


## GAM --------------------------------------------------------------------

# Estimate logistic GAM for the probability of infection given the
# baseline covariates and the surrogate. Note that we missing predictors: the
# antibody titers are missing according to the case-cohort sampling mechanism. 
gam_fit_baseline = gam(
  formula = Y ~ s(risk_score) + BMI_stratum + s(Age),
  data = ipd_tbl, 
  weights = weight,
  subset = (ipd_tbl$Delta == 1) | (ipd_tbl$vax == 0),
  family = quasibinomial()
)
gam_fit_spike = gam(
  formula = Y ~ detected_bindSpike + s(bindSpike) +
    s(risk_score) + BMI_stratum + s(Age),
  data = ipd_tbl, 
  weights = weight,
  subset = (ipd_tbl$Delta == 1) | (ipd_tbl$vax == 0),
  family = quasibinomial()
)
gam_fit_neut = gam(
  formula = Y ~ detected_pseudoneut50 + s(pseudoneutid50) +
    s(risk_score) + BMI_stratum + s(Age),
  data = ipd_tbl, 
  weights = weight,
  subset = (ipd_tbl$Delta == 1) | (ipd_tbl$vax == 0),
  family = quasibinomial()
)
gam_fit_both = gam(
  formula = Y ~ detected_bindSpike + s(bindSpike) +
    detected_pseudoneut50 + s(pseudoneutid50) +
    s(risk_score) + BMI_stratum + s(Age),
  data = ipd_tbl, 
  weights = weight,
  subset = (ipd_tbl$Delta == 1) | (ipd_tbl$vax == 0),
  family = quasibinomial()
)


surrogate_index_models_tbl = surrogate_index_models_tbl %>%
  bind_rows(tibble(
    fitted_model = list(gam_fit_spike, gam_fit_neut, gam_fit_both, gam_fit_baseline),
    method = rep("gam", 4),
    surrogate = c("bindSpike", "pseudoneutid50", "bindSpike + pseudoneutid50", "none")
  ))

## HAL ------------------------------------------------------------------


glm_baseline = glm(
  Y ~ risk_score + Age,
  ipd_tbl,
  family = gaussian(),
  weights = weight,
  subset = (ipd_tbl$Delta == 1) |
    (ipd_tbl$vax == 0)
)
hal_fit_baseline = fit_hal(
  X = model.matrix(glm_baseline)[, -1],
  Y = glm_baseline$model$Y,
  family = "binomial",
  max_degree = 2,
  num_knots = num_knots,
  weights = glm_baseline$prior.weights
)

glm_spike = glm(
  Y ~ detected_bindSpike + bindSpike + risk_score + Age,
  ipd_tbl,
  family = gaussian(),
  weights = weight,
  subset = (ipd_tbl$Delta == 1) |
    (ipd_tbl$vax == 0)
)
hal_fit_spike = fit_hal(
  X = model.matrix(glm_spike)[, -1],
  Y = glm_spike$model$Y,
  family = "binomial",
  max_degree = 2,
  num_knots = num_knots,
  weights = glm_spike$prior.weights
)

glm_neut = glm(
  Y ~ detected_pseudoneut50 + pseudoneutid50 + risk_score + Age,
  ipd_tbl,
  family = gaussian(),
  weights = weight,
  subset = (ipd_tbl$Delta == 1) |
    (ipd_tbl$vax == 0),
  
)
hal_fit_neut = fit_hal(
  X = model.matrix(glm_neut)[, -1],
  Y = glm_neut$model$Y,
  family = "binomial",
  max_degree = 2,
  num_knots = num_knots,
  weights = glm_neut$prior.weights
)

glm_both = glm(
  Y ~ detected_pseudoneut50 + pseudoneutid50 + detected_bindSpike + bindSpike + risk_score + Age,
  ipd_tbl,
  family = gaussian(),
  weights = weight,
  subset = (ipd_tbl$Delta == 1) |
    (ipd_tbl$vax == 0)
)
hal_fit_both = fit_hal(
  X = model.matrix(glm_both)[, -1],
  Y = glm_both$model$Y,
  family = "binomial",
  max_degree = 2,
  num_knots = num_knots,
  weights = glm_both$prior.weights
)

surrogate_index_models_tbl = surrogate_index_models_tbl %>%
  bind_rows(tibble(
    fitted_model = list(hal_fit_spike, hal_fit_neut, hal_fit_both, hal_fit_baseline),
    method = rep("hal", 4),
    surrogate = c(
      "bindSpike",
      "pseudoneutid50",
      "bindSpike + pseudoneutid50",
      "none"
    )
  ))

## Prediction Accuracies of all Prediction Models --------------------------

# Add predictions to a new data set. We first add two artificial data sets where
# the surrogate is the original surrogate, either bindSpike or pseudoneutid50.
ipd_surr_indices_tbl = bind_rows(
  ipd_tbl %>%
    mutate(method = "untransformed surrogate", surrogate = "pseudoneutid50") %>%
    rename(surrogate_index = "pseudoneutid50"),
  ipd_tbl %>%
    mutate(method = "untransformed surrogate", surrogate = "bindSpike") %>%
    rename(surrogate_index = "bindSpike")
)
ipd_surr_indices_tbl = bind_rows(
  ipd_surr_indices_tbl,
  surrogate_index_models_tbl %>%
    rowwise(method, surrogate) %>%
    filter(method != "hal") %>%
    reframe(tibble(
      surrogate_index = predict(fitted_model, newdata = ipd_tbl, type = "response"),
      ipd_tbl
    ))
)

# For doing predictions based on HAL, we have to reconstruct design matrices for
# the entire data set. Missing predictor values are replaced by zeroes. The
# corresponding predictions have to be replaced with NAs in a second step. 
ipd_tbl_imputed = ipd_tbl %>%
  mutate(detected_bindSpike = ifelse(is.na(detected_bindSpike), "Detected", detected_bindSpike),
         detected_pseudoneut50 = ifelse(is.na(detected_pseudoneut50), "Detected", detected_pseudoneut50),
         bindSpike = ifelse(is.na(bindSpike), 0, bindSpike),
         pseudoneutid50 = ifelse(is.na(pseudoneutid50), 0, pseudoneutid50))
glm_baseline = glm(Y ~ risk_score + Age, ipd_tbl_imputed, family = gaussian())
glm_spike = glm(
  Y ~ detected_bindSpike + bindSpike + risk_score + Age,
  ipd_tbl_imputed,
  family = gaussian()
)
glm_neut = glm(
  Y ~ detected_pseudoneut50 + pseudoneutid50 + risk_score + Age,
  ipd_tbl_imputed,
  family = gaussian()
)
glm_both = glm(
  Y ~ detected_pseudoneut50 + pseudoneutid50 + detected_bindSpike + bindSpike + risk_score + Age,
  ipd_tbl_imputed,
  family = gaussian()
)

ipd_surr_indices_tbl = bind_rows(
  ipd_surr_indices_tbl,
  surrogate_index_models_tbl %>%
    rowwise(method, surrogate) %>%
    filter(method == "hal", surrogate == "none") %>%
    reframe(tibble(
      surrogate_index = predict(fitted_model, new_data = model.matrix(glm_baseline)[, -1], type = "response"),
      ipd_tbl
    )),
  surrogate_index_models_tbl %>%
    rowwise(method, surrogate) %>%
    filter(method == "hal", surrogate == "bindSpike") %>%
    reframe(tibble(
      surrogate_index = predict(fitted_model, new_data = model.matrix(glm_spike)[, -1], type = "response"),
      ipd_tbl
    )),
  surrogate_index_models_tbl %>%
    rowwise(method, surrogate) %>%
    filter(method == "hal", surrogate == "pseudoneutid50") %>%
    reframe(tibble(
      surrogate_index = predict(fitted_model, new_data = model.matrix(glm_neut)[, -1], type = "response"),
      ipd_tbl
    )),  
  surrogate_index_models_tbl %>%
    rowwise(method, surrogate) %>%
    filter(method == "hal", surrogate == "bindSpike + pseudoneutid50") %>%
    reframe(tibble(
      surrogate_index = predict(fitted_model, new_data = model.matrix(glm_both)[, -1], type = "response"),
      ipd_tbl
    ))
)

# Replace predictions based on imputed predictors with NAs.
ipd_surr_indices_tbl = bind_rows(
  ipd_surr_indices_tbl %>% filter(surrogate == "none"),
  ipd_surr_indices_tbl %>%
    filter(surrogate != "none") %>%
    mutate(surrogate_index = ifelse((Delta == 0) &
                                      (vax == 1), NA, surrogate_index))
)


# Plot ROC and compute AUC to assess predictive performance of the estimated
# logistic regression model.

roc_tbl = ipd_surr_indices_tbl %>%
  group_by(method, surrogate) %>%
  filter(!is.na(surrogate_index), weight != 0) %>%
  reframe(WeightedROC(
    guess = surrogate_index,
    label = Y, 
    weight = weight
  ))

roc_tbl %>%
  ggplot(aes(color = surrogate, linetype = method)) +
  geom_path(aes(FPR, TPR)) +
  coord_equal()

ggsave(filename = "results/figures/application/roc-surrogate-indices.pdf", 
       height = double_height,
       width = double_width,
       device = "pdf",
       units = "cm")

roc_tbl %>% group_by(method, surrogate) %>%
summarise(AUC = WeightedAUC(pick(c(TPR, FPR, FP, FN, threshold)))) %>%
  pivot_wider(names_from = "method", values_from = "AUC") %>%
  write.csv("results/tables/application/auc-surrogate-indices.csv")

# Meta-Analysis -----------------------------------------------------------

## Helper Functions -------------------------------------------------------

estimate_treatment_effect_surrogate_index = function(data, VE_surr) {
  # Fit linear regression model adjusted for baseline covariates and weighted by
  # the case-cohort-sampling weights.
  lm_fit = lm(
    formula = surrogate_index ~ vax + vax*bs(risk_score) + BMI_stratum + bs(Age),
    data = data,
    subset = (data$Delta == 1) | (data$vax == 0),
    weights = weight
  )
  
  # Predict outcome probabilities given the observed baseline covariates
  # combining both treatment groups and setting the treatment to placebo and
  # vaccine. We quantify the treatment effect on the same scale as VE.
  trt_effect_tbl = data %>%
    mutate(
      predicted_prob_placebo = predict(lm_fit, newdata = pick(everything()) %>%
                                         mutate(vax = 0)),
      predicted_prob_vax = predict(lm_fit, newdata = pick(everything()) %>%
                                     mutate(vax = 1))
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

estimate_treatment_effect_clinical = function(data) {
  # Fit logistic regression model adjusted for baseline covariates.
  glm_fit = glm(
    formula = Y ~ vax + vax*risk_score + BMI_stratum + Age,
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
          mutate(vax = 0),
        type = "response"
      ),
      predicted_prob_vax = predict(
        glm_fit,
        newdata = pick(everything()) %>%
          mutate(vax = 1),
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
  VE_est = estimate_treatment_effect_clinical(data)
  
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
        # Re-estimate the case-cohort sampling weights.
        data_bs = data_bs %>%
          select(-weight) %>%
          left_join(
            data_bs %>%
              group_by(CC_stratum) %>%
              summarize(weight = 1 / mean(Delta)) %>%
              ungroup(),
            by = "CC_stratum"
          )
        
        # The weights should be set to one for placebo patients. Their titers were set
        # to the lower limits of the assays.
        data_bs = data_bs %>%
          mutate(weight = ifelse(vax == 0, 1, weight))
        
        # Set weights to zero for patients in the vaccine groups for whom the titers
        # were not measured.
        data_bs = data_bs %>%
          mutate(weight = ifelse((vax == 1) &
                                   (Delta == 0), 0, weight))
        
        # Estimate treatment effects on surrogate index and clinical endpoint.
        trt_effect_surrogate_index_est = estimate_treatment_effect_surrogate_index(data_bs, VE_surr)
        VE_est = estimate_treatment_effect_clinical(data_bs)
        
        if (log_RR_alpha) {
          trt_effect_surrogate_index_est = log(1 - trt_effect_surrogate_index_est)
        }
        if (log_RR_beta) {
          VE_est = log(1 - VE_est)
        }
        
        return(c(trt_effect_surrogate_index_est, VE_est))
      }
    )
    vcov_est = var(as.matrix(do.call(rbind, estimates_list)))
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
  group_by(method, surrogate, trial.lbl) %>%
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
  ungroup()
  

## Summarize Trial-level Results ------------------------------------------

# Summarize the estimated trial-level treatment effects in meta-analytic plots.
ma_trt_effects_tbl %>% 
  filter(method == "untransformed surrogate") %>%
  ggplot(aes(x = trt_effect_surrogate_index_est, y = -1 * log_RR_est, color = trial.lbl)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = -1 * (log_RR_est - 1.96 * log_RR_est_se),
    ymax = -1 * (log_RR_est + 1.96 * log_RR_est_se)
  ),
  width = 0.01) +
  geom_errorbarh(
    aes(
      xmin = trt_effect_surrogate_index_est - 1.96 * trt_effect_surrogate_index_est_se,
      xmax = trt_effect_surrogate_index_est + 1.96 * trt_effect_surrogate_index_est_se
    ),
    height = 0.01
  ) +
  facet_grid(~surrogate) +
  xlab("Estimated Treatment Effect on Titer") +
  ylab("Estimated -log(1 - VE)")

ggsave(
  filename = "results/figures/application/ma-standard.pdf",
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

ma_trt_effects_tbl %>% filter(method != "untransformed surrogate") %>%
  ggplot(aes(x = -1 * trt_effect_surrogate_index_est, y = -1 * log_RR_est, color = trial.lbl)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = -1 * (log_RR_est - 1.96 * log_RR_est_se),
    ymax = -1 * (log_RR_est + 1.96 * log_RR_est_se)
  ),
  width = 0.01) +
  geom_errorbarh(
    aes(
      xmin = -1 * (trt_effect_surrogate_index_est - 1.96 * trt_effect_surrogate_index_est_se),
      xmax = -1 * (trt_effect_surrogate_index_est + 1.96 * trt_effect_surrogate_index_est_se)
    ),
    height = 0.01
  ) +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(method~surrogate) +
  xlab("Estimated -log(1 - VE) on Surr. Index") +
  ylab("Estimated -log(1 - VE)")

ggsave(
  filename = "results/figures/application/ma-surrogate-indices.pdf",
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

## Estimate Surrogacy Parameters ------------------------------------------

source("R/helper-functions/moment-based-estimator.R")
source("R/helper-functions/multiplier-bootstrap.R")
source("R/helper-functions/delta-method-rho-trial.R")

# Helper function that returns the trial-level correlation estimate given a
# weighted data set.
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

# Estimate the surrogacy parameters on each data set of trial-level treatment
# effect estimates.
surrogate_results_tbl = ma_trt_effects_tbl %>%
  group_by(surrogate, method) %>%
  summarise(
    moment_estimate = list(
      moment_estimator(
        alpha_hat = trt_effect_surrogate_index_est,
        beta_hat = log_RR_est,
        vcov_list = vcov,
        estimator_adjustment = "N - 1",
        sandwich_adjustment = "N - 1"
      )
    ),
    bootstrap_ci = list(
      multiplier_bootstrap_ci(
        data = pick(everything()) %>%
          rename(trt_surr = "trt_effect_surrogate_index_est", trt_clin = 'log_RR_est'),
        statistic = statistic_f,
        B = B_multiplier,
        type = "BCa"
      )
    )
  ) %>%
  ungroup() 

surrogate_results_tbl = surrogate_results_tbl %>%
  mutate(
    d_alpha = map_dbl(moment_estimate, function(x) x$coefs[3]),
    d_beta = map_dbl(moment_estimate, function(x) x$coefs[4]),
    d_alphabeta = map_dbl(moment_estimate, function(x) x$coefs[5]),
    rho_trial = d_alphabeta / sqrt(d_alpha * d_beta),
    residual_var = map_dbl(moment_estimate, function(x) x$residual_var)
  )


  


# Summarize inferences in a table.
surrogate_results_tbl %>%
  mutate(
    CI_lower_bs = purrr::map_dbl(bootstrap_ci, "ci_lower"),
    CI_upper_bs = purrr::map_dbl(bootstrap_ci, "ci_upper")
  ) %>%
  select(-moment_estimate, -bootstrap_ci) %>%
  write.csv(file = "results/tables/application/surrogacy-inferences.csv")








