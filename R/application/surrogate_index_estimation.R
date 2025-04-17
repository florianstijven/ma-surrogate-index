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

## Analysis Parameters -------------------------------------------------- 

# Number of bootstrap replications for computing within-trial covariance
# matrices.
time_cumulative_incidence = 120

## Data Preparation ----------------------------------------------------


# Load data. We immediately start working with tibbles instead of classic data
# frames.
ipd_tbl = read.csv("data/processed_data.csv") %>%
  tibble() %>%
  # Compute the binary endpoint as observed infection before 120 days. Note that
  # we're ignoring censoring to some extent by doing so because subjects
  # censored before 120 days cannot have the event.
  mutate(infection_120d = ifelse(event == 1 &
                                   time_to_event < time_cumulative_incidence, 1L, 0L)) %>%
  # Code the BMI dummy variables into a single factor variable.
  mutate(
    BMI_stratum = case_when(
      BMI_underweight == 1 ~ "Underweight",
      BMI_normal == 1 ~ "Normal",
      BMI_overweight == 1 ~ "Overweight",
      BMI_obese == 1 ~ "Obese",
      .default = NA
    ),
    BMI_stratum = as.factor(BMI_stratum)
  )

# Center the risk scores by trial. 
ipd_tbl = ipd_tbl %>%
  group_by(trial) %>%
  mutate(risk_score_centered = risk_score - mean(risk_score))

# Split the Sanofi trials according to serological status.
ipd_tbl = ipd_tbl %>%
  mutate(trial = ifelse(
    trial == "Sanofi-1",
    ifelse(Bserostatus == 1, "Sanofi-1 non-naive", "Sanofi-1 naive"),
    ifelse(
      trial == "Sanofi-2",
      ifelse(Bserostatus == 1, "Sanofi-2 non-naive", "Sanofi-2 naive"),
      trial
    )
  ))

# Add an indicator variable for each of the three analysis sets.
ipd_tbl = ipd_tbl %>%
  mutate(
    first_four = trial %in% c("AstraZeneca", "Moderna", "Janssen", "Novavax"),
    naive_only = !(trial %in% c(
      "Sanofi-1 non-naive", "Sanofi-2 non-naive"
    )),
    mixed = TRUE
  )


# Make sure that each trial-treatment group combination gets the same total
# weight. We don't want one trial to have an undue influence on the predictions.
ipd_tbl = ipd_tbl %>%
  left_join(ipd_tbl %>%
              group_by(trial, treatment) %>%
              summarize(total_weight = sum(case_cohort_weight_nAb * !is.na(pseudoneutid50))) %>%
              ungroup()) %>%
  mutate(weight_prediction = case_cohort_weight_nAb / total_weight)

# Every trial now receives the same total weight.
ipd_tbl %>%
  group_by(trial, treatment) %>%
  summarise(sum(weight_prediction * !is.na(pseudoneutid50)))

# Compute pseudo-values. These are computed for each trial separately. We first
# fit a separate KM curve by trial and treatment group.
surv_fit_all = survfit(Surv(time_to_event, event) ~ strata(treatment, trial), data = ipd_tbl)
ipd_tbl$pseudo_value = 1 - pseudo(surv_fit_all, times = time_cumulative_incidence)

# Compute placebo cumulative incidence rate per trial. This will be used later
# on on some plots to standardize some things across trials.
cumulative_incidence_control_tbl = ipd_tbl %>%
  group_by(trial, treatment) %>%
  summarise(prob_infection_free = mean(pseudo_value)) %>%
  ungroup() %>%
  filter(treatment == 0) %>%
  select(-treatment)

## Add placebo cumulative incidence rate as a covariate to the trial. We can
## adjust for this variable later to indirectly adjust for different forces of
## infection across the trials.
ipd_tbl = ipd_tbl %>%
  left_join(cumulative_incidence_control_tbl %>%
              mutate(logit_prob_infection_free = log(
                prob_infection_free / (1 - prob_infection_free)
              )))

# Drop variables that are not needed further on.
ipd_tbl = ipd_tbl %>%
  select(-X,
         -Ptid,
         -USUBJID,
         -(BMI_underweight:abrogation_coefficient),
         -total_weight, 
         -Wstratum, 
         -age.geq.65,
         -Bserostatus, 
         -HighRiskInd)

# Convert Character variables to factors. This is more efficient in terms of
# memory.
ipd_tbl = ipd_tbl %>%
  mutate(trial = as.factor(trial), protocol = as.factor(protocol))


# Prediction Models -------------------------------------------------------

# We first construct a tibble that contains on each row information about a
# particular prediction model. 
prediction_model_settings = tibble(
  surrogate = c("pseudoneutid50", "bindSpike", "pseudoneutid50_adjusted"),
  weights_chr = c(
    "case_cohort_weight_nAb",
    "case_cohort_weight_bAb",
    "case_cohort_weight_nAb"
  ),
  case_cohort_ind_chr = c("Delta_nAb", "Delta_bAb", "Delta_nAb"),
  weighting = rep("unnormalized", 3)
)

# For the current analysies, we only consider the unnormalized weights.
# prediction_model_settings = prediction_model_settings %>%
#   bind_rows(tibble(
#     surrogate = c("pseudoneutid50", "bindSpike", "pseudoneutid50_adjusted"),
#     weights_chr = rep("weight_prediction", 3),
#     case_cohort_ind_chr = c("Delta_nAb", "Delta_bAb", "Delta_nAb"),
#     weighting = rep("Unnormalized", 3)
#   ))

prediction_model_settings = prediction_model_settings %>%
  cross_join(tibble(analysis_set = c("first_four", "naive_only", "mixed")))

## Parametric Prediction Model --------------------------------------------

# Estimate logistic regression model for the probability of infection given the
# baseline covariates and the surrogate. Note that we missing predictors: the
# antibody titers are missing according to the case-cohort sampling mechanism. 


# Function that fits a logistic regression model with the given predictors.
glm_fitter = function(predictors_chr,
                      weights_chr,
                      case_cohort_ind_chr,
                      analysis_set) {
  data_temp = ipd_tbl %>%
    filter(.data[[analysis_set]])
  # Redefine the predictors as smooth functions.
  predictors_chr = paste0("bs(", predictors_chr, ")")
  # Define formula
  string_formula = paste0(
    "infection_120d ~ bs(risk_score_centered) + BMI_stratum + bs(Age) + Sex + logit_prob_infection_free +",
    paste(predictors_chr, collapse = " + ")
  )
  formula_final = as.formula(string_formula)
  
  # Fit logistic regression model.
  glm_fit = glm(
    formula = formula_final,
    data = data_temp,
    weights = data_temp %>%
      pull(any_of(weights_chr)),
    family = quasibinomial(),
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  return(glm_fit)
}

glm_models_tbl = prediction_model_settings %>%
  rowwise(surrogate, weighting, analysis_set) %>%
  summarise(fitted_model = list(glm_fitter(surrogate, weights_chr, case_cohort_ind_chr, analysis_set)))

surrogate_index_models_tbl = glm_models_tbl %>%
  mutate(method = "glm")

rm("glm_models_tbl")


## GAM --------------------------------------------------------------------

# Estimate logistic GAM for the probability of infection given the
# baseline covariates and the surrogate. Note that we missing predictors: the
# antibody titers are missing according to the case-cohort sampling mechanism. 

# Function that fits a logistic regression model with the given predictors.
gam_fitter = function(predictors_chr, weights_chr, case_cohort_ind_chr, analysis_set) {
  data_temp = ipd_tbl %>%
    filter(.data[[analysis_set]])
  # Redefine the predictors as smooth functions.
  predictors_chr = paste0("s(", predictors_chr, ")")
  # Define formula
  string_formula = paste0(
    "infection_120d ~ s(risk_score_centered) + BMI_stratum + s(Age) + Sex + logit_prob_infection_free +",
    paste(predictors_chr, collapse = " + ")
  )
  formula_final = as.formula(string_formula)
  
  # Fit logistic regression model.
  gam_fit = gam(
    formula = formula_final,
    data = data_temp %>% 
      as.data.frame(), 
    weights = data_temp %>% 
      pull(any_of(weights_chr)),
    # subset = (data_temp[, case_cohort_ind_chr] %>% pull(any_of(case_cohort_ind_chr)) == 1) | (data_temp$treatment == 0), 
    family = quasibinomial() 
  )
  return(gam_fit)
}

gam_models_tbl = prediction_model_settings %>%
  rowwise(surrogate, weighting, analysis_set) %>%
  summarise(fitted_model = list(gam_fitter(surrogate, weights_chr, case_cohort_ind_chr, analysis_set)))

surrogate_index_models_tbl = surrogate_index_models_tbl %>%
  bind_rows(gam_models_tbl %>%
              mutate(method = "gam"))

rm("gam_models_tbl")

## Cox PH model  ---------------------------------------------------------

# Estimate Cox PH models using the time-to-event endpoint directly. Note that we
# will be stratifying be trial in these models. We therefore do not include
# `logit_prob_infection_free` as predictor in these models. 

# Function that fits a Cox PH model with the given predictors.
cox_fitter = function(predictors_chr, weights_chr, case_cohort_ind_chr, analysis_set) {
  data_temp = ipd_tbl %>%
    filter(.data[[analysis_set]])
  # We truncate follow-up at 120 days to match the definition of the clinical
  # endpoint.
  data_temp = data_temp %>%
    mutate(event = ifelse(time_to_event > time_cumulative_incidence, 0, event))
  
  # Redefine the predictors as smooth functions.
  predictors_chr = paste0("bs(", predictors_chr, ")")
  # Define formula
  string_formula = paste0("Surv(time_to_event, event) ~ bs(risk_score_centered) + BMI_stratum + bs(Age) + Sex + strata(trial) +",
                          paste(predictors_chr, collapse = " + "))
  formula_final = as.formula(string_formula)
  
  # Fit logistic regression model.
  cox_fit = coxph(
    formula = formula_final,
    data = ipd_tbl %>% as.data.frame(), 
    weights = ipd_tbl %>% pull(any_of(weights_chr)),
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  return(cox_fit)
}

cox_models_tbl = prediction_model_settings %>%
  rowwise(surrogate, weighting, analysis_set) %>%
  summarise(fitted_model = list(cox_fitter(surrogate, weights_chr, case_cohort_ind_chr, analysis_set)))

surrogate_index_models_tbl = surrogate_index_models_tbl %>%
  bind_rows(cox_models_tbl %>%
              mutate(method = "cox"))

rm("cox_models_tbl")

## Prediction Accuracies of all Prediction Models --------------------------

# Add predictions to a new data set. We first add two artificial data sets where
# the surrogate is the original surrogate, either bindSpike or pseudoneutid50.
ipd_surr_indices_tbl = bind_rows(
  ipd_tbl %>%
    mutate(
      method = "untransformed surrogate",
      surrogate = "pseudoneutid50",
      weighting = "NA",
      include_Bserostatus = FALSE
    ) %>%
    rename(surrogate_index = "pseudoneutid50"),
  ipd_tbl %>%
    mutate(
      method = "untransformed surrogate",
      surrogate = "bindSpike",
      weighting = "NA",
      include_Bserostatus = FALSE
    ) %>%
    rename(surrogate_index = "bindSpike"),
  ipd_tbl %>%
    mutate(
      method = "untransformed surrogate",
      surrogate = "pseudoneutid50_adjusted",
      weighting = "NA",
      include_Bserostatus = FALSE
    ) %>%
    rename(surrogate_index = "pseudoneutid50_adjusted")
)
ipd_surr_indices_tbl = bind_rows(
  ipd_surr_indices_tbl,
  surrogate_index_models_tbl %>%
    ungroup() %>%
    filter(method %in% c("gam", "glm")) %>%
    rowwise(method, surrogate, weighting, analysis_set) %>%
    reframe(tibble(
      surrogate_index = predict(fitted_model, newdata = ipd_tbl, type = "response"),
      ipd_tbl
    )) %>%
    ungroup(),
  surrogate_index_models_tbl %>%
    ungroup() %>%
    filter(method %in% c("cox")) %>%
    rowwise(method, surrogate, weighting, analysis_set) %>%
    reframe(tibble(
      surrogate_index = 1 - predict(fitted_model, newdata = ipd_tbl %>%
                                      mutate(time_to_event = time_cumulative_incidence), type = "survival"),
      ipd_tbl
    )) %>%
    ungroup()
)
# Add variable containing the "correct" weight for the corresponding surrogate
# index.
ipd_surr_indices_tbl = ipd_surr_indices_tbl %>%
  mutate(
    sample_weight = ifelse(
      surrogate == "bindSpike",
      case_cohort_weight_bAb,
      case_cohort_weight_nAb
    )
  ) %>%
  # Only keep variables that are needed further on.
  select(
    method,
    weighting,
    analysis_set,
    trial,
    surrogate,
    sample_weight,
    surrogate_index,
    treatment,
    infection_120d,
    event,
    time_to_event,
    risk_score
  )

## Prediction Model Performance -------------------------------------------

roc_tbl = ipd_surr_indices_tbl %>%
  group_by(method, surrogate, weighting, trial, analysis_set) %>%
  filter(!is.na(surrogate_index)) %>%
  reframe(WeightedROC(
    guess = surrogate_index,
    label = infection_120d, 
    weight = sample_weight
  )) 

# Make and save the plots for all combinations of `surrogate` and `weighting`.
roc_ggplots = roc_tbl %>%
  filter(method != "untransformed surrogate") %>%
  group_by(analysis_set, surrogate) %>%
  summarise(data = list(pick(everything()))) %>%
  ungroup() %>%
  mutate(ggplot_object = purrr::pmap(
    .l = list(
      data = data,
      analysis_set = analysis_set, 
      surrogate = surrogate
    ),
    .f = function(data, analysis_set, surrogate) {
      data %>%
        ggplot(aes(color = method)) +
        geom_path(aes(FPR, TPR)) +
        coord_equal() +
        facet_wrap(~ trial) +
        theme(legend.position = "bottom", legend.box = "vertical") +
        scale_color_discrete(name = "Method") +
        geom_abline(intercept = 0, slope = 1) +
        ggtitle(paste0(analysis_set, " - ", surrogate)) 
    }
  ))

roc_ggplots %>%
  rowwise(analysis_set, surrogate) %>%
  summarise(
    ggsave(
      plot = ggplot_object,
      paste0("results/figures/application/roc-", surrogate, "-", analysis_set, ".pdf"),
      height = double_height,
      width = double_width,
      device = "pdf",
      units = "cm"
    )
  )


roc_tbl %>%
  group_by(method, surrogate, weighting, trial, analysis_set) %>%
  filter(method != "untransformed surrogate") %>%
  summarise(AUC = WeightedAUC(pick(c(
    TPR, FPR, FP, FN, threshold
  )))) %>%
  pivot_wider(names_from = "method", values_from = "AUC") %>%
  write.csv("results/tables/application/auc-surrogate-indices.csv")

rm("roc_tbl")

# Saving Results ----------------------------------------------------------

# Save data with estimated surrogate index to file.
saveRDS(ipd_surr_indices_tbl, file = "R/application/ipd_surr_indices_tbl.rds")
