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
library(sl3)
library(origami)
library(splines)

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
  data_location = "data/processed_data.csv"
  out_file = "results/raw-results/application/ipd_surr_indices_tbl.rds"
  
  # Specify options for saving the plots to files
  figures_dir = "results/figures/application/surrogate-index"
  tables_dir = "results/tables/application/surrogate-index"
} else if (data_set == "synthetic") {
  data_location = "data/processed_data_synthetic.csv"
  out_file = "results/raw-results/application-synthetic/ipd_surr_indices_tbl.rds"
  
  # Specify options for saving the plots to files
  figures_dir = "results/figures/application-synthetic/surrogate-index"
  tables_dir = "results/tables/application-synthetic/surrogate-index"
}

## Analysis Parameters --------------------------------------------------

# Number of bootstrap replications for computing within-trial covariance
# matrices.
time_cumulative_incidence = 80

## Data Preparation ----------------------------------------------------


# Load data. We immediately start working with tibbles instead of classic data
# frames.
ipd_tbl = read.csv(data_location) %>%
  tibble() %>%
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


# Compute pseudo-values. These are computed for each trial separately. We first
# fit a separate KM curve by trial and treatment group.
surv_fit_all = survfit(Surv(time_to_event, event) ~ strata(treatment, trial), data = ipd_tbl)
ipd_tbl$pseudo_value = 1 - pseudo(surv_fit_all, times = time_cumulative_incidence)

# Compute placebo cumulative incidence rate per trial. This will be used later
# on on some plots to standardize some things across trials.
cumulative_incidence_control_tbl = ipd_tbl %>%
  group_by(trial, treatment) %>%
  summarise(prob_infection = mean(pseudo_value)) %>%
  ungroup() %>%
  filter(treatment == 0) %>%
  select(-treatment)

ipcw_estimator  = function(time_to_event, event) {
  survfit_object = survfit(Surv(time_to_event, 1 - event) ~ 1)
  
  time =  unique(pmin(time_to_event, time_cumulative_incidence))
  surv_probs_tbl = summary(survfit_object, times =  time)[c("surv", "time")] %>%
    as_tibble()
  
  censoring_probs = tibble(time = pmin(time_to_event, time_cumulative_incidence))  %>%
    left_join(surv_probs_tbl) %>%
    pull(surv)

  return(1 / censoring_probs)
}

ipd_tbl = ipd_tbl %>%
  group_by(trial, treatment) %>%
  mutate(ipcw = ipcw_estimator(time_to_event, event),
  ipcw = ifelse(time_to_event < time_cumulative_incidence & event == 0, 0, ipcw)) %>%
  ungroup()


ipd_tbl = ipd_tbl %>%
  mutate(# If the patient has an infection or censoring event before 80 days, we set
    # the infection event to NA if the patients was censored. If a patient is
    # under observation for more than 80 days, the infection event is set to
    # zero.
    infection_120d = ifelse(
      time_to_event > time_cumulative_incidence,
      0,
      ifelse(
        time_to_event == time_cumulative_incidence,
        event,
        ifelse(event == 1, 1, NA)
      )
    ))

ipd_tbl = ipd_tbl %>%
  mutate(weight_default_nAb = case_cohort_weight_nAb * ipcw,
         weight_default_bAb = case_cohort_weight_bAb * ipcw)

# Add placebo cumulative incidence rate as a covariate to the trial. We can
# adjust for this variable later to indirectly adjust for different forces of
# infection across the trials.
ipd_tbl = ipd_tbl %>%
  left_join(cumulative_incidence_control_tbl %>%
              mutate(logit_prob_infection = log(
                prob_infection / (1 - prob_infection)
              )))


# Make sure that each trial-treatment group combination gets the same total
# weight. We don't want one trial to have an undue influence on the predictions.
ipd_tbl = ipd_tbl %>%
  left_join(
    ipd_tbl %>%
      group_by(trial, treatment) %>%
      summarize(
        # Complete case weights for models that use infection at day x as outcome.
        total_weight_nAb = sum(
          weight_default_nAb * !is.na(pseudoneutid50) * !is.na(infection_120d)
        ),
        total_weight_bAb = sum(
          weight_default_bAb * !is.na(bindSpike) * !is.na(infection_120d)
        ),
        # Case-cohort weights (i.e., for the cox models)
        total_weight_nAb_cox = sum(case_cohort_weight_nAb * !is.na(pseudoneutid50)),
        total_weight_bAb_cox = sum(case_cohort_weight_bAb * !is.na(bindSpike))
      ) %>%
      ungroup()
  ) %>%
  mutate(
    weight_normalized_nAb = (weight_default_nAb) / total_weight_nAb,
    weight_normalized_bAb = (weight_default_bAb) / total_weight_bAb,
    weight_normalized_nAb_cox = case_cohort_weight_nAb / total_weight_nAb_cox,
    weight_normalized_bAb_cox = case_cohort_weight_bAb / total_weight_bAb_cox
  )

# Every trial now receives the same total weight.
ipd_tbl %>%
  group_by(trial, treatment) %>%
  summarise(
    sum(
      weight_normalized_nAb * !is.na(pseudoneutid50) * !is.na(infection_120d)
    ),
    sum(
      weight_normalized_bAb * !is.na(bindSpike) * !is.na(infection_120d)
    ),
    sum(weight_normalized_nAb_cox * !is.na(pseudoneutid50)),
    sum(weight_normalized_bAb_cox * !is.na(bindSpike))
  )

# # Drop variables that are not needed further on.
# ipd_tbl = ipd_tbl %>%
#   select(
#     -X,-Ptid,-USUBJID,-(BMI_underweight:abrogation_coefficient),-total_weight_nAb,-total_weight_bAb,-Wstratum,-Bserostatus
#   )

# Convert Character variables to factors. This is more efficient in terms of
# memory.
ipd_tbl = ipd_tbl %>%
  mutate(trial = as.factor(trial))


# Prediction Models -------------------------------------------------------

# We first construct a tibble that contains on each row information about a
# particular prediction model.
prediction_model_settings = tibble(
  surrogate = c("pseudoneutid50", "bindSpike", "pseudoneutid50_adjusted"),
  weights_chr = c(
    "weight_default_nAb",
    "weight_default_bAb",
    "weight_default_nAb"
  ),
  weights_chr_cox = c(
    "case_cohort_weight_nAb",
    "case_cohort_weight_bAb",
    "case_cohort_weight_nAb"
  ),
  case_cohort_ind_chr = c("Delta_nAb", "Delta_bAb", "Delta_nAb"),
  weighting = rep("unnormalized", 3)
)

# # For the current analyses, we only consider the unnormalized weights.
# prediction_model_settings = prediction_model_settings %>%
#   bind_rows(
#     tibble(
#       surrogate = c("pseudoneutid50", "bindSpike", "pseudoneutid50_adjusted"),
#       weights_chr = c(
#         "weight_normalized_nAb",
#         "weight_normalized_bAb",
#         "weight_normalized_nAb"
#       ),
#       weights_chr_cox = c(
#         "weight_normalized_nAb_cox",
#         "weight_normalized_bAb_cox",
#         "weight_normalized_nAb_cox"
#       ),
#       case_cohort_ind_chr = c("Delta_nAb", "Delta_bAb", "Delta_nAb"),
#       weighting = rep("normalized", 3)
#     )
#   )

prediction_model_settings = prediction_model_settings %>%
  cross_join(tibble(analysis_set = c("naive_only"))) %>%
  cross_join(tibble(include_risk_score = c(FALSE)))

## Parametric Prediction Model --------------------------------------------

# We're currently not using the GLM because it's very similar to the GAM.

# # Estimate logistic regression model for the probability of infection given the
# # baseline covariates and the surrogate. Note that we missing predictors: the
# # antibody titers are missing according to the case-cohort sampling mechanism.
#
#
# # Function that fits a logistic regression model with the given predictors.
# glm_fitter = function(predictors_chr,
#                       weights_chr,
#                       case_cohort_ind_chr,
#                       analysis_set) {
#   data_temp = ipd_tbl %>%
#     filter(.data[[analysis_set]])
#   # Redefine the predictors as smooth functions.
#   predictors_chr = paste0("bs(", predictors_chr, ")")
#   # Define formula
#   string_formula = paste0(
#     "infection_120d ~ bs(risk_score_centered) + BMI_stratum + bs(Age) + Sex + logit_prob_infection_free +",
#     paste(predictors_chr, collapse = " + ")
#   )
#   formula_final = as.formula(string_formula)
#
#   # Fit logistic regression model.
#   glm_fit = glm(
#     formula = formula_final,
#     data = data_temp,
#     weights = data_temp %>%
#       pull(any_of(weights_chr)),
#     family = quasibinomial(),
#     model = FALSE,
#     x = FALSE,
#     y = FALSE
#   )
#   return(glm_fit)
# }
#
# glm_models_tbl = prediction_model_settings %>%
#   rowwise(surrogate, weighting, analysis_set) %>%
#   summarise(fitted_model = list(glm_fitter(surrogate, weights_chr, case_cohort_ind_chr, analysis_set)))
#
# surrogate_index_models_tbl = glm_models_tbl %>%
#   mutate(method = "glm")
#
# rm("glm_models_tbl")




# ## GAM --------------------------------------------------------------------
# 
# # Estimate logistic GAM for the probability of infection given the
# # baseline covariates and the surrogate. Note that we missing predictors: the
# # antibody titers are missing according to the case-cohort sampling mechanism.
# 
# # Function that fits a logistic regression model with the given predictors.
# gam_fitter = function(predictors_chr,
#                       weights_chr,
#                       case_cohort_ind_chr,
#                       analysis_set, 
#                       include_risk_score) {
#   # Select the subset of the data corresponding to `analysis_set`
#   data_temp = ipd_tbl %>%
#     filter(.data[[analysis_set]])
#   # Compute the inverse probability weights as the predict of the inverse
#   # probability of censoring and case-cohort weights.
#   weights = data_temp %>%
#     pull(any_of(weights_chr))
#   
#   # Redefine the predictors as smooth functions.
#   predictors_chr = paste0("s(", predictors_chr, ", k = 4)")
#   # Define formula
#   string_formula = paste0(
#     "infection_120d ~ Sex + HighRiskInd + BMI_stratum + s(Age) + logit_prob_infection +",
#     paste(predictors_chr, collapse = " + ")
#   )
#   if (include_risk_score) {
#     string_formula = paste0(
#       string_formula,
#       ' + risk_score_centered'
#     )
#   }
#   formula_final = as.formula(string_formula)
#   
#   # Fit logistic regression model.
#   gam_fit = gam(
#     formula = formula_final,
#     data = data_temp %>%
#       as.data.frame(),
#     weights = data_temp %>%
#       pull(any_of(weights_chr)),
#     family = quasibinomial()
#   )
#   return(gam_fit)
# }
# 
# gam_models_tbl = prediction_model_settings %>%
#   rowwise(surrogate, weighting, analysis_set, include_risk_score) %>%
#   summarise(fitted_model = list(
#     gam_fitter(surrogate, weights_chr, case_cohort_ind_chr, analysis_set, include_risk_score)
#   ))
# 
# surrogate_index_models_tbl = gam_models_tbl %>%
#   mutate(method = "gam")
# 
# rm("gam_models_tbl")

## Cox PH model  ---------------------------------------------------------

# Estimate Cox PH models using the time-to-event endpoint directly. Note that we
# will be stratifying be trial in these models. We therefore do not include
# `logit_prob_infection_free` as predictor in these models.

# Function that fits a Cox PH model with the given predictors.
cox_fitter = function(predictors_chr,
                      weights_chr,
                      case_cohort_ind_chr,
                      analysis_set, 
                      include_risk_score) {
  data_temp = ipd_tbl %>%
    filter(.data[[analysis_set]])
  # Compute the inverse probability weights as the predict of the inverse
  # probability of censoring and case-cohort weights.
  weights = data_temp %>%
    pull(any_of(weights_chr))
  # We truncate follow-up at 120 days to match the definition of the clinical
  # endpoint.
  data_temp = data_temp %>%
    mutate(
      event = ifelse(time_to_event > time_cumulative_incidence + 40, 0, event),
      time_to_event = ifelse(
        time_to_event > time_cumulative_incidence,
        time_cumulative_incidence,
        time_to_event
      )
    )
  
  # Redefine the predictors as smooth functions.
  predictors_chr = paste0("ns(", predictors_chr, ", knots = c(1.5, 2.5), Boundary.knots = c(0, 4))")
  # Define formula
  string_formula = paste0(
    "Surv(time_to_event, event) ~ Sex + HighRiskInd + BMI_stratum + bs(Age) + strata(trial) +",
    paste(predictors_chr, collapse = " + ")
  )
  if (include_risk_score) {
    string_formula = paste0(
      string_formula,
      ' + risk_score_centered'
    )
  }
  formula_final = as.formula(string_formula)
  
  # Fit logistic regression model.
  cox_fit = coxph(
    formula = formula_final,
    data = data_temp %>% as.data.frame(),
    weights = weights,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  return(cox_fit)
}

cox_models_tbl = prediction_model_settings %>%
  rowwise(surrogate, weighting, analysis_set, include_risk_score) %>%
  summarise(fitted_model = list(
    cox_fitter(surrogate, weights_chr_cox, case_cohort_ind_chr, analysis_set, include_risk_score)
  ))


# x = seq(from = 0, to = 5, length.out = 100)
# ns_matrix = ns(x = x, knots = c(1.5, 2.5), Boundary.knots = c(0, 4)) 
# coefs = coef(cox_models_tbl$fitted_model[[8]])[10:12]
# plot(x, ns_matrix %*% coefs)

surrogate_index_models_tbl = cox_models_tbl %>%
  mutate(method = "cox")

rm("cox_models_tbl")

## Superlearner -----------------------------------------------------------

# Predict the infection outcome using a Superlearner with a modified
# leave-one-trial-out CV procedure.
sl_fitter = function(predictors_chr,
                     weights_chr,
                     case_cohort_ind_chr,
                     analysis_set,
                     include_risk_score) {
  # Define the covariates that will be included as predictors.
  covariates = c("Sex", "HighRiskInd", "BMI_stratum", "Age", "logit_prob_infection", predictors_chr)
  if (include_risk_score) {
    covariates = c(covariates, "risk_score_centered")
  }
  
  # Select the subset of the data corresponding to `analysis_set`
  data_temp = ipd_tbl %>%
    filter(.data[[analysis_set]]) %>%
    # Drop the missing observations.
    filter(!is.na(infection_120d)) %>%
    filter(if_all(all_of(covariates), ~ !is.na(.)))
  # Compute the inverse probability weights as the predict of the inverse
  # probability of censoring and case-cohort weights.
  weights = data_temp %>%
    pull(any_of(weights_chr))
  data_temp$weights = weights
  
  # Instantiate a set of learners.
  lrn_glm = Lrnr_glm$new()
  lrn_glm_surr_only = Lrnr_glm$new(formula = paste0("~ logit_prob_infection + ", predictors_chr))
  lrn_glm_bs2 = Lrnr_glm$new(
    formula = paste0(
      "~ .",
      " + ns(",
      predictors_chr,
      ", knots = c(1.5, 2.5), Boundary.knots = c(0, 4))"
    )
  )
  lrn_glm_interactions = Lrnr_glm$new(formula = paste0("~.", " + ", predictors_chr, ":Sex"))
  lrn_glm_interactions_bs2 = Lrnr_glm$new(
    formula = paste0(
      "~ .",
      " + ",
      predictors_chr,
      ":Sex",
      " + ns(",
      predictors_chr,
      ", knots = c(1.5, 2.5), Boundary.knots = c(0, 4))"
    )
  )
  lrn_gam = Lrnr_gam$new()

  
  stack = Stack$new(lrn_glm, lrn_glm_surr_only, lrn_glm_bs2, lrn_glm_interactions, lrn_glm_interactions_bs2, lrn_gam)

  task = make_sl3_Task(
    data = data_temp,
    outcome = "infection_120d",
    covariates = covariates,
    weights = "weights",
    outcome_type = variable_type("binomial"),
    id = "trial",
    folds = make_folds(
      n = nrow(data_temp),
      cluster_ids = data_temp$trial,
      fold_fun = folds_loo
    )
  )

  sl = Lrnr_sl$new(
    learners = stack,
    keep_extra = FALSE,
    metalearner = Lrnr_solnp$new(eval_function = loss_loglik_binomial)
  )

  sl_fit = sl$train(task = task)

  return(sl_fit)
}

sl_prediction_f = function(sl_fit, newdata, predictors_chr, include_risk_score) {
  # Define the covariates that will be included as predictors.
  covariates = c("Sex", "HighRiskInd", "BMI_stratum", "Age", "logit_prob_infection", predictors_chr)
  if (include_risk_score) {
    covariates = c(covariates, "risk_score_centered")
  }

  # Add row number, which will be need later on.
  newdata = newdata %>%
    mutate(row_number = row_number())

  # Split the data set into one with one with no missing value and one with missing values.
  newdata_no_missing = newdata %>%
    filter(if_all(all_of(covariates), ~ !is.na(.)))
  newdata_missing = newdata %>%
    filter(if_any(all_of(covariates), is.na))

  # Do predictions for the rows with no missing values.
  prediction_task = make_sl3_Task(
    data = newdata_no_missing,
    covariates = covariates,
    outcome = "infection_120d"
  )

  newdata_no_missing$pred = sl_fit$predict(task = prediction_task)
  newdata_missing$pred = NA
  # Return the predictions in the original data's order.
  return(
    bind_rows(
      newdata_missing,
      newdata_no_missing
    ) %>% arrange(row_number) %>%
      pull(pred)
  )
}


sl_models_tbl = prediction_model_settings %>%
  mutate(
    fitted_model = future_pmap(
      .l = list(
        predictors_chr = surrogate,
        weights_chr = weights_chr,
        case_cohort_ind_chr = case_cohort_ind_chr,
        analysis_set = analysis_set,
        include_risk_score = include_risk_score
      ),
      .f = sl_fitter, 
      .options = furrr_options(
        packages = "splines",
        seed = TRUE
      )
    )
  )

sl_models_tbl = sl_models_tbl %>%
  select(fitted_model, surrogate, weighting, analysis_set, include_risk_score)

surrogate_index_models_tbl = surrogate_index_models_tbl %>%
  bind_rows(sl_models_tbl %>%
              mutate(method = "sl"))

rm("sl_models_tbl")

plan(sequential)

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
  # surrogate_index_models_tbl %>%
  #   ungroup() %>%
  #   filter(method %in% c("gam", "glm")) %>%
  #   rowwise(method, surrogate, weighting, analysis_set, include_risk_score) %>%
  #   reframe(tibble(
  #     surrogate_index = predict(fitted_model, newdata = ipd_tbl, type = "response"),
  #     ipd_tbl
  #   )) %>%
  #   ungroup(),
  surrogate_index_models_tbl %>%
    ungroup() %>%
    filter(method %in% c("cox")) %>%
    rowwise(method, surrogate, weighting, analysis_set, include_risk_score) %>%
    reframe(tibble(
      surrogate_index = 1 - predict(
        fitted_model,
        newdata = ipd_tbl %>%
          mutate(time_to_event = time_cumulative_incidence, trial = "Janssen"),
        type = "survival"
      ),
      ipd_tbl
    )) %>%
    ungroup(),
  surrogate_index_models_tbl %>%
    ungroup() %>%
    filter(method %in% c("sl")) %>%
    rowwise(method, surrogate, weighting, analysis_set, include_risk_score) %>%
    reframe(tibble(
      surrogate_index = sl_prediction_f(
        sl_fit = fitted_model,
        newdata = ipd_tbl,
        include_risk_score = include_risk_score,
        predictors_chr = surrogate
      ),
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
    include_risk_score,
    trial,
    surrogate,
    sample_weight,
    surrogate_index,
    treatment,
    infection_120d,
    event,
    time_to_event,
    risk_score,
    ipcw,
    Sex,
    HighRiskInd,
    BMI_stratum,
    Age
  )

## Prediction Model Performance -------------------------------------------


roc_tbl = ipd_surr_indices_tbl %>%
  # Reorder trial factor.
  mutate(
    trial = forcats::fct_recode(
      trial,
      "Moderna (naive)" = "Moderna",
      "AstraZeneca (naive)" = "AstraZeneca",
      "Janssen (naive)" = "Janssen",
      "Novavax (naive)" = "Novavax",
      "Sanofi 1 (naive)" = "Sanofi-1 naive",
      "Sanofi 1 (non-naive)" = "Sanofi-1 non-naive",
      "Sanofi 2 (naive)" = "Sanofi-2 naive",
      "Sanofi 2 (non-naive)" = "Sanofi-2 non-naive"
    ),
    trial = fct_relevel(
      trial,
      "Moderna (naive)",
      "AstraZeneca (naive)",
      "Janssen (naive)",
      "Novavax (naive)",
      "Sanofi 1 (naive)",
      "Sanofi 1 (non-naive)",
      "Sanofi 2 (naive)",
      "Sanofi 2 (non-naive)"
    )
  ) %>%
  filter(!is.na(surrogate_index) & !(is.na(infection_120d))) %>%
  group_by(method, surrogate, weighting, trial, analysis_set, include_risk_score) %>%
  reframe(
    WeightedROC(
      guess = surrogate_index,
      label = infection_120d,
      weight = sample_weight * ipcw
    )
  )

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
      surrogate_chr = switch(
        surrogate,
        bindSpike = "IgG Spike",
        pseudoneutid50 = "nAb ID50",
        pseudoneutid50_adjusted = "adjusted nAb ID50"
      )
      analysis_set_chr = switch(
        analysis_set,
        first_four = "first four trials",
        naive_only = "all trials, naive subjects only",
        mixed = "all trials"
      )
      subtitle = paste0("Surrogate index for ",
                        surrogate_chr,
                        " using ",
                        analysis_set_chr)
      data %>%
        mutate(
          weighting = ifelse(weighting == "normalized", "Normalized", "Unnormalized"),
          method = ifelse(method == "gam", "GAM logistic regression", ifelse(method == "cox", "Cox Model", "SuperLearner"))
        ) %>%
        ggplot(aes(color = method)) +
        geom_path(aes(FPR, TPR)) +
        facet_wrap(~ trial) +
        theme(legend.position = "bottom", legend.box = "vertical") +
        scale_color_discrete(name = "Method") +
        geom_abline(intercept = 0, slope = 1) +
        scale_color_discrete(name = "Estimated Surrogate Index") +
        ggtitle("ROC for estimated surrogate index") +
        labs(subtitle = subtitle) +
        theme(legend.spacing.y = unit(0, "cm"),
              legend.box.spacing = unit(0, "cm"))
    }
  ))

roc_ggplots %>%
  rowwise(analysis_set, surrogate) %>%
  summarise(
    ggsave(
      plot = ggplot_object,
      paste0("roc-", surrogate, "-", analysis_set, ".pdf"),
      path = figures_dir,
      height = double_height,
      width = double_width,
      device = "pdf",
      units = "cm"
    )
  )


roc_tbl %>%
  group_by(method, surrogate, weighting, trial, analysis_set, include_risk_score) %>%
  filter(method != "untransformed surrogate") %>%
  summarise(AUC = WeightedAUC(pick(c(
    TPR, FPR, FP, FN, threshold
  )))) %>%
  pivot_wider(names_from = "method", values_from = "AUC") %>%
  write.csv(paste0(tables_dir, "/auc-surrogate-indices.csv"))

rm("roc_tbl")

# Saving Results ----------------------------------------------------------

# Save data with estimated surrogate index to file.
saveRDS(ipd_surr_indices_tbl, file = out_file)
