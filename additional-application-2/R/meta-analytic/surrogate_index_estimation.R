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
library(dbarts)
library(glmnet)
library(mgcv)
library(speedglm)
library(SuperLearner)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}

data_location = "additional-application-2/data/processed_data_synthetic_high_d.rds"
out_file = "additional-application-2/results/raw-results/ipd_surr_indices_tbl.rds"

# Specify options for saving the plots to files
figures_dir = "additional-application-2/results/figures/"
tables_dir = "additional-application-2/results/tables/"


## Analysis Parameters --------------------------------------------------

# Number of bootstrap replications for computing within-trial covariance
# matrices.
time_cumulative_incidence = 80

## Data Preparation ----------------------------------------------------

# Load data. 
ipd_tbl = readRDS(data_location) %>%
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
  dplyr::select(-treatment)

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

# Convert Character variables to factors. This is more efficient in terms of
# memory.
ipd_tbl = ipd_tbl %>%
  mutate(trial = as.factor(trial))


# Prediction Models -------------------------------------------------------

# We first construct a tibble that contains on each row information about a
# particular prediction model.
prediction_model_settings = tibble(
  surrogate = list("pseudoneutid50_adjusted"),
  high_d_surrogate = FALSE,
  weights_chr = c("weight_default_nAb"),
  weights_chr_cox = c("case_cohort_weight_nAb"),
  case_cohort_ind_chr = c("Delta_nAb"),
  weighting = rep("unnormalized", 1)
)

# Add settings with high-dimensional surrogates included.
colnames_high_d_surrogates = ipd_tbl %>%
  colnames() %>%
  str_subset(pattern = "surrogate_marker")
prediction_model_settings = bind_rows(
  prediction_model_settings,
  prediction_model_settings %>%
    rowwise() %>%
    mutate(
      surrogate = list(c(surrogate, colnames_high_d_surrogates)),
      high_d_surrogate = TRUE
    )
)


prediction_model_settings = prediction_model_settings %>%
  cross_join(tibble(analysis_set = c("naive_only"))) %>%
  cross_join(tibble(include_risk_score = c(FALSE)))

# To save memory, we define ipd_tbl_cc which only contains the complete cases
# for either surrogates. The other rows are anyhow not used for fitting the
# prediction models.
ipd_tbl_cc = ipd_tbl %>%
  filter(weight_normalized_nAb != 0 | weight_normalized_bAb != 0)
  

## Superlearner -----------------------------------------------------------

# Predict the infection outcome using a Superlearner with a modified
# leave-one-trial-out CV procedure.
sl_fitter = function(predictors_chr,
                     weights_chr,
                     case_cohort_ind_chr,
                     analysis_set,
                     include_risk_score) {
  # Define the covariates that will be included as predictors.
  covariates = c("Sex", "HighRiskInd", "Age", "logit_prob_infection", predictors_chr)
  if (include_risk_score) {
    covariates = c(covariates, "risk_score_centered")
  }
  
  # Select the subset of the data corresponding to `analysis_set`
  data_temp = ipd_tbl_cc %>%
    filter(.data[[analysis_set]]) %>%
    # Drop the missing observations.
    filter(!is.na(infection_120d)) %>%
    filter(if_all(all_of(covariates), ~ !is.na(.)))
  # Compute the inverse probability weights as the predict of the inverse
  # probability of censoring and case-cohort weights.
  weights = data_temp %>%
    pull(any_of(weights_chr))
  data_temp$weights = weights
  
  # Define formula that excludes high-dimensional markers
  simple_formula_chr = paste0("infection_120d ~ ", paste(covariates %>% str_subset(pattern = "surrogate_marker", negate = TRUE), collapse = " + "))
  
  # Instantiate a set of learners.
  
  # Generalized linear model that does not use the high-dimensional marker.
  lrn_glm = Lrnr_glm_fast$new(formula = as.formula(simple_formula_chr))
  # Lasso GLM with all covariates and markers.
  lrn_glmnet = Lrnr_glmnet$new(nlambda = 30)
  
  # Alternative approach that first screens predictors, and then fits a glm on
  # the selected predictors.
  slscreener <- Lrnr_pkg_SuperLearner_screener$new("screen.glmnet")
  glm_learner <- Lrnr_glm$new()
  screen_and_glm <- Pipeline$new(slscreener, glm_learner)

  
  stack = Stack$new(
    screen_and_glm,
    lrn_glm,
    lrn_glmnet
  )
  

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
  covariates = c("Sex", "HighRiskInd", "Age", "logit_prob_infection", predictors_chr)
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
      .f = sl_fitter
    )
  )


plan(sequential)

## Prediction Accuracies of all Prediction Models --------------------------

# Add predictions to a new data set.
ipd_surr_indices_tbl = bind_rows(
  sl_models_tbl %>%
    rowwise(surrogate, weighting, analysis_set, include_risk_score, high_d_surrogate) %>%
    reframe(tibble(
      surrogate_index = sl_prediction_f(
        sl_fit = fitted_model,
        newdata = ipd_tbl,
        include_risk_score = include_risk_score,
        predictors_chr = surrogate
      ),
      ipd_tbl %>%
        dplyr::select(-str_subset(colnames(ipd_tbl), pattern = "surrogate_marker"))
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
  dplyr::select(
    high_d_surrogate,
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
  group_by(weighting, analysis_set, include_risk_score, high_d_surrogate, trial) %>%
  reframe(
    WeightedROC(
      guess = surrogate_index,
      label = infection_120d,
      weight = sample_weight * ipcw
    )
  )

# Make and save the plots for all combinations of `surrogate` and `weighting`.
roc_ggplots = roc_tbl %>%
  group_by(high_d_surrogate) %>%
  summarise(data = list(pick(everything()))) %>%
  ungroup() %>%
  mutate(ggplot_object = purrr::pmap(
    .l = list(
      data = data,
      high_d_surrogate = high_d_surrogate
    ),
    .f = function(data, high_d_surrogate) {
      data %>%
        filter(high_d_surrogate == .env$high_d_surrogate) %>%
        mutate(
          weighting = ifelse(weighting == "normalized", "Normalized", "Unnormalized")
        ) %>%
        ggplot() +
        geom_path(aes(FPR, TPR)) +
        facet_wrap(~ trial) +
        theme(legend.position = "bottom", legend.box = "vertical") +
        geom_abline(intercept = 0, slope = 1) +
        scale_color_discrete(name = "Estimated Surrogate Index") +
        ggtitle("ROC for estimated surrogate index") +
        # labs(subtitle = subtitle) +
        theme(legend.spacing.y = unit(0, "cm"),
              legend.box.spacing = unit(0, "cm"))
    }
  ))

roc_ggplots %>%
  rowwise(high_d_surrogate) %>%
  summarise(
    ggsave(
      plot = ggplot_object,
      paste0("roc-", high_d_surrogate, ".pdf"),
      path = figures_dir,
      height = double_height,
      width = double_width,
      device = "pdf",
      units = "cm"
    )
  )


roc_tbl %>%
  group_by(weighting, trial, analysis_set, include_risk_score, high_d_surrogate) %>%
  summarise(AUC = WeightedAUC(pick(c(
    TPR, FPR, FP, FN, threshold
  )))) %>%
  write.csv(paste0(tables_dir, "auc-surrogate-indices.csv"))

rm("roc_tbl")

# Saving Results ----------------------------------------------------------

# Save data with estimated surrogate index to file.
saveRDS(ipd_surr_indices_tbl, file = out_file)
