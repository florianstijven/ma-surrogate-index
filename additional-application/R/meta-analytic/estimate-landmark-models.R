args = commandArgs(trailingOnly=TRUE)

# Load packages.
library(tidyverse)
library(survival)
library(GGally)
library(mgcv)
library(future)
library(furrr)
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
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}

# Options ----------------------------------------------------------------

landmark_times = 365.25 * c(1 / 12, 2 / 12, 3 / 12, 4 / 12, 5 / 12, 0.5, 1, 1.5, 2)
time_cumulative_incidence = 2.5 * 365.25

full_data_new_endpoints_long = readRDS(file = "additional-application/R/data/intermediate-objects/full_data_tbl_long.rds")

endpoints = c("CVMB_DTH")

# Data Preparation --------------------------------------------------------

## Landmark data format -----------------------------------------------

# Convert the data to a format suitable for estimating landmark prediction
# models. For each landmark time, we find the three most recent measurements of
# the surrogates for each patient and construct a data set in which each row
# corresponds to a patient and the columns include the three most recent
# measurements of the surrogates before the landmark time, the baseline
# measurements of the surrogates, and the time-to-event endpoints.

full_data_new_endpoints_landmark_tbl = full_data_new_endpoints_long %>%
  filter(FALSE) %>%
  mutate(timing = numeric(), landmark_time = numeric()) %>%
  dplyr::select(-VISNAM1A, days, months)
for (landmark_time in landmark_times) {
  # Given the landmark time, find the times of the three most recent measurements.
  days_measurements = unique(full_data_new_endpoints_long$days) %>%
    na.omit()
  days = sort(days_measurements[days_measurements <= landmark_time])
  latest_days = tail(days, 3)
  # Only post-randomization time points should be included. The baseline
  # measurements of blood pressure and sitting pulse will be added later on.
  latest_days = latest_days[latest_days > 0]
  # Construct the correct data set for the given landmark time.
  full_data_new_endpoints_landmark_tbl = bind_rows(
    full_data_new_endpoints_landmark_tbl,
    full_data_new_endpoints_long %>%
      filter(days %in% latest_days) %>%
      group_by(SID1A) %>%
      arrange(-days, .by_group = TRUE) %>%
      mutate(timing = row_number()) %>%
      dplyr::select(-days, -months, -VISNAM1A) %>%
      pivot_wider(
        values_from = c(SBPAVG3N, DBPAVG3N, STNPLS1N),
        names_from = timing
      ) %>% 
      mutate(landmark_time = landmark_time)
  )
}

# Add the baseline measurements (at day -14) to the landmark data set.
full_data_new_endpoints_landmark_tbl = full_data_new_endpoints_landmark_tbl %>%
  left_join(
    full_data_new_endpoints_long %>%
      filter(days == -14) %>%
      dplyr::select(
        SID1A,
        SBPAVG3N_baseline = SBPAVG3N,
        DBPAVG3N_baseline = DBPAVG3N,
        STNPLS1N_baseline = STNPLS1N
      ),
    by = "SID1A"
  )

# Restructure the data set to have one row per patient and endpoint, with
# separate columns for the time and censoring status for each endpoint.
full_data_new_endpoints_landmark_tbl = full_data_new_endpoints_landmark_tbl %>%
  pivot_longer(cols = any_of(contains(endpoints))) %>%
  mutate(
    variable_type = ifelse(stringr::str_starts(name, "C4"), "censored", "time"),
    endpoint = stringr::str_sub(name, start = 3)
  ) %>%
  dplyr::select(-name) %>%
  pivot_wider(values_from = "value", names_from = "variable_type")

## IPCW -------------------------------------------------------------------

# For general prediction models, we will use the event status at month 36 as
# outcome, denoted by Y. The prediction algorithms estimate the following
# regression function: E[Y | X, S, T > landmark_time], where X are the baseline
# covariates, S are the surrogates, and T is the time to event. We have to deal
# with two types of censoring: (i) patients that are censored before the
# landmark time, and (ii) patients that are censored after the landmark time but
# before month 36.

# For the first type of censoring, we simply exclude these patients from the
# analysis, as they do not contribute to the estimation of the regression
# function. This is valid under the assumption that the censoring is independent
# of the event time conditionally on the covariates and surrogates included in
# the model.

# For the second type of censoring, we can include these patients in the
# analysis, but we have to account for the fact that their event status at month
# 36 is not observed. We will do this by using inverse probability of censoring
# weights (IPCW). Let D be the censoring status at month 36 (1 if censored, 0 if
# event observed). The IPCW for a patient is defined as w_i = 1 / P(C > min(y_i,
# t_0) | C > t, T > t) where C is the censoring time, y_i is the observed event
# time, t is the landmark time and t_0 = 36 is the time of interest. We can
# estimate the censoring distribution using the Kaplan-Meier estimator, where we
# treat events as censored and censoring times as events.


# Function to estimate inverse probability of censoring weights. It takes a
# vector of event times and event indicators (1 if event, 0 if censored), and
# the landmark time. It returns the IPCWs for each subject.
ipcw_estimator  = function(time_to_event, event, landmark_time) {
  # We only use patients that are still under observation at the landmark time for
  # estimating the censoring distribution. We further set the IPCW to zero for
  # patients that are censored before the landmark time, as they do not contribute
  # to the estimation of the regression function.
  included = time_to_event > landmark_time
  included[is.na(included)] = FALSE
  
  time_to_event_landmark = time_to_event[included] - landmark_time
  event_landmark = event[included]
  
  survfit_object = survfit(Surv(time_to_event_landmark, event_landmark) ~ 1)
  
  time = unique(pmin(
    time_to_event_landmark,
    time_cumulative_incidence - landmark_time
  ))
  surv_probs_tbl = summary(survfit_object, times =  time)[c("surv", "time")] %>%
    as_tibble()
  
  censoring_probs = tibble(time = pmin(
    time_to_event_landmark,
    time_cumulative_incidence - landmark_time
  ))  %>%
    left_join(surv_probs_tbl) %>%
    pull(surv)
  
  weights = rep(NA, length(time_to_event))
  weights[included] = 1 / censoring_probs
  
  if (any(weights == Inf |
          weights == 0, na.rm = TRUE))
    stop("invalid weights computed")
  
  return(weights)
}

# Compute the IPCW for each patient. These weights are based on the KM estimate
# by country and treatment arm. We immediately also do this for each of the new
# composite endpoints separately to take into account the possibility that
# the censoring time is not the same across all endpoints.
full_data_new_endpoints_landmark_tbl = full_data_new_endpoints_landmark_tbl %>%
  group_by(COU1A, TRTREG1C, endpoint, landmark_time) %>%
  mutate(
    ipcw = ipcw_estimator(time, censored, landmark_time[1]),
    ipcw = ifelse(time < time_cumulative_incidence &
                    censored == 1, 0, ipcw)
  ) %>%
  ungroup()

# Landmark Prediction Models ----------------------------------------------

## Cox PH Models ----------------------------------------------------------

# For patients that remain under observations up to a certain time point, we
# estimate a landmark prediction models. This model is based on the patients at
# under observation at X months. We further use penalized cox models.

# Function to estimate landmark model.

landmark_cox_model = function(data, type = "full") {
  formula_landmark = formula(
    Surv(time, I(1 - censored)) ~
      PCI_B + SEX1C + CURSMK1C + CVSDIS5C + CVSDIS7C +
      SBPAVG3N_baseline + DBPAVG3N_baseline + STNPLS1N_baseline +
      splines::ns(
        BMI_1N,
        knots = c(20, 25, 30),
        Boundary.knots = c(18, 40)
      ) +
      splines::ns(
        AGE_1N,
        knots = c(67.5, 75),
        Boundary.knots = c(60, 85)
      )
  )
  if (type == "full") {
    # Add the intermediate endpoints to the formula.
    formula_landmark = update.formula(
      formula_landmark,
      . ~ . + splines::ns(
        SBPAVG3N_1,
        knots = c(120, 140),
        Boundary.knots = c(110, 180)
      ) +
        splines::ns(
          SBPAVG3N_1 - DBPAVG3N_1,
          knots = c(45, 55, 65),
          Boundary.knots = c(40, 80)
        ) +
        splines::ns(
          STNPLS1N_1,
          knots = c(60, 100),
          Boundary.knots = c(40, 140)
        ) +
        splines::ns(
          SBPAVG3N_1 - SBPAVG3N_2,
          knots = c(-15, 15),
          Boundary.knots = c(-40, 40)
        ) +
        splines::ns(
          (SBPAVG3N_1 - DBPAVG3N_1) - (SBPAVG3N_2 - DBPAVG3N_2),
          knots = c(-15, 15),
          Boundary.knots = c(-35, 35)
        ) +
        splines::ns(
          STNPLS1N_1 - STNPLS1N_2,
          knots = c(-8, 8),
          Boundary.knots = c(-15, 15)
        )
    )
  }
  
  cox_fit = coxph(formula = formula_landmark,
                  data = data,
                  model = TRUE)
  return(cox_fit)
}

# For each time-to-event endpoint,
landmark_cox_models_tbl = expand_grid(landmark_time_model = landmark_times, endpoint_model = endpoints) %>%
  group_by(landmark_time_model, endpoint_model) %>%
  summarize(
    cox_fit = list(
      landmark_cox_model(
        full_data_new_endpoints_landmark_tbl %>%
          filter(
            landmark_time == landmark_time_model[1],
            endpoint == endpoint_model[1]
          ) %>%
          # The landmark model only uses data on patients that are still under
          # observation at the landmark time.
          filter(time > landmark_time_model[1])
      )
    ),
    cox_fit_base = list(
      landmark_cox_model(
        full_data_new_endpoints_landmark_tbl %>%
          filter(
            landmark_time == landmark_time_model[1],
            endpoint == endpoint_model[1]
          ) %>%
          # The landmark model only uses data on patients that are still under
          # observation at the landmark time.
          filter(time > landmark_time_model[1]),
        type = "base"
      )
    )
  ) %>%
  rename(landmark_time = landmark_time_model, endpoint = endpoint_model)


# Restructure to have one model per row.
landmark_cox_models_tbl = landmark_cox_models_tbl %>%
  pivot_longer(
    cols = c("cox_fit", "cox_fit_base"),
    values_to = "cox_model",
    names_to = "model_type"
  ) %>%
  mutate(model_type = ifelse(model_type == "cox_fit", "full", "base"))

## SuperLearner Models ----------------------------------------------------

# Predict the infection outcome using a SuperLearner with a modified
# leave-one-trial-out CV procedure.
sl_fitter = function(predictors_chr, endpoint, landmark) {
  # Define the covariates that will be included as predictors.
  covariates = c(
    "PCI_B",
    "SEX1C",
    "DBT1",
    "REN1",
    "CURSMK1C",
    "CVSDIS5C",
    "CVSDIS7C",
    "BMI_1N",
    "AGE_1N",
    "WAICIR2N",
    "B_CREA",
    "B_GLUC",
    "B_POT",
    "B_CHOL",
    "B_HDL",
    "SBPAVG3N_baseline",
    "DBPAVG3N_baseline",
    "STNPLS1N_baseline",
    predictors_chr
  )
  
  # Base formula for the ML methods.
  base_formula_chr = paste0("event_status ~ ", paste(covariates, collapse = " + "))
  # Formula for a simple learner that only includes the surrogates.
  simple_formula_chr = paste0("event_status ~ ", paste(c("STNPLS1N_1", "REN1", "AGE_1N"), collapse = " + "))
  
  
  # Filter the data set to the relevant clinical endpoint.
  data_temp = full_data_new_endpoints_landmark_tbl %>%
    filter(endpoint == .env$endpoint) %>%
    # Filter only rows corresponding to the selected landmark time.
    filter(landmark_time == .env$landmark) %>%
    # The landmark model only uses data on patients that are still under
    # observation at the landmark time.
    filter(time > .env$landmark) %>%
    # Drop the missing observations.
    filter(!is.na(time)) %>%
    filter(if_all(all_of(covariates), ~ !is.na(.))) %>%
    # Compute the binary endpoint (i.e., event status at month 36).
    mutate(event_status = ifelse(time <= time_cumulative_incidence &
                                   censored == 0, 1, 0)) %>%
    # Only keep rows where the event_status was actually observed.
    filter(!(time < time_cumulative_incidence & censored == 1)) %>%
    # keep only variables needed further on.
    dplyr::select(any_of(c(
      "COU1A", "event_status", "ipcw", covariates
    )))
  
  # Check for missing values.
  data_temp %>%
    filter(!if_all(all_of(covariates), ~ !is.na(.))) %>%
    nrow()
  
  # Instantiate a set of learners.
  lrn_mean = Lrnr_mean$new()
  lrn_glm = Lrnr_glm_fast$new(formula = base_formula_chr)
  lrn_glmnet = Lrnr_glmnet$new()
  lrn_dbarts = Lrnr_dbarts$new()
  # lrn_xgboost = Lrnr_xgboost$new()
  lrn_nnet = Lrnr_nnet$new(trace = FALSE)
  
  slscreener <- Lrnr_pkg_SuperLearner_screener$new("screen.glmnet")
  glm_learner <- Lrnr_glm$new()
  screen_and_glm <- Pipeline$new(slscreener, glm_learner)
  
  
  stack = Stack$new(
    lrn_glm,
    lrn_glmnet,
    lrn_dbarts,
    lrn_nnet,
    screen_and_glm,
    lrn_mean
  )
  
  task = make_sl3_Task(
    data = data_temp,
    outcome = "event_status",
    covariates = covariates,
    weights = "ipcw",
    outcome_type = variable_type("binomial"),
    id = "COU1A",
    folds = make_folds(
      n = nrow(data_temp),
      cluster_ids = data_temp$COU1A,
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

# In the current version of the code, we use the last two measurements of the
# surrogates for estimating the surrogate index.
surrogate_1 = c(
  "SBPAVG3N_1",
  "SBPAVG3N_2",
  "SBPAVG3N_3",
  "DBPAVG3N_1",
  "DBPAVG3N_2",
  "DBPAVG3N_3",
  "STNPLS1N_1",
  "STNPLS1N_2",
  "STNPLS1N_3"
)
surrogate_2 = c("SBPAVG3N_1",
                "SBPAVG3N_2",
                "DBPAVG3N_1",
                "DBPAVG3N_2",
                "STNPLS1N_1",
                "STNPLS1N_2")

sl_models_tbl = expand_grid(
  landmark_time_model = landmark_times[-1],
  endpoint_model = endpoints,
  surrogates = list(surrogate_1)
) %>%
  bind_rows(
    expand_grid(
      landmark_time_model = landmark_times[1],
      endpoint_model = endpoints,
      surrogates = list(surrogate_2)
    )
  ) 

sl_models_tbl$fitted_model = future_pmap(
  .l = list(
    endpoint = sl_models_tbl$endpoint_model,
    landmark = sl_models_tbl$landmark_time_model,
    predictors_chr = sl_models_tbl$surrogates
  ),
  .f = sl_fitter,
  .options = furrr_options(packages = "splines", seed = TRUE)
)

# Save Results -----------------------------------------------------------

# Save the data in landmark-prediction format and the fitted landmark cox
# models.
saveRDS(full_data_new_endpoints_landmark_tbl, file =
          "additional-application/R/meta-analytic/intermediate-objects/full_data_new_endpoints_landmark_tbl.rds")
saveRDS(landmark_cox_models_tbl, file =
          "additional-application/R/meta-analytic/intermediate-objects/landmark_cox_models_tbl.rds")
saveRDS(sl_models_tbl,
        file = paste0("additional-application/R/meta-analytic/intermediate-objects/sl_models_tbl.rds"))
