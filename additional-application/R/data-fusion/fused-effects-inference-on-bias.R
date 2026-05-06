args = commandArgs(trailingOnly=TRUE)
# This R script takes only one argument: the variable used as clustering unit.
clustering_variable_chr = args[1]

# Load packages.
library(tidyverse)
library(survival)
library(GGally)
library(mgcv)
library(boot)
library(parallel)

# Data Preparation --------------------------------------------------------
# Load intermediate results.
intermediate_object_location = paste0("R/data-fusion/intermediate-objects/",
                                      clustering_variable_chr,
                                      "/")
source("R/data-fusion/helper-functions.R")
landmark_cox_models_tbl = readRDS(file = paste0(intermediate_object_location, "landmark_cox_models_tbl.rds"))
full_data_new_endpoints_long = readRDS(file = "R/data/intermediate-objects/full_data_tbl_long.rds") %>%
  rename(clustering_variable = any_of(clustering_variable_chr))
  

# Number of bootstrap replications.
R = 500

# The treatment effect is estimated on the time-to-event endpoint. We look at
# the 1, 2, 3, and 4-year survival probabilities. Other treatment effect
# measures (e.g., RMST) can be added later.
times = 365.25 * 1:4 

# Bootstrap Helper Functions ----------------------------------------------

# To do inference about the bias in the predicted treatment effect, we first
# need the within-trial correlation between the estimated and predicted
# treatment effects. This represents the correlation due to sampling
# variability. We compute it through non-parametric bootstrapping and ignore the
# uncertainty in the prediction rule (i.e., the prediction rule is treated as
# fixed).

# Function to extract estimated survival probability and corresponding SE for a
# given treatment group.
survival_prob = function(surv_fit, times) {
  surv_fit_summary = summary(surv_fit, times = times)
  surv_fit_summary[c("time", "surv", "strata", "std.err")] %>%
    as_tibble() %>%
    mutate(TRTREG1C = stringr::str_sub(strata, start = -1, end = -1)) %>%
    select(-strata)
}

# Function that estimates the difference in survival probabilities at various
# time points.
surv_diff = function(data, pred_time, endpoint) {
  data = data %>% 
    group_by(SID1A) %>%
    slice_min(days) %>%
    ungroup()
  if (endpoint == "DTH") {
    data = data %>%
      rename(time = T2DTH,
             censored = C4DTH)
  }
  else if (endpoint == "CVMB_DTH") {
    data = data %>%
      rename(time = T2CVMB_DTH,
             censored = C4CVMB_DTH)
  }
  # In the data argument, the observed time is time and the censoring indicator
  # is censored.
  cox_fit = survfit(Surv(time, I(1 - censored))~strata(TRTREG1C), data)
  surv_diff_estimates = survival_prob(cox_fit, pred_time) %>%
    pivot_wider(
      names_from = c("TRTREG1C"),
      values_from = c("surv", "std.err")) %>%
    mutate(surv_diff = surv_A - surv_B) %>%
    pull(surv_diff)
  names(surv_diff_estimates) = paste0("surv_diff-", pred_time)
  return(surv_diff_estimates)
}

# Function that estimates the predicted treatment effect at various time points
# given a prediction rule.
pred_surv_diff = function(data, pred_time, pred_rule, landmark_time, endpoint) {
  if (endpoint == "DTH") {
    data = data %>%
      rename(time = T2DTH,
             censored = C4DTH)
  }
  else if (endpoint == "CVMB_DTH") {
    data = data %>%
      rename(time = T2CVMB_DTH,
             censored = C4CVMB_DTH)
  }
  # Create new data set in the landmark-prediction format.
  full_data_new_endpoints_landmark_tbl = data %>%
    filter(FALSE) %>%
    mutate(timing = numeric(), 
           landmark_time = numeric()) %>%
    select(-VISNAM1A, days, months)
  # Given the landmark time, find the times of the two most recent measurements.
  days_measurements = unique(data$days) %>%
    na.omit()
  latest_days = sort(days_measurements[days_measurements <= landmark_time], decreasing = TRUE)[1:2]
  full_data_new_endpoints_landmark_tbl = bind_rows(
    full_data_new_endpoints_landmark_tbl,
    data %>%
      filter(days %in% latest_days) %>%
      group_by(SID1A) %>%
      arrange(-days, .by_group = TRUE) %>%
      mutate(timing = row_number()) %>%
      select(-days, -months, -VISNAM1A) %>%
      pivot_wider(values_from = c(SBPAVG3N, DBPAVG3N, STNPLS1N),
                  names_from = timing) %>%
      mutate(landmark_time = landmark_time)
  )
  pred_surv_diff_estimates = tibble(pred_time) %>%
    cross_join(full_data_new_endpoints_landmark_tbl) %>%
    group_by(.data$pred_time) %>%
    reframe(tibble(
      pred_surv_prob = pred_rule(
        newdata = pick(everything()),
        pred_time = pred_time[1]
      ),
      TRTREG1C, 
      SID1A
    )) %>%
    ungroup() %>%
    group_by(.data$pred_time, TRTREG1C) %>%
    summarise(mean_pred_prob = mean(pred_surv_prob, na.rm = TRUE)) %>%
    pivot_wider(names_from = c("TRTREG1C"),
                values_from = c("mean_pred_prob")) %>%
    mutate(pred_surv_diff = A - B) %>%
    pull(pred_surv_diff)
  names(pred_surv_diff_estimates) = paste0("pred-surv_diff-", pred_time, "-", landmark_time)
  return(pred_surv_diff_estimates)
}

pred_rule_constructor = function(cox_model, landmark_time){
  pred_rule = function(newdata, pred_time) {
    landmark_prediction(newdata, cox_model, landmark_time, pred_time)
  }
}

estimate_surv_effects = function(data, pred_time, landmark_time, pred_rule, endpoint){
  est_surv_diff = NA
  try(expr = {est_surv_diff = surv_diff(data = data, pred_time = pred_time, endpoint = endpoint)})
  est_pred_surv_diff = pred_surv_diff(
    data = data,
    pred_time = pred_time,
    pred_rule = pred_rule,
    landmark_time = landmark_time,
    endpoint =  endpoint
  )
  return(c(est_surv_diff, est_pred_surv_diff))
}

# Function to generate bootstrap resamples of the original data.
ran.gen_function = function(data, mle) {
  # Vector with all subject ids
  patient_ids = unique(data$SID1A)
  # Resample patient ids with replacement.
  patient_ids_sampled = sample(patient_ids, size = length(patient_ids), replace = TRUE)
  return(tibble(SID1A = patient_ids_sampled) %>%
           left_join(data, relationship = "many-to-many"))
}

# Function that performs a non-parametric bootstrap for the estimated and
# predicted treatment effects, given a (i) data set, (ii) a landmark time, (iii)
# prediction window, and (iv) a prediction rule..
np_bootstrap_surv_effects = function(data, landmark_time, pred_time, pred_rule, endpoint, R) {
  np_boot_object = boot(
    data = data,
    statistic = estimate_surv_effects,
    sim = "parametric",
    ran.gen = ran.gen_function,
    R = R,
    pred_time = pred_time,
    landmark_time = landmark_time,
    pred_rule = pred_rule,
    endpoint = endpoint
  )
  return(np_boot_object)
}



# Bootstrap ---------------------------------------------------------------

# Start from big data set and for each cluster, compute the required statistics. 
landmark_prediction_settings_tbl = expand_grid(pred_time = times, endpoint = c("DTH", "CVMB_DTH")) %>%
  left_join(landmark_cox_models_tbl %>%
              filter(model_type == "full"),
            relationship = "many-to-many")  %>%
  filter(landmark_time < pred_time) %>%
  left_join(
    full_data_new_endpoints_long %>%
      group_by(clustering_variable) %>%
      summarize(data = list(pick(everything(
        
      ))))
  )

cl = makePSOCKcluster(detectCores() - 2)
clusterExport(
  cl,
  varlist = c(
    "survival_prob",
    "surv_diff",
    "pred_surv_diff",
    "pred_rule_constructor",
    "estimate_surv_effects",
    "ran.gen_function",
    "np_bootstrap_surv_effects",
    "landmark_prediction"
  )
)
clusterExport(
  cl,
  varlist = c(
    "R"
  )
)
clusterEvalQ(cl, {library("boot"); library("survival"); library("tidyverse")})
np_boot_object_list = clusterMap(
  cl = cl,
  data = landmark_prediction_settings_tbl$data,
  cox_model = landmark_prediction_settings_tbl$cox_model,
  landmark_time = landmark_prediction_settings_tbl$landmark_time,
  pred_time = landmark_prediction_settings_tbl$pred_time,
  endpoint = landmark_prediction_settings_tbl$endpoint
  ,
  fun = function(data,
                 cox_model,
                 landmark_time,
                 pred_time,
                 endpoint) {
    pred_rule = pred_rule_constructor(cox_model, landmark_time)
    np_boot_object = np_bootstrap_surv_effects(
      data = data,
      landmark_time = landmark_time,
      pred_time = pred_time,
      pred_rule = pred_rule,
      endpoint = endpoint,
      R = R
    )
    return(np_boot_object)
  }
)
stopCluster(cl)

bias_inference_tbl = landmark_prediction_settings_tbl %>%
  select(-cox_model, -data) %>%
  mutate(
    np_boot_object = np_boot_object_list,
    vcov_matrix = lapply(
      X = np_boot_object,
      FUN = function(x) {
        var(x$t, na.rm = TRUE)
      }
    ),
    est_surv_diff = lapply(
      X = np_boot_object,
      FUN = function(x) {
        x$t0[1]
      }
    ) %>% as.numeric(),
    pred_surv_diff = lapply(
      X = np_boot_object,
      FUN = function(x) {
        x$t0[2]
      }
    ) %>% as.numeric(),
    est_diff_se = lapply(
      X = vcov_matrix,
      FUN = function(x) {
        sqrt(x[1, 1] + x[2, 2] - 2 * x[1, 2])
      }
    ) %>% as.numeric(),
    perc_ci = lapply(
      X = np_boot_object,
      FUN = function(x) {
        diff_est = x$t[, 1] - x$t[, 2]
        ci = quantile(diff_est, probs = c(0.025, 0.975), na.rm = TRUE)
        return(ci)
      }
    )
  )

# Save the data set with inferences on the bias for the predicted treatment
# effects.
saveRDS(
  bias_inference_tbl,
  file = paste0(
    "R/data-fusion/intermediate-objects/",
    clustering_variable_chr,
    "/bias_inference_tbl.rds"
  )
)