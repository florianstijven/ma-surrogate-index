#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# This R script takes only one argument: the variable used as clustering unit.
clustering_variable_chr = args[1]

# Load packages.
library(tidyverse)
library(survival)

pred_time = 2.5 * 365.25

# Data Preparation --------------------------------------------------------

full_data_tbl_wide = readRDS("R/data/intermediate-objects/full_data_tbl_wide.rds") %>%
  rename(clustering_variable = any_of(clustering_variable_chr))
landmark_predictions_tbl = readRDS(file = "R/meta-analytic/intermediate-objects/landmark_predictions_tbl.rds")
full_data_tbl_long = readRDS("R/data/intermediate-objects/full_data_tbl_long.rds") %>%
  rename(clustering_variable = any_of(clustering_variable_chr))

endpoints = c("CORON", "CVMB_DTH", "CVMB", "CVMM")

# Convert full long data set to long format with one row per potential surrogate
# endpoint. This data set is then joined to the data set with the landmark
# predictions.
surrogate_endpoints_long_tbl = full_data_tbl_long %>%
  pivot_longer(
    cols = c("SBPAVG3N", "DBPAVG3N", "STNPLS1N"),
    names_to = "surrogate",
    values_to = "value"
  ) %>%
  mutate(
    surrogate = case_when(
      stringr::str_detect(surrogate, "SBP") ~ "Systolic Blood Pressure",
      stringr::str_detect(surrogate, "DBP") ~ "Diastolic Blood Pressure",
      stringr::str_detect(surrogate, "STNPLS") ~ "Sitting Pulse"
    )
  ) %>%
  select(SID1A, surrogate, value, months, clustering_variable, TRTREG1C, SEX1C, REN1) %>%
  rename(
    landmark_time = months
  ) %>%
  filter(landmark_time >= 1) %>%
  bind_rows(
    landmark_predictions_tbl %>%
      pivot_longer(
        cols = c("predicted_prob_cox", "predicted_prob_sl"),
        names_to = "model_type",
        values_to = "predicted_prob"
      ) %>%
      mutate(
        surrogate = ifelse(model_type == "predicted_prob_cox", "Pred. S (Cox)", "Pred. S (SL)"),
        landmark_time = 12 * landmark_time / 365.25
      ) %>%
      rename(value = predicted_prob) %>%
      select(SID1A, surrogate, value, landmark_time, endpoint) %>%
      # Add clustering variable and treatment indicator for every patient.
      left_join(full_data_tbl_wide %>%
                  select(SID1A, clustering_variable, TRTREG1C))
  ) %>%
  mutate(landmark_time = as.integer(landmark_time))

# We only look at measurement times up to 24 months.
surrogate_endpoints_long_tbl = surrogate_endpoints_long_tbl %>%
  filter(landmark_time <= 24)

# Convert the wide data set into a long format with one row per time-to-event
# endpoint.
event_time_data_tbl_long = full_data_tbl_wide %>%
  pivot_longer(cols = any_of(contains(endpoints))) %>%
  mutate(
    variable_type = ifelse(stringr::str_starts(name, "C4"), "censored", "time"),
    endpoint = stringr::str_sub(name, start = 3)
  ) %>%
  select(-name) %>%
  pivot_wider(values_from = "value", names_from = "variable_type")

# Determine the number of events for every cluster-endpoint combination.
n_events_tbl = event_time_data_tbl_long %>%
  group_by(clustering_variable, endpoint) %>%
  summarize(n_events = sum(1 - censored, na.rm = TRUE))

# Analysis ----------------------------------------------------------------

## Center-Level Treatment Effects -----------------------------------------

# Function to estimate cox model

# Estimate the treatment effect on the time-to-event endpoints with the most
# events.
cluster_models_tbl_long <- event_time_data_tbl_long %>%
  left_join(n_events_tbl) %>%
  filter(n_events >= 5) %>%
  group_by(endpoint, clustering_variable) %>%
  summarize(
    surv_diff = {
      # Fit KM curves by treatment
      fit <- survfit(Surv(time, I(1 - censored)) ~ TRTREG1C, 
                     data = pick(everything()))
      
      # Extract estimates at time t0
      s <- summary(fit, times = pred_time, extend = TRUE)
      survs <- s$surv
      ses   <- s$std.err
      
      # Make sure order matches treatment factor levels
      # Here: difference = treatment (level 2) - control (level 1)
      delta <- survs[2] - survs[1]
      se_delta <- sqrt(ses[1]^2 + ses[2]^2)  # assuming independent arms
      
      list(
        surv_treat = survs[2],
        surv_control = survs[1],
        delta = delta,
        se_delta = se_delta
      )
    } %>%
      list(),
    cox_fit = list(coxph(
      formula = Surv(time, I(1 - censored)) ~ TRTREG1C,
      data = pick(everything())
    )),
    .groups = "drop"
  ) %>%
  tidyr::unnest_wider(surv_diff)


# Extract estimated treatment effects and SEs in the Cox models.
cluster_models_tbl_long = cluster_models_tbl_long %>%
  ungroup() %>%
  mutate(
    cox_coef = purrr::map_dbl(
      .x = cox_fit,
      .f = function(fit) {
        if (is.na(fit[[1]]))
          return(NA)
        else
          return(coef(fit)[1])
      }
    ),
    cox_coef_se = purrr::map_dbl(
      .x = cox_fit,
      .f = function(fit) {
        if (is.na(fit[[1]]))
          return(NA)
        else
          return(sqrt(vcov(fit)[1, 1]))
      }
    )
  )

# Function to fit linear model that handles special cases (i.e., when a cluster
# contains only one treatment level).
fit_lm = function(value, trt) {
  if (length(unique(trt[!is.na(value)])) <= 1) {
    return(NA)
  }
  else {
    return(lm(value~trt))
  }
}
# Estimate the treatment effect on the blood-pressure measurements at 1, 2, 3,
# and 6 months. Also for the sitting pulse at the same moments.
cluster_models_tbl_full = cluster_models_tbl_long %>%
  left_join(
    # For the surrogate endpoints, we can ignore the type of endpoint for joining.
    surrogate_endpoints_long_tbl %>%
      filter(is.na(endpoint)) %>%
      select(-endpoint) %>%
      group_by(clustering_variable, surrogate, landmark_time) %>%
      summarize(lm_fit = list(fit_lm(value, TRTREG1C))) %>%
      ungroup(),
    relationship = "many-to-many"
  ) %>%
  bind_rows(
    cluster_models_tbl_long %>%
      left_join(
        # For the predicted surrogates, we should take the type of endpoint into
        # account for joining (predictions are endpoint specific).
        surrogate_endpoints_long_tbl %>%
          filter(!is.na(endpoint)) %>%
          group_by(clustering_variable, surrogate, landmark_time, endpoint) %>%
          summarize(lm_fit = list(fit_lm(value, TRTREG1C))) %>%
          ungroup(),
        by = c("clustering_variable", "endpoint"),
        relationship = "one-to-many"
      )
  )

# Extract estimated treatment effects from the fitted linear models.
cluster_models_tbl_full = cluster_models_tbl_full %>%
  mutate(lm_coef = purrr::map_dbl(
    .x = lm_fit,
    .f = ~ {# If the linear model couldn't be fit, return NA.
      if (is.na(.x)[1]) return(NA)
      # Otherwise, return the estimated treatment effect.
      coef(.x)[2]}
  ),
  lm_coef_se = purrr::map_dbl(
    .x = lm_fit,
    .f = ~ { # If the linear model couldn't be fit, return NA.
      if (is.na(.x)[1]) return(NA)
      # Otherwise, return the standard error of the estimated treatment effect.
      sqrt(vcov(.x)[2, 2])}
  ))
  
# Drop variables that are no longer needed.
cluster_models_tbl_full = cluster_models_tbl_full %>%
  select(-cox_fit, -lm_fit)

## Plots ------------------------------------------------------------------
# Recode month variable to a character variable for nicer plots.
cluster_models_tbl_full = cluster_models_tbl_full %>%
  mutate(landmark_time_chr = paste0("Month ", landmark_time))

# Define function to make desired plots.
ma_plot_function = function(plot_type, endpoint, add_predicted_surrogate, months) {
  temp_data = cluster_models_tbl_full %>%
    filter(abs(cox_coef) < 5, landmark_time %in% months, 
           endpoint == .env$endpoint) %>% 
    mutate(surrogate = factor(surrogate))
  if (!add_predicted_surrogate) {
    temp_data = temp_data %>%
      filter(!(surrogate %in% c("Pred. S (Cox)", "Pred. S (SL)")))
  }
  p = temp_data %>%
    # Reorder the factor levels for the surrogate variable.
    mutate(
      surrogate = forcats::lvls_reorder(surrogate, c(1, 5, 4, 3, 2)),
      surrogate = forcats::fct_recode(
        .f = surrogate,
        "Systolic BP" = "Systolic Blood Pressure",
        "Diastolic BP" = "Diastolic Blood Pressure"
      )
    ) %>%
    ggplot(aes(x = lm_coef, y = cox_coef)) +
    geom_point() + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
  if (plot_type == "scatter") {
    p = p + geom_smooth(method = "lm")
  }
  else if (plot_type == "error") {
    p = p +
      geom_linerange(aes(xmin = lm_coef - lm_coef_se, xmax = lm_coef + lm_coef_se),
                     alpha = 0.2) +
      geom_linerange(
        aes(
          ymin = cox_coef - cox_coef_se,
          ymax = cox_coef + cox_coef_se
        ),
        alpha = 0.2
      )
  }
  p + facet_grid(landmark_time_chr ~ surrogate, scales = "free_x") +
    xlab("Estimated Mean Difference") +
    ylab("Estimated log HR") +
    ggtitle(endpoint)
  
  if (add_predicted_surrogate) {
    plot_type = paste0(plot_type, "-pred_surr")
  }
  
  ggsave(
    filename = paste0(
      "figures/meta-analytic/",
      clustering_variable_chr,
      "/",
      endpoint,
      "-trial-level-",
      plot_type,
      "-plots-",
      months[1],
      "_",
      months[length(months)],
      ".pdf"
    ),
    device = "pdf",
    width = double_width,
    height = double_height,
    units = "cm",
    dpi = res
  )
  return(1)
}

expand_grid(
  plot_type = c("scatter", "error"),
  endpoint = endpoints,
  add_predicted_surrogate = c(FALSE, TRUE),
  months = list(c(1, 2, 3, 6), c(12, 18, 24))
) %>%
  rowwise() %>%
  summarise(temp = ma_plot_function(plot_type[1], endpoint[1], add_predicted_surrogate[1], months))

ma_plot_function_surv_diff = function(plot_type, endpoint, add_predicted_surrogate, months) {
  temp_data = cluster_models_tbl_full %>%
    filter(landmark_time %in% months, 
           endpoint == .env$endpoint) %>% 
    mutate(surrogate = factor(surrogate))
  if (!add_predicted_surrogate) {
    temp_data = temp_data %>%
      filter(!(surrogate %in% c("Pred. S (Cox)", "Pred. S (SL)")))
  }
  p = temp_data %>%
    # Reorder the factor levels for the surrogate variable.
    mutate(
      surrogate = forcats::lvls_reorder(surrogate, c(1, 5, 4, 3, 2)),
      surrogate = forcats::fct_recode(
        .f = surrogate,
        "Systolic BP" = "Systolic Blood Pressure",
        "Diastolic BP" = "Diastolic Blood Pressure"
      )
    ) %>%
    ggplot(aes(x = lm_coef, y = delta)) +
    geom_point() + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
  if (plot_type == "scatter") {
    p = p + geom_smooth(method = "lm")
  }
  else if (plot_type == "error") {
    p = p +
      geom_linerange(aes(xmin = lm_coef - lm_coef_se, xmax = lm_coef + lm_coef_se),
                     alpha = 0.2) +
      geom_linerange(
        aes(
          ymin = delta - se_delta,
          ymax = delta + se_delta
        ),
        alpha = 0.2
      )
  }
  p + facet_grid(landmark_time_chr ~ surrogate, scales = "free_x") +
    xlab("Estimated Mean Difference") +
    ylab("Estimated Diff in Surv Prob") +
    ggtitle(endpoint)
  
  if (add_predicted_surrogate) {
    plot_type = paste0(plot_type, "-pred_surr")
  }
  
  ggsave(
    filename = paste0(
      "figures/meta-analytic/",
      clustering_variable_chr,
      "/",
      endpoint,
      "-survdiff",
      "-trial-level-",
      plot_type,
      "-plots-",
      months[1],
      "_",
      months[length(months)],
      ".pdf"
    ),
    device = "pdf",
    width = double_width,
    height = double_height,
    units = "cm",
    dpi = res
  )
  return(1)
}

expand_grid(
  plot_type = c("scatter", "error"),
  endpoint = endpoints,
  add_predicted_surrogate = c(FALSE, TRUE),
  months = list(c(1, 2, 3, 6), c(12, 18, 24))
) %>%
  rowwise() %>%
  summarise(temp = ma_plot_function_surv_diff(plot_type[1], endpoint[1], add_predicted_surrogate[1], months))


# # Fit linear regression models for the three surrogates at each time point.
# ma_models_per_month = cluster_models_tbl_full %>%
#   filter(surrogate != "Prediction Surrogate") %>%
#   select(-lm_coef_se) %>%
#   pivot_wider(
#     names_from = surrogate,
#     values_from = lm_coef
#   ) %>%
#   group_by(landmark_time, endpoint) %>%
#   summarize(
#     lm_fit = list(
#       lm(
#         cox_coef ~ `Systolic Blood Pressure` + `Diastolic Blood Pressure` + `Sitting Pulse`,
#         data = pick(everything())
#       )
#     )
#   )
# 
# sink(file = paste0("tables/meta-analytic/", clustering_variable_chr, "/ma-models-per-month-rsquared.txt"))
# cat(
# "Estimated R-squared values for linear models with the estimated log-hazard ratios
# for CVMM, CVMB, CORON, and CVMB_DTH as outcome and the treatment effects on diastolic and 
# systolic blood pressure, and sitting pulse as predictors. Each row corresponds to 
# distinct months when these variables were measured. \n"
# )
# ma_models_per_month %>%
#   pivot_wider(
#     names_from = "endpoint",
#     values_from = lm_fit
#   ) %>%
#   ungroup() %>%
#   filter(landmark_time %in% c(1, 2, 3, 6)) %>%
#   rowwise(landmark_time) %>%
#   summarize(across(
#     .cols = 1:4,
#     .fns = ~ summary(.x)$r.squared
#   )) %>%
#   ungroup()
# sink()

# # Fit linear regression models for the three surrogates at all time points up to
# # six months.
# ma_models_up_to_six_months = cluster_models_tbl_long %>%
#   pivot_wider(
#     names_from = c(endpoint, month),
#     values_from = lm_coef,
#     id_cols = c(any_of(clustering_variable_chr), contains("coxph"))
#   ) %>%
#   rename(CVMM = coxph_model_CVMM_coef, CORON = coxph_model_CORON_coef, CVMB = coxph_model_CVMB_coef, CVMB_DTH = coxph_model_CVMB_DTH_coef) %>%
#   summarize(CVMM_fit = list(
#     lm(
#       CVMM ~ `Systolic Blood Pressure_1` + `Diastolic Blood Pressure_1` + `Sitting Pulse_1` +
#         `Systolic Blood Pressure_2` + `Diastolic Blood Pressure_2` + `Sitting Pulse_2` +
#         `Systolic Blood Pressure_3` + `Diastolic Blood Pressure_3` + `Sitting Pulse_3` +
#         `Systolic Blood Pressure_6` + `Diastolic Blood Pressure_6` + `Sitting Pulse_6`,
#       data = pick(everything())
#     )
#   ),
#   CORON_fit = list(
#     lm(
#       CORON ~ `Systolic Blood Pressure_1` + `Diastolic Blood Pressure_1` + `Sitting Pulse_1` +
#         `Systolic Blood Pressure_2` + `Diastolic Blood Pressure_2` + `Sitting Pulse_2` +
#         `Systolic Blood Pressure_3` + `Diastolic Blood Pressure_3` + `Sitting Pulse_3` +
#         `Systolic Blood Pressure_6` + `Diastolic Blood Pressure_6` + `Sitting Pulse_6`,
#       data = pick(everything())
#     )
#   ),
#   CVMB_fit = list(
#     lm(
#       CVMB ~ `Systolic Blood Pressure_1` + `Diastolic Blood Pressure_1` + `Sitting Pulse_1` +
#         `Systolic Blood Pressure_2` + `Diastolic Blood Pressure_2` + `Sitting Pulse_2` +
#         `Systolic Blood Pressure_3` + `Diastolic Blood Pressure_3` + `Sitting Pulse_3` +
#         `Systolic Blood Pressure_6` + `Diastolic Blood Pressure_6` + `Sitting Pulse_6`,
#       data = pick(everything())
#     )
#   ),
#   CVMB_DTH_fit = list(
#     lm(
#       CVMB_DTH ~ `Systolic Blood Pressure_1` + `Diastolic Blood Pressure_1` + `Sitting Pulse_1` +
#         `Systolic Blood Pressure_2` + `Diastolic Blood Pressure_2` + `Sitting Pulse_2` +
#         `Systolic Blood Pressure_3` + `Diastolic Blood Pressure_3` + `Sitting Pulse_3` +
#         `Systolic Blood Pressure_6` + `Diastolic Blood Pressure_6` + `Sitting Pulse_6`,
#       data = pick(everything())
#     )
#   ))
# 
# sink(file = paste0("tables/meta-analytic/", clustering_variable_chr, "/ma-models-up-to-six-months.txt"))
# cat(
# "Estimated R-squared values for linear models with the estimated log-hazard ratios
# for CVMM, CVMB, CORON, and CVMB_DTH as outcome and the treatment effects (for month 1, 2, 3,
# and 6) on diastolic and systolic blood pressure, and sitting pulse as predictors. \n"
# )
# ma_models_up_to_six_months %>%
#   rowwise() %>%
#   summarize(across(
#     .cols = 1:3,
#     .fns = ~ summary(.x)$r.squared
#   ))
# sink()


# Tests for Heterogeneity -------------------------------------------------
library(meta)
# Fit univariate meta-analytic models.
fit_ma = function(x, y) {
  metagen(x, y)
}

univariate_ma_models = cluster_models_tbl_full %>%
  # The estimated log-hazard ratios are replicated across rows (e.g., due to
  # multiple landmark times). We need to select the unique log-hazard ratio
  # estimates (i.e., one estimate per endpoint-cluster combination).
  group_by(endpoint, clustering_variable) %>%
  slice_head() %>%
  ungroup() %>%
  group_by(endpoint) %>%
  summarize(fit_ma = list(
    fit_ma(cox_coef, cox_coef_se)
  ))

# Print the results for the tests of heterogeneity.
sink(file = paste0("tables/meta-analytic/", clustering_variable_chr, "/heterogeneity-tests.txt"))
cat(
"Q-test for heterogeneity in the center-level log HR estimates on CVMM,
CORON, CVMB, and CVMB_DTH. \n"
    )
univariate_ma_models %>%
  rowwise(endpoint) %>%
  summarise(pval.Q = fit_ma[["pval.Q"]])
cat(
"
Estimated between-center variance in the center-level log HR estimated on CVMM,
CORON, CVMB, and CVMB_DTH. \n"
)
univariate_ma_models %>%
  rowwise(endpoint) %>%
  summarise(tau2 = fit_ma[["tau2"]])
sink()

