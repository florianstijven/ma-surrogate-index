#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load packages.
library(tidyverse)
library(survival)
library(sl3)

# Data Preparation --------------------------------------------------------

source("R/data-fusion/helper-functions.R")
landmark_cox_models_tbl = readRDS(file = "R/meta-analytic/intermediate-objects/landmark_cox_models_tbl.rds")
sl_models_tbl = readRDS(file = "R/meta-analytic/intermediate-objects/sl_models_tbl.rds")
full_data_landmark_tbl = readRDS(file = "R/meta-analytic/intermediate-objects/full_data_new_endpoints_landmark_tbl.rds")

landmark_times = 365.25 * c(1 / 12, 2 / 12, 3 / 12, 4 / 12, 5 / 12, 0.5, 1, 1.5, 2)


# Predicted Surrogates ----------------------------------------------------

# For patients that are still under observations at the landmark time, predict
# the probability of X-year survival. 

# For each landmark_time-endpoint combination, we let the data be a single
# element in the tibble. This condensed data set is then joined with the tibble
# of models.
landmark_predictions_tbl = full_data_landmark_tbl %>%
  group_by(landmark_time, endpoint) %>%
  summarize(data_set = list(pick(everything()))) %>%
  left_join(
    landmark_cox_models_tbl %>%
      filter(model_type == "full"),
    by = c("landmark_time", "endpoint")
  ) %>%
  left_join(
    sl_models_tbl %>%
      rename(
        "landmark_time" = landmark_time_model,
        "endpoint" = endpoint_model,
        "sl_model" = fitted_model
      ) %>%
      select(-surrogates),
    by = c("landmark_time", "endpoint")
  ) %>%
  # We only look at the 30 months survival probability.
  mutate(times = 2.5 * 365.25) %>%
  # Only keep entries where the landmark time is smaller than the X-year time
  # point.
  filter(landmark_time < times) %>%
  # For each row, which contains a full data set and the fitted gam, the
  # predictions are computed and we return to a normal tibble with one row per
  # data entry. To reduce memory usage, we only keep the variables that are
  # needed further on.
  group_by(landmark_time, endpoint) %>%
  reframe(bind_cols(
    tibble(
      predicted_prob_cox = landmark_prediction(data_set[[1]], cox_model[[1]], landmark_time[1], times[1], model_type = "cox"),
      predicted_prob_sl = landmark_prediction(data_set[[1]], sl_model[[1]], landmark_time[1], times[1], model_type = "sl")
    ),
    data_set[[1]] %>%
      select(SID1A, time, censored, TRTREG1C)
  )) %>%
  mutate(censored_before_landmark = (time <= landmark_time) &
           censored == 1)



# Saving Results -----------------------------------------------------------

saveRDS(landmark_predictions_tbl,
        file = "additional-application/R/meta-analytic/intermediate-objects/landmark_predictions_tbl.rds")
