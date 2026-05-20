# Function to predict X-Year survival probability given landmark cox model or
# superlearner, landmark time, and a new data set.
landmark_prediction = function(newdata, model, landmark_time, time, model_type = "cox") {
  if (model_type == "cox") {
    landmark_predictions_tbl = bind_cols(
      tibble(predicted_prob_landmark = predict_surv_prob_cox(newdata, model, time)),
      newdata %>%
        select(SID1A, time, censored, TRTREG1C)
    )
  } else if (model_type == "sl") {
    landmark_predictions_tbl = bind_cols(
      tibble(predicted_prob_landmark = 1 - sl_prediction_f(newdata, model, time)),
      newdata %>%
        select(SID1A, time, censored, TRTREG1C)
    )
  } else {
    stop("The model type is not recognized.")
  }
  
  # For each landmark prediction, we replace the predicted survival probability
  # with a zero if the event was observed by the landmark time.
  landmark_predictions_tbl = landmark_predictions_tbl %>%
    mutate(
      predicted_prob = ifelse((time <= landmark_time) &
                                censored == 0,
                              0,
                              predicted_prob_landmark
      )
    )
  # If the patient was censored before the landmark time, then we replace the
  # predicted probability with a NA. These entries are also indicated by a
  # separate variable.
  landmark_predictions_tbl = landmark_predictions_tbl %>%
    mutate(
      predicted_prob = ifelse((time <= landmark_time) &
                                censored == 1, NA, predicted_prob),
      censored_before_landmark = (time <= landmark_time) & censored == 1
    )
  # Return the complete landmark predictions.
  return(landmark_predictions_tbl %>% pull(predicted_prob))
}

# Function to predict the survival probability at a given time point for a
# fitted Cox PH model.
predict_surv_prob_cox = function(newdata, cox_model, time) {
  exp(
    -1 * predict(
      cox_model,
      newdata = newdata %>%
        mutate(time = .env$time, censored = 1),
      type = "expected"
    )
  )
}

# Functions that does predictions on new data based on a fitted SuperLearner.
# The `time` argument has no effect but is used for consistency with the other
# functions.
sl_prediction_f = function(newdata, sl_fit, time) {
  # Define the covariates that will be included as predictors.
  covariates = sl_fit$training_task$nodes$covariates
  
  # Add row number, which will be need later on.
  newdata = newdata %>%
    mutate(row_number = row_number()) %>%
    mutate(event_status = 1,
           ipcw = 1)
  
  # Split the data set into one with one with no missing value and one with missing values.
  newdata_no_missing = newdata %>%
    filter(if_all(all_of(covariates), ~ !is.na(.)))
  newdata_missing = newdata %>%
    filter(if_any(all_of(covariates), is.na))
  
  # Do predictions for the rows with no missing values.
  prediction_task = make_sl3_Task(
    data = newdata_no_missing,
    covariates = covariates,
    outcome = "event_status",
    weights = "ipcw",
    id = "COU1A"
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