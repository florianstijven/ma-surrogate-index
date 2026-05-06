# Load packages.
library(tidyverse)
library(survival)

# Data Preparation --------------------------------------------------------

full_data_tbl_wide = readRDS("R/data/intermediate-objects/full_data_tbl_wide.rds")
full_data_tbl_long = readRDS("R/data/intermediate-objects/full_data_tbl_long.rds")

# Vector of the censoring indicator names for all time-to-event variables
censoring_indicator_names = c("C4DTH",
                              "C4CVDTH",
                              "C4MI",
                              "C4STRK",
                              "C4USA",
                              "C4CORON",
                              "C4CVMB",
                              "C4CVMM", 
                              "C4CVMB_DTH")
# Vector of time-to-event variable names for all time-to-event variables
event_time_names = c("T2DTH",
                     "T2CVDTH",
                     "T2MI",
                     "T2STRK",
                     "T2USA",
                     "T2CORON",
                     "T2CVMB",
                     "T2CVMM",
                     "T2CVMB_DTH")

# Helper Functions --------------------------------------------------------

## Cox Models with Time-Varying Covariates --------------------------------

### Fit Cox Models --------------------------------------------------------

# Association between blood pressure measurements and the time-to-event
# endpoints.
fit_cox_model = function(censoring_indicator, event_time) {
  # Select time-to-event information for the desired variable. The generic names
  # censoring_indicator and event_time are used within this function.
  temp = full_data_tbl_wide %>%
    rename(censoring_indicator = {{censoring_indicator}},
           event_time = {{event_time}}) %>%
    mutate(event_indicator = 1 - censoring_indicator) %>%
    select(SID1A, event_indicator, event_time)
  # Remove all missing observations from the trimmed data set.
  temp = temp %>%
    na.omit()
  # Convert the selected data into the counting process format.
  temp1 = tmerge(temp,
                 temp, 
                 id = SID1A,
                 event = event(event_time, event_indicator)
  )
  temp1 = tmerge(
    temp1,
    full_data_tbl_long %>%
      filter(!is.na(months)),
    id = SID1A,
    systolic_bp = tdc(days,SBPAVG3N),
    diastolic_bp = tdc(days, DBPAVG3N),
    stnpls = tdc(days, STNPLS1N)
  )
  cox_fit = coxph(
    formula = Surv(tstart, tstop, event) ~ splines::ns(
      systolic_bp,
      knots = c(120, 140),
      Boundary.knots = c(110, 180)
    ) + splines::ns(
      diastolic_bp,
      knots = c(80, 90),
      Boundary.knots = c(70, 100)
    ) + splines::ns(
      stnpls,
      knots = c(60, 100),
      Boundary.knots = c(40, 140)
    ),
    data = temp1,
    model = TRUE, 
    cluster = SID1A
  )
  return(cox_fit)
}

### Plot Cox Models -------------------------------------------------------

plot_effect_covariate = function(cox_fit, covariate, grid) {
  n_grid = length(grid)
  # Create a tible that sets the covariate not in [covariate] to arbitrary
  # values and sets [covariate] to the grid.
  grid_tbl = tibble(
    systolic_bp = rep(120, n_grid),
    diastolic_bp = rep(80, n_grid),
    stnpls = rep(70, n_grid)
  ) %>%
    mutate(covariate = grid) %>%
    select(!all_of(covariate))
  colnames(grid_tbl)[colnames(grid_tbl) == "covariate"] = covariate
  predict_lp = predict(
    cox_fit,
    newdata = grid_tbl,
    type = "lp", 
    se.fit = TRUE
  )
  final_plot = tibble(
    linear_predictor = predict_lp$fit,
    se = predict_lp$se.fit,
    grid_tbl
  ) %>% 
    mutate(linear_predictor = linear_predictor - min(predict_lp$fit)) %>%
    ggplot(aes(x = .data[[covariate]], y = exp(linear_predictor))) + 
    geom_line() +
    scale_y_continuous(transform = "log2") +
    geom_line(aes(x = .data[[covariate]], y = exp(linear_predictor + 1.96 * se)), linetype = "dashed") +
    geom_line(aes(x = .data[[covariate]], y = exp(linear_predictor - 1.96 * se)), linetype = "dashed") 
  return(final_plot)
}

plot_effect_facet = function(cox_fit, covariates, grids, references) {
  tibble_covariates = list()
  for (i in seq_along(covariates)) {
    grid = grids[[i]]
    covariate = covariates[[i]]
    reference = references[[i]]
    n_grid = length(grid)
    # Create a tible that sets the covariate not in [covariate] to arbitrary
    # values and sets [covariate] to the grid.
    grid_tbl = tibble(
      systolic_bp = rep(120, n_grid),
      diastolic_bp = rep(80, n_grid),
      stnpls = rep(70, n_grid)
    ) %>%
      mutate(covariate = grid) %>%
      select(!all_of(covariate))
    colnames(grid_tbl)[colnames(grid_tbl) == "covariate"] = covariate
    predict_lp = predict(
      cox_fit,
      newdata = grid_tbl,
      type = "lp", 
      se.fit = TRUE
    )
    # Prediction for th reference value of the current covariate.
    reference_tbl = grid_tbl[1, ] %>%
      select(!all_of(covariate)) %>%
      mutate(covariate = {{reference}})
    colnames(reference_tbl)[colnames(reference_tbl) == "covariate"] = covariate
    predict_lp_reference = predict(
      cox_fit,
      newdata = reference_tbl,
      type = "lp", 
      se.fit = TRUE
    )
    # For further processing, we need to rename {{ covariate }} to covariate.
    colnames(grid_tbl)[colnames(grid_tbl) == covariate] = "value"
    # Add tibble for the current covariate with the linear predictor, the other
    # covariate set to arbitrary values, and the current covariate over the
    # specified grid. We also add an endpoint variable to indicator which
    # covariate is the current covariate.
    tibble_covariates[[i]] = tibble(
      linear_predictor = predict_lp$fit,
      se = predict_lp$se.fit,
      grid_tbl %>%
        select(value),
      endpoint = {{ covariate }},
      reference_lp = predict_lp_reference$fit
    ) 
  }
  
  final_plot = bind_rows(tibble_covariates) %>% 
    mutate(endpoint = case_match(
      endpoint,
      "diastolic_bp" ~ "Diastolic Blood Pressure",
      "systolic_bp" ~ "Systolic Blood Pressure",
      "stnpls" ~ "Sitting Pulse"
    )) %>%
    group_by(endpoint) %>%
    mutate(linear_predictor = linear_predictor - reference_lp) %>%
    ungroup() %>%
    ggplot(aes(x = value, y = exp(linear_predictor))) + 
    geom_line() +
    scale_y_continuous(transform = "log2") +
    geom_line(aes(x = value, y = exp(linear_predictor + 1.96 * se)), linetype = "dashed") +
    geom_line(aes(x = value, y = exp(linear_predictor - 1.96 * se)), linetype = "dashed") +
    facet_grid(~endpoint, scales = "free")
  return(final_plot)
}


## Landmark Cox Models ----------------------------------------------------



# Analyses ----------------------------------------------------------------

# Fit Cox-models with time-dependent covariates for all eight time-to-event
# outcomes. A tibble containing the fitted models is produced. This facilities
# further parallel processing.
fitted_cox_models = tibble(censoring_indicator_names, event_time_names) %>%
  mutate(
    fitted_cox_models = purrr::map2(.x = censoring_indicator_names, .y = event_time_names, .f = fit_cox_model)
  )
# Add plots that described the functional forms of the estimated effects for
# systolic bp, diastolic bp, and sitting pulse.
fitted_cox_models = fitted_cox_models %>%
  mutate(
    effect_plots_systolic_bp = purrr::map(
      .x = fitted_cox_models,
      .f = plot_effect_covariate, 
      covariate = "systolic_bp",
      grid = 90:190
    ),
    effect_plots_diastolic_bp = purrr::map(
      .x = fitted_cox_models,
      .f = plot_effect_covariate, 
      covariate = "diastolic_bp",
      grid = 50:120
    ),
    effect_plots_stnpls = purrr::map(
      .x = fitted_cox_models,
      .f = plot_effect_covariate, 
      covariate = "stnpls",
      grid = 30:140
    ),
    effect_plots = purrr::map(
      .x = fitted_cox_models,
      .f = plot_effect_facet,
      covariates = c("systolic_bp", "diastolic_bp", "stnpls"),
      grids = list(90:190, 50:120, 30:140),
      references = c(120, 80, 65)
    )
  )
# Save the effect plots to disk.
fitted_cox_models %>%
  pivot_longer(
    cols = effect_plots_systolic_bp:effect_plots_stnpls,
    values_to = "effect_plot",
    names_to = "covariate"
  ) %>%
  mutate(purrr::pmap(
    .l = list(
      event_time_name = event_time_names,
      name = covariate,
      effect_plot = effect_plot
    ),
    .f = function(event_time_name, name, effect_plot) {
      x_axis_name = case_when(
        stringr::str_detect(name, "systolic") ~ "Systolic Blood Pressure",
        stringr::str_detect(name, "diastolic") ~ "Diastolic Blood Pressure",
        stringr::str_detect(name, "stnpls") ~ "Sitting Pulse",
      )
      ggsave(
        filename = paste0("figures/associational-analyses/", event_time_name, "-", name, ".pdf"),
        plot = effect_plot + ylab("Estimated HR") + xlab(x_axis_name),
        device = "pdf",
        width = single_width,
        height = single_height,
        units = "cm",
        dpi = res
      )
    }
  ))
fitted_cox_models %>%
  mutate(purrr::map2(
    .x = effect_plots,
    .y = event_time_names,
    .f = function(effect_plot, event_time_name) {
      ggsave(
        filename = paste0(
          "figures/associational-analyses/",
          event_time_name,
          "-effect-plots",
          ".pdf"
        ),
        plot = effect_plot + 
          xlab("Value") + 
          ylab("Estimated HR") + 
          ggtitle(event_time_name) + 
          geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.5) +
          coord_cartesian(ylim = c(0.25, 16)),
        device = "pdf",
        width = double_width,
        height = double_height,
        units = "cm",
        dpi = res
      )
    }
  ))