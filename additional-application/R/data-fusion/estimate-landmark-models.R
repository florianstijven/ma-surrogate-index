args = commandArgs(trailingOnly=TRUE)
# This R script takes only one argument: the variable used as clustering unit.
clustering_variable_chr = args[1]

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

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}

# Data Preparation --------------------------------------------------------

# Add user-defined clustering variable to the data set. Currently, two options
# are implemented: (i) regions (US and Nordic) and (ii) country. 
full_data_new_endpoints = readRDS(file = "R/data/intermediate-objects/full_data_tbl_wide.rds") %>%
  rename(clustering_variable = any_of(clustering_variable_chr))

# The original time-to-event data are competing risks data. This introduces
# additional difficulties. We will therefore use alternative composite endpoints
# that are not subject to competing risks. Specifically, we will consider time
# to death, and time to cardiovascular event or death. 
full_data_new_endpoints = full_data_new_endpoints %>%
  select(any_of(
    c(
      "SID1A",
      "TRTREG1C",
      "T2DTH",
      "C4DTH",
      "T2CVMB_DTH",
      "C4CVMB_DTH",
      "clustering_variable"
    )
  )) 

full_data_new_endpoints_long = readRDS(file = "R/data/intermediate-objects/full_data_tbl_long.rds") %>%
  rename(clustering_variable = any_of(clustering_variable_chr))



# Estimate the Kaplan-Meier functions for the alternative endpoints for each
# endpoint and clustering unit separately.
cox_models_new_endpoints_tbl = full_data_new_endpoints  %>%
  pivot_longer(cols = contains("DTH")) %>%
  mutate(
    variable_type = ifelse(stringr::str_starts(name, "C4"), "censored", "time"),
    endpoint = stringr::str_sub(name, start = 3)
  ) %>%
  pivot_wider(
    id_cols = c("SID1A", "endpoint", "TRTREG1C", "clustering_variable"),
    values_from = "value",
    names_from = "variable_type"
  ) %>%
  filter(!is.na(clustering_variable)) %>%
  # For each type of time to event, the treatment effect is estimated.
  group_by(endpoint, clustering_variable) %>%
  summarise(cox_fit = list(survfit(Surv(time, I(1 - censored))~strata(TRTREG1C))))

# Landmark Prediction Models ----------------------------------------------

# For patients that remain under observations up to a certain time point, we
# estimate a landmark prediction models. This model is based on the patients at
# under observation at X months. We further use penalized cox models. 
landmark_times = 365.25 * c(1 / 12, 2 / 12, 3 / 12, 4 / 12, 5 / 12, 0.5, 1:2)

# Function to estimate landmark model.

landmark_cox_model = function(data, type = "full") {
  formula_landmark = formula(
    Surv(time, I(1 - censored)) ~
      TRTREG1C + PCI_B + SEX1C + CURSMK1C + CVSDIS5C + CVSDIS7C +
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

full_data_new_endpoints_landmark_tbl = full_data_new_endpoints_long %>%
  filter(FALSE) %>%
  mutate(timing = numeric(), 
         landmark_time = numeric()) %>%
  select(-VISNAM1A, days, months)
for (landmark_time in landmark_times) {
  # Given the landmark time, find the times of the two most recent measurements.
  days_measurements = unique(full_data_new_endpoints_long$days) %>%
    na.omit()
  latest_days = sort(days_measurements[days_measurements <= landmark_time], decreasing = TRUE)[1:2]
  full_data_new_endpoints_landmark_tbl = bind_rows(
    full_data_new_endpoints_landmark_tbl,
    full_data_new_endpoints_long %>%
      filter(days %in% latest_days) %>%
      group_by(SID1A) %>%
      arrange(-days, .by_group = TRUE) %>%
      mutate(timing = row_number()) %>%
      select(-days, -months, -VISNAM1A) %>%
      pivot_wider(values_from = c(SBPAVG3N, DBPAVG3N, STNPLS1N),
                  names_from = timing) %>%
      mutate(landmark_time = landmark_time)
  )
}

full_data_new_endpoints_landmark_tbl = full_data_new_endpoints_landmark_tbl %>%
  pivot_longer(cols = any_of(c(
    "T2CVMB_DTH", "T2DTH", "C4CVMB_DTH", "C4DTH"
  ))) %>%
  mutate(
    variable_type = ifelse(stringr::str_starts(name, "C4"), "censored", "time"),
    endpoint = stringr::str_sub(name, start = 3)
  ) %>%
  select(-name) %>%
  pivot_wider(values_from = "value", names_from = "variable_type")


# For each level of the clustering variable, we fit a cox model excluding the
# data from that level. 
landmark_cox_models_tbl = expand_grid(
  landmark_time_model = landmark_times,
  endpoint_model = c("CVMB_DTH", "DTH"),
  clustering_variable_model = levels(full_data_new_endpoints$clustering_variable)
) %>%
  group_by(landmark_time_model,
           endpoint_model,
           clustering_variable_model) %>%
  summarize(
    cox_fit = list(
      landmark_cox_model(
        full_data_new_endpoints_landmark_tbl %>%
          filter(
            landmark_time == landmark_time_model[1],
            endpoint == endpoint_model[1],
            clustering_variable != clustering_variable_model[1]
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
            endpoint == endpoint_model[1],
            clustering_variable != clustering_variable_model[1]
          ) %>%
          # The landmark model only uses data on patients that are still under
          # observation at the landmark time.
          filter(time > landmark_time_model[1]),
        type = "base"
      )
    )
  ) %>%
  rename(
    landmark_time = landmark_time_model,
    endpoint = endpoint_model,
    clustering_variable = clustering_variable_model
  )


# Restructure to have one model per row.
landmark_cox_models_tbl = landmark_cox_models_tbl %>%
  pivot_longer(cols = c("cox_fit", "cox_fit_base"), 
               values_to = "cox_model", 
               names_to = "model_type") %>%
  mutate(model_type = ifelse(model_type == "cox_fit", "full", "base"))

# Save the data in landmark-prediction format and the fitted landmark cox
# models.
saveRDS(
  cox_models_new_endpoints_tbl,
  file = paste0(
    "R/data-fusion/intermediate-objects/",
    clustering_variable_chr,
    "/cox_models_new_endpoints_tbl.rds"
  )
)

saveRDS(
  full_data_new_endpoints_landmark_tbl,
  file = paste0(
    "R/data-fusion/intermediate-objects/",
    clustering_variable_chr,
    "/full_data_new_endpoints_landmark_tbl.rds"
  )
)
saveRDS(
  landmark_cox_models_tbl,
  file = paste0(
    "R/data-fusion/intermediate-objects/",
    clustering_variable_chr,
    "/landmark_cox_models_tbl.rds"
  )
)

