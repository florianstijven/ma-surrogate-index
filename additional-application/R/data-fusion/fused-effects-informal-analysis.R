args = commandArgs(trailingOnly=TRUE)
# This R script takes only one argument: the variable used as clustering unit.
clustering_variable_chr = args[1]

# Load packages.
library(tidyverse)
library(survival)
library(GGally)
library(mgcv)
library(boot)

# Data Preparation --------------------------------------------------------
# Load intermediate results.
intermediate_object_location = paste0("R/data-fusion/intermediate-objects/",
                                      clustering_variable_chr,
                                      "/")
source("R/data-fusion/helper-functions.R")
cox_models_new_endpoints_tbl = readRDS(file = paste0(
  intermediate_object_location,
  "cox_models_new_endpoints_tbl.rds"
))
landmark_cox_models_tbl = readRDS(file = paste0(intermediate_object_location, "landmark_cox_models_tbl.rds"))
full_data_new_endpoints_landmark_tbl = readRDS(file = paste0(
  intermediate_object_location,
  "full_data_new_endpoints_landmark_tbl.rds"
))

landmark_times = 365.25 * c(1 / 12, 2 / 12, 3 / 12, 4 / 12, 5 / 12, 0.5, 1:2)

# Overall Treatment Effects -----------------------------------------------

## Kaplan-Meier Plots -----------------------------------------------------

cox_models_new_endpoints_tbl %>%
  group_by(endpoint, clustering_variable) %>%
  reframe(summary(cox_fit[[1]], times = seq(from = 0, to = 365.25*5, by = 10))[c("time", "surv", "lower", "upper", "strata")] %>%
            as_tibble()) %>%
  mutate(TRTREG1C = stringr::str_sub(strata, start = -1, end = -1)) %>%
  ggplot(aes(x = time / 365.25)) +
  geom_line(aes(y = surv, color = TRTREG1C)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = TRTREG1C), alpha = 0.2) +
  facet_grid(endpoint ~ clustering_variable, scales = "free_y") +
  scale_x_continuous(name = "Time (years)") +
  scale_y_continuous(name = "Survival Probability") +
  theme(legend.position = "bottom")
ggsave(
  filename = paste0("figures/data-fusion/", clustering_variable_chr, "/KM-composite-endpoints.pdf"),
  device = "pdf",
  width = double_width,
  height = double_height,
  units = "cm",
  dpi = res
)


## Effects on Survival Probabilities --------------------------------------

# The treatment effect is estimated on the time-to-event endpoint. We look at
# the 1, 2, 3, and 4-year survival probabilities. Other treatment effect
# measures (e.g., RMST) can be added later.
times = 365.25 * 1:4 
# Function to extract estimated survival probability and corresponding SE for a
# given treatment group.
survival_prob = function(surv_fit, times) {
  surv_fit_summary = summary(surv_fit, times = times)
  surv_fit_summary[c("time", "surv", "strata", "std.err")] %>%
    as_tibble() %>%
    mutate(TRTREG1C = stringr::str_sub(strata, start = -1, end = -1)) %>%
    select(-strata)
}
# Estimate the treatment effects on the X-year survival probabilities.
surv_probs_tbl = cox_models_new_endpoints_tbl %>%
  group_by(endpoint, clustering_variable) %>%
  reframe(survival_prob(cox_fit[[1]], times)) %>%
  pivot_wider(
    names_from = c("TRTREG1C"),
    values_from = c("surv", "std.err")) %>%
  mutate(surv_diff = surv_A - surv_B,
         surv_diff_se = sqrt(std.err_A ** 2 + std.err_B ** 2),
         surv_diff_z = surv_diff / surv_diff_se)

# Plot the estimated survival probability differences + CIs.
surv_probs_tbl %>%
  ggplot(aes(x = time/365.25, color = clustering_variable, group = clustering_variable)) +
  geom_errorbar(aes(ymin = surv_diff - 1.96 * surv_diff_se,
                    ymax = surv_diff + 1.96 * surv_diff_se),
                position = position_dodge(width = 0.4, preserve = "total"), 
                width = 0.4
  ) +
  geom_point(aes(y = surv_diff),
             position = position_dodge(width = 0.4)) +
  scale_x_continuous(name = "Time (years)") +
  scale_y_continuous(name = "Difference in Survival Probability") +
  facet_wrap(~endpoint) +
  theme(legend.position = "bottom")
ggsave(
  filename = paste0("figures/data-fusion/", clustering_variable_chr, "/trt-effects-survival-probability.pdf"),
  device = "pdf",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)


# Predicted Treatment Effects ---------------------------------------------

## Computation ------------------------------------------------------------

# For patients that are still under observations at the landmark time, predict
# the probability of X-year survival. 

# For each landmark_time-region combination, we let the data be a single element
# in the tibble. This condensed data set is then joined with the tibble of
# models.
landmark_predictions_tbl = full_data_new_endpoints_landmark_tbl %>%
  group_by(landmark_time, clustering_variable, endpoint) %>%
  summarize(data_set = list(pick(everything()))) %>%
  left_join(landmark_cox_models_tbl,
            by = c("landmark_time", "endpoint", "clustering_variable")) %>%
  # We expand the tibble to include all X-year survival probabilities.
  cross_join(tibble(times = times)) %>%
  # Only keep entries where the landmark time is smaller than the X-year time
  # point.
  filter(landmark_time < times) %>%
  # filter(landmark_time == 1095.75, clustering_variable == "FIN", endpoint == "DTH", times ==
  #          1461, model_type == "full") %>%
  # For each row, which contains a full data set and the fitted gam, the
  # predictions are computed and we return to a normal tibble with one row per
  # data entry. To reduce memory usage, we only keep the variables that are
  # needed further on.
  group_by(landmark_time,
           clustering_variable,
           endpoint,
           times,
           model_type) %>%
  reframe(bind_cols(
    tibble(predicted_prob = landmark_prediction(data_set[[1]], cox_model[[1]], landmark_time[1], times[1])),
    data_set[[1]] %>%
      select(SID1A, time, censored, TRTREG1C)
  )) %>%
  mutate(censored_before_landmark = (time <= landmark_time) & censored == 1)

# We don't need the data set in landmark format anymore.
rm(full_data_new_endpoints_landmark_tbl)

# Proportion of patients censored before each landmark time. This proportion is
# very small (except at 3 years), so we further ignored these patients. The
# amount of NAs is a bit larger, this is due to missing intermediate endpoints.
landmark_predictions_tbl %>%
  group_by(landmark_time, clustering_variable, endpoint, model_type) %>%
  summarize(
    prop_censored_before_landmark = mean(censored_before_landmark, na.rm = TRUE),
    prop_NA = mean(is.na(predicted_prob))
  )

# Estimate the mean of the predicted survival probabilities for each landmark time,
# X-year survival, treatment group, and region.
predicted_surv_probs_tbl = landmark_predictions_tbl %>%
  group_by(landmark_time, times, TRTREG1C, clustering_variable, endpoint, model_type) %>%
  summarize(
    predicted_surv_prob = mean(predicted_prob, na.rm = TRUE)
  ) %>%
  pivot_wider(
    values_from = c("predicted_surv_prob"),
    names_from = c("TRTREG1C"),
    names_prefix = "predicted_surv_prob_"
  ) %>%
  mutate(
    predicted_surv_diff = predicted_surv_prob_A - predicted_surv_prob_B
  ) 

# Combine direct treatment effect estimates tibble with the predicted treatment
# effect estimates tibble.
fused_effects_tbl = surv_probs_tbl %>%
  rename(times = time) %>%
  right_join(predicted_surv_probs_tbl, by = c("endpoint", "clustering_variable", "times")) %>%
  mutate(
    landmark_month = 12 * landmark_time / 365.25,
    landmark_month_chr = factor(
      x = landmark_month,
      levels = 12 * landmark_times / 365.25,
      labels = paste0("Month ", 12 * landmark_times / 365.25)
    )
  )

## Plots: Predictions Versus Estimates -------------------------------------

# Plot direct estimates of treatment effects against those based on landmark
# predictions. We use a function to produce a set of similar plots.
fused_predictions_plot = function(endpoint_plot, model_type_plot, include_year_four) {
  if (model_type_plot == "both") {
    model_types_plot = c("base", "full")
  }
  else {
    model_types_plot = model_type_plot
  }
  if (include_year_four) {
    landmark_upper = 37
    times_upper = 365.25 * 4 + 1
  }
  else {
    landmark_upper = 36
    times_upper = 365.25 * 4
  }
  p = fused_effects_tbl %>%
    filter(endpoint == endpoint_plot, model_type %in% model_types_plot, landmark_month < landmark_upper, times < times_upper)  %>%
    ggplot(aes(
      x = times / 365.25,
      color = clustering_variable,
      group = clustering_variable
    )) +
    geom_errorbar(
      aes(
        ymin = surv_diff - 1.96 * surv_diff_se,
        ymax = surv_diff + 1.96 * surv_diff_se
      ),
      position = position_dodge(width = 0.80),
      width = 0.4,
      alpha = 0.35,
      linewidth = 0.35
    ) +
    geom_point(aes(y = surv_diff),
               position = position_dodge(width = 0.80),
               alpha = 0.35) +
    geom_point(
      aes(y = predicted_surv_diff, shape = model_type),
      position = position_dodge(width = 0.80),
      data = fused_effects_tbl %>%
        filter(endpoint == endpoint_plot, model_type %in% model_types_plot, landmark_month < landmark_upper, times < times_upper)
    ) +
    geom_hline(aes(yintercept = 0), alpha = 0.5) +
    geom_vline(aes(xintercept = landmark_month / 12), linetype = "dashed") +
    scale_x_continuous(name = "Time (years)") +
    scale_y_continuous(name = "Difference in Survival Probability") +
    facet_wrap( ~ landmark_month_chr, ) +
    coord_cartesian(ylim = c(-0.03, 0.075)) +
    scale_color_discrete(name = "Clustering Unit") +
    scale_shape_manual(values = c(0, 2), name = "Prediction Model") +
    theme(legend.position = "bottom", legend.box = "vertical")
  if (model_type_plot != "both") {
    p = p +
      ggtitle(paste0(endpoint_plot, " - ", model_type_plot, " prediction model")) +
      guides(shape = FALSE)
  }
  else {
    p = p +
      ggtitle("CVMB_DTH")
  }
  
  include_year_four_name = ifelse(include_year_four, "-4-years", "")
  filename = paste0(
    "figures/data-fusion/",
    clustering_variable_chr,
    "/trt-effects-predicted-estimated-",
    endpoint_plot,
    include_year_four_name,
    "-",
    model_type_plot,
    ".pdf"
  )
  ggsave(
    filename = filename,
    device = "pdf",
    width = double_width,
    height = double_height,
    units = "cm",
    dpi = res
  )
  return(1)
}

# Enumerate all plotting options in a tibble. Call plotting function for each
# combinations of plotting options.
expand_grid(
  endpoint_plot = c("CVMB_DTH", "DTH"),
  model_type = c("full", "both"),
  include_year_four = c(TRUE, FALSE)
) %>%
  group_by(endpoint_plot, model_type, include_year_four) %>%
  summarise(fused_predictions_plot(endpoint_plot, model_type, include_year_four))
