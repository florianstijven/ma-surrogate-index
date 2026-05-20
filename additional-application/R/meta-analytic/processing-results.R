# Setup ------------------------------------------------------------------
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
library(patchwork)
library(ggpp)

# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

ma_trt_effects_tbl_location = "additional-application/results/raw-results/ma_trt_effects_tbl.rds"
# Specify options for saving the plots to files
figures_dir = "additional-application/results/figures"
tables_dir = "additional-application/results/tables"


# Read in nonparametric results
surrogate_results_tbl = read.csv(file = paste0(tables_dir, "/", "surrogacy-inferences.csv")) %>%
  # Add Number of months as variable.
  mutate(
    landmark_time_months = 12 * landmark_time / 365.25,
    landmark_time_months_chr = paste0(landmark_time_months, " months"),
    landmark_time_months_fct = forcats::fct_reorder(landmark_time_months_chr, landmark_time_months)
  )

# Load data with trial-level treatment effects. 
ma_trt_effects_tbl = readRDS(ma_trt_effects_tbl_location)

# Add better trial names.
ma_trt_effects_tbl = ma_trt_effects_tbl %>%
  # Reorder trial factor.
  mutate(
    trial_fct = as.factor(COU1A)
  ) %>%
  # Add Number of months as variable.
  mutate(
    landmark_time_months = 12 * landmark_time / 365.25,
    landmark_time_months_chr = paste0(landmark_time_months, " months"),
    landmark_time_months_fct = forcats::fct_reorder(landmark_time_months_chr, landmark_time_months)
  )

# Duplicate original data where the vcov matrices are divided by 2, to get
# results that better illustrate the methods.
ma_trt_effects_tbl = bind_rows(
  ma_trt_effects_tbl %>%
    mutate(vcov_multiplier = 1),
  ma_trt_effects_tbl %>%
    mutate(
      trt_effect_surrogate_index_est_se = purrr::map_dbl(trt_effect_surrogate_index_est_se, ~ .x / sqrt(2)),
      trt_effect_clinical_est_se = purrr::map_dbl(trt_effect_clinical_est_se, ~ .x / sqrt(2)),
      vcov_multiplier = 0.5
    )
)

# Tables ----------------------------------------------------------------

# Save tables for the trial-level correlation estimates with corresponding
# confidence intervals. Tables are saves separately by type of CI and whether the
# original data or pseudo-real data were used.
format_CI = function(estimate, ci_lower, ci_upper, digits = 2) {
  paste0(round(estimate, digits), " (", round(ci_lower, digits), ", ", round(ci_upper, digits), ")")
}

# Plots -----------------------------------------------------------------

## MA plots -------------------------------------------------------------

# Function to create MA plots for given time period, SI estimation method.
plot_ma_f = function(time_start, time_end, method, vcov_multiplier) {
  ma_trt_effects_tbl %>%
    filter(
      landmark_time >= time_start,
      landmark_time < time_end,
      method == .env$method,
      vcov_multiplier == .env$vcov_multiplier
    ) %>%
    ggplot(
      aes(
        x = trt_effect_surrogate_index_est,
        y = trt_effect_clinical_est,
        color = trial_fct,
        shape = trial_fct
      )
    ) +
    geom_point() +
    geom_errorbar(
      aes(
        ymin = trt_effect_clinical_est - 1.96 * trt_effect_clinical_est_se,
        ymax = trt_effect_clinical_est + 1.96 * trt_effect_clinical_est_se
      ),
      width = 0.01
    ) +
    geom_errorbarh(
      aes(
        xmin = trt_effect_surrogate_index_est - 1.96 * trt_effect_surrogate_index_est_se,
        xmax = trt_effect_surrogate_index_est + 1.96 * trt_effect_surrogate_index_est_se
      ),
      height = 0.01
    ) +
    facet_grid(endpoint ~ landmark_time_months_fct) +
    xlab("Treatment effect on Est. SI") +
    ylab("Diff in Surv Prob") +
    theme(legend.position = "bottom", legend.title = element_blank())
}

# Create MA plots for the different time periods and SI estimation methods.
plotting_pms = expand_grid(
  time_start = 0,
  time_end = 365.25/2,
  method = c("predicted_prob_cox", "predicted_prob_sl"),
  vcov_multiplier = c(1, 0.5)
) %>%
  bind_rows(
    expand_grid(
      time_start = 365.25/2,
      time_end = 1e4,
      method = c("predicted_prob_cox", "predicted_prob_sl"),
      vcov_multiplier = c(1, 0.5)
    )
  ) 

# Save plots.
plotting_pms %>%
  rowwise() %>%
  summarise(
    ggsave(
      filename = paste0("ma-", method, "-", time_start, "-", time_end, "-", vcov_multiplier, ".pdf"),
      plot = plot_ma_f(time_start, time_end, method, vcov_multiplier),
      device = "pdf",
      width = double_width,
      height = double_height,
      units = "cm",
      dpi = res, 
      path = figures_dir
    )
  )

## Plot of estimated trial-level correlations + CIs -------------------------

# Function to create plot of trial-level correlation estimates with confidence
# intervals.
plot_rho_f = function() {
  surrogate_results_tbl %>%
    ggplot(aes(x = landmark_time_months, y = rho_trial, color = method, shape = method)) +
    geom_point() +
    geom_errorbar(aes(ymin = CI_lower_sandwich, ymax = CI_upper_sandwich), width = 0.01) +
    facet_grid(endpoint~.) +
    coord_cartesian(ylim = c(-1, 1)) +
    xlab("Time") +
    ylab("Trial-level correlation") +
    theme(legend.position = "bottom", legend.title = element_blank())
}
