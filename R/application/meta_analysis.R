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

# # Set up parallel computing
# if (parallelly::supportsMulticore()) {
#   plan("multicore")
# } else {
#   plan(multisession)
# }

## Analysis Parameters -------------------------------------------------- 

# Number of bootstrap replications for the multiplier bootstrap for the
# meta-analytic parameters.
B_multiplier = 1e3

time_cumulative_incidence = 120

# Set of trials corresponding to each analysis set.

trials_first_four = c("AstraZeneca", "Moderna", "Janssen", "Novavax")
trials_naive_only = c("AstraZeneca",
                      "Moderna",
                      "Janssen",
                      "Novavax",
                      "Sanofi-1 naive",
                      "Sanofi-2 naive")
trials_mixed = c(
  "AstraZeneca",
  "Moderna",
  "Janssen",
  "Novavax",
  "Sanofi-1 naive",
  "Sanofi-2 naive",
  "Sanofi-1 non-naive",
  "Sanofi-2 non-naive"
)


## Intermediate Results  --------------------------------------------------

# Load data with trial-level treatment effects. 
ma_trt_effects_tbl = readRDS("R/application/ma_trt_effects_tbl.rds")

# Add poper name of the surrogates.
ma_trt_effects_tbl = ma_trt_effects_tbl %>%
  mutate(surrogate_name = case_when(
    surrogate == "bindSpike" ~ "IgG Spike",
    surrogate == "pseudoneutid50" ~ "nAb ID50",
    surrogate == "pseudoneutid50_adjusted" ~ "adjusted nAb ID50"
  ))

# Add indicator for whether only naive individuals are in a given trial. 
ma_trt_effects_tbl = ma_trt_effects_tbl %>%
  mutate(
    only_naive = !(trial %in% c("Sanofi-1 non-naive", "Sanofi-2 non-naive"))
  )

# Vector of trial names. We exclude 
trial_names = levels(ma_trt_effects_tbl %>% pull(trial))



# Meta-Analytic Plots -----------------------------------------------------

# Define transformation for VE scale.
transform_VE = new_transform(
  name = "VE",
  transform = function(x) -1 * log(1 - x),
  inverse = function(x) 1 - exp(-x)
)

# Summarize the estimated trial-level treatment effects in meta-analytic plots.
ma_trt_effects_tbl %>% 
  filter(method == "untransformed surrogate") %>%
  ggplot(aes(x = trt_effect_surrogate_index_est, y = 1 - exp(log_RR_est), color = trial, shape = trial)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = 1 - exp(log_RR_est - 1.96 * log_RR_est_se),
    ymax = 1 - exp(log_RR_est + 1.96 * log_RR_est_se)
  ),
  width = 0.01) +
  geom_errorbarh(
    aes(
      xmin = trt_effect_surrogate_index_est - 1.96 * trt_effect_surrogate_index_est_se,
      xmax = trt_effect_surrogate_index_est + 1.96 * trt_effect_surrogate_index_est_se
    ),
    height = 0.01
  ) +
  scale_y_continuous(transform = transform_VE, breaks = c(0, 0.5, 0.75, 0.90, 0.95)) +
  scale_shape_manual(labels = trial_names, values = c(19, 19, 19, 19, 19, 17, 19, 17), drop = FALSE) +
  coord_cartesian(ylim = c(-0.2, 0.95)) +
  facet_wrap(~surrogate_name, scales = "free_x", nrow = 2) +
  xlab("Estimated Treatment Effect on Titer") +
  ylab("Estimated VE") +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(
  filename = "results/figures/application/ma-standard.pdf",
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

ma_trt_effects_tbl %>% 
  filter(method == "untransformed surrogate") %>%
  ggplot(aes(x = trt_effect_surrogate_index_est, y = 1 - exp(log_RR_est), color = trial)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = 1 - exp(log_RR_est - 1.96 * log_RR_est_se),
    ymax = 1 - exp(log_RR_est + 1.96 * log_RR_est_se)
  ),
  width = 0.01) +
  geom_errorbarh(
    aes(
      xmin = trt_effect_surrogate_index_est - 1.96 * trt_effect_surrogate_index_est_se,
      xmax = trt_effect_surrogate_index_est + 1.96 * trt_effect_surrogate_index_est_se
    ),
    height = 0.01
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 0.75, 0.90, 0.95)) +
  scale_shape_manual(labels = trial_names, values = c(19, 19, 19, 19, 19, 17, 19, 17), drop = FALSE) +
  coord_cartesian(ylim = c(-0.2, 0.95)) +
  facet_wrap(~surrogate_name, scales = "free_x", nrow = 2) +
  xlab("Estimated Treatment Effect on Titer") +
  ylab("Estimated VE") +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(
  filename = "results/figures/application/ma-standard-ve-scale.pdf",
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

# Helper function for MA plots for the surrogate indices by weighting and
# include_Bserostatus.
ma_plot_helper = function(weighting, analysis_set, VE_scale) {
  # Set of included trials according to analysis_set
  if (analysis_set == "first_four") {
    trials_included = trials_first_four
  } else if (analysis_set == "naive_only") {
    trials_included = trials_naive_only
  } else if (analysis_set == "mixed") {
    trials_included = trials_mixed
  }
  ggplot_object = ma_trt_effects_tbl %>% filter(method != "untransformed surrogate", 
                                analysis_set == .env$analysis_set, 
                                trial %in% trials_included) %>%
    ggplot(aes(
      x = 1 - exp(trt_effect_surrogate_index_est),
      y = 1 - exp(log_RR_est),
      color = trial
    )) +
    geom_point() +
    geom_errorbar(aes(
      ymin = 1 - exp(log_RR_est - 1.96 * log_RR_est_se),
      ymax = 1 - exp(log_RR_est + 1.96 * log_RR_est_se)
    ), width = 0.01) +
    geom_errorbarh(aes(
      xmin = 1 - exp(
        trt_effect_surrogate_index_est - 1.96 * trt_effect_surrogate_index_est_se
      ),
      xmax = 1 - exp(
        trt_effect_surrogate_index_est + 1.96 * trt_effect_surrogate_index_est_se
      )
    ), height = 0.01) +
    geom_abline(intercept = 0, slope = 1) +
    scale_shape_manual(labels = trial_names, values = c(19, 19, 19, 19, 19, 17, 19, 17), drop = FALSE) +
    coord_cartesian(ylim = c(-0.2, 0.95)) +
    facet_grid(method ~ surrogate_name) +
    xlab("Estimated VE on Surr. Index") +
    ylab("Estimated VE") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  if (VE_scale) {
    ggplot_object = ggplot_object +
      scale_y_continuous(breaks = c(0, 0.5, 0.75, 0.90, 0.95)) +
      scale_x_continuous(breaks = c(0, 0.5, 0.75, 0.90, 0.95))
  } else{
    ggplot_object = ggplot_object +
      scale_y_continuous(transform = transform_VE, breaks = c(0, 0.5, 0.75, 0.90, 0.95)) +
      scale_x_continuous(transform = transform_VE, breaks = c(0, 0.5, 0.75, 0.90, 0.95))
  }
  
  weighting_chr = ifelse(weighting == "Normalized By Trial", "normalized", "unnormalized")
  VE_scale_chr = ifelse(VE_scale, "ve-scale", "transformed-scale")
  
  outfile = paste0(
    "results/figures/application/ma-standard",
    "-",
    weighting_chr,
    "-",
    analysis_set,
    "-",
    VE_scale_chr,
    ".pdf"
    
  )
  
  ggsave(
    plot = ggplot_object,
    filename = outfile,
    height = double_height,
    width = double_width,
    dpi = res,
    device = "pdf",
    units = "cm"
  )
}

expand_grid(
  weighting = c("Unnormalized"),
  analysis_set = c("first_four", "naive_only", "mixed"),
  VE_scale = c(FALSE, TRUE)
) %>%
  rowwise(everything()) %>%
  summarise(ma_plot_helper(weighting, analysis_set, VE_scale))

# Create extra plots for the four-trial analysis set where the treatment effects
# for the excluded trials are shown.
ma_trt_effects_tbl %>% filter(
  method != "untransformed surrogate",
  analysis_set == "first_four",!(trial %in% trials_first_four)
) %>%
  ggplot(aes(
    x = 1 - exp(trt_effect_surrogate_index_est),
    y = 1 - exp(log_RR_est),
    color = trial
  )) +
  geom_point() +
  geom_errorbar(aes(
    ymin = 1 - exp(log_RR_est - 1.96 * log_RR_est_se),
    ymax = 1 - exp(log_RR_est + 1.96 * log_RR_est_se)
  ), width = 0.01) +
  geom_errorbarh(aes(
    xmin = 1 - exp(
      trt_effect_surrogate_index_est - 1.96 * trt_effect_surrogate_index_est_se
    ),
    xmax = 1 - exp(
      trt_effect_surrogate_index_est + 1.96 * trt_effect_surrogate_index_est_se
    )
  ), height = 0.01) +
  geom_abline(intercept = 0, slope = 1) +
  scale_shape_manual(
    labels = trial_names,
    values = c(19, 19, 19, 19, 19, 17, 19, 17),
    drop = FALSE
  ) +
  facet_grid(method ~ surrogate_name) +
  xlab("Estimated VE on Surr. Index") +
  ylab("Estimated VE") +
  theme(legend.position = "bottom", legend.title = element_blank()) 


ggsave(
  filename = "results/figures/application/ma-predictions.pdf",
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)




# Formal Meta-Analysis -----------------------------------------------------

source("R/helper-functions/moment-based-estimator.R")
source("R/helper-functions/multiplier-bootstrap.R")
source("R/helper-functions/delta-method-rho-trial.R")

# Helper function that returns the trial-level correlation estimate given a
# weighted data set.
statistic_f = function(data, weights) {
  moment_estimate = moment_estimator(
    alpha_hat = data$treatment_effect_surr,
    beta_hat = data$treatment_effect_clin,
    vcov_list = data$covariance_matrix,
    estimator_adjustment = "N - 1",
    weights = weights,
    nearest_PD = TRUE
  )
  rho = rho_delta_method(
    coefs = moment_estimate$coefs,
    vcov = moment_estimate$vcov,
    method = "t-adjustment",
    # N is only used for the t-adjustment, it doesn't matter for the estimate
    # or SE.
    N = 5
  )
  
  return(list(estimate = rho$rho, se = rho$se))
}

# Construct data set for the treatment effect estimates on the untransformed
# surrogates. We duplicate rows such that we have a set of rows corresponding to
# each analysis set.
ma_trt_effects_tbl_untransformed = bind_rows(
  ma_trt_effects_tbl %>%
    filter(method == "untransformed surrogate") %>%
    filter(trial %in% trials_first_four) %>%
    mutate(analysis_set = "first_four"),
  ma_trt_effects_tbl %>%
    filter(method == "untransformed surrogate") %>%
    filter(trial %in% trials_naive_only) %>%
    mutate(analysis_set = "naive_only"),
  ma_trt_effects_tbl %>%
    filter(method == "untransformed surrogate") %>%
    filter(trial %in% trials_mixed) %>%
    mutate(analysis_set = "mixed")
)


# Estimate the surrogacy parameters on each data set of trial-level treatment
# effect estimates.
surrogate_results_tbl = ma_trt_effects_tbl %>%
  filter(method != "untransformed surrogate") %>%
  bind_rows(ma_trt_effects_tbl_untransformed) %>%
  # Make sure that only the correct trials are included for each analysis set.
  filter(ifelse(
    analysis_set == "first_four",
    trial %in% trials_first_four,
    ifelse(analysis_set == "naive_only", trial %in% trials_naive_only, TRUE)
  )) %>%
  group_by(surrogate, method, weighting, analysis_set) %>%
  summarise(
    moment_estimate = list(
      moment_estimator(
        alpha_hat = trt_effect_surrogate_index_est,
        beta_hat = log_RR_est,
        vcov_list = vcov,
        estimator_adjustment = "N - 1",
        sandwich_adjustment = "N - 1",
        nearest_PD = TRUE
      )
    ),
    bootstrap_ci = list(
      multiplier_bootstrap_ci(
        data = pick(everything()) %>%
          rename(
            treatment_effect_surr = "trt_effect_surrogate_index_est",
            treatment_effect_clin = 'log_RR_est',
            covariance_matrix = "vcov"
          ),
        statistic = statistic_f,
        B = B_multiplier,
        type = "percentile"
      )
    )
  ) %>%
  ungroup() 

surrogate_results_tbl = surrogate_results_tbl %>%
  mutate(
    d_alpha = map_dbl(moment_estimate, function(x) x$coefs[3]),
    d_beta = map_dbl(moment_estimate, function(x) x$coefs[4]),
    d_alphabeta = map_dbl(moment_estimate, function(x) x$coefs[5]),
    rho_trial = d_alphabeta / sqrt(d_alpha * d_beta),
    residual_var = map_dbl(moment_estimate, function(x) x$residual_var)
  )

# Summarize inferences in a table.
surrogate_results_tbl = surrogate_results_tbl %>%
  mutate(
    CI_lower_bs = purrr::map_dbl(bootstrap_ci, "ci_lower"),
    CI_upper_bs = purrr::map_dbl(bootstrap_ci, "ci_upper")
  ) %>%
  select(-moment_estimate, -bootstrap_ci) 

surrogate_results_tbl %>%
  write.csv(file = "results/tables/application/surrogacy-inferences.csv")

# Summarize the results in a plot. 
surrogate_results_tbl %>%
  # Multiply the correlations for the untransformed surrogate with -1 to get
  # positive correlations which are comparable with those for the surrogate
  # indices.
  mutate(rho_trial = ifelse(method == "untransformed surrogate", -1 * rho_trial, rho_trial),
         CI_lower_bs = ifelse(method == "untransformed surrogate", -1 * CI_lower_bs, CI_lower_bs),
         CI_upper_bs = ifelse(method == "untransformed surrogate", -1 * CI_upper_bs, CI_upper_bs)) %>%
  ggplot(aes(x = analysis_set, color = method)) +
  geom_point(aes(y = rho_trial), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_lower_bs, ymax = CI_upper_bs), position = position_dodge(width = 0.5), width = 0.2) +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_grid(surrogate~.) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin())
ggsave(
  "results/figures/application/surrogacy-measures-summary.pdf",
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

