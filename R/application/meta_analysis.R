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
library(RColorBrewer)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}

# Specify options for saving the plots to files
figures_dir = "results/figures/application/meta-analysis"
tables_dir = "results/tables/application/meta-analysis"

## Analysis Parameters -------------------------------------------------- 

# Number of bootstrap replications for the multiplier bootstrap for the
# meta-analytic parameters.
B_multiplier = 1e5

time_cumulative_incidence = 80

# Set of trials corresponding to each analysis set.

trials_first_four = c("Moderna", "AstraZeneca", "Janssen", "Novavax")
trials_first_four_fct = trials_first_four
trials_naive_only = c("Moderna",
                      "AstraZeneca",
                      "Janssen",
                      "Novavax",
                      "Sanofi-1 naive",
                      "Sanofi-2 naive")
trials_naive_only_fct = c("Moderna",
                          "AstraZeneca",
                          "Janssen",
                          "Novavax",
                          "Sanofi 1 (naive)",
                          "Sanofi 2 (naive)")
trials_mixed = c(
  "Moderna",
  "AstraZeneca",
  "Janssen",
  "Novavax",
  "Sanofi-1 naive",
  "Sanofi-1 non-naive",
  "Sanofi-2 naive",
  "Sanofi-2 non-naive"
)
trials_mixed_fct = c(
  "Moderna",
  "AstraZeneca",
  "Janssen",
  "Novavax",
  "Sanofi 1 (naive)",
  "Sanofi 1 (non-naive)",
  "Sanofi 2 (naive)",
  "Sanofi 2 (non-naive)"
)


## Intermediate Results  --------------------------------------------------

# Load data with trial-level treatment effects. 
ma_trt_effects_tbl = readRDS("R/application/ma_trt_effects_tbl.rds")

# Add proper name of the surrogates.
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

# Add better trial names.
ma_trt_effects_tbl = ma_trt_effects_tbl %>%
  # Reorder trial factor.
  mutate(
    trial_fct = forcats::fct_recode(
      trial,
      "Moderna" = "Moderna",
      "AstraZeneca" = "AstraZeneca",
      "Janssen" = "Janssen",
      "Novavax" = "Novavax",
      "Sanofi 1 (naive)" = "Sanofi-1 naive",
      "Sanofi 1 (non-naive)" = "Sanofi-1 non-naive",
      "Sanofi 2 (naive)" = "Sanofi-2 naive",
      "Sanofi 2 (non-naive)" = "Sanofi-2 non-naive"
    ),
    trial_fct = fct_relevel(
      trial_fct,
      "Moderna",
      "AstraZeneca",
      "Janssen",
      "Novavax",
      "Sanofi 1 (naive)",
      "Sanofi 1 (non-naive)",
      "Sanofi 2 (naive)",
      "Sanofi 2 (non-naive)"
    )
  )

# Create a named vector as the palette
class <- levels(ma_trt_effects_tbl$trial_fct)  # All factor levels across figures
colors <- brewer.pal(length(class), "Dark2")  # Colors
my_palette <- set_names(colors, class)  # Named vector
# Create corresponding set of shapes that distinguish beween naive and non-naive
# trials.
my_shapes = set_names(c(19, 19, 19, 19, 19, 17, 19, 17), class)

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
# ma_trt_effects_tbl %>% 
#   filter(method == "untransformed surrogate") %>%
#   ggplot(aes(x = trt_effect_surrogate_index_est, y = 1 - exp(log_RR_est), color = trial_fct, shape = trial_fct)) +
#   geom_point() +
#   geom_errorbar(aes(
#     ymin = 1 - exp(log_RR_est - 1.96 * log_RR_est_se),
#     ymax = 1 - exp(log_RR_est + 1.96 * log_RR_est_se)
#   ),
#   width = 0.01) +
#   geom_errorbarh(
#     aes(
#       xmin = trt_effect_surrogate_index_est - 1.96 * trt_effect_surrogate_index_est_se,
#       xmax = trt_effect_surrogate_index_est + 1.96 * trt_effect_surrogate_index_est_se
#     ),
#     height = 0.01
#   ) +
#   scale_y_continuous(transform = transform_VE, breaks = c(0, 0.5, 0.75, 0.90, 0.95)) +
#   scale_color_manual(values = my_palette, limits = trials_mixed_fct) +
#   scale_shape_manual(values = my_shapes, limits = trials_mixed_fct) +
#   coord_cartesian(ylim = c(-0.2, 0.95)) +
#   facet_wrap(~surrogate_name, scales = "free_x", nrow = 2) +
#   xlab("Estimated Treatment Effect on Antibody Marker") +
#   ylab("Estimated VE") +
#   theme(legend.position = "bottom", legend.title = element_blank())
# 
# ggsave(
#   filename = "ma-standard-all.pdf",
#   path = figures_dir,
#   height = double_height,
#   width = double_width,
#   dpi = res,
#   device = "pdf",
#   units = "cm"
# )

# Summarize the estimated trial-level treatment effects in meta-analytic plots.
ma_trt_effects_tbl %>% 
  filter(method == "untransformed surrogate", trial %in% trials_naive_only) %>%
  ggplot(aes(x = trt_effect_surrogate_index_est, y = 1 - exp(log_RR_est), color = trial_fct, shape = trial_fct)) +
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
  scale_color_manual(values = my_palette, limits = trials_naive_only_fct) +
  scale_shape_manual(values = my_shapes, limits = trials_naive_only_fct) +
  coord_cartesian(ylim = c(-0.2, 0.95)) +
  facet_wrap(~surrogate_name, scales = "free_x", nrow = 2) +
  xlab("Estimated Treatment Effect on Antibody Marker") +
  ylab("Estimated VE") +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(
  filename = "ma-standard-naive-only.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)




# Helper function for MA plots for the surrogate indices by weighting and
# analysis_set.
ma_plot_helper = function(weighting, analysis_set, VE_scale, include_risk_score) {
  # Set of included trials according to analysis_set
  if (analysis_set == "first_four") {
    trials_included = trials_first_four_fct
  } else if (analysis_set == "naive_only") {
    trials_included = trials_naive_only_fct
  } else if (analysis_set == "mixed") {
    trials_included = trials_mixed_fct
  }
  ggplot_object = ma_trt_effects_tbl %>% filter(method != "untransformed surrogate", 
                                analysis_set == .env$analysis_set, 
                                trial_fct %in% trials_included,
                                weighting == .env$weighting, 
                                include_risk_score == .env$include_risk_score) %>%
    ggplot(aes(
      x = 1 - exp(trt_effect_surrogate_index_est),
      y = 1 - exp(log_RR_est),
      color = trial_fct,
      shape = trial_fct
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
    scale_color_manual(values = my_palette, limits = trials_included) +
    scale_shape_manual(values = my_shapes, limits = trials_included) +
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
  
  weighting_chr = ifelse(weighting == "normalized", "normalized", "unnormalized")
  VE_scale_chr = ifelse(VE_scale, "ve-scale", "transformed-scale")
  include_risk_score_chr = ifelse(include_risk_score, "w-riskscore", "wo-riskscore")
  
  outfile = paste0(
    "ma",
    "-",
    weighting_chr,
    "-",
    analysis_set,
    "-",
    include_risk_score_chr, 
    "-",
    VE_scale_chr,
    ".pdf"
    
  )
  
  ggsave(
    plot = ggplot_object,
    filename = outfile,
    path = figures_dir,
    height = double_height,
    width = double_width,
    dpi = res,
    device = "pdf",
    units = "cm"
  )
}

expand_grid(
  weighting = c("unnormalized"),
  analysis_set = c("naive_only"),
  VE_scale = c(FALSE),
  include_risk_score = c(FALSE)
) %>%
  rowwise(everything()) %>%
  summarise(ma_plot_helper(weighting, analysis_set, VE_scale, include_risk_score))

# Create extra plots for the four-trial analysis set where the treatment effects
# for the excluded trials are shown.
ma_trt_effects_tbl %>% filter(
  method != "untransformed surrogate",
  analysis_set == "naive_only",
  !(trial %in% trials_naive_only),
  weighting == "unnormalized",
  include_risk_score == FALSE
) %>%
  ggplot(aes(
    x = 1 - exp(trt_effect_surrogate_index_est),
    y = 1 - exp(log_RR_est),
    color = trial_fct,
    shape = trial_fct
  )) +
  geom_point(show.legend = TRUE) +
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
  scale_color_manual(values = my_palette, limits = trials_mixed_fct[!(trials_mixed_fct %in% trials_first_four_fct)]) +
  scale_shape_manual(values = my_shapes, limits = trials_mixed_fct[!(trials_mixed_fct %in% trials_first_four_fct)]) +
  scale_y_continuous(transform = transform_VE, breaks = c(0, 0.5, 0.75, 0.90, 0.95)) +
  scale_x_continuous(transform = transform_VE, breaks = c(0, 0.5, 0.75, 0.90, 0.95)) +
  coord_cartesian(ylim = c(-0.2, 0.95), xlim = c(-0.2, 0.95)) +
  facet_grid(method ~ surrogate_name) +
  xlab("Predicted VE/Estimated VE on Surr. Index") +
  ylab("Estimated VE") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2))

ggsave(
  filename = "ma-predictions-unnormalized-wo-riskscore.pdf",
  path = figures_dir,
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
statistic_f_rho = function(data, weights) {
  moment_estimate = moment_estimator(
    alpha_hat = data$treatment_effect_surr,
    beta_hat = data$treatment_effect_clin,
    vcov_list = data$covariance_matrix,
    estimator_adjustment = "N - 1",
    weights = weights,
    nearest_PD = FALSE
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

statistic_f_residual_var = function(data, weights) {
  moment_estimate = moment_estimator(
    alpha_hat = data$treatment_effect_surr,
    beta_hat = data$treatment_effect_clin,
    vcov_list = data$covariance_matrix,
    estimator_adjustment = "N - 1",
    weights = weights,
    nearest_PD = FALSE,
    SE = FALSE
  )
  # Residual variance
  residual_var = moment_estimate$residual_var
  
  
  return(list(estimate = residual_var, se = NA))
}

statistic_f_residual_var_prop = function(data, weights) {
  moment_estimate = moment_estimator(
    alpha_hat = data$treatment_effect_surr,
    beta_hat = data$treatment_effect_clin,
    vcov_list = data$covariance_matrix,
    estimator_adjustment = "N - 1",
    weights = weights,
    nearest_PD = TRUE,
    SE = FALSE
  )
  # Residual variance
  residual_var = max(moment_estimate$residual_var, 1e-5)
  var_beta = max(moment_estimate$coefs[4], 1e-5)
  
  # Proportion of variance in beta explained by the identity line
  prop_explained = 1 - (residual_var / var_beta)
  
  return(list(estimate = prop_explained, se = NA))
}

# Construct data set for the treatment effect estimates on the untransformed
# surrogates. We duplicate rows such that we have a set of rows corresponding to
# each analysis set.
ma_trt_effects_tbl_untransformed = bind_rows(
  ma_trt_effects_tbl %>%
    filter(method == "untransformed surrogate") %>%
    filter(trial %in% trials_naive_only) %>%
    mutate(analysis_set = "naive_only")
  # ma_trt_effects_tbl %>%
  #   filter(method == "untransformed surrogate") %>%
  #   filter(trial %in% trials_mixed) %>%
  #   mutate(analysis_set = "mixed")
)

# Construct data set with a set of rows corresponding to each analysis set.
ma_trt_effects_tbl_modified = ma_trt_effects_tbl %>%
  filter(method != "untransformed surrogate") %>%
  bind_rows(ma_trt_effects_tbl_untransformed) %>%
  # Make sure that only the correct trials are included for each analysis set.
  filter(ifelse(
    analysis_set == "first_four",
    trial %in% trials_first_four,
    ifelse(analysis_set == "naive_only", trial %in% trials_naive_only, TRUE)
  ))

# Add pseudo-real data by (i) cloning each trial four times and (ii) dividing
# the within-trial variance matrices by four.
ma_trt_effects_tbl_modified = bind_rows(
  ma_trt_effects_tbl_modified %>%
    mutate(scenario = "real data"),
  bind_rows(
    ma_trt_effects_tbl_modified,
    ma_trt_effects_tbl_modified,
    ma_trt_effects_tbl_modified,
    ma_trt_effects_tbl_modified
  ) %>%
    mutate(scenario = "four clones", analysis_set == "naive_only"),
  ma_trt_effects_tbl_modified %>%
    filter(analysis_set == "naive_only") %>%
    mutate(
      vcov = purrr::map(
        .x = vcov,
        .f = function(x) {
          x / 4
        }
      ),
      scenario = "precise trials"
    )
)


# Estimate the surrogacy parameters on each data set of trial-level treatment
# effect estimates.
surrogate_results_tbl = ma_trt_effects_tbl_modified %>%
  group_by(surrogate, method, weighting, analysis_set, include_risk_score, scenario) %>%
  summarise(data_tbl = list(pick(everything())), N = nrow(data_tbl[[1]])) %>%
  ungroup() %>%
  mutate(
    moment_estimate = purrr::map(
      .x = data_tbl,
      .f = function(data_tbl) {
        moment_estimator(
          alpha_hat = data_tbl$trt_effect_surrogate_index_est,
          beta_hat = data_tbl$log_RR_est,
          vcov_list = data_tbl$vcov,
          estimator_adjustment = "N - 1",
          sandwich_adjustment = "N - 1",
          nearest_PD = TRUE
        )
      }
    ),
    bootstrap_ci = future_map(
      .x = data_tbl,
      .f = function(data_tbl) {
        multiplier_bootstrap_ci(
          data = data_tbl %>%
            rename(
              treatment_effect_surr = "trt_effect_surrogate_index_est",
              treatment_effect_clin = 'log_RR_est',
              covariance_matrix = "vcov"
            ),
          statistic = statistic_f_rho,
          B = B_multiplier,
          type = "BCa",
          alpha = 0.05
        )
      }
    ),
    bootstrap_ci_residual_var = future_map(
      .x = data_tbl,
      .f = function(data_tbl) {
        multiplier_bootstrap_ci(
          data = data_tbl %>%
            rename(
              treatment_effect_surr = "trt_effect_surrogate_index_est",
              treatment_effect_clin = 'log_RR_est',
              covariance_matrix = "vcov"
            ),
          statistic = statistic_f_residual_var,
          B = B_multiplier,
          type = "BCa",
          alpha = 0.05
        )
      }
    ),
    # bootstrap_ci_residual_var_prop = future_map(
    #   .x = data_tbl,
    #   .f = function(data_tbl) {
    #     multiplier_bootstrap_ci(
    #       data = data_tbl %>%
    #         rename(
    #           treatment_effect_surr = "trt_effect_surrogate_index_est",
    #           treatment_effect_clin = 'log_RR_est',
    #           covariance_matrix = "vcov"
    #         ),
    #       statistic = statistic_f_residual_var_prop,
    #       B = B_multiplier,
    #       type = "BCa",
    #       alpha = 0.05
    #     )
    #   }
    # )
  ) %>%
  mutate(rho_sandwich_inference = purrr::map2(
    .x = moment_estimate,
    .y = N,
    .f = function(moment_estimate, N) {
      rho_delta_method(
        coefs = moment_estimate$coefs,
        vcov = moment_estimate$vcov,
        method = "t-adjustment",
        N = N
      )
    }
  ))


surrogate_results_tbl = surrogate_results_tbl %>%
  mutate(
    d_alpha = map_dbl(moment_estimate, function(x)
      x$coefs[3]),
    d_beta = map_dbl(moment_estimate, function(x)
      x$coefs[4]),
    d_alphabeta = map_dbl(moment_estimate, function(x)
      x$coefs[5]),
    rho_trial = d_alphabeta / sqrt(d_alpha * d_beta),
    residual_var = map_dbl(moment_estimate, function(x)
      x$residual_var),
    CI_lower_bs = purrr::map_dbl(bootstrap_ci, "ci_lower"),
    CI_upper_bs = purrr::map_dbl(bootstrap_ci, "ci_upper"),
    CI_lower_bs_residual_var = purrr::map_dbl(bootstrap_ci_residual_var, "ci_lower"),
    CI_upper_bs_residual_var = purrr::map_dbl(bootstrap_ci_residual_var, "ci_upper"),
    # CI_lower_bs_residual_var_prop = purrr::map_dbl(bootstrap_ci_residual_var_prop, "ci_lower"),
    # CI_upper_bs_residual_var_prop = purrr::map_dbl(bootstrap_ci_residual_var_prop, "ci_upper"),
    CI_lower_sandwich = map_dbl(rho_sandwich_inference, function(x)
      x$ci[[1]]),
    CI_upper_sandwich = map_dbl(rho_sandwich_inference, function(x)
      x$ci[[2]]),
    rho_se = map_dbl(rho_sandwich_inference, function(x)
      x$se)
  )

# Summarize inferences in a table.
surrogate_results_tbl = surrogate_results_tbl %>%
  select(-moment_estimate, -bootstrap_ci, -bootstrap_ci_residual_var, -rho_sandwich_inference, -data_tbl) 

surrogate_results_tbl %>%
  write.csv(file = paste0(tables_dir, "/surrogacy-inferences.csv"))

