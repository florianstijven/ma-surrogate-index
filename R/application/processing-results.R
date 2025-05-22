# Setup ------------------------------------------------------------------
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)

# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

# The first argument indicates whether the analysis should be conducted on the
# original data or on the synthetic data.
data_set = args[1]
if (data_set == "real") {
  ma_trt_effects_tbl_location = "results/raw-results/application/ma_trt_effects_tbl.rds"
  rho_long_tbl_location = "results/raw-results/application/rho_long_tbl.rds"
  # Specify options for saving the plots to files
  figures_dir = "results/figures/application/meta-analysis"
  tables_dir = "results/tables/application/meta-analysis"
} else if (data_set == "synthetic") {
  ma_trt_effects_tbl_location = "results/raw-results/application-synthetic/ma_trt_effects_tbl.rds"
  rho_long_tbl_location = "results/raw-results/application-synthetic/rho_long_tbl.rds"
  # Specify options for saving the plots to files
  figures_dir = "results/figures/application-synthetic/meta-analysis"
  tables_dir = "results/tables/application-synthetic/meta-analysis"
}

# Read in nonparametric results
surrogate_results_tbl = read.csv(file = paste0(tables_dir, "/", "surrogacy-inferences.csv"))
# Change some variables to facilite plotting.
surrogate_results_tbl = surrogate_results_tbl %>%
  mutate(`Estimated Surrogate Index` = case_when(
    method == "cox" ~ "Cox Model",
    method == "sl" ~ "SuperLearner",
    method == "untransformed surrogate" ~ "Ab Marker"
  ),
  surrogate_name = case_when(
    surrogate == "bindSpike" ~ "IgG Spike",
    surrogate == "pseudoneutid50" ~ "nAb ID50",
    surrogate == "pseudoneutid50_adjusted" ~ "adjusted nAb ID50"
  )) %>%
  mutate(
    surrogate_name = factor(surrogate_name, levels = c("IgG Spike", "nAb ID50", "adjusted nAb ID50"), ordered = TRUE)
  )
# Read in Bayesian results
rho_long_tbl = readRDS("results/raw-results/application/rho_long_tbl.rds") %>%
  mutate(`Estimated Surrogate Index` = case_when(
    method == "cox" ~ "Cox Model",
    method == "sl" ~ "SuperLearner",
    method == "untransformed surrogate" ~ "Ab Marker"
  ),
  surrogate_name = case_when(
    surrogate == "bindSpike" ~ "IgG Spike",
    surrogate == "pseudoneutid50" ~ "nAb ID50",
    surrogate == "pseudoneutid50_adjusted" ~ "adjusted nAb ID50"
  )) %>%
  mutate(
    surrogate_name = factor(surrogate_name, levels = c("IgG Spike", "nAb ID50", "adjusted nAb ID50"), ordered = TRUE)
  )

# Load data with trial-level treatment effects. 
ma_trt_effects_tbl = readRDS(ma_trt_effects_tbl_location)


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

# Tables ----------------------------------------------------------------

# Save tables for the trial-level correlation estimates with corresponding
# confidence intervals. Tables are saves separately by type of CI and whether the
# original data or pseudo-real data were used.
format_CI = function(estimate, ci_lower, ci_upper, digits = 2) {
  paste0(round(estimate, digits), " (", round(ci_lower, digits), ", ", round(ci_upper, digits), ")")
}

save_inferences_table = function(CI_type, data_type) {
  if (CI_type == "BCa") {
    table_temp = surrogate_results_tbl %>%
      rename(CI_lower = CI_lower_bs, CI_upper = CI_upper_bs)
  } else if (CI_type == "sandwich") {
    table_temp = surrogate_results_tbl %>%
      rename(CI_lower = CI_lower_sandwich, CI_upper = CI_upper_sandwich)
  }
  
  if (data_type == "real-data") {
    table_temp = table_temp %>% filter(scenario == "real data")
  } else {
    table_temp = table_temp %>%
      filter(scenario != "real data")
  }
  
  outfile = paste0(tables_dir, "/table-inferences-", CI_type, "-", data_type, ".csv")

  table_temp %>%
    filter(analysis_set == "naive_only") %>%
    rowwise(everything()) %>%
    mutate(rho_inference = format_CI(rho_trial, CI_lower, CI_upper)) %>%
    ungroup() %>%
    select(surrogate, method, rho_inference, scenario) %>%
    pivot_wider(values_from = "rho_inference", names_from = "method") %>%
    write.csv(file = outfile)
}

save_inferences_table("BCa", "real-data")
save_inferences_table("BCa", "pseudo-data")
save_inferences_table("sandwich", "real-data")
save_inferences_table("sandwich", "pseudo-data")
 


# Plots -----------------------------------------------------------------

## MA plots -------------------------------------------------------------

# Add proper name of the surrogates.
ma_trt_effects_tbl = ma_trt_effects_tbl %>%
  mutate(`Estimated Surrogate Index` = case_when(
    method == "cox" ~ "Cox Model",
    method == "sl" ~ "SuperLearner",
    method == "untransformed surrogate" ~ "Ab Marker"
  ),
  surrogate_name = case_when(
    surrogate == "bindSpike" ~ "IgG Spike",
    surrogate == "pseudoneutid50" ~ "nAb ID50",
    surrogate == "pseudoneutid50_adjusted" ~ "adjusted nAb ID50"
  )) %>%
  mutate(
    surrogate_name = factor(surrogate_name, levels = c("IgG Spike", "nAb ID50", "adjusted nAb ID50"), ordered = TRUE)
  )

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

# Define transformation for VE scale.
transform_VE = new_transform(
  name = "VE",
  transform = function(x) -1 * log(1 - x),
  inverse = function(x) 1 - exp(-x)
)

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
    facet_grid(`Estimated Surrogate Index` ~ surrogate_name) +
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
  scale_color_manual(values = my_palette, limits = trials_mixed_fct[!(trials_mixed_fct %in% trials_naive_only_fct)]) +
  scale_shape_manual(values = my_shapes, limits = trials_mixed_fct[!(trials_mixed_fct %in% trials_naive_only_fct)]) +
  scale_y_continuous(transform = transform_VE, breaks = c(0, 0.5, 0.75, 0.90, 0.95)) +
  scale_x_continuous(transform = transform_VE, breaks = c(0, 0.5, 0.75, 0.90, 0.95)) +
  coord_cartesian(ylim = c(-0.2, 0.95), xlim = c(-0.2, 0.95)) +
  facet_grid(`Estimated Surrogate Index` ~ surrogate_name) +
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



## Non-Parametric MA ----------------------------------------------------

# Helper function to make plots.
conf_int_plot_f = function(include_risk_score, type, res_var_prop, scenario) {
  plotting_data = surrogate_results_tbl %>%
    filter(include_risk_score == .env$include_risk_score |
             method == "untransformed surrogate") %>%
    filter(scenario == .env$scenario, analysis_set == "naive_only")
  if (res_var_prop) {
    plotting_data = plotting_data %>%
      rename(CI_lower = CI_lower_bs_residual_var_prop, CI_upper = CI_upper_bs_residual_var_prop) %>%
      filter(method != "untransformed surrogate")
  } else {
    if (type == "bs") {
      plotting_data = plotting_data %>%
        rename(CI_lower = CI_lower_bs, CI_upper = CI_upper_bs)
    } else {
      plotting_data = plotting_data %>%
        rename(CI_lower = CI_lower_sandwich, CI_upper = CI_upper_sandwich)
    }
    plotting_data = plotting_data %>%
      mutate(
        rho_trial = ifelse(method == "untransformed surrogate", -1 * rho_trial, rho_trial),
        CI_lower = ifelse(method == "untransformed surrogate", -1 * CI_lower, CI_lower),
        CI_upper = ifelse(method == "untransformed surrogate", -1 * CI_upper, CI_upper)
      )
  }


  conf_int_plot = plotting_data %>%
    ggplot(aes(x = surrogate_name, color = `Estimated Surrogate Index`)) +
    geom_point(aes(y = rho_trial), position = position_dodge(width = 0.5)) +
    geom_errorbar(
      aes(ymin = CI_lower, ymax = CI_upper),
      position = position_dodge(width = 0.5),
      width = 0.2
    ) +
    coord_cartesian(ylim = c(-1, 1)) +
    scale_y_continuous(name = expr(rho[trial])) +
    scale_x_discrete(name = NULL) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(),
      legend.text = element_text(margin = margin(l = 0))
    ) +
    guides(color = guide_legend(direction = "vertical", nrow = 1))
  
  risk_score_chr = ifelse(include_risk_score, "w-riskscore", "wo-riskscore")
  type_chr = ifelse(type == "bs", "bs", "sandwich")
  res_var_prop_chr = ifelse(res_var_prop, "res_var_prop", "cor")
  outfile = paste0("/surrogacy-measures-summary-", risk_score_chr, "-", type_chr, "-", res_var_prop_chr, ".pdf")
  
  ggsave(
    outfile,
    path = figures_dir,
    height = single_height,
    width = single_width,
    dpi = res,
    device = "pdf",
    units = "cm"
  )
}

# conf_int_plot_f(TRUE, "bs", FALSE)
# conf_int_plot_f(TRUE, "sandwich", FALSE)
conf_int_plot_f(FALSE, "bs", FALSE, "real data")
conf_int_plot_f(FALSE, "sandwich", FALSE, "real data")

# conf_int_plot_f(TRUE, "bs", TRUE)
# conf_int_plot_f(FALSE, "bs", TRUE, "real data")

## Bayesian MA -------------------------------------------------------------

# Define helper function to make and save plots.
posterior_plots_f = function(assume_proportional_line, include_risk_score) {
  posterior_plot = rho_long_tbl %>% filter(
    assume_proportional_line == .env$assume_proportional_line,
    include_risk_score == .env$include_risk_score | is.na(include_risk_score),
    analysis_set == "naive_only"
  ) %>%
    # The trial-level correlation is "reversed" for the untransformed Ab Marker.
    # To have a better comparison, we multiply the corresponding trial-level
    # correlation with -1.
    mutate(
      rho = ifelse(`Estimated Surrogate Index` == "Ab Marker", -1 * rho, rho) 
    ) %>%
    ggplot(aes(x = rho, y = `Estimated Surrogate Index`, fill = `Estimated Surrogate Index`)) +
    geom_density_ridges(quantile_lines = TRUE,
                        quantiles = c(0.025, 0.5, 0.975)) +
    facet_grid(surrogate_name~.) +
    scale_y_discrete(labels = NULL) +
    labs(title = "Posterior Distributions of Trial-Level Correlation",
         x = expr(rho[trial]),
         y = NULL) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin()) +
    guides(fill = guide_legend(direction = "vertical", nrow = 1))
  
  assume_proportional_line_chr = ifelse(assume_proportional_line, "prop-line", "default")
  risk_score_chr = ifelse(include_risk_score, "w-riskscore", "wo-riskscore")
  outfile = paste0("/posterior-distributions-", assume_proportional_line_chr, "-", risk_score_chr, ".pdf")

  ggsave(
    outfile,
    path = figures_dir,
    height = double_height,
    width = double_width,
    dpi = res,
    device = "pdf",
    units = "cm"
  )
}

# Make plots for specified scenarios
# posterior_plots_f(TRUE, TRUE)
posterior_plots_f(TRUE, FALSE)
# posterior_plots_f(FALSE, TRUE)
posterior_plots_f(FALSE, FALSE)

# Compute posterior mean, median, and quantiles.
rho_long_tbl %>%
  group_by(
    surrogate,
    method,
    weighting,
    analysis_set,
    include_risk_score,
    assume_proportional_line
  ) %>%
  summarise(
    mean = mean(rho),
    median = median(rho),
    quantile_2.5 = quantile(rho, 0.025),
    quantile_5 = quantile(rho, 0.05),
    quantile_95 = quantile(rho, 0.95),
    quantile_97.5 = quantile(rho, 0.975)
  ) %>%
  write.csv(paste0(tables_dir, "/posterior-summaries.csv"))

