# Setup ------------------------------------------------------------------
library(tidyverse)
library(ggridges)

# Specify options for saving the plots to files
figures_dir = "results/figures/application/meta-analysis"
tables_dir = "results/tables/application/meta-analysis"

# Read in nonparametric results
surrogate_results_tbl = read.csv(file = "results/tables/application/meta-analysis/surrogacy-inferences.csv")
# Read in Bayesian results
surrogate_results_bayesian_tbl = readRDS("R/application/bayesian_ma_results.rds")

# Plots -----------------------------------------------------------------

## Non-Parametric MA ----------------------------------------------------

# Helper function to make plots.
conf_int_plot_f = function(include_risk_score, type) {
  plotting_data = surrogate_results_tbl %>%
    filter(include_risk_score == .env$include_risk_score |
             method == "untransformed surrogate") 
  if (type == "bs") {
    plotting_data = plotting_data %>%
      rename(CI_lower = CI_lower_bs, 
             CI_upper = CI_upper_bs)
  } else {
    plotting_data = plotting_data %>%
      rename(CI_lower = CI_lower_sandwich, 
             CI_upper = CI_upper_sandwich)
  }

  plotting_data = plotting_data %>%
    mutate(rho_trial = ifelse(method == "untransformed surrogate", -1 * rho_trial, rho_trial),
           CI_lower = ifelse(method == "untransformed surrogate", -1 * CI_lower, CI_lower),
           CI_upper = ifelse(method == "untransformed surrogate", -1 * CI_upper, CI_upper))
  
  
  conf_int_plot  = plotting_data %>%
    ggplot(aes(x = analysis_set, color = method)) +
    geom_point(aes(y = rho_trial), position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), position = position_dodge(width = 0.5), width = 0.2) +
    coord_cartesian(ylim = c(-1, 1)) +
    facet_grid(surrogate~.) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin())
  
  risk_score_chr = ifelse(include_risk_score, "w-riskscore", "wo-riskscore")
  type_chr = ifelse(type == "bs", "bs", "sandwich")
  outfile = paste0("/surrogacy-measures-summary-", risk_score_chr, "-", type_chr, ".pdf")
  
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

conf_int_plot_f(TRUE, "bs")
conf_int_plot_f(TRUE, "sandwich")
conf_int_plot_f(FALSE, "bs")
conf_int_plot_f(FALSE, "sandwich")

## Bayesian MA -------------------------------------------------------------

# Extract and reshape posterior samples
rho_long_tbl <- surrogate_results_bayesian_tbl %>%
  mutate(rho_samples = map(stan_fit, ~ as.data.frame(.x)$rho)) %>%
  unnest(rho_samples) %>%
  rename(rho = rho_samples)

# Define helper function to make and save plots.
posterior_plots_f = function(assume_proportional_line, include_risk_score) {
  posterior_plot = rho_long_tbl %>% filter(
    method != "untransformed surrogate",
    assume_proportional_line == .env$assume_proportional_line,
    include_risk_score == .env$include_risk_score
  ) %>%
    ggplot(aes(x = rho, y = method, fill = method)) +
    geom_density_ridges(quantile_lines = TRUE,
                        quantiles = c(0.025, 0.5, 0.975)) +
    facet_grid(surrogate ~ analysis_set) +
    labs(title = "Posterior Distributions of Trial-Level Correlation (rho)",
         x = "rho",
         y = NULL,
         fill = "Method") +
    theme(legend.position = "bottom")
  
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
posterior_plots_f(TRUE, TRUE)
posterior_plots_f(TRUE, FALSE)
posterior_plots_f(FALSE, TRUE)
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
