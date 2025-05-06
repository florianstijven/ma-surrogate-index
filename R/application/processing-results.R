# Setup ------------------------------------------------------------------
library(tidyverse)
library(ggridges)

# Specify options for saving the plots to files
figures_dir = "results/figures/application/meta-analysis"
tables_dir = "results/tables/application/meta-analysis"

# Read in Bayesian results
surrogate_results_tbl = readRDS("R/application/bayesian_ma_results.rds")

# Plots -----------------------------------------------------------------

# Extract and reshape posterior samples
rho_long_tbl <- surrogate_results_tbl %>%
  mutate(rho_samples = map(stan_fit, ~ as.data.frame(.x)$rho)) %>%
  unnest(rho_samples) %>%
  rename(rho = rho_samples)

# Define helper function to make and save plots.
posterior_plots_f = function(assume_proportional_line) {
  posterior_plot = rho_long_tbl %>% filter(method != "untransformed surrogate", assume_proportional_line == .env$assume_proportional_line) %>%
    ggplot(aes(x = rho, y = method, fill = method)) +
    geom_density_ridges(
      quantile_lines = TRUE,
      quantiles = c(0.025, 0.5, 0.975)
    ) +
    facet_grid(surrogate ~ analysis_set) +
    labs(title = "Posterior Distributions of Trial-Level Correlation (rho)",
         x = "rho",
         y = NULL,
         fill = "Method") +
    theme(legend.position = "bottom")
  
  if (assume_proportional_line) {
    outfile = "/posterior-distributions-prop-line.pdf"
  } else {
    outfile = "/posterior-distributions.pdf"
  }
  
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
posterior_plots_f(TRUE)
posterior_plots_f(FALSE)
