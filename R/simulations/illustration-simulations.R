# Preliminaries -----------------------------------------------------------
library(tidyverse)
library(ggpubr)
# Source helper functions.
source("R/helper-functions/simulation-functions.R")
source("R/helper-functions/train-clinical-prediction-models.R")
# Specify options for saving the plots to files
figures_dir = "results/figures/simulations/"
tables_dir = "results/tables/simulations/"

set.seed(1)

# Number of Independent trials
N = 6

# Tibble with information about the data to simulate
simulation_parameters_tbl = tibble(
  N = N,
  n = c(2e3, 2e3, 5e3, 5e3),
  scenario = c("proof-of-concept", "proof-of-concept", "vaccine", "vaccine"),
  sd_beta = list(c(0.05, 0.05), c(0.1, 0.1), c(0.125, 0.125), c(0.25, 0.25)),
  violations = c("slight", "moderate", "slight", "moderate")
)


# Simulate Data -----------------------------------------------------------

# Generate IPD for the different simulation scenarios
simulated_ipd_tbl = simulation_parameters_tbl %>%
  rowwise(everything()) %>%
  summarize(ipd_data = list(simulate_trials_with_random_coefficients(N, n, sd_beta, scenario))) %>%
  ungroup()

# Estimate surrogate indices for the simulated data sets.
simulated_ipd_tbl = simulated_ipd_tbl %>%
  rowwise(everything()) %>%
  summarize(surrogate_index_f = list(train_clinical_prediction_model(
    ipd_data, ifelse(scenario == "proof-of-concept", "lm", "gam")
  ))) %>%
  ungroup()

simulated_ipd_tbl = simulated_ipd_tbl %>%
  rowwise(everything()) %>%
  # Apply surrogate index to the simulated data sets.
  summarise(
    ipd_data = ipd_data %>%
      mutate(surrogate_index = surrogate_index_f(ipd_data)) %>%
      list()
  ) %>%
  ungroup() %>%
  select(-sd_beta, -surrogate_index_f)

# Expand `simulated_ipd_tbl` to tidy format with one row per observation.
simulated_ipd_tbl = simulated_ipd_tbl %>%
  group_by(scenario, violations) %>%
  reframe(ipd_data[[1]])

# Plots ------------------------------------------------------------------

## Proof of Concept ------------------------------------------------------

# Plot to illustrate data-generating mechanism for the proof-of-concept
# scenario. We plot the surrogate against the clinical endpoint, coloring by the
# value for X1. This illustrates that the direction of the relation between the
# surrogate and clinical endpoint depends on X1.
by_trial_plots = simulated_ipd_tbl %>%
  filter(scenario == "proof-of-concept") %>%
  ggplot(aes(x = surrogate, y = clinical, color = as.factor(X1))) +
  geom_point(alpha = 0.05) +
  geom_smooth(se = FALSE) +
  facet_grid(violations~trial) +
  theme(legend.position = "none") + 
  scale_x_continuous(name = "Surrogate") +
  scale_y_continuous(name = "Clinical Endpoint")

agrregated_plots = simulated_ipd_tbl %>%
  filter(scenario == "proof-of-concept") %>%
  ggplot(aes(
    x = surrogate,
    y = clinical,
    color = as.factor(X1) ,
    group = interaction(trial, X1)
  )) +
  geom_smooth(se = FALSE) +
  facet_grid(violations ~ .) +
  theme(legend.position = "none") +
  scale_x_continuous(name = "Surrogate") +
  scale_y_continuous(name = "Clinical Endpoint")


# Combine the plots into a single one.
ggarrange(by_trial_plots, agrregated_plots, labels = "auto", nrow = 2)

ggsave(
  filename = "illustration-proof-of-concept.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)


## Vaccine ---------------------------------------------------------------

by_trial_plots = simulated_ipd_tbl %>%
  filter(scenario == "vaccine") %>%
  ggplot(aes(x = surrogate_index, y = clinical)) +
  geom_jitter(alpha = 0.01, width = 0, height = 0.05) +
  geom_smooth(se = FALSE) +
  facet_grid(violations~trial) +
  theme(legend.position = "none") +
  scale_x_continuous(name = "Surrogate Index", limits = c(0, 0.35)) + 
  scale_y_continuous(name = "Clinical Endpoint")

agrregated_plots = simulated_ipd_tbl %>%
  filter(scenario == "vaccine") %>%
  ggplot(aes(
    x = surrogate_index,
    y = clinical,
    group = interaction(trial)
  )) +
  geom_smooth(se = FALSE) +
  facet_grid(violations ~ .) +
  theme(legend.position = "none") +
  scale_x_continuous(name = "Surrogate Index") +
  scale_x_continuous(name = "Surrogate Index", limits = c(0, 0.35)) + 
  scale_y_continuous(name = "Clinical Endpoint")


# Combine the plots into a single one.
ggarrange(by_trial_plots, agrregated_plots, labels = "auto", nrow = 2)

ggsave(
  filename = "illustration-vaccine.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)















