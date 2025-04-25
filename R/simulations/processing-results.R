# Preliminaries ----------------------------------------------------
# load packages
library(tidyverse)
library(ggpubr)

# Specify options for saving the plots to files
figures_dir = "results/figures/simulations/"
tables_dir = "results/tables/simulations/"

# Read in simulation results
ma_sim_results <- bind_rows(
  readRDS(
    "results/raw-results/simulations/ma-sim-results-proof-of-concept-large.rds"
  ) %>%
    mutate(scenario = "proof-of-concept", setting = "large N, small n"),
  readRDS(
    "results/raw-results/simulations/ma-sim-results-proof-of-concept-small.rds"
  ) %>%
    mutate(scenario = "proof-of-concept", setting = "small N, large n"),
  readRDS(
    "results/raw-results/simulations/ma-sim-results-vaccine-small.rds"
  ) %>%
    mutate(scenario = "vaccine", setting = "small N, large n")
)

# The sandwich CIs are missing whenever the trial-level correlation estimate not
# in the [-1, 1] interval. Furthermore, the sandwich CI is just the entire
# interval if the correlation estimate is -1 or 1. These "issues" are caused by
# computing the CI on the Fisher's Z scale and then backtransforming. In these
# causes, we computed the CI using the sandwich SE on the original correlation
# scale. It is always possible to compute this CI, but it will have limits
# outside the [-1, 1] interval.
ma_sim_results = bind_rows(
  ma_sim_results %>%
    filter(CI_type != "sandwich" | is.na(CI_type)),
  ma_sim_results %>%
    filter(CI_type == "sandwich" & !(abs(rho_est) > 0.9999 | is.nan(rho_ci_lower))),
  ma_sim_results %>%
    filter(CI_type == "sandwich" & (abs(rho_est) > 0.9999 | is.nan(rho_ci_lower))) %>%
    mutate(rho_ci_lower = rho_est - qt(0.975, df = N - 1) * rho_se,
           rho_ci_upper = rho_est + qt(0.975, df = N - 1))
)

# Code character columns as factors.
ma_sim_results = ma_sim_results %>%
  mutate(
    surrogate_index_estimator = as.factor(surrogate_index_estimator),
    SI_violation = as.factor(SI_violation),
    estimator_adjustment = as.factor(estimator_adjustment),
    sandwich_adjustment = as.factor(sandwich_adjustment),
    CI_type = as.factor(CI_type),
    setting = as.factor(setting),
    scenario = as.factor(scenario),
  )


# Compute summaries from the simulation results
ma_sim_summary = ma_sim_results %>%
  group_by(
    surrogate_index_estimator,
    SI_violation,
    estimator_adjustment,
    sandwich_adjustment,
    CI_type,
    setting,
    scenario,
    N,
    n, 
    nearest_PD
  ) %>%
  summarise(
    coverage = mean((rho_true >= rho_ci_lower) &
                      (rho_true <= rho_ci_upper), na.rm = TRUE),
    mean_bias = mean(rho_est - rho_true, na.rm = TRUE),
    median_bias = median(rho_est - rho_true, na.rm = TRUE),
    emp_SD = sqrt(mean((rho_est - rho_true)^2, na.rm = TRUE)),
    mean_SE = mean(rho_se, na.rm = TRUE),
    MSE = mean((rho_true - rho_est) ^ 2, na.rm = TRUE),
    power_0.75 = mean(rho_ci_lower >= 0.75)
  )

# Summary of the Simulation Results ---------------------------------------

## Distribution of the estimands ------------------------------------------
estimand_plot_1 = ma_sim_results %>%
  filter(
    surrogate_index_estimator != "surrogate",
    setting == "small N, large n",
    CI_type == "sandwich",
    scenario == "proof-of-concept"
  ) %>%
  mutate(surrogate_index_estimator = fct_drop(surrogate_index_estimator, only = c("surrogate"))) %>%
  ggplot(aes(x = rho_true, color = surrogate_index_estimator)) +
  geom_density(show.legend = TRUE) +
  geom_vline(
    aes(xintercept = rho_true),
    data = ma_sim_results %>%
      filter(
        surrogate_index_estimator == "surrogate",
        setting == "small N, large n"
      ) %>%
      mutate(rho_true = abs(rho_true)) %>%
      filter(scenario == "proof-of-concept") %>%
      group_by(scenario, N, SI_violation) %>%
      slice_head(n = 1) %>%
      ungroup()
  ) +
  scale_x_continuous(lim = c(0.4, 1), name = expr(rho[trial])) +
  scale_color_discrete(name = "Surr. Index Estimator", drop = FALSE) +
  facet_grid(SI_violation ~ N, scales = "free")

estimand_plot_1
ggsave(
  filename = "distribution-estimands-proof-of-concept.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)



estimand_plot_2 = ma_sim_results %>%
  filter(
    surrogate_index_estimator != "surrogate",
    setting == "small N, large n",
    CI_type == "sandwich",
    scenario == "vaccine"
  ) %>%
  ggplot(aes(x = rho_true, color = surrogate_index_estimator)) +
  geom_density(show.legend = TRUE) +
  geom_vline(
    aes(xintercept = rho_true),
    data = ma_sim_results %>%
      filter(
        surrogate_index_estimator == "surrogate",
        setting == "small N, large n"
      ) %>%
      mutate(rho_true = abs(rho_true)) %>%
      filter(scenario == "vaccine") %>%
      group_by(scenario, N, SI_violation) %>%
      slice_head(n = 1) %>%
      ungroup()
  ) +
  scale_x_continuous(lim = c(0.82, 1), name = expr(rho[trial])) +
  scale_color_discrete(name = "Surr. Index Estimator", drop = FALSE) +
  facet_grid(SI_violation ~ N, scales = "free")

estimand_plot_2
ggsave(
  filename = "distribution-estimands-vaccine.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

# ggarrange(estimand_plot_1, estimand_plot_2, common.legend = TRUE, legend = "bottom", labels = "auto")

# ggsave(
#   filename = "distribution-estimands.pdf",
#   path = figures_dir,
#   height = double_height,
#   width = double_width,
#   dpi = res,
#   device = "pdf",
#   units = "cm"
# )


## Estimation Accuracies --------------------------------------------------

### Bias ------------------------------------------------------------------

ma_sim_summary %>%
  filter(setting == "small N, large n", CI_type == "sandwich") %>%
  ggplot(aes(
    x = N,
    y = mean_bias,
    color = surrogate_index_estimator,
    linetype = nearest_PD
  )) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0) +
  scale_x_continuous(breaks = c(6, 12, 24)) +
  scale_y_continuous(name = "Mean Bias") +
  scale_color_discrete(name = "Surr. Index Estimator") +
  facet_grid(SI_violation ~ scenario) + 
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(
  filename = "mean-bias.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

ma_sim_summary %>%
  filter(setting == "small N, large n", CI_type == "sandwich") %>%
  ggplot(aes(
    x = N,
    y = median_bias,
    color = surrogate_index_estimator,
    linetype = nearest_PD
  )) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0) +
  scale_x_continuous(breaks = c(6, 12, 24)) +
  scale_y_continuous(name = "Median Bias") +
  scale_color_discrete(name = "Surr. Index Estimator") +
  facet_grid(SI_violation ~ scenario) + 
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(
  filename = "median-bias.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

### MSE -------------------------------------------------------------------

ma_sim_summary %>%
  filter(setting == "small N, large n", CI_type == "sandwich") %>%
  ggplot(aes(
    x = N,
    y = MSE,
    color = surrogate_index_estimator,
    linetype = nearest_PD
  )) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 0) +
  scale_x_continuous(breaks = c(6, 12, 24)) +
  scale_y_continuous(name = "MSE", transform = "log10") +
  scale_color_discrete(name = "Surr. Index Estimator") +
  facet_grid(SI_violation ~ scenario) + 
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(
  filename = "mse.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)


## Performance of Inferences ----------------------------------------------

### Coverage --------------------------------------------------------------

ma_sim_summary %>%
  filter(setting == "small N, large n", CI_type == "sandwich") %>%
  ggplot(aes(
    x = N,
    y = coverage,
    color = surrogate_index_estimator,
    linetype = nearest_PD
  )) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_line() +
  geom_abline(intercept = 0.95, slope = 0) +
  scale_x_continuous(breaks = c(6, 12, 24)) +
  scale_y_continuous(name = "Coverage") +
  scale_color_discrete(name = "Surr. Index Estimator") +
  facet_grid(SI_violation ~ scenario) + 
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(
  filename = "coverage-sandwich.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

ma_sim_summary %>%
  filter(setting == "small N, large n", CI_type == "multiplier bootstrap") %>%
  ggplot(aes(
    x = N,
    y = coverage,
    color = surrogate_index_estimator,
    linetype = nearest_PD
  )) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_line() +
  geom_abline(intercept = 0.95, slope = 0) +
  scale_x_continuous(breaks = c(6, 12, 24)) +
  scale_y_continuous(name = "Coverage") +
  scale_color_discrete(name = "Surr. Index Estimator") +
  facet_grid(SI_violation ~ scenario) + 
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(
  filename = "coverage-bootstrap.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

### Distribution of lower limits ------------------------------------------

ma_sim_results %>%
  filter(
    surrogate_index_estimator != "surrogate",
    setting == "small N, large n",
    scenario == "proof-of-concept",
    nearest_PD == TRUE,
    CI_type == "multiplier bootstrap"
  ) %>%
  mutate(surrogate_index_estimator = fct_drop(surrogate_index_estimator, only = c("surrogate"))) %>%
  ggplot(aes(x = rho_ci_lower), fill = "gray") +
  geom_histogram(alpha = 0.5, position = "identity", color = "black") +
  scale_x_continuous(lim = c(-1, 1), name = expr(rho[trial])) +
  scale_color_discrete(name = "Surr. Index Estimator", drop = FALSE) +
  facet_grid(SI_violation ~ N, scales = "free") +
  theme(legend.position = "bottom")

ggsave(
  filename = "lower-limit-proof-of-concept-nearest-PD-bootstrap.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

ma_sim_results %>%
  filter(
    setting == "small N, large n",
    scenario == "vaccine",
    nearest_PD == TRUE,
    CI_type == "multiplier bootstrap"
  ) %>%
  mutate(surrogate_index_estimator = fct_drop(surrogate_index_estimator, only = c("surrogate"))) %>%
  # There is a negative trial-level correlation for the untransformed surrogate.
  # We therefore use -1 * the upper limit for this setting to be comparable to
  # with the estimated surrogate index.
  mutate(
    rho_ci_lower = ifelse(
      surrogate_index_estimator == "surrogate",
      -1 * rho_ci_upper,
      rho_ci_lower
    )
  ) %>%
  ggplot(aes(x = rho_ci_lower, fill = SI_violation)) +
  geom_histogram(alpha = 0.5,
                 position = "identity",
                 color = "black",
                 bins = 20) +
  scale_x_continuous( name = expr(rho[trial])) +
  scale_color_discrete(name = "Surr. Index Estimator", drop = FALSE) +
  facet_grid(surrogate_index_estimator ~ N, scales = "free") +
  theme(legend.position = "bottom")

ggsave(
  filename = "lower-limit-vaccine-nearest-PD-bootstrap.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)

### SE estimator ----------------------------------------------------------


## Failure of Asymptotics -------------------------------------------------

ma_sim_results %>% filter(setting == "large N, small n") %>%
  mutate(
    Z = (rho_est - rho_true) / rho_se,
    N_f = factor(
      as.character(N),
      levels = c("5000", "2000", "1000", "500", "100"),
      ordered = TRUE
    )
  ) %>%
  ggplot(aes(x = Z, fill = N_f)) +
  geom_histogram(aes(y = after_stat(density)), color = "black") +
  # Add standard normal density.
  geom_line(aes(x = x, y = d), 
            data = tibble(
    x = seq(
      from = -4,
      to = 4,
      length.out = 1e3
    ),
    d = dnorm(seq(
      from = -4,
      to = 4,
      length.out = 1e3
    ))
  ) %>% cross_join(tibble(
    N_f = factor(
      as.character(c(5000, 2000, 1000, 500, 100)),
      levels = c("5000", "2000", "1000", "500", "100"),
      ordered = TRUE
    )
  ))) +
  xlim(c(-4, 4)) +
  facet_grid(N_f ~ n)

ggsave(
  filename = "failure-asymptotics.pdf",
  path = figures_dir,
  height = double_height,
  width = double_width,
  dpi = res,
  device = "pdf",
  units = "cm"
)
