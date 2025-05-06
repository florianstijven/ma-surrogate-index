library(rstan)
library(bayesplot)
library(tidyverse)
library(tidybayes)
library(posterior)

# Specify options for saving the plots to files
figures_dir = "results/figures/application/meta-analysis"
tables_dir = "results/tables/application/meta-analysis"

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

fit_surrogacy_model <- function(data, iter = 10000, warmup = 5000, chains = 4, seed = 123) {
  
  stan_code <- "
  data {
    int<lower=1> N;              // number of studies
    vector[2] y[N];              // observed effects: [S, T]
    cov_matrix[2] Sigma[N];      // known within-trial covariance matrices
  }
  parameters {
    real mu;                     // mean treatment effect on surrogate
    vector<lower=0>[2] tau;      // between-trial SDs
    real<lower=-1, upper=1> rho; // trial-level correlation
    vector[2] theta[N];          // latent true effects per trial
  }
  model {
    mu ~ normal(-2, 3);           
    tau ~ normal(0, 2);        
    rho ~ uniform(0, 1);        

    matrix[2,2] Omega;
    Omega[1,1] = 1;
    Omega[1,2] = rho;
    Omega[2,1] = rho;
    Omega[2,2] = 1;

    real mu_beta;
    mu_beta = rho * (tau[2] / tau[1]) * mu;

    vector[2] mu_full;
    mu_full[1] = mu;
    mu_full[2] = mu_beta;

    for (i in 1:N)
      theta[i] ~ multi_normal(mu_full, quad_form_diag(Omega, tau));

    for (i in 1:N)
      y[i] ~ multi_normal(theta[i], Sigma[i]);
  }
  "
  
  # Build input data for Stan
  N <- nrow(data)
  y_matrix <- cbind(data$trt_effect_surrogate_index_est, data$log_RR_est)
  
  Sigma_array <- array(NA, dim = c(N, 2, 2))
  for (i in 1:N) {
    Sigma_array[i,,] <- data$vcov[[i]]
  }
  
  stan_data <- list(
    N = N,
    y = y_matrix,
    Sigma = Sigma_array
  )
  
  # Fit model
  fit <- stan(
    model_code = stan_code,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed,
    control = list(adapt_delta = 0.95)
  )
  
  return(fit)
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
  filter(method != "untransformed surrogate", !include_risk_score) %>%
  bind_rows(ma_trt_effects_tbl_untransformed) %>%
  # Make sure that only the correct trials are included for each analysis set.
  filter(ifelse(
    analysis_set == "first_four",
    trial %in% trials_first_four,
    ifelse(analysis_set == "naive_only", trial %in% trials_naive_only, TRUE)
  )) %>%
  group_by(surrogate, method, weighting, analysis_set, include_risk_score) %>%
  summarise(stan_fit = list(fit_surrogacy_model(data = pick(everything()), chains = 1, iter = 2e3, warmup = 1e3)))


# Extract rho samples and prepare for plotting
rho_posterior_tbl <- surrogate_results_tbl %>%
  mutate(draws = purrr::map(stan_fit, ~ as_draws_df(.x))) %>%
  mutate(rho_draws = purrr::map(draws, ~ .x$rho)) %>%
  select(surrogate, method, analysis_set, rho_draws, weighting) %>%
  unnest(cols = c(rho_draws)) %>%
  rename(rho = rho_draws)

# Plotting posterior distributions with 90% and 95% intervals
ggplot(rho_posterior_tbl, aes(x = rho, fill = method, color = method)) +
  stat_halfeye(
    .width = c(0.95, 0.90),
    point_interval = mean_qi,
    slab_alpha = 0.7
  ) +
  facet_grid(surrogate ~ .) +
  labs(
    title = "Posterior Distributions of ρ (rho)",
    x = expression(rho),
    y = NULL,
    fill = "Method",
    color = "Method"
  )

library(tidyverse)
library(bayesplot)

# Extract and summarize posterior for rho from each model
rho_summary_tbl <- surrogate_results_tbl %>%
  mutate(rho_samples = map(stan_fit, ~ as.data.frame(.x)$rho)) %>%
  unnest(rho_samples) %>%
  group_by(surrogate, method, analysis_set) %>%
  summarise(
    median = median(rho_samples),
    ci_90_low = quantile(rho_samples, 0.05),
    ci_90_high = quantile(rho_samples, 0.95),
    ci_95_low = quantile(rho_samples, 0.025),
    ci_95_high = quantile(rho_samples, 0.975),
    .groups = "drop"
  )

# Plot using ggplot2
ggplot(rho_summary_tbl, aes(x = median, y = method, color = method)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = ci_95_low, xmax = ci_95_high), width = 0.2, size = 0.8) +
  geom_errorbar(aes(xmin = ci_90_low, xmax = ci_90_high), width = 0.5, size = 1.2) +
  facet_grid(surrogate ~ analysis_set) +
  labs(
    title = "Posterior Estimates of Trial-Level Correlation (rho)",
    x = "Posterior Median and Credible Intervals for rho",
    y = NULL,
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

library(tidyverse)
library(ggridges)

# Extract and reshape posterior samples
rho_long_tbl <- surrogate_results_tbl %>%
  mutate(rho_samples = map(stan_fit, ~ as.data.frame(.x)$rho)) %>%
  unnest(rho_samples) %>%
  rename(rho = rho_samples)

# Plot full posterior distributions
rho_long_tbl %>% filter(method != "untransformed surrogate") %>%
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




