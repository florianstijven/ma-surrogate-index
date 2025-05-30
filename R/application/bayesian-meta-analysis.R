library(rstan)
library(bayesplot)
library(tidyverse)
library(tidybayes)
library(posterior)
library(furrr)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}

# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

# The first argument indicates whether the analysis should be conducted on the
# original data or on the synthetic data.
data_set = args[1]
if (data_set == "real") {
  ma_trt_effects_tbl_location = "results/raw-results/application/ma_trt_effects_tbl.rds"
  out_file1 = "results/raw-results/application/bayesian_ma_results.rds"
  out_file2 = "results/raw-results/application/rho_long_tbl.rds"
  
  # Specify options for saving the plots to files
  figures_dir = "results/figures/application/meta-analysis"
  tables_dir = "results/tables/application/meta-analysis"
} else if (data_set == "synthetic") {
  ma_trt_effects_tbl_location = "results/raw-results/application-synthetic/ma_trt_effects_tbl.rds"
  out_file1 = "results/raw-results/application-synthetic/bayesian_ma_results.rds"
  out_file2 = "results/raw-results/application-synthetic/rho_long_tbl.rds"
  
  # Specify options for saving the plots to files
  figures_dir = "results/figures/application-synthetic/meta-analysis"
  tables_dir = "results/tables/application-synthetic/meta-analysis"
}

# Source Helper function to fit Bayesian Model. 
source("R/helper-functions/bayesian-model.R")

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
ma_trt_effects_tbl = readRDS(ma_trt_effects_tbl_location)

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
  bind_rows(ma_trt_effects_tbl_untransformed) %>%
  # Make sure that only the correct trials are included for each analysis set.
  filter(ifelse(
    analysis_set == "first_four",
    trial %in% trials_first_four,
    ifelse(analysis_set == "naive_only", trial %in% trials_naive_only, TRUE)
  )) %>%
  group_by(surrogate, method, weighting, analysis_set, include_risk_score) %>%
  summarise(data = list(pick(everything()))) %>%
  cross_join(tibble(assume_proportional_line = c(FALSE, TRUE))) %>%
  mutate(
    stan_fit = future_map2(
      .x = data,
      .y = assume_proportional_line,
      .f = fit_surrogacy_model,
      chains = 4,
      iter = 2e4,
      warmup = 5e3,
      .options = furrr_options(seed = TRUE)
    )
  )

# Saving Results ----------------------------------------------------------

# Save fitted Bayesian models. 
saveRDS(surrogate_results_tbl, file = out_file1)

# The fitted Bayesian models are very large objects. So, we also save a subset
# of this information next, which is more convenient to work with.
rho_long_tbl <- surrogate_results_tbl %>%
  mutate(rho_samples = map(stan_fit, ~ as.data.frame(.x)$rho)) %>%
  unnest(rho_samples) %>%
  rename(rho = rho_samples) %>%
  select(-stan_fit, -data)

saveRDS(rho_long_tbl, file = out_file2)



