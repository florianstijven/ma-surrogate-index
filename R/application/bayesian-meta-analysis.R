library(rstan)
library(bayesplot)
library(tidyverse)
library(tidybayes)
library(posterior)
library(furrr)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}

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

stan_code_prop_line = 
  "
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
      mu ~ normal(log(0.5), 1);           
      tau ~ normal(0, 2);        
      rho ~ uniform(-1, 1);        

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
stan_code_default = 
  "
    data {
      int<lower=1> N;              // number of studies
      vector[2] y[N];              // observed effects: [S, T]
      cov_matrix[2] Sigma[N];      // known within-trial covariance matrices
    }
    parameters {
      vector[2] mu;                     // mean treatment effect on surrogate and clinical endpoint
      vector<lower=0>[2] tau;      // between-trial SDs
      real<lower=-1, upper=1> rho; // trial-level correlation
      vector[2] theta[N];          // latent true effects per trial
    }
    model {
      mu ~ normal(log(0.5), 1);           
      tau ~ normal(0, 2);        
      rho ~ uniform(-1, 1);        

      matrix[2,2] Omega;
      Omega[1,1] = 1;
      Omega[1,2] = rho;
      Omega[2,1] = rho;
      Omega[2,2] = 1;

      for (i in 1:N)
        theta[i] ~ multi_normal(mu, quad_form_diag(Omega, tau));

      for (i in 1:N)
        y[i] ~ multi_normal(theta[i], Sigma[i]);
    }
    "

stan_model_prop_line = stan_model(model_code = stan_code_prop_line)
stan_model_default = stan_model(model_code = stan_code_default)

fit_surrogacy_model <- function(data, assume_proportional_line, iter = 10000, warmup = 5000, chains = 4, seed = 123) {
  # Define Stan code dynamically based on the assumption of an identity regression line
  stan_model <- if (assume_proportional_line) {
    stan_model_prop_line
  } else {
    stan_model_default
  }
  
  # Prepare input data for Stan
  N <- nrow(data)
  y_matrix <- cbind(data$trt_effect_surrogate_index_est, data$log_RR_est)
  
  Sigma_array <- array(NA, dim = c(N, 2, 2))
  for (i in 1:N) {
    Sigma_array[i,,] <- as.matrix(Matrix::nearPD(data$vcov[[i]])$mat)
  }
  
  stan_data <- list(
    N = N,
    y = y_matrix,
    Sigma = Sigma_array
  )
  
  # Fit model
  fit <- sampling(
    object = stan_model,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed
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
      iter = 5e3,
      warmup = 2e3,
      .options = furrr_options(seed = TRUE)
    )
  )

# Saving Results ----------------------------------------------------------

# Save fitted Bayesian models. 
saveRDS(surrogate_results_tbl, file = "results/raw-results/application/bayesian_ma_results.rds")

# The fitted Bayesian models are very large objects. So, we also save a subset
# of this information next, which is more convenient to work with.
rho_long_tbl <- surrogate_results_bayesian_tbl %>%
  mutate(rho_samples = map(stan_fit, ~ as.data.frame(.x)$rho)) %>%
  unnest(rho_samples) %>%
  rename(rho = rho_samples)

saveRDS(rho_long_tbl, file = "results/raw-results/application/rho_long_tbl.rds")



