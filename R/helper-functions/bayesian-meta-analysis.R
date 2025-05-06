library(rstan)
library(bayesplot)

# Specify options for saving the plots to files
figures_dir = "results/figures/application/meta-analysis"
tables_dir = "results/tables/application/meta-analysis"

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
  summarise(stan_fit = list(fit_surrogacy_model(data = pick(everything()))))


