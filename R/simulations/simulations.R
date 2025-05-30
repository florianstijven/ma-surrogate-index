# Preliminaries -----------------------------------------------------------

a = Sys.time()
# Load all R packages
library(tidyverse)
library(future)
library(furrr)
library(mgcv)
library(splines)
library(sl3)
library(origami)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}
# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

# Under which scenario should the simulations be conducted? ("vaccine" or
# "proof-of-concept")
scenario = args[1]
# Are the simulations small or large sample? Large sample is just to illustrate
# problems under some asymptotic regimes. ("small" or "large")
regime = args[2]
# Number of Monte Carlo simulations.
N_MC = as.numeric(args[3])

# This script simulates data. We set a seed to ensure reproducibility.
set.seed(123)

# Set different kinds of sample-size parameters depending on the scenario.
if (scenario == "proof-of-concept") {
  # Within-trial sample size
  n = 2e3
  # Number of independent trials
  N = c(6, 12, 24)
  
  # Number of Monte Carlo trial replications for the approximation of the true
  # trial-level correlation.
  N_approximation_MC = 2e4
  # Number of patients in each trial for the approximation of the true
  # trial-level correlation. This can be small for the proof-of-concept scenario
  # because that is based on differences in sample means for which the standard
  # covariance estimator is unbiased.
  n_approximation_MC = 2e2
} else if (scenario == "vaccine") {
  # Within-trial sample size
  n = 5e3
  # Number of independent trials
  N = c(6, 12, 24)
  
  
  # Number of Monte Carlo trial replications for the approximation of the true
  # trial-level correlation.
  N_approximation_MC = 5e3
  # Number of patients in each trial for the approximation of the true
  # trial-level correlation. This should be large in the vaccine scenario
  # because the standard covariance estimator (based on the delta method) for
  # log RR estimators is consistent but biased.
  n_approximation_MC = 3e3
}


if (regime == "small") {
  # sd_beta represents the variance of important trial-level parameters that
  # control the surrogacy and comparability assumptions. A non-zero variance
  # implies some violation of these assumptions.
  if (scenario == "vaccine") {
    # Which SDs correspond to slight and moderate violations of surrogayc and
    # comparability depend on the scenario.
    sd_beta = list(c(0.15, 0.15), c(0.30, 0.30))
  } else if (scenario == "proof-of-concept") {
    sd_beta = list(c(0.05, 0.05), c(0.1, 0.1))
  }
  # The elements of sd_beta correspond to "slight" and "moderate" violations.
  SI_violation = c("slight", "moderate")
  
  # Define the number of bootstrap replications depending on the number of
  # independent trials.
  B_N_lookup = tibble(N = c(6, 12, 24), B = c(1e5, 5e4, 2.5e4))
  
  # The types of bootstrap CIs studied in the simulations. 
  bootstrap_types = c("BCa", "BC percentile")
  
  # Indicator of whether Bayesian inferences are included in the simulations.
  bayesian = TRUE
} else if (regime == "large") {
  if (scenario == "vaccine") {
    stop("The large sample regime is not suitable for the vaccine scenario.")
  }
  # This regime corresponds to the large N small n asymptotic regime. This
  # scenario is artificial and just meant to show when methods can break.
  N = c(100, 500, 1e3, 2e3, 5e3)  # Number of trials in each meta-analytic data set
  # The within-trial sample size is set to something small.
  n = c(100, 500, 1e3, 2e3, 5e3)
  
  # The approximation accuracy for the true rho is increased for the large N
  # setting. In this setting, the SE of the estimators will be much smaller; so,
  # the MC error in the true rho is relatively more important.
  N_approximation_MC = 2e5
  n_approximation_MC = 5e2
  
  # Set the number of bootstrap replications to something lower.
  # Define the number of bootstrap replications depending on the number of
  # independent trials.
  B_N_lookup = tibble(N = N, B = 1e3)
  bootstrap_types = "BCa"
  
  bayesian = FALSE
  
  sd_beta = list(c(0.10, 0.10))
  SI_violation = c("moderate")
}




# Helper Functions --------------------------------------------------------

# Source helper functions.
source("R/helper-functions/simulation-functions.R")
source("R/helper-functions/moment-based-estimator.R")
source("R/helper-functions/delta-method-rho-trial.R")
source("R/helper-functions/multiplier-bootstrap.R")
source("R/helper-functions/train-clinical-prediction-models.R")
source("R/helper-functions/bayesian-model.R")

# Set the surrogate index estimation methods.
if (scenario == "proof-of-concept") {
  # The clinical endpoint for the proof-of-concept scenario is continuous and
  # the underlying DGM is quite simple; so, we use linear regression here.
  surrogate_index_estimator = c("surrogate", "lm")
  # We only look at the original surrogate for the large regime. This saves some
  # computational time because the true correlation parameter is
  # non-data-adaptive and thus only has to be computed once.
  if (regime == "large") {
    surrogate_index_estimator = c("surrogate")
  }
} else if (scenario == "vaccine") {
  # The clinical endpoint for the vaccine scenario is binary. So, we use
  # logistic regression or a GAM logistic regression model to estimate the
  # surrogate index here.
  surrogate_index_estimator = c("surrogate", "superlearner")
}

# Function factory that returns "statistic" functions which take the MA data and
# corresponding weights as input, and return the estimated trial-level
# correlation and SE. This function factory allows for different versions of the
# correlation etimator: (i) finite-sample estimator adjustment and (ii)
# truncating the estima
statistic_function_factory = function(estimator_adjustment, nearest_PD) {
  statistic_f = function(data, weights) {
    moment_estimate = moment_estimator(
      alpha_hat = data$treatment_effect_surr,
      beta_hat = data$treatment_effect_clin,
      vcov_list = data$covariance_matrix,
      estimator_adjustment = estimator_adjustment,
      weights = weights,
      nearest_PD = nearest_PD
    )
    
    coefs_est = moment_estimate$coefs
    # We're not computing the SE because we're not using the studentized
    # interval in these simulations. If we were to use the studentized interval,
    # we would have to modify this part of the code.
    return(list(estimate = coefs_est[5] / sqrt(coefs_est[4] * coefs_est[3]), se = NA))
  }
  return(statistic_f)
}

# Helper Compute CI based on Bivariate Bayesian Model.
compute_bayesian_ci = function(data) {
  # Temporarily disable parallelization
  old_plan <- future::plan("sequential")
  on.exit(future::plan(old_plan), add = TRUE) # Restore original plan after execution
  
  data = data %>%
    rename(
      trt_effect_surrogate_index_est = treatment_effect_surr,
      log_RR_est = treatment_effect_clin,
      vcov = covariance_matrix
    )
  
  fit = fit_surrogacy_model(
    data,
    assume_proportional_line = FALSE,
    iter = 1e4,
    warmup = 5e3,
    chains = 1,
    seed = 1
  )
  # Extract 95% credible interval
  summary_fit <- rstan::summary(fit, probs = c(0.025, 0.975))$summary
  rm("fit")
  gc()
  
  # Extract the row corresponding to 'rho'
  rho_summary <- summary_fit["rho", ]
  
  # View the 95% credible interval
  rho_ci <- rho_summary[c("2.5%", "97.5%")]
  
  return(rho_ci)
}


## Simulation Function --------------------------------------------------

# Function to analyze one data set.
analyze = function(alpha_hat,
                   beta_hat,
                   vcov_list,
                   estimator_adjustment,
                   sandwich_adjustment,
                   nearest_PD,
                   bootstraps,
                   B,
                   bayesian,
                   N_approximation_MC,
                   n_approximation_MC) {
  # Construct tibble with a cross-tabulation of all analysis options.
  analysis_options = expand_grid(estimator_adjustment, sandwich_adjustment, nearest_PD)
  
  N = length(alpha_hat)
  
  # Estimate MA model parameters.
  inferences_tbl = analysis_options %>%
    rowwise(everything()) %>%
    summarise(
      # Estimate the meta-analytic parameters using the moment-based estimator.
      moment_estimates = list(
        moment_estimator(
          alpha_hat,
          beta_hat,
          vcov_list,
          estimator_adjustment = estimator_adjustment,
          sandwich_adjustment = sandwich_adjustment,
          nearest_PD = nearest_PD
        )
      ),
      # Estimate the trial-level Pearson correlation using the delta method.
      rho_delta_method = list(
        rho_delta_method(
          moment_estimates$coefs,
          moment_estimates$vcov,
          method = "t-adjustment",
          N = N
        )
      ),
      rho_est = rho_delta_method$rho,
      rho_se = rho_delta_method$se,
      rho_ci_lower = rho_delta_method$ci[1],
      rho_ci_upper = rho_delta_method$ci[2]
    ) %>%
    ungroup() %>%
    select(-moment_estimates, -rho_delta_method) %>%
    mutate(CI_type = "sandwich")
  
  # Compute boostrap confidence intervals if required.
  if (length(bootstraps) >= 1) {
    inferences_bootstrap_tbl = analysis_options %>%
      right_join(inferences_tbl %>% select(-CI_type, -rho_ci_upper, -rho_ci_lower)) %>%
      cross_join(tibble(CI_type = bootstraps)) %>%
      rowwise(everything()) %>%
      summarise(
        bootstrap_ci =
          multiplier_bootstrap_ci(
            data = tibble(
              treatment_effect_surr = alpha_hat,
              treatment_effect_clin = beta_hat,
              covariance_matrix = vcov_list
            ),
            statistic = statistic_function_factory(estimator_adjustment, nearest_PD),
            B = B,
            alpha = 0.05,
            type = CI_type
          ),
        rho_ci_lower = bootstrap_ci$ci_lower,
        rho_ci_upper = bootstrap_ci$ci_upper,
      ) %>%
      ungroup() %>%
      select(-bootstrap_ci)
    
    
    # Join to inferences_tbl
    inferences_tbl = inferences_tbl %>%
      bind_rows(inferences_bootstrap_tbl)
  }
  
  
  # Compute Bayesian CIs if required
  if (bayesian) {
    credible_interval = compute_bayesian_ci(
      tibble(
        treatment_effect_surr = alpha_hat,
        treatment_effect_clin = beta_hat,
        covariance_matrix = vcov_list
      )
    )
    inferences_tbl = inferences_tbl %>%
      bind_rows(
        tibble(
          rho_ci_lower = credible_interval[1],
          rho_ci_upper = credible_interval[2],
          CI_type = "Bayesian"
        )
      )
  }
  return(inferences_tbl)
}

# Function to simulate and analyze one data set.
simulate_and_analyze = function(sd_beta,
                                n,
                                N,
                                surrogate_index_estimators,
                                scenario,
                                regime,
                                bootstraps,
                                B,
                                bayesian,
                                N_approximation_MC,
                                n_approximation_MC) {
  # Temporarily disable parallelization
  old_plan <- future::plan("sequential")
  on.exit(future::plan(old_plan), add = TRUE) # Restore original plan after execution
  
  # Generate the actual MA data
  meta_analytic_data = generate_meta_analytic_data(
    sd_beta = sd_beta,
    n = n,
    N = N,
    surrogate_index_estimators = surrogate_index_estimators,
    scenario = scenario,
    regime = regime
  )
  
  
  
  # Compute/approximate the true trial-level correlation with the given
  # estimated surrogate index.
  rho_approximated = rep(NA, length(surrogate_index_estimators))
  inferences_tbl = tibble(
    estimator_adjustment = character(),
    sandwich_adjustment = character(),
    nearest_PD = logical(),
    rho_est = numeric(),
    rho_se = numeric(),
    rho_ci_lower = numeric(),
    rho_ci_upper = numeric(),
    CI_type = character(),
    surrogate_index_estimator = character()
  )
  for (i in seq_along(1:length(surrogate_index_estimators))) {
    # Estimate the meta-analytic parameters.
    inferences_tbl = inferences_tbl %>%
      bind_rows(
        analyze(
          alpha_hat = meta_analytic_data[[1]][[i]]$treatment_effect_surr,
          beta_hat = meta_analytic_data[[1]][[i]]$treatment_effect_clin,
          vcov_list = meta_analytic_data[[1]][[i]]$covariance_matrix,
          estimator_adjustment = "N - 1",
          sandwich_adjustment = "N - 1",
          nearest_PD = c(TRUE, FALSE),
          bootstraps = bootstraps,
          B = B,
          bayesian = bayesian,
          n_approximation_MC = n_approximation_MC,
          N_approximation_MC = N_approximation_MC
        ) %>%
          mutate(surrogate_index_estimator = surrogate_index_estimators[i])
      )
    
    if (surrogate_index_estimators[i] == "surrogate") {
      # We will not compute the true rho for the surrogate because this is the
      # same true rho across simulations. This value is computed later on just
      # once per setting, so avoiding numerically computing the same thing many
      # times.
      rho_approximated[i] = NA
    } else {
      # Compute the true rho given the surrogate index f.
      rho_approximated[i] = rho_MC_approximation(
        sd_beta = sd_beta,
        # We are assuming here that there is only one non-trivial
        # surrogate index estimator.
        f = meta_analytic_data[[2]][[i]],
        N_approximation_MC = N_approximation_MC,
        n_approximation_MC = n_approximation_MC,
        scenario = scenario,
        regime = regime
      )
    }
  }
  # Add approximated rhos to tibble with inferences.
  inferences_tbl = inferences_tbl %>%
    left_join(
      tibble(rho_true = rho_approximated, surrogate_index_estimator = surrogate_index_estimators)
    )
  
  return(inferences_tbl)
}


# Simulation Study ---------------------------------------------------------

# Tibble with simulation parameters. Each row corresponds to a data-generating
# mechanism.
dgm_param_tbl = expand_grid(tibble(sd_beta, SI_violation), N, n) %>%
  mutate(surrogate_index_estimators = list(surrogate_index_estimator)) %>%
  left_join(B_N_lookup)

## Simulate IPD data and compute trial-level estimates  --------------------

# Generate N meta-analytic data sets and compute the trial-level Pearson
# correlation for (i) the treatment effects on the surrogate and clinical
# endpoints and (ii) the treatment effects on the surrogate index and clinical
# endpoint.
data_set_indicator = 1:N_MC

simulations_results_tbl <- expand_grid(data_set_indicator, dgm_param_tbl)

simulations_results_tbl$inferences_tbl = future_pmap(
  .l = list(
    sd_beta = simulations_results_tbl$sd_beta,
    n = simulations_results_tbl$n,
    N = simulations_results_tbl$N,
    surrogate_index_estimators = simulations_results_tbl$surrogate_index_estimators,
    B = simulations_results_tbl$B
  ),
  .f = simulate_and_analyze,
  bayesian = bayesian,
  N_approximation_MC = N_approximation_MC,
  n_approximation_MC = n_approximation_MC,
  bootstraps = bootstrap_types,
  scenario = scenario,
  regime = regime,
  .options = furrr_options(
    seed = TRUE,
    stdout = FALSE,
    conditions = character()
  )
)

simulations_results_tbl = simulations_results_tbl %>%
  rowwise(!(c(inferences_tbl, surrogate_index_estimators, sd_beta))) %>%
  reframe(inferences_tbl)

print(Sys.time() - a)

## Approximate the true trial-level correlation parameter ------------------

# Approximate the true trial-level correlation using the surrogate endpoint.
rho_true_surrogate_tbl = tibble(sd_beta, SI_violation) %>%
  mutate(surrogate_index_estimator = "surrogate") 

rho_true_surrogate_tbl$rho_true = future_map_dbl(
  .x = rho_true_surrogate_tbl$sd_beta,
  .f = rho_MC_approximation,
  N_approximation_MC = N_approximation_MC,
  n_approximation_MC = n_approximation_MC,
  f = NULL,
  scenario = scenario,
  regime = regime,
  .options = furrr_options(seed = TRUE)
)

print(Sys.time() - a)

if (regime == "large") {
  plan(sequential)
}


# Add the true trial-level rho parameters to the tibble with the simulated data
# sets.
simulations_results_tbl = simulations_results_tbl %>%
  left_join(
    rho_true_surrogate_tbl %>%
      select(-sd_beta),
    by = c("SI_violation", "surrogate_index_estimator"),
    relationship = "many-to-one"
  ) %>%
  # rho_true is present in both tibbles before joining. Hence, we get rho_true.x
  # and .y in the joined data set. We keep the non-missing value.
  mutate(rho_true = coalesce(rho_true.x, rho_true.y)) %>%
  select(-rho_true.x, -rho_true.y)


print(Sys.time() - a)

# Save Results -------------------------------------------------------------

# Save the simulated MA data sets with the corresponding rho estimates and
# confidence intervals.
saveRDS(
  simulations_results_tbl,
  paste0(
    "results/raw-results/simulations/ma-sim-results-",
    scenario,
    "-",
    regime,
    ".rds"
  )
)

print(Sys.time() - a)


