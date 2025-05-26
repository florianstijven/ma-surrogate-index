# Preliminaries -----------------------------------------------------------

a = Sys.time()
# Load all R packages
library(tidyverse)
library(sandwich)
library(lmtest)
library(future)
library(furrr)
library(mgcv)
library(splines)
library(sl3)
library(origami)
library(rstan)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
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
    sd_beta = list(c(0.15, 0.15), c(0.30, 0.30))
  } else if (scenario == "proof-of-concept") {
    sd_beta = list(c(0.05, 0.05), c(0.1, 0.1))
  }
  SI_violation = c("slight", "moderate")
  
  # Define the number of bootstrap replications depending on the number of
  # independent trials.
  B_N_lookup = tibble(N = c(6, 12, 24), B = c(1e5, 5e4, 2.5e4))
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
  B_N_lookup = tibble(N = N, B = 2e3)
  
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

meta_analytic_data <- expand_grid(data_set_indicator, dgm_param_tbl) %>%
  mutate(
    list_of_ma_data_objects = future_pmap(
      .l = list(
        sd_beta = sd_beta,
        n = n,
        N = N,
        surrogate_index_estimators = surrogate_index_estimators
      ),
      .f = function(sd_beta, n, N, surrogate_index_estimators, scenario, regime, N_approximation_MC, n_approximation_MC) {
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
        rho_approximated = rep(NA, length(surrogate_index_estimator))
        for (i in seq_along(1:length(surrogate_index_estimator))) {
          if (surrogate_index_estimator[i] == "surrogate") {
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
              f = meta_analytic_data[[2]][[2]],
              N_approximation_MC = N_approximation_MC,
              n_approximation_MC = n_approximation_MC,
              scenario = scenario,
              regime = regime
            )
          }
        }

        
        return(
          list(
            meta_analytic_data = meta_analytic_data,
            rho_approximated = rho_approximated
          )
        )
      },
      N_approximation_MC = N_approximation_MC,
      n_approximation_MC = n_approximation_MC,
      scenario = scenario,
      regime = regime,
      .options = furrr_options(seed = TRUE)
    )
  ) %>%
  rowwise(everything()) %>%
  reframe(
    tibble(
      treatment_effects = list_of_ma_data_objects$meta_analytic_data[[1]],
      rho_true = list_of_ma_data_objects$rho_approximated
    )
  ) %>%
  ungroup() %>%
  select(-list_of_ma_data_objects, -surrogate_index_estimators)

print(Sys.time() - a)

## Approximate the true trial-level correlation parameter ------------------

print(Sys.time() - a)

# Approximate the true trial-level correlation using the surrogate endpoint.
rho_true_surrogate_tbl = tibble(sd_beta, SI_violation) %>%
  mutate(surrogate_index_estimator = "surrogate") %>%
  mutate(
    rho_true = future_map_dbl(
      .x = sd_beta,
      .f = rho_MC_approximation,
      N_approximation_MC = N_approximation_MC,
      n_approximation_MC = n_approximation_MC,
      f = NULL,
      scenario = scenario,
      regime = regime,
      .options = furrr_options(seed = TRUE)
    )
  )

print(Sys.time() - a)


# Add the true trial-level rho parameters to the tibble with the simulated data
# sets.
meta_analytic_data = meta_analytic_data %>%
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

## Moment-based estimator --------------------------------------------------

# Estimate trial-level Pearson correlation for the generated meta-analytic data
# sets of treatment effects and add everything to a tibble.
meta_analytic_data_simulated = meta_analytic_data %>%
  cross_join(expand_grid(
    estimator_adjustment = c("N - 1"),
    sandwich_adjustment = c("N - 1"),
    nearest_PD = c(TRUE, FALSE)
  )) %>%
  rowwise(everything()) %>%
  summarise(
    # Estimate the meta-analytic parameters using the moment-based estimator.
    moment_estimates = list(
      moment_estimator(
        treatment_effects$treatment_effect_surr,
        treatment_effects$treatment_effect_clin,
        treatment_effects$covariance_matrix,
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
  # Drop superfluous variables
  select(-moment_estimates, -rho_delta_method)

print(Sys.time() - a)

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
    rho = rho_delta_method(
      coefs = moment_estimate$coefs,
      vcov = moment_estimate$vcov,
      method = "t-adjustment",
      # N is only used for the t-adjustment, it doesn't matter for the estimate
      # or SE.
      N = 5
    )
    
    return(list(estimate = rho$rho, se = rho$se))
  }
  return(statistic_f)
}

# Helper: Compute bootstrap CIs of specified type for a batch of rows
compute_bootstrap_cis <- function(data, type) {
  furrr::future_pmap(
    .l = list(
      treatment_effects = data$treatment_effects,
      estimator_adjustment = data$estimator_adjustment,
      nearest_PD = data$nearest_PD,
      B = data$B
    ),
    .f = function(treatment_effects, estimator_adjustment, nearest_PD, B) {
      multiplier_bootstrap_ci(
        data = treatment_effects,
        statistic = statistic_function_factory(estimator_adjustment, nearest_PD),
        B = B,
        alpha = 0.05,
        type = type
      )
    },
    .options = furrr_options(seed = TRUE)
  )
}

# Source helper function to fit Bayesian model. 
source("R/helper-functions/bayesian-model.R")

# Helper Compute CI based on Bivariate Bayesian Model.
compute_bayesian_ci = function(data) {
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
    chains = 3,
    seed = 1
  )
  # Extract 95% credible interval
  summary_fit <- summary(fit, probs = c(0.025, 0.975))$summary
  
  # Extract the row corresponding to 'rho'
  rho_summary <- summary_fit["rho", ]
  
  # View the 95% credible interval
  rho_ci <- rho_summary[c("2.5%", "97.5%")]
  
  return(rho_ci)
}

compute_bayesian_cis = function(data) {
  furrr::future_map(
    .x = data$treatment_effects,
    .f = compute_bayesian_ci,
    .options = furrr_options(seed = TRUE)
  )
}

if (regime == "small") {
  # 1. "Sandwich" CIs (no bootstrap, just copy as-is)
  sandwich_tbl <- meta_analytic_data_simulated %>%
    mutate(CI_type = "sandwich")
  
  # 2. For nearest_PD == FALSE, compute BCa, studentized, and BC percentile CIs
  pd_false <- meta_analytic_data_simulated %>% filter(nearest_PD == FALSE)
  
  # Compute all bootstrap CIs OUTSIDE mutate!
  bca_cis <- compute_bootstrap_cis(pd_false, "BCa")
  studentized_cis <- compute_bootstrap_cis(pd_false, "studentized")
  bcperc_cis <- compute_bootstrap_cis(pd_false, "BC percentile")
  
  bayesian_cis = compute_bayesian_cis(pd_false)
  
  # Add CIs as columns (using map_dbl to extract bounds)
  pd_false_bca <- pd_false %>%
    mutate(
      rho_ci_lower = purrr::map_dbl(bca_cis, 1),
      rho_ci_upper = purrr::map_dbl(bca_cis, 2),
      CI_type = "BCa"
    )
  
  pd_false_studentized <- pd_false %>%
    mutate(
      rho_ci_lower = purrr::map_dbl(studentized_cis, 1),
      rho_ci_upper = purrr::map_dbl(studentized_cis, 2),
      CI_type = "studentized"
    )
  
  pd_false_bcperc <- pd_false %>%
    mutate(
      rho_ci_lower = purrr::map_dbl(bcperc_cis, 1),
      rho_ci_upper = purrr::map_dbl(bcperc_cis, 2),
      CI_type = "BC percentile"
    )
  
  pd_false_bayesian = pd_false %>%
    mutate(
      rho_ci_lower = purrr::map_dbl(bayesian_cis, 1),
      rho_ci_upper = purrr::map_dbl(bayesian_cis, 2),
      CI_type = "Bayesian"
    )
  
  # Combine all results
  meta_analytic_data_simulated <- bind_rows(
    sandwich_tbl,
    pd_false_bca,
    pd_false_studentized,
    pd_false_bcperc,
    pd_false_bayesian
  )
  
} else {
  # For the large N regime, only BCa CIs for nearest_PD == FALSE
  sandwich_tbl <- meta_analytic_data_simulated %>%
    mutate(CI_type = "sandwich")
  
  pd_false <- meta_analytic_data_simulated %>% filter(nearest_PD == FALSE)
  bca_cis <- compute_bootstrap_cis(pd_false, "BCa")
  
  # bayesian_cis = compute_bayesian_cis(pd_false)
  
  pd_false_bca <- pd_false %>%
    mutate(
      rho_ci_lower = purrr::map_dbl(bca_cis, 1),
      rho_ci_upper = purrr::map_dbl(bca_cis, 2),
      CI_type = "BCa"
    )
  
  # pd_false_bayesian = pd_false %>%
  #   mutate(
  #     rho_ci_lower = purrr::map_dbl(bayesian_cis, 1),
  #     rho_ci_upper = purrr::map_dbl(bayesian_cis, 2),
  #     CI_type = "Bayesian"
  #   )
  
  meta_analytic_data_simulated <- bind_rows(
    sandwich_tbl,
    pd_false_bca
    # pd_false_bayesian
  )
}




# Drop columns that are not needed.
meta_analytic_data_simulated = meta_analytic_data_simulated %>%
  select(-treatment_effects)

print(Sys.time() - a)

# Save Results -------------------------------------------------------------

# Save the simulated MA data sets with the corresponding rho estimates and
# confidence intervals.
saveRDS(
  meta_analytic_data_simulated,
  paste0(
    "results/raw-results/simulations/ma-sim-results-",
    scenario,
    "-",
    regime,
    ".rds"
  )
)

print(Sys.time() - a)


