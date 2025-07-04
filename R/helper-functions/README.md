# Helper Functions

This folder contains R scripts with helper functions that support the simulation
studies and analyses in this project. These files are sourced by scripts in
`R/simulations/` and `R/application/` such that the functions become available.
Each script in this folder is described below:

### `bayesian-model.R`

Provides a function (`fit_surrogacy_model()`) that fits the Bayesian
meta-analytic surrogate model. This model is fitted with the `rstan` package.
This function allows one to fit the standard meta-analytic surrogate model as
well as the corresponding model that assumes a proportional regression function
for regressing the surrogate (index) treatment effects on the clinical endpoint
treatment effects.

When running this script, the stan model is compiled (which takes some time) and 
retained in the R environment, so it does not have to be recompiled each time 
`fit_surrogacy_model()` is called. 


### `delta-method-rho-trial.R`

Contains the `rho_delta_method()` function that computes the confidence interval
for the trial-level Pearson correlation parameter based on the estimated means
and covariance parameters and the corresponding sandwich variance estimate. This
confidence interval is based on the delta method, but some finite-sample
adjustments are available.

### `moment-based-estimator.R`

Contains functions to estimate the mean and covariance parameters of the
trial-level treatment effects adjusted for the measurement error/sampling
variability in the trial-level treatment effect estimates. The sandwich estimate
for the variance of these estimators is also computed. The main function is
`moment_based_estimator()`, which takes the trial-level treatment effect
estimates and their estimated variance matrices as input.

### `multiplier-bootstrap.R`

Implements the Bayesian bootstrap. Several bootstrap CI types are implemented,
including BCa, percentile, and studentized intervals. Only the BCa intervals are 
used in the simulations and application, however. 

### `simulation-functions.R`

Contains a set of helper functions that are used to simulate individual-level
data from multiple trials. This script also contains a function to numerically 
approximate the trial-level correlation for a given data-generating mechanism 
and (estimated) surrogate index. 

### `train-clinical-prediction-models.R`

Provides functions to train various clinical endpoint prediction models on the
simulated data. Supported model types includes parametric regression models,
random forests, SuperLearner, and the highly adaptive lasso (HAL). Only a subset
of these functions are actually used in the simulations.