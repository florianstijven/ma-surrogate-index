# Application Functions

This folder contains R scripts with functions specifically supporting the *application* 
of the meta-analytic surrogate index methods to the COVID-19 vaccine trials.

All R scripts in this folder take the following command-line argument:

1. Should the analysis be run on the `"real"` or `"synthetic"` data`?
   - If `"real"`, it uses the actual COVID-19 vaccine trial data. This requires 
   the original data to be present in `data/processed_data.csv`. Note that these
   data are not included in the repository. 
   - If `"synthetic"`, it uses the synthetic data generated in `data-generation.R`.
   These data are included in the repository in `data/processed_data_synthetic.csv`.

### `surrogate_index_estimation.R`

This script estimates the surrogate index using various methods.

### `trial-level-effects.R`

This script estimates the trial-level treatment effects and corresponding variance
matrices (based on the bootstrap). These results are saved to `results/raw-results/application/`
for further use in `meta_analysis.R`. 

### `meta_analysis.R`

This script performs the non-parametric frequentist meta-analysis using the
trial-level treatment effect estimates (and variance matrices) computed in
`trial-level-effects.R`.

### `bayesian-meta-analysis.R`

This script performs the Bayesian meta-analysis using the trial-level treatment
effect estimates (and variance matrices) computed in `trial-level-effects.R`.

### `processing-results.R`

This script processes the results from the meta-analyses performed in
`meta_analysis.R` and `bayesian-meta-analysis.R`. The tables and figures are
saved to `results/figures/application/meta-analysis/` and
`results/tables/application/meta-analysis/` for the original data, and to
`results/figures/application-synthetic/meta-analysis/` and
`results/tables/application-synthetic/meta-analysis/` for the synthetic data.


### `data-exploration.R`

This script performs exploratory data analysis on the COVID-19 vaccine trial data.

