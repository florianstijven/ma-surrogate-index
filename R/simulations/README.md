# Simulations

This folder contains the R scripts used for conducting, illustrating, and 
analyzing the simulations. The following R scripts are present:

### `simulations.R`

This script is used for actually conducting the simulation studies. It takes the
following command-line arguments in the presented order:

1. The simulation scenario, "proof-of-concept" or "vaccine".
2. The sample size regime, "small" or "large".
3. The number of simulation replicates to run.

This script saves the raw results to `results/raw-results/simulations/` for 
further processing (e.g., tables and plots).

### `processing-results.R`

This script processes the results from completed simulation runs. It reads raw 
results and summarizes them in more digestible tables and plots which are saved 
to `results/tables/simulations/` and `results/figures/simulations/` for reporting
and further analysis.

### `illustration-simulations.R`

This script generates data and uses these data to make plots that illustrate 
the different data-generating mechanisms used in the simulations. These plots 
are saved to `results/figures/simulations/`.


