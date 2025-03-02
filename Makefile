.PHONY: all
all: results/raw-results/simple-simulation/meta_analytic_data_simulated.rds

results/raw-results/simple-simulation/meta_analytic_data_simulated.rds: R/simple-simulation.R
	Rscript R/simple-simulation.R