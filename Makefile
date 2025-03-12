B = 5000
analysishelpers = R/helper-functions/delta-method-rho-trial.R R/helper-functions/moment-based-estimator.R R/helper-functions/multiplier-bootstrap.R
simulationhelpers = R/helper-functions/simulation-functions.R
data = data/processed_data.csv

.PHONY: all
all: results/raw-results/simple-simulation/ma-sim-results-proof-of-concept.rds \
results/raw-results/simple-simulation/ma-sim-results-vaccine.rds
	

results/raw-results/simple-simulation/ma-sim-results-proof-of-concept.rds: R/simple-simulation.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simple-simulation.R proof-of-concept
	
results/raw-results/simple-simulation/ma-sim-results-vaccine.rds: R/simple-simulation.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simple-simulation.R vaccine
	
# R/application.Rout: R/application.R $(analysishelpers) $(data)
#	 Rscript --verbose R/application.R  > $@ 2> $@