B = 5000
analysishelpers = R/helper-functions/delta-method-rho-trial.R R/helper-functions/moment-based-estimator.R R/helper-functions/multiplier-bootstrap.R R/helper-functions/train-clinical-prediction-models.R
simulationhelpers = R/helper-functions/simulation-functions.R
data = data/processed_data.csv

.PHONY: all
all: results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-small.rds \
	results/raw-results/simple-simulation/ma-sim-results-vaccine-small.rds \
	results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-large.rds
	

results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-small.rds: R/simple-simulation.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simple-simulation.R proof-of-concept small 5
	
results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-large.rds: R/simple-simulation.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simple-simulation.R proof-of-concept large 5
	
results/raw-results/simple-simulation/ma-sim-results-vaccine-small.rds: R/simple-simulation.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simple-simulation.R vaccine small 5
	
# R/application.Rout: R/application.R $(analysishelpers) $(data)
#	 Rscript --verbose R/application.R  > $@ 2> $@