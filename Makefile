analysishelpers = R/helper-functions/delta-method-rho-trial.R R/helper-functions/moment-based-estimator.R R/helper-functions/multiplier-bootstrap.R R/helper-functions/train-clinical-prediction-models.R
simulationhelpers = R/helper-functions/simulation-functions.R
data = data/processed_data.csv

.PHONY: all
all: results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-small.rds \
	results/raw-results/simple-simulation/ma-sim-results-vaccine-small.rds \
	results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-large.rds \
	R/application/data-exploration.Rout \
	R/application/meta_analysis.Rout
	

results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-small.rds: R/simulations/simulations.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simple-simulation.R proof-of-concept small 5
	
results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-large.rds: R/simulations/simulations.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simple-simulation.R proof-of-concept large 5
	
results/raw-results/simple-simulation/ma-sim-results-vaccine-small.rds: R/simulations/simulations.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simple-simulation.R vaccine small 5
	
R/application/data-exploration.Rout: R/application/data-exploration.R $(data)
	Rscript --verbose R/application/data-exploration.R  > $@ 2> $@
	
R/application/ipd_surr_indices_tbl.rds: R/application/surrogate_index_estimation.R $(data)
	Rscript --verbose R/application/surrogate_index_estimation.R  > $@ 2> $@
	
R/application/ma_trt_effects_tbl.rds: R/application/trial-level-effects.R R/application/ipd_surr_indices_tbl.rds
	Rscript --verbose R/application/trial-level-effects.R  > $@ 2> $@
	
R/application/meta_analysis.Rout: R/application/ma_trt_effects_tbl.rds $(analysishelpers)
	Rscript --verbose R/application/meta_analysis.R  > $@ 2> $@