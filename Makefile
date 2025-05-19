analysishelpers = R/helper-functions/delta-method-rho-trial.R R/helper-functions/moment-based-estimator.R R/helper-functions/multiplier-bootstrap.R R/helper-functions/train-clinical-prediction-models.R
simulationhelpers = R/helper-functions/simulation-functions.R
data = data/processed_data.csv

.PHONY: all application simulation
all: simulation \
	application
	
simulation: results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-small.rds \
	results/raw-results/simple-simulation/ma-sim-results-vaccine-small.rds \
	results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-large.rds \
	R/simulations/processing-results.Rout 
	
application: R/application/data-exploration.Rout \
	R/application/ipd_surr_indices_tbl.rds \
	R/application/ma_trt_effects_tbl.rds \
	R/application/meta_analysis.Rout \
	R/results/raw-results/application/bayesian_ma_results.rds \
	R/application/processing-results.Rout
	

results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-small.rds: R/simulations/simulations.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simulations/simulations.R proof-of-concept small 5
	
results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-large.rds: R/simulations/simulations.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simulations/simulations.R proof-of-concept large 5
	
results/raw-results/simple-simulation/ma-sim-results-vaccine-small.rds: R/simulations/simulations.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simulations/simulations.R vaccine small 5
	
R/simulations/processing-results.Rout: results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-small.rds results/raw-results/simple-simulation/ma-sim-results-proof-of-concept-large.rds results/raw-results/simple-simulation/ma-sim-results-vaccine-small.rds
	Rscript --verbose R/simulations/processing-results.R  > $@ 2> $@
	
R/application/data-exploration.Rout: R/application/data-exploration.R $(data)
	Rscript --verbose R/application/data-exploration.R  > $@ 2> $@
	
R/application/ipd_surr_indices_tbl.rds: R/application/surrogate_index_estimation.R $(data)
	Rscript --verbose R/application/surrogate_index_estimation.R  > R/application/surrogate_index_estimation.Rout 2> R/application/surrogate_index_estimation.Rout
	
R/application/ma_trt_effects_tbl.rds: R/application/trial-level-effects.R R/application/ipd_surr_indices_tbl.rds
	Rscript --verbose R/application/trial-level-effects.R  > R/application/trial-level-effects.Rout 2> R/application/trial-level-effects.Rout
	
R/application/meta_analysis.Rout: R/application/meta_analysis.R R/application/ma_trt_effects_tbl.rds $(analysishelpers)
	Rscript --verbose R/application/meta_analysis.R  > $@ 2> $@
	
R/results/raw-results/application/bayesian_ma_results.rds: R/application/bayesian-meta-analysis.R R/application/ma_trt_effects_tbl.rds
	Rscript --verbose R/application/bayesian-meta-analysis.R  > R/application/bayesian-meta-analysis.Rout 2> R/application/bayesian-meta-analysis.Rout
	
R/application/processing-results.Rout: R/application/processing-results.R R/application/bayesian_ma_results.rds R/application/meta_analysis.Rout
	Rscript --verbose R/application/processing-results.R  > $@ 2> $@