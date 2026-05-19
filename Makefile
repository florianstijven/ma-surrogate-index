analysishelpers = R/helper-functions/delta-method-rho-trial.R R/helper-functions/moment-based-estimator.R R/helper-functions/multiplier-bootstrap.R R/helper-functions/train-clinical-prediction-models.R R/helper-functions/bayesian-model.R
simulationhelpers = R/helper-functions/simulation-functions.R
data = data/processed_data.csv
data-synthetic = data/processed_data_synthetic.csv

.PHONY: all application simulation additional-application
all: simulation \
	application \
	application-synthetic \
	additional-application
	
simulation: results/raw-results/simulations/ma-sim-results-proof-of-concept-small.rds \
	results/raw-results/simulations/ma-sim-results-vaccine-small.rds \
	results/raw-results/simulations/ma-sim-results-proof-of-concept-large.rds \
	R/simulations/processing-results.Rout 
	
application: R/application/data-exploration.Rout \
	results/raw-results/application/ipd_surr_indices_tbl.rds \
	results/raw-results/application/ma_trt_effects_tbl.rds \
	R/application/meta_analysis.Rout \
	results/raw-results/application/bayesian_ma_results.rds \
	R/application/processing-results.Rout

application-synthetic: R/application/data-exploration-synthetic.Rout \
	results/raw-results/application-synthetic/ipd_surr_indices_tbl.rds \
	results/raw-results/application-synthetic/ma_trt_effects_tbl.rds \
	R/application/meta_analysis-synthetic.Rout \
	results/raw-results/application-synthetic/bayesian_ma_results.rds \
	R/application/processing-results-synthetic.Rout
	
additional-application: additional-application/R/meta-analytic/estimate-landmark-models.Rout \
	R/meta-analytic/intermediate-objects/landmark_predictions_tbl.rds \
	additional-application/results/raw-results/ma_trt_effects_tbl.rds \
	additional-application/results/tables/surrogacy-inferences.csv \
	additional-application/R/meta-analytic/processing-results.Rout

	

results/raw-results/simulations/ma-sim-results-proof-of-concept-small.rds: R/simulations/simulations.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simulations/simulations.R proof-of-concept small 500
	
results/raw-results/simulations/ma-sim-results-proof-of-concept-large.rds: R/simulations/simulations.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simulations/simulations.R proof-of-concept large 100
	
results/raw-results/simulations/ma-sim-results-vaccine-small.rds: R/simulations/simulations.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simulations/simulations.R vaccine small 500
	
R/simulations/processing-results.Rout: results/raw-results/simulations/ma-sim-results-proof-of-concept-small.rds results/raw-results/simulations/ma-sim-results-proof-of-concept-large.rds results/raw-results/simulations/ma-sim-results-vaccine-small.rds
	Rscript --verbose R/simulations/processing-results.R  > $@ 2> $@
	
	

	
R/application/data-exploration.Rout: R/application/data-exploration.R $(data)
	Rscript --verbose R/application/data-exploration.R real > $@ 2> $@
	
results/raw-results/application/ipd_surr_indices_tbl.rds: R/application/surrogate_index_estimation.R $(data)
	Rscript --verbose R/application/surrogate_index_estimation.R real > R/application/surrogate_index_estimation.Rout 2> R/application/surrogate_index_estimation.Rout
	
results/raw-results/application/ma_trt_effects_tbl.rds: R/application/trial-level-effects.R results/raw-results/application/ipd_surr_indices_tbl.rds
	Rscript --verbose R/application/trial-level-effects.R real > R/application/trial-level-effects.Rout 2> R/application/trial-level-effects.Rout
	
R/application/meta_analysis.Rout: R/application/meta_analysis.R results/raw-results/application/ma_trt_effects_tbl.rds $(analysishelpers)
	Rscript --verbose R/application/meta_analysis.R real > $@ 2> $@
	
results/raw-results/application/bayesian_ma_results.rds: R/application/bayesian-meta-analysis.R results/raw-results/application/ma_trt_effects_tbl.rds
	Rscript --verbose R/application/bayesian-meta-analysis.R  real > R/application/bayesian-meta-analysis.Rout 2> R/application/bayesian-meta-analysis.Rout
	
R/application/processing-results.Rout: R/application/processing-results.R results/raw-results/application/bayesian_ma_results.rds R/application/meta_analysis.Rout
	Rscript --verbose R/application/processing-results.R real > $@ 2> $@
	
	
	
R/application/data-exploration-synthetic.Rout: R/application/data-exploration.R $(data-synthetic)
	Rscript --verbose R/application/data-exploration.R synthetic > R/application/data-exploration-synthetic.Rout 2> R/application/data-exploration-synthetic.Rout
	
results/raw-results/application-synthetic/ipd_surr_indices_tbl.rds: R/application/surrogate_index_estimation.R $(data-synthetic)
	Rscript --verbose R/application/surrogate_index_estimation.R synthetic > R/application/surrogate_index_estimation-synthetic.Rout 2> R/application/surrogate_index_estimation-synthetic.Rout
	
results/raw-results/application-synthetic/ma_trt_effects_tbl.rds: R/application/trial-level-effects.R results/raw-results/application-synthetic/ipd_surr_indices_tbl.rds
	Rscript --verbose R/application/trial-level-effects.R synthetic > R/application/trial-level-effects.Rout 2> R/application/trial-level-effects.Rout
	
R/application/meta_analysis-synthetic.Rout: R/application/meta_analysis.R results/raw-results/application-synthetic/ma_trt_effects_tbl.rds $(analysishelpers)
	Rscript --verbose R/application/meta_analysis.R synthetic > R/application/meta_analysis-synthetic.Rout 2> R/application/meta_analysis-synthetic.Rout
	
results/raw-results/application-synthetic/bayesian_ma_results.rds: R/application/bayesian-meta-analysis.R results/raw-results/application-synthetic/ma_trt_effects_tbl.rds
	Rscript --verbose R/application/bayesian-meta-analysis.R synthetic > R/application/bayesian-meta-analysis.Rout 2> R/application/bayesian-meta-analysis.Rout
	
R/application/processing-results-synthetic.Rout: R/application/processing-results.R results/raw-results/application-synthetic/bayesian_ma_results.rds R/application/meta_analysis-synthetic.Rout
	Rscript --verbose R/application/processing-results.R synthetic > R/application/processing-results-synthetic.Rout 2> R/application/processing-results-synthetic.Rout
	
	
	
additional-application/R/meta-analytic/estimate-landmark-models.Rout: additional-application/R/meta-analytic/estimate-landmark-models.R
	Rscript --verbose additional-application/R/meta-analytic/estimate-landmark-models.R > $@ 2> $@
	
R/meta-analytic/intermediate-objects/landmark_predictions_tbl.rds: additional-application/R/meta-analytic/estimate-landmark-models.Rout additional-application/R/meta-analytic/predicted-surrogates.R
	Rscript --verbose additional-application/R/meta-analytic/predicted-surrogates.R > $@ 2> $@
	
additional-application/results/raw-results/ma_trt_effects_tbl.rds: additional-application/R/meta-analytic/trial-level-effects.R R/meta-analytic/intermediate-objects/landmark_predictions_tbl.rds
	Rscript --verbose additional-application/R/meta-analytic/trial-level-effects.R > $@ 2> $@
	
additional-application/results/tables/surrogacy-inferences.csv: additional-application/R/meta-analytic/meta-analysis.R additional-application/results/raw-results/ma_trt_effects_tbl.rds
	Rscript --verbose additional-application/R/meta-analytic/meta-analysis.R > $@ 2> $@

additional-application/R/meta-analytic/processing-results.Rout: additional-application/R/meta-analytic/processing-results.R additional-application/results/raw-results/ma_trt_effects_tbl.rds additional-application/results/tables/surrogacy-inferences.csv
	Rscript --verbose additional-application/R/meta-analytic/processing-results.R > $@ 2> $@







