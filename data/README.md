# R scripts

This directory contains the following three R scripts. 

* `read_data.R` reads in and combines the data sets from the different COVID-19
trials. The combined data set is saved as `CrossProtocolData.csv` (which is not
included in the repo).
* `data-preparation.R` prepares `CrossProtocolData.csv` for analysis.
The pre-processed data are saved as `processed_data.csv`, which are again not
included in this repo.
* `generate-synthetic-data.R` reads in `processed_data.csv` and generates a
synthetic data set `processed_data_synthetic.csv` with the same structure as
`processed_data.csv`. This script also generates figures and tables that compare
the original with the synthetic data. These are saved in
`original-vs-synthetic-data/`.

The above R scripts use data sets as input that cannot be shared. Hence, none of
these R scripts can be rerun by other researchers. The data file output by
`generate-synthetic-data.R` (i.e., `processed_data_synthetic.csv`) is included
in this repo.

# Data

`processed_data_synthetic.csv` (and `processed_data.csv`) contains the following
variables:

* `treatment`. Treatment received. `treatment = 0` for placebo patients, `treatment = 1` for
patients who received the vaccine treatment. Note that the exact vaccine a
patient received depends on the trial.
* `bindSpike`. Titer in log_10 BAU/ml for the binding antibody. This value is set to the
universal lower limit of detection (defined below) for all placebo patients. For 
vaccine patients, this value is missing unless `Delta = 1`. 
* `pseudoneutID50`. Titer in log_10 ID50 for the neutralizing antibody. This value is set to the
universal lower limit of detection (defined below) for all placebo patients. For 
vaccine patients, this value is missing unless `Delta = 1`. 
* `case_cohort_weight_bAb`. Case-cohort sampling weight for the binding antibody. 
This is the inverse of the probability of being sampled to have the marker measured.
* `HighRiskInd`. Indicator for whether the participant has a high risk of severe
COVID-19. `HighRiskInd = 1` if the participant has a high risk, `HighRiskInd =
0` otherwise.
* `risk_score`. The COVID-19 baseline behavioral risk score defined for each trial 
based on an ensemble statistical learning algorithm trained on the placebo arm 
USG COVID-19 Response Team / Coronavirus Prevention Network (CoVPN) Biostatistics Team
et al. (2022). The risk score is defined as the log of the predicted risk.
* `event`. Event indicator for the SARS-CoV-2 infection outcome. `event = 1` if the 
participant had a virologically confirmed SARS-CoV-2 infection, `event = 0` otherwise.
* `Delta_bAb`. Case-cohort sampling indicator for the binding antibody. `Delta_bAb = 1` if the
participant was sampled to have the binding antibody titer measured, `Delta_bAb
= 0` otherwise.
* `time_to_event`. Time to infection event or censoring in days post marker measurement
visit.
* `Sex`. Biological sex.
* `Age`. Age in years.
* `Bserostatus`. Serostatus at baseline. `Bserostatus = 1` if the participant was seropositive
for SARS-CoV-2 at baseline, `Bserostatus = 0` otherwise.
* `case_cohort_weight_nAb`. Case-cohort sampling weight for the neutralizing antibody.
This is the inverse of the probability of being sampled to have the marker
measured.
* `Delta_nAb`. Case-cohort sampling indicator for the neutralizing antibody. `Delta_nAb = 1` if the
participant was sampled to have the neutralizing antibody titer measured,
`Delta_nAb = 0` otherwise.
* `trial`. This variables indicates the trial (or trial subunit for the Sanofi trial).
* `BMI_underweight`, `BMI_normal`, `BMI_overweight`, and `BMI_obese`. Indicator variables
for the BMI categories defined as follows:
  - Underweight: BMI < 18.5
  - Normal: 18.5 <= BMI < 25
  - Overweight: 25 <= BMI < 30
  - Obese: BMI > = 30
* `abrogation_coefficient`. This is the abrogation coefficient used for matching
the neutralizing titer to the circulating strains. The matched titer is
`pseudoneutid50_adjusted`.
* `pseudoneutid50_adjusted`. This is the adjusted neutralizing titer. See manuscript
for details on how this is computed.

# Others

The `VariantNeutralizationScores.docx` file contains details about the GMT estimates that 
were used for adjusting neutralization titers against the circulating variants. 



## Universal Lower Limit of Detection

The "universal lower limit of detection" is defined as the maximum of the
trial-specific lower limits of detection. The latter are defined as the lowest
antibody marker observed in any vaccine patient in the given trial.

Placebo patients get assigned this universal lower limit of detection to
maintain consistency across treatment groups. A patient from the vaccine group
with no measurable antibody marker has a value equal to the lower limit of detection. A
placebo patient, who by definition should have no measurable antibodies, should have
the same value. This would not be true if placebo patients get assigned a antibody marker 
value of zero.


# References

USG COVID-19 Response Team / Coronavirus Prevention Network (CoVPN) Biostatistics
Team, Gilbert, P. B., Fong, Y., Benkeser, D., Andriesen, J., Borate, B., Carone, M., Carpp,
L. N., Diaz, I., Fay, M. P., Fiore-Gartland, A., Hejazi, N. S., Huang, Y., Huang, Y.,
Hyrien, O., Janes, H. E., Juraska, M., Li, K., Luedtke, A., Nason, M., Randhawa, A. K.,
van der Laan, L., Williamson, B., Zhang, W., and Follmann, D. (2022). USG COVID-
19 Response Team / CoVPN vaccine efficacy trial immune correlates Statistical Analysis Plan.
https://figshare.com/articles/online_resource/CoVPN_OWS_COVID-19_
Vaccine_Efficacy_Trial_Immune_Correlates_SAP/13198595. Version 0.4.
