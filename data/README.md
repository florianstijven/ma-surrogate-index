# Rscripts

This directory contains the two Rscripts. 

* `read_data.R` reads in and combines the data sets from the different COVID-19 trials. The combined data set is saved as
`CrossProtocolData.csv` (which is not included in the repo).
* `data-preparation.R` processes `CrossProtocolData.csv` for analysis.
The processed data are saved as `processed_data.csv`, which are again not
included in this repo.
* `generate-synthetic-data.R` reads in `processed_data.csv` and generates a
synthetic data set `processed_data_synthetic.csv` with the same structure as
`processed_data.csv`. This script also generates figures and tables that compare
the original with the synthetic data. These are saved in
`original-vs-synthetic-data/`

# Data

`processed_data_synthetic.csv` contains the following variables:

* `A`. For placebo patients from any trial, we have `A = 0`. For vaccine patients,
`A` varies from `1` to `8` where each integer corresponds to a different trial (or 
  trial subunit for J&J). 
* `bindSpike`. Titer in log_10 BAU/ml for the binding antibody. This value is set to the
universal lower limit of detection (defined below) for all placebo patients. For 
vaccine patients, this value is missing unless `Delta = 1`. 
* `pseudoneutID50`. Titer in log_10 ID50 for the neutralizing antibody. This value is set to the
universal lower limit of detection (defined below) for all placebo patients. For 
vaccine patients, this value is missing unless `Delta = 1`. 
* `age.geq.65`. Age category indicator. `age.geq.65 = 1` for patients with age 65 years 
or higher; `age.geq.65 = 0` otherwise. 
* `risk_score`. The COVID-19 baseline behavioral risk score defined for each trial 
based on an ensemble statistical learning algorithm trained on the placebo arm 
USG COVID-19 Response Team / Coronavirus Prevention Network (CoVPN) Biostatistics Team
et al. (2022). The risk score is defined as the log of the predicted risk.
* `Y`. Infection outcome. `Y = 1` if the infection outcome ocurred, `Y = 0` otherwise.
* `Delta`. Case-cohort sampling indicator. `Delta = 1` if the patient was sampled to have
the titers measured, `Delta = 0` otherwise.
* `Age`. Age in years.
* `trial`. This variables indicates the trial as a integer between `1` and `8` where each integer corresponds to a different trial (or 
  trial subunit for J&J). This variable equals `A` for vaccine patients. 
* `vax`. Treatment received. `vax = 0` for placebo patients, `vax = 1` for 
patients who received the vaccine treatment. Note that the exact vaccine a patient 
received depends on the trial.
* `CC_stratum`. Case-cohort stratum. Note that all placebo patients are assigned to 
the `Placebo` stratum because their probability of having the titer measured is zero.
* `trial.lbl`. Character vector indicating the name of the trial to which a patients belongs.
This variables "corresponds" with the integers in `trial`.
* `BMI_underweight`, `BMI_normal`, `BMI_overweight`, and `BMI_obese`. Indicator variables
for the BMI categories defined as follows:
  - Underweight: BMI < 18.5
  - Normal: 18.5 <= BMI < 25
  - Overweight: 25 <= BMI < 30
  - Obese: BMI > = 30


## Universal Lower Limit of Detection

The "universal lower limit of detection" is defined as the maximum of the
trial-specific lower limits of detection. The latter are defined as the lowest
titer observed in any vaccine patients in the given trial.

Placebo patients get assigned this universal lower limit of detection to
mainting consistency across treatment groups. A patient from the vaccine group
with no measurable Abs has a titer equal to the lower limit of detection. A
placebo patient, who by definition should have no measurable Abs, should have
the same value. This would not be true if placebo patients get assigned a titer 
value of zero (which seems reasonable at first sight).


# References

USG COVID-19 Response Team / Coronavirus Prevention Network (CoVPN) Biostatistics
Team, Gilbert, P. B., Fong, Y., Benkeser, D., Andriesen, J., Borate, B., Carone, M., Carpp,
L. N., Diaz, I., Fay, M. P., Fiore-Gartland, A., Hejazi, N. S., Huang, Y., Huang, Y.,
Hyrien, O., Janes, H. E., Juraska, M., Li, K., Luedtke, A., Nason, M., Randhawa, A. K.,
van der Laan, L., Williamson, B., Zhang, W., and Follmann, D. (2022). USG COVID-
19 Response Team / CoVPN vaccine efficacy trial immune correlates Statistical Analysis Plan.
https://figshare.com/articles/online_resource/CoVPN_OWS_COVID-19_
Vaccine_Efficacy_Trial_Immune_Correlates_SAP/13198595. Version 0.4.



