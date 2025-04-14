library(tidyverse)


p3001_file = "S:/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_20230919.csv"
p3002_file = "S:/p3002/analysis/correlates/Part_A_Blinded_Phase_Data/adata/azd1222_data_processed_with_riskscore.csv"
p3003_file = "S:/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/janssen_pooled_partA_data_processed_with_riskscore_20240305.csv"
p3004_file = "S:/p3004/analysis/correlates/Part_A_Blinded_Phase_Data/adata/prevent19_data_processed_20250325.csv"
p3005_file = "S:/p3005/analysis/correlates/Part_A_Blinded_Phase_Data/adata/vat08_combined_data_processed_20250321.csv"


# Read the data from the individual trials.

# Moderna
p3001 = read.csv(p3001_file) %>% filter(Bserostatus == 0) %>%
  dplyr::select(
    Ptid,
    Trt,
    Day57bindSpike,
    Day57pseudoneutid50,
    wt.D57,
    age.geq.65,
    HighRiskInd,
    risk_score,
    CalendarDateEnrollment,
    Wstratum,
    EventIndPrimaryD57,
    ph1.D57,
    ph2.D57,
    EventTimePrimaryD57,
    WhiteNonHispanic,
    Sex,
    Age,
    BMI
  ) %>%
  rename(
    bindSpike = Day57bindSpike,
    pseudoneutid50 = Day57pseudoneutid50,
    case_cohort_weight_bAb = wt.D57,
    phase1_ind = ph1.D57,
    Delta_bAb = ph2.D57,
    event = EventIndPrimaryD57,
    time_to_event = EventTimePrimaryD57
  ) %>%
  mutate(
    USUBJID = paste0("mRNA-1273-P301-", stringr::str_sub(Ptid, 1, 5), "-", stringr::str_sub(Ptid, 6, 9)),
    case_cohort_weight_nAb = case_cohort_weight_bAb,
    Delta_nAb = Delta_bAb,
    protocol = "p3001",
    trial = "Moderna"
  )

# AstraZeneca
p3002 = read.csv(p3002_file) %>% filter(Bserostatus == 0) %>%
  dplyr::select(
    Ptid,
    Trt,
    Day57bindSpike,
    Day57pseudoneutid50,
    wt.D57,
    age.geq.65,
    HighRiskInd,
    risk_score,
    CalendarDateEnrollment,
    Wstratum,
    EventIndPrimaryD57,
    ph1.D57,
    ph2.D57,
    EventTimePrimaryD57,
    Country,
    WhiteNonHispanic,
    Sex,
    Age,
    BMI
  ) %>%
  rename(
    bindSpike = Day57bindSpike,
    pseudoneutid50 = Day57pseudoneutid50,
    case_cohort_weight_bAb = wt.D57,
    phase1_ind = ph1.D57,
    Delta_bAb = ph2.D57,
    event = EventIndPrimaryD57,
    time_to_event = EventTimePrimaryD57
  ) %>%
  mutate(
    USUBJID = Ptid,
    case_cohort_weight_nAb = case_cohort_weight_bAb,
    Delta_nAb = Delta_bAb,
    protocol = "p3002",
    trial = "AstraZeneca"
  )

# J&J
p3003 = read.csv(p3003_file) %>% filter(Bserostatus == 0) %>%
  dplyr::select(
    Ptid,
    Trt,
    Day29bindSpike,
    Day29pseudoneutid50,
    wt.D29,
    age.geq.65,
    HighRiskInd,
    risk_score,
    CalendarDateEnrollment,
    Country,
    Wstratum,
    EventIndPrimaryD29,
    ph1.D29,
    ph2.D29,
    EventTimePrimaryD29,
    Country,
    WhiteNonHispanic,
    Sex,
    Age,
    BMI
  ) %>%
  rename(
    bindSpike = Day29bindSpike,
    pseudoneutid50 = Day29pseudoneutid50,
    case_cohort_weight_bAb = wt.D29,
    phase1_ind = ph1.D29,
    Delta_bAb = ph2.D29,
    event = EventIndPrimaryD29,
    time_to_event = EventTimePrimaryD29
  ) %>%
  mutate(
    USUBJID = Ptid,
    case_cohort_weight_nAb = case_cohort_weight_bAb,
    Delta_nAb = Delta_bAb,
    protocol = "p3003",
    trial = "Janssen"
  )

# Novavax
p3004 = read.csv(p3004_file) %>% filter(Bserostatus == 0) %>%
  dplyr::select(
    Ptid,
    Trt,
    Day35bindSpike,
    Day35pseudoneutid50,
    wt.D35,
    age.geq.65,
    HighRiskInd,
    risk_score2,
    CalendarDateEnrollment,
    Wstratum,
    EventIndPrimaryD35,
    ph1.D35,
    ph2.D35,
    EventTimePrimaryD35,
    Country,
    WhiteNonHispanic,
    Sex,
    Age,
    BMI
  ) %>%
  rename(
    bindSpike = Day35bindSpike,
    pseudoneutid50 = Day35pseudoneutid50,
    case_cohort_weight_bAb = wt.D35,
    phase1_ind = ph1.D35,
    Delta_bAb = ph2.D35,
    event = EventIndPrimaryD35,
    time_to_event = EventTimePrimaryD35,
    risk_score = risk_score2
  ) %>%
  mutate(
    USUBJID = paste0("2019nCoV-", stringr::str_sub(Ptid, 9, 22)),
    case_cohort_weight_nAb = case_cohort_weight_bAb,
    Delta_nAb = Delta_bAb,
    protocol = "p3004",
    trial = "Novavax"
  )

# Sanofi
p3005 = read.csv(p3005_file) %>% 
  filter(Bserostatus == 0) %>%
  dplyr::select(
    Ptid,
    Trt,
    Day43pseudoneutid50,
    Day43bindSpike,
    wt.D43.bAb,
    wt.D43.nAb,
    age.geq.65,
    HighRiskInd, 
    risk_score, 
    CalendarDateEnrollment,
    Wstratum, 
    EventIndPrimaryD43,
    ph1.D43, 
    ph2.D43.bAb,
    ph2.D43.nAb,
    EventTimePrimaryD43, 
    Country,
    WhiteNonHispanic,
    Sex,
    Age,
    BMI
  ) %>%
  rename(
    bindSpike = Day43bindSpike,
    pseudoneutid50 = Day43pseudoneutid50,
    case_cohort_weight_bAb = wt.D43.bAb,
    case_cohort_weight_nAb = wt.D43.nAb,
    phase1_ind = ph1.D43,
    Delta_bAb = ph2.D43.bAb,
    Delta_nAb = ph2.D43.nAb,
    event = EventIndPrimaryD43,
    time_to_event = EventTimePrimaryD43,
    risk_score = risk_score
  ) %>%
  mutate(
    USUBJID = Ptid,
    protocol = "p3005",
    trial = "Sanofi"
  )


# combine trials into one dataframe
data_all = bind_rows(
  p3001,
  p3002,
  p3003,
  p3004,
  p3005
) %>%
  rename(treatment = Trt) 

# Record number of rows before filtering the data. 
n1 = nrow(data_all)

# Drop patients that were not included in the phase 1 set. 
data_all = data_all %>%
  filter(phase1_ind == TRUE) 

n2 = nrow(data_all)
n1 - n2
# 9899 rows were dropped.


# Drop patients with missing covariates. 
data_all = data_all %>%
  filter(if_all(
    c(
      age.geq.65,
      HighRiskInd,
      CalendarDateEnrollment,
      Sex,
      WhiteNonHispanic,
      BMI
    ),
    complete.cases
  )) 

n3 = nrow(data_all)
n2 - n3 
# 461 rows were dropped.


write.csv(data_all, "data/CrossProtocolData.csv")

