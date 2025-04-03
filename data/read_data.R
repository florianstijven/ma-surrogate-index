library(tidyverse)


p3001_file = "S:/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/moderna_real_data_processed_20230919.csv"
p3002_file = "S:/p3002/analysis/correlates/Part_A_Blinded_Phase_Data/adata/azd1222_data_processed_with_riskscore.csv"
p3003_file = "S:/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/janssen_pooled_partA_data_processed_with_riskscore_20240305.csv"
p3004_file = "S:/p3004/analysis/correlates/Part_A_Blinded_Phase_Data/adata/prevent19_data_processed_20250325.csv"



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
    wt.est = wt.D57,
    ph1 = ph1.D57,
    Delta = ph2.D57,
    Y = EventIndPrimaryD57,
    TTY = EventTimePrimaryD57
  ) %>%
  mutate(
    trial = 1,
    Y. = ifelse(Y == 1 & TTY <= 120, 1, 0),
    TTY. = ifelse(Y. == 1, TTY, ifelse(TTY <= 120, TTY, 120)),
    trial.lbl = "Moderna",
    protocol = "p3001"
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
    wt.est = wt.D57,
    ph1 = ph1.D57,
    Delta = ph2.D57,
    Y = EventIndPrimaryD57,
    TTY = EventTimePrimaryD57
  ) %>%
  mutate(
    trial = 2,
    Trt = ifelse(Trt == 1, 2, 0),
    Y. = ifelse(Y == 1 & TTY <= 120, 1, 0),
    TTY. = ifelse(Y. == 1, TTY, ifelse(TTY <= 120, TTY, 120)),
    trial.lbl = "AstraZeneca",
    protocol = "p3002"
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
    wt.est = wt.D29,
    ph1 = ph1.D29,
    Delta = ph2.D29,
    Y = EventIndPrimaryD29,
    TTY = EventTimePrimaryD29
  ) %>%
  mutate(
    Y. = ifelse(Y == 1 & TTY <= 120, 1, 0),
    TTY. = ifelse(Y. == 1, TTY, ifelse(TTY <= 120, TTY, 120)),
    trial.lbl = "Novavax",
    protocol = "p3003"
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
    wt.est = wt.D35,
    ph1 = ph1.D35,
    Delta = ph2.D35,
    Y = EventIndPrimaryD35,
    TTY = EventTimePrimaryD35,
    risk_score = risk_score2
  ) %>%
  mutate(
    trial = 3,
    Trt = ifelse(Trt == 1, 3, 0),
    Y. = ifelse(Y == 1 & TTY <= 120, 1, 0),
    TTY. = ifelse(Y. == 1, TTY, ifelse(TTY <= 120, TTY, 120)),
    protocol = "p3004"
  )


# separate p3003 into 5 separate trial units
p3003_brazil = p3003 %>% filter(Country == 2) %>%
  mutate(trial = 4,
         Trt = ifelse(Trt == 1, 4, 0),
         trial.lbl = "J&J (Brazil)")

p3003_colombia = p3003 %>% filter(Country == 4) %>%
  mutate(trial = 5,
         Trt = ifelse(Trt == 1, 5, 0),
         trial.lbl = "J&J (Colombia)")

p3003_southam = p3003 %>% filter(Country %in% c(1, 3, 5, 6)) %>%
  mutate(trial = 6,
         Trt = ifelse(Trt == 1, 6, 0),
         trial.lbl = "J&J (S. America)")

p3003_zaf = p3003 %>% filter(Country == 7) %>%
  mutate(trial = 7,
         Trt = ifelse(Trt == 1, 7, 0),
         trial.lbl = "J&J (S. Africa)")

p3003_usa = p3003 %>% filter(Country == 0) %>%
  mutate(trial = 8,
         Trt = ifelse(Trt == 1, 8, 0),
         trial.lbl = "J&J (USA)") 

# combine trials into one dataframe
dat = bind_rows(
  p3001,
  p3002,
  p3004,
  p3003_brazil,
  p3003_colombia,
  p3003_southam,
  p3003_zaf,
  p3003_usa
) %>%
  rename(A = Trt) %>%
  filter(ph1 == TRUE) %>% # drops 9570 pts
  mutate(Wstratum = as.factor(Wstratum)) %>%
  filter(if_all(
    c(
      age.geq.65,
      HighRiskInd,
      CalendarDateEnrollment,
      Sex,
      WhiteNonHispanic
    ),
    complete.cases
  )) # filters out 2 missing values of HighRiskInd 


dat = dat %>% rename(w = wt.est)

## Replace relative calendar date of enrollment for J&J trials
dat = dat %>% group_by(trial) %>% mutate(CalendarDate2 = ifelse(
  trial %in% c(4, 5, 6, 7, 8),
  CalendarDateEnrollment -
    min(CalendarDateEnrollment),
  CalendarDateEnrollment
))

## EARLYIND - was the pt enrolled during the first or second half of the enrollment period?
dat <- dat %>% group_by(trial) %>% mutate(EarlyInd = ifelse(CalendarDate2 < (max(CalendarDate2 /
                                                                                   2)), 1, 0))


write.csv(dat, "data/CrossProtocolData.csv")

