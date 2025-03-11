# Load required packages
library(tidyr)
library(dplyr)

# Load the data set. Wstratum refers to the strata that determine the
# probability of being sampled for measuring S. 
df = read.csv("data/CrossProtocolData_moderna_tgt_FINAL.csv") %>% 
  mutate(Wstratum = as.factor(Wstratum)) %>%
  droplevels() 

# Add variable indicating whether the patient was vaccinated with any vaccine or
# received placebo.
df = df %>%
  mutate(vax = ifelse(df$A == 0, 0, 1))

# Compute LODs for each trial and type of surrogate.
LOD <- df %>% group_by(trial) %>% dplyr::summarise(
  LLOD.Spike = min(bindSpike, na.rm = T),
  LLOD.Neut = min(pseudoneutid50, na.rm = T)
)

# Maximum LOQs across all trials. 
maxLLODSpike <- max(LOD$LLOD.Spike)
maxLLODNeut <- max(LOD$LLOD.Neut)

# Set all values of placebo to the max of the LLODs
df <- df %>% mutate(bindSpike = ifelse(A == 0, maxLLODSpike, bindSpike))
df <- df %>% mutate(pseudoneutid50 = ifelse(A == 0, maxLLODNeut, pseudoneutid50))

# Truncate surrogate measurements to the LOD across trials. 
df <- df %>% mutate(bindSpike = ifelse(Delta == T &
                                         bindSpike < maxLLODSpike, maxLLODSpike, bindSpike))
df <- df %>% mutate(
  pseudoneutid50 = ifelse(
    Delta == T & pseudoneutid50 < maxLLODNeut,
    maxLLODNeut,
    pseudoneutid50
  )
)

# Rowwise deletetion for missing BMI values. 
df <- df %>% filter(!is.na(BMI))

# Compute strata for the case-cohort sampling. 
df = df %>%
  mutate(
    CC_stratum = paste(trial, Y, A, Wstratum, sep = "_"),
    CC_stratum = ifelse(A == 0, "Placebo", CC_stratum),
    CC_stratum = as.factor(CC_stratum)
  )

# Drop variables which are currently not used in any analyses.
df = df %>%
  select(
    -X,
    -Ptid,
    -HighRiskInd,
    -CalendarDateEnrollment,
    -ph1,
    -WhiteNonHispanic,
    -STUDYID,
    -USUBJID, 
    -protocol, 
    -hd_mean,
    -Country, 
    -SUBJID, 
    -EarlyInd,
    -TTY,
    -Sex,
    -Wstratum
  )

# Add label for the trials.
df = df %>%
  mutate(trial.lbl = sapply(
    X = trial,
    FUN = function(x) {
      switch (
        as.character(x),
        "1" = "Moderna",
        "2" = "AstraZeneca",
        "3" = "Novavax",
        "4" = "J&J (Brazil)",
        "5" = "J&J (Colombia)",
        "6" = "J&J (S. America)",
        "7" = "J&J (S. Africa)",
        "8" = "J&J (USA)"
      )
    }
  ))

# For the non-JnJ trials, BMI should be categorized into <18.5, [18.5, 25), [25,
# 30), >= 30.
df = df %>%
  mutate(
    BMI_underweight = as.integer(ifelse(
      stringr::str_detect(trial.lbl, "J&J"), BMI == 1, BMI < 18.5
    )),
    BMI_normal = as.integer(ifelse(
      stringr::str_detect(trial.lbl, "J&J"),
      BMI == 2,
      BMI >=  18.5 &
        BMI < 25
    )),
    BMI_overweight = as.integer(ifelse(
      stringr::str_detect(trial.lbl, "J&J"), BMI == 3, BMI >= 25 &
        BMI < 30
    )),
    BMI_obese = as.integer(ifelse(
      stringr::str_detect(trial.lbl, "J&J"), BMI == 4, BMI >= 30
    )),
    BMI_underweight_normal = as.integer(ifelse(
      stringr::str_detect(trial.lbl, "J&J"), BMI <= 2, BMI < 25
    ))
  )

df = df %>%
  filter(!is.na(risk_score)) 

# Remove variables that are used in any of the current analyses.
df = df %>%
  select(-wt.est, -BMI)


# Save processed data set.
write.csv(df, "data/processed_data.csv")
