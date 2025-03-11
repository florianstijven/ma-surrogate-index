library(synthpop)
library(dplyr)

set.seed(1)

df = read.csv("data/processed_data.csv") %>%
  select(-X)

# Recode variables into the correct format for the the synthpop package. Factors
# should be factors, not 0/1 dummy variables. 
df = df %>%
  mutate(CC_stratum = as.factor(CC_stratum),
         BMI_stratum = case_when(
           BMI_underweight == 1 ~ "Underweight",
           BMI_normal == 1 ~ "Normal",
           BMI_overweight == 1 ~ "Overweight",
           BMI_obese == 1 ~ "Obese",
           .default = NA
         ),
         BMI_stratum = as.factor(BMI_stratum),
         trial.lbl = as.factor(trial.lbl))

# Some patients have titer values event though they were not sampled in the
# case-cohort sampling. We set the titer values for these patients to NA.
df = df %>%
  mutate(
    bindSpike = ifelse((Delta == 0) & (vax == 1), NA, bindSpike),
    pseudoneutid50 = ifelse((Delta == 0) & (vax == 1), NA, pseudoneutid50)
  )


# Check whether the synthpop recognizes the variables correctly. 
codebook.syn(df)

# Order in which variables are predicted. Variables which are not used for
# predictions (i.e., as predictors) and which are not predicted themselves,
# should not be in this character vector. The value of such variables is just
# kept unchanged.
visit.sequence.ini <- c(
  "Delta",
  "vax",
  "age.geq.65",
  "BMI_stratum",
  "CC_stratum",
  "risk_score",
  "Age",
  "bindSpike",
  "pseudoneutid50",
  "Y"
)

# Method used to predict each variable in the data set. Note that we have on
# element for each variable in the data set and the order of the character
# vector below corresponds to the order of the columns in the data set. An empty
# string means that that variable is not predicted (i.e., the original values
# are kept in the  synthetic data set).
method.ini <- c("",
                "ctree",
                "ctree",
                "",
                "ctree",
                "ctree",
                "",
                "ctree",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "")

# Generate synthetic data set. The variables are generated in each trial
# separately. We don't use the syn.strata() function because that function
# samples from the strata; consequently, the trial-specific sample sizes are not
# the same in the synthetic data as in the original data. 
sds_list = lapply(
  X = levels(df$trial.lbl),
  FUN = function(trial) {
    data_temp = df %>%
      filter(trial.lbl == .env$trial)
    sds <- syn(
      data = data_temp,
      visit.sequence = visit.sequence.ini,
      method = method.ini,
      m = 1,
      # Variables that are not predicted are kept in the synthetic data set.
      drop.not.used = FALSE,
      drop.pred.only = FALSE,
      # The CC_stratum variable has 58 categories, so we have to ensure that that
      # variable is kept as a categorical variable.
      maxfaclevels = 60,
      minnumlevels = 20,
      # The titers should be missing for patients with Delta == 0 in the vaccine
      # group.
      rules = list(bindSpike = "(Delta == 0) & vax == 1", 
                   pseudoneutid50 = "(Delta == 0) & vax == 1",
                   A = "vax == 0"),
      rvalues = list(bindSpike = NA, pseudoneutid50 = NA, A = 0)
    )
  }
)

# Combine the synthetic data sets contains in the sds objects in the list
# generated above into a single data set.
synthetic_df = lapply(
  X = sds_list,
  FUN = function(x) x$syn
) %>%
  bind_rows() %>%
  # Remove helper variables that were not present in the original data from
  # which we started in this script.
  select(-BMI_stratum)

# Save synthetic data set.
write.csv(synthetic_df, file = "data/processed_data_synthetic.csv")


