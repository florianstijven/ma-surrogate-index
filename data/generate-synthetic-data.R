library(synthpop)
library(tidyverse)

# Set seed for reproducibility.
set.seed(2)

# Location to save figures/tables comparing the original and synthetic data.
results_dir = "data/original-vs-synthetic-data/"

df = read.csv("data/processed_data.csv") %>%
  select(-X)

# Recode variables into the correct format for the the synthpop package. Factors
# should be factors, not 0/1 dummy variables.
df = df %>%
  mutate(
    CC_stratum = as.factor(CC_stratum),
    BMI_stratum = case_when(
      BMI_underweight == 1 ~ "Underweight",
      BMI_normal == 1 ~ "Normal",
      BMI_overweight == 1 ~ "Overweight",
      BMI_obese == 1 ~ "Obese",
      .default = NA
    ),
    BMI_stratum = as.factor(BMI_stratum),
    trial.lbl = as.factor(trial.lbl)
  )

# Some patients have titer values event though they were not sampled in the
# case-cohort sampling. We set the titer values for these patients to NA.
df = df %>%
  mutate(
    bindSpike = ifelse((Delta == 0) & (vax == 1), NA, bindSpike),
    pseudoneutid50 = ifelse((Delta == 0) &
                              (vax == 1), NA, pseudoneutid50),
    Delta = as.integer(Delta)
  )


# Check whether the synthpop recognizes the variables correctly.
codebook.syn(df)

# Order in which variables are predicted. Variables which are not used for
# predictions (i.e., as predictors) and which are not predicted themselves,
# should not be in this character vector. The value of such variables is just
# kept unchanged.
visit.sequence.ini <- c(
  "Delta",
  "BMI_stratum",
  "age.geq.65",
  "risk_score",
  "Age",
  "bindSpike",
  "pseudoneutid50"
)

# Method used to predict each variable in the data set. Note that we have one
# element for each variable in the data set and the order of the character
# vector below corresponds to the order of the columns in the data set. An empty
# string means that that variable is not predicted (i.e., the original values
# are kept in the synthetic data set).
method.ini <- c(
  "",
  "parametric",
  "parametric",
  "satcat",
  "satcat",
  "",
  "",
  "parametric",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "",
  "satcat"
)

# Generate synthetic data set. The variables are generated in each trial
# separately. We don't use the syn.strata() function because that function
# samples from the strata; consequently, the trial-specific sample sizes are not
# the same in the synthetic data as in the original data.
strata_tbl = expand_grid(trial = levels(df$trial.lbl),
                         vax = 0:1,
                         Y = 0:1)
sds_list = purrr::pmap(
  .l = list(
    trial = strata_tbl$trial,
    vax = strata_tbl$vax,
    Y = strata_tbl$Y
  ),
  .f = function(trial, vax, Y) {
    # In the placebo group, we should not generate Ab titers since they are
    # constant by definition.
    if (vax == 0) {
      visit.sequence.ini_temp = visit.sequence.ini[-c(6, 7)]
      method.ini_temp = method.ini
      method.ini_temp[2:3] = ""
      
    } else {
      visit.sequence.ini_temp = visit.sequence.ini
      method.ini_temp = method.ini
    }
    data_temp = df %>%
      filter(trial.lbl == .env$trial, vax == .env$vax, Y == .env$Y)
    sds <- syn(
      data = data_temp,
      visit.sequence = visit.sequence.ini_temp,
      method = method.ini_temp,
      m = 1,
      # Variables that are not predicted are kept in the synthetic data set.
      drop.not.used = FALSE,
      drop.pred.only = FALSE,
      # The CC_stratum variable has 58 categories, so we have to ensure that that
      # variable is kept as a categorical variable.
      maxfaclevels = 60,
      minnumlevels = 20
    )
  }
)

summary(sds_list[[1]])

# Combine the synthetic data sets contains in the sds objects in the list
# generated above into a single data set.
synthetic_df = lapply(
  X = sds_list,
  FUN = function(x)
    x$syn
) %>%
  bind_rows() %>%
  # Make sure that the BMI dummy variables are consistent with BMI_stratum.
  mutate(
    BMI_underweight = ifelse(BMI_stratum == "Underweight", 1, 0),
    BMI_normal = ifelse(BMI_stratum == "Normal", 1, 0),
    BMI_overweight = ifelse(BMI_stratum == "Overweight", 1, 0),
    BMI_obese = ifelse(BMI_stratum == "Obese", 1, 0),
    BMI_underweight_normal = BMI_underweight + BMI_normal
  ) %>%
  # Remove helper variables that were not present in the original data from
  # which we started in this script.
  select(-BMI_stratum) %>%
  mutate(Delta = Delta == 1)

# Apply the missing data rules to the synthetic data. In principle, the syn()
# function can handle this. However, it samples the variables involved in the
# rules, while we want to keep those variables fixed.
synthetic_df = synthetic_df %>%
  mutate(
    bindSpike = ifelse((Delta == 0) & (vax == 1), NA, bindSpike),
    pseudoneutid50 = ifelse((Delta == 0) &
                              (vax == 1), NA, pseudoneutid50),
    A = ifelse(vax == 0, 0, A)
  )

# For some reason, NAs are generated for the titers in the synthetic data. 
synthetic_df %>%
  filter(Delta == 1, is.na(bindSpike))

# No such NAs are present in the original data. 
df %>%
  filter(Delta == 1, is.na(bindSpike))

# We solve this by replacing the NAs with a randomly sampled variable. 
synthetic_df = synthetic_df %>%
  mutate(
    bindSpike = ifelse(
      is.na(bindSpike) &
        Delta == 1 &
        vax == 1,
      runif(min = 0.75, max = 1, n = 1),
      bindSpike
    ),
    pseudoneutid50 = ifelse(
      is.na(pseudoneutid50) &
        Delta == 1 &
        vax == 1,
      runif(min = 0.5, max = 1, n = 1),
      pseudoneutid50
    )
  )

# Set the universal lower limits of detection for binding and neutralizing
# antibody titers.
llod_spike = log(10.84, base = 10)
llod_neut = log(2.61, base = 10)
# Observed values below the LLOD will be truncated the the LLOD divided by 2.
llod_spike_truncated = log(10.84 / 2, base = 10)
llod_neut_truncated = log(2.61 / 2, base = 10)


# Set all values of placebo to the LLOD divided by 2.
synthetic_df <- synthetic_df %>% mutate(bindSpike = ifelse(A == 0, llod_spike_truncated, bindSpike))
synthetic_df <- synthetic_df %>% mutate(pseudoneutid50 = ifelse(A == 0, llod_neut_truncated, pseudoneutid50))

# Patients with measured titers, who have a titer below the LLOD, will have
# their measurement truncated to the LLOD divided by 2.
synthetic_df <- synthetic_df %>% mutate(bindSpike = ifelse(
  Delta == T &
    bindSpike < llod_spike,
  llod_spike_truncated,
  bindSpike
))
synthetic_df <- synthetic_df %>% mutate(
  pseudoneutid50 = ifelse(
    Delta == T & pseudoneutid50 < llod_neut,
    llod_neut_truncated,
    pseudoneutid50
  )
)

# Compare the synthetic and original data in terms of case-cohort sampling and
# the number of events.
sink(file = paste0(results_dir, "n-events-delta-comparison.txt"))
bind_rows(df %>%
            mutate(data = "original"),
          synthetic_df %>%
            mutate(data = "synthetic")) %>%
  group_by(trial.lbl, vax, data, CC_stratum) %>%
  summarise(n = n(),
            n_events = sum(Y),
            n_Delta = sum(Delta)) %>%
  ungroup() %>%
  arrange(CC_stratum) %>%
  print(n = 500)
sink()

bind_rows(df %>%
            mutate(Data = "Original"),
          synthetic_df %>%
            mutate(Data = "Synthetic")) %>%
  mutate(Treatment = ifelse(vax == 1, "Vaccine", "Placebo"),
         Age = ifelse(age.geq.65, "Age >= 65", "Age < 65")) %>%
  ggplot(aes(x = risk_score, y = trial.lbl, fill = Data)) +
  geom_violin() +
  facet_grid(. ~ Age) +
  ylab(NULL) +
  xlab("Risk Score")

ggsave(
  filename = paste0(results_dir, "risk-score-by-trial.pdf"),
  device = "pdf",
  width = double_width,
  height = double_height,
  units = unit
)

bind_rows(df %>%
            mutate(Data = "Original"),
          synthetic_df %>%
            mutate(Data = "Synthetic")) %>%
  filter(vax == 1) %>%
  ggplot(aes(x = risk_score, y = bindSpike, color = Data)) +
  geom_point(alpha = 0.05) +
  geom_smooth(se = FALSE) +
  facet_wrap("trial.lbl") +
  xlab("Risk Score") +
  ylab(expression(paste("Spike Protein IgG (", log[10], "BAU/ml)")))
ggsave(
  filename = paste0(results_dir, "bindSpike-vs-risk-score.pdf"),
  device = "pdf",
  width = double_width,
  height = double_height,
  units = unit
)

bind_rows(df %>%
            mutate(Data = "Original"),
          synthetic_df %>%
            mutate(Data = "Synthetic")) %>%
  filter(vax == 1) %>%
  ggplot(aes(x = risk_score, y = pseudoneutid50, color = Data)) +
  geom_point(alpha = 0.05) +
  geom_smooth(se = FALSE) +
  facet_wrap("trial.lbl") +
  xlab("Risk Score") +
  ylab(expression(paste("Neutralizing Antibody (", log[10], "ID50)")))
ggsave(
  filename = paste0(results_dir, "pseudoneutid50-vs-risk-score.pdf"),
  device = "pdf",
  width = double_width,
  height = double_height,
  units = unit
)

# Save synthetic data set.
write.csv(synthetic_df, file = "data/processed_data_synthetic.csv")


