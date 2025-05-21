library(synthpop)
library(tidyverse)

# Set seed for reproducibility.
set.seed(2)

# Location to save figures/tables comparing the original and synthetic data.
results_dir = "data/original-vs-synthetic-data/"

df = read.csv("data/processed_data.csv") %>%
  select(-X)

# Drop redundant variables.
df = df %>%
  select(-Ptid, -age.geq.65, -Wstratum, -BMI, -USUBJID, -protocol, -total_proportion_row,, -total_proportion_row_new)

# Recode variables into the correct format for the the synthpop package. Factors
# should be factors, not 0/1 dummy variables.
df = df %>%
  mutate(
    BMI_stratum = case_when(
      BMI_underweight == 1 ~ "Underweight",
      BMI_normal == 1 ~ "Normal",
      BMI_overweight == 1 ~ "Overweight",
      BMI_obese == 1 ~ "Obese",
      .default = NA
    ),
    BMI_stratum = as.factor(BMI_stratum),
    Delta_bAb = as.integer(Delta_bAb),
    Delta_nAb = as.integer(Delta_nAb)
  )

# Check whether the synthpop recognizes the variables correctly.
codebook.syn(df)

# Order in which variables are predicted. Variables which are not used for
# predictions (i.e., as predictors) and which are not predicted themselves,
# should not be in this character vector. The value of such variables is just
# kept unchanged.
visit.sequence.ini <- c(
  "Delta_nAb",
  "Delta_bAb",
  "BMI_stratum",
  "Age",
  "Sex",
  "risk_score",
  "bindSpike",
  "pseudoneutid50",
  "abrogation_coefficient"
)

# Method used to predict each variable in the data set. Note that we have one
# element for each variable in the data set and the order of the character
# vector below corresponds to the order of the columns in the data set. An empty
# string means that that variable is not predicted (i.e., the original values
# are kept in the synthetic data set).
method.ini <- c(
  "",
  "norm",
  "norm",
  "",
  "logreg",
  "norm",
  "survctree",
  "",
  "survctree",
  "logreg",
  "norm",
  "logreg",
  rep("", 8),
  "cart",
  "",
  "polyreg"
)

# Generate synthetic data set. The variables are generated in each trial
# separately. We don't use the syn.strata() function because that function
# samples from the strata; consequently, the trial-specific sample sizes are not
# the same in the synthetic data as in the original data.
strata_tbl = expand_grid(
  trial = unique(df$trial),
  treatment = 0:1,
  event = 0:1,
  Bserostatus = 0:1
) %>%
  # Only in the Sanofi trials are there subjects with Bserostatus == 1.
  filter((trial %in% c("Sanofi-1", "Sanofi-2") |
            (!(
              trial %in% c("Sanofi-1", "Sanofi-2")
            ) & Bserostatus == 0)))
sds_list = purrr::pmap(
  .l = list(
    trial = strata_tbl$trial,
    treatment = strata_tbl$treatment,
    event = strata_tbl$event,
    Bserostatus = strata_tbl$Bserostatus
  ),
  .f = function(trial, treatment, event, Bserostatus) {
    # In the placebo group, we should not generate Ab titers since they are
    # constant by definition.
    if (treatment == 0 & !(trial %in% c("Sanofi-1", "Sanofi-2"))) {
      visit.sequence.ini_temp = visit.sequence.ini[-c(7, 8)]
      method.ini_temp = method.ini
      method.ini_temp[2:3] = ""
      
    } else {
      visit.sequence.ini_temp = visit.sequence.ini
      method.ini_temp = method.ini
    }
    data_temp = df %>%
      filter(trial == .env$trial, treatment == .env$treatment, event == .env$event, Bserostatus == .env$Bserostatus)
    sds <- syn(
      data = data_temp,
      visit.sequence = visit.sequence.ini_temp,
      method = method.ini_temp,
      m = 1,
      # Variables that are not predicted are kept in the synthetic data set.
      drop.not.used = FALSE,
      drop.pred.only = FALSE,
      minnumlevels = 5,
      event = list(time_to_event = "event")
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
  select(-BMI_stratum)

# Apply the missing data rules to the synthetic data. In principle, the syn()
# function can handle this. However, it samples the variables involved in the
# rules, while we want to keep those variables fixed.
synthetic_df = synthetic_df %>%
  mutate(
    bindSpike = ifelse((Delta_bAb == 0) & (treatment == 1), NA, bindSpike),
    pseudoneutid50 = ifelse((Delta_nAb == 0) &
                              (treatment == 1), NA, pseudoneutid50)
  )


# For some reason, NAs are generated for the titers in the synthetic data. 
synthetic_df %>%
  filter(Delta_bAb == 1, is.na(bindSpike))

synthetic_df %>%
  filter(Delta_nAb == 1, is.na(pseudoneutid50))

# No such NAs are present in the original data. 
df %>%
  filter(Delta_bAb == 1, is.na(bindSpike))

df %>%
  filter(Delta_nAb == 1, is.na(pseudoneutid50))

# We solve this by replacing the NAs with a randomly sampled variable. 
synthetic_df = synthetic_df %>%
  mutate(
    bindSpike = ifelse(
      is.na(bindSpike) &
        Delta_bAb == 1 &
        treatment == 1,
      runif(min = 0.75, max = 1, n = 1),
      bindSpike
    ),
    pseudoneutid50 = ifelse(
      is.na(pseudoneutid50) &
        Delta_nAb == 1 &
        treatment == 1,
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
synthetic_df <- synthetic_df %>% mutate(bindSpike = ifelse(treatment == 0 & !(trial %in% c("Sanofi-1", "Sanofi-2")), llod_spike_truncated, bindSpike))
synthetic_df <- synthetic_df %>% mutate(pseudoneutid50 = ifelse(treatment == 0 & !(trial %in% c("Sanofi-1", "Sanofi-2")), llod_neut_truncated, pseudoneutid50))

# Patients with measured titers, who have a titer below the LLOD, will have
# their measurement truncated to the LLOD divided by 2.
synthetic_df <- synthetic_df %>% mutate(bindSpike = ifelse(
  Delta_bAb == T &
    bindSpike < llod_spike,
  llod_spike_truncated,
  bindSpike
))
synthetic_df <- synthetic_df %>% mutate(
  pseudoneutid50 = ifelse(
    Delta_nAb == T & pseudoneutid50 < llod_neut,
    llod_neut_truncated,
    pseudoneutid50
  )
)

# Compute the adjusted neutralizing titer from the sampled titer and abrogation
# coefficient.
synthetic_df = synthetic_df %>%
  mutate(pseudoneutid50_adjusted = pseudoneutid50 - log(abrogation_coefficient, base = 10))

# Because of the adjustment, some patients may not have an adjusted titer below
# the LLOD. We truncated these values as before.
# Patients with measured titers, who have a titer below the LLOD, will have
# their measurement truncated to the LLOD divided by 2.
synthetic_df <- synthetic_df %>% mutate(
  pseudoneutid50_adjusted = ifelse(
    pseudoneutid50_adjusted < llod_neut,
    llod_neut_truncated,
    pseudoneutid50_adjusted
  )
)

# Compare the synthetic and original data in terms of case-cohort sampling and
# the number of events.
sink(file = paste0(results_dir, "n-events-delta-comparison.txt"))
bind_rows(df %>%
            mutate(data = "original"),
          synthetic_df %>%
            mutate(data = "synthetic")) %>%
  group_by(trial, treatment, data) %>%
  summarise(n = n(),
            n_events = sum(event),
            n_Delta_nAb = sum(Delta_nAb),
            n_Delta_bAb = sum(Delta_bAb)) %>%
  ungroup() %>%
  print(n = 500)
sink()

bind_rows(df %>%
            mutate(Data = "Original"),
          synthetic_df %>%
            mutate(Data = "Synthetic")) %>%
  mutate(Treatment = ifelse(treatment == 1, "Vaccine", "Placebo")) %>%
  ggplot(aes(x = risk_score, y = trial, fill = Data)) +
  geom_violin() +
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
  filter(treatment == 1) %>%
  ggplot(aes(x = risk_score, y = bindSpike, color = Data)) +
  geom_point(alpha = 0.05) +
  geom_smooth(se = FALSE) +
  facet_wrap("trial") +
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
  filter(treatment == 1) %>%
  ggplot(aes(x = risk_score, y = pseudoneutid50, color = Data)) +
  geom_point(alpha = 0.05) +
  geom_smooth(se = FALSE) +
  facet_wrap("trial") +
  xlab("Risk Score") +
  ylab(expression(paste("Neutralizing Antibody (", log[10], "ID50)")))
ggsave(
  filename = paste0(results_dir, "pseudoneutid50-vs-risk-score.pdf"),
  device = "pdf",
  width = double_width,
  height = double_height,
  units = unit
)

bind_rows(df %>%
            mutate(Data = "Original"),
          synthetic_df %>%
            mutate(Data = "Synthetic")) %>%
  filter(treatment == 1) %>%
  ggplot(aes(x = risk_score, y = pseudoneutid50_adjusted, color = Data)) +
  geom_point(alpha = 0.05) +
  geom_smooth(se = FALSE) +
  facet_wrap("trial") +
  xlab("Risk Score") +
  ylab(expression(paste("Neutralizing Antibody (", log[10], "ID50)")))
ggsave(
  filename = paste0(results_dir, "pseudoneutid50-adjusted-vs-risk-score.pdf"),
  device = "pdf",
  width = double_width,
  height = double_height,
  units = unit
)

# Save synthetic data set.
write.csv(synthetic_df, file = "data/processed_data_synthetic.csv")


