#!/usr/bin/env Rscript

## Setup ----------------------------------------------------------------------------

# Specify options for saving the plots to files
figures_dir = "results/figures/application/data-exploration"
tables_dir = "results/tables/application/data-exploration"

# Load required packages. 
library(tidyverse)
library(patchwork)
library(survival)
library(ggsurvfit)
library(tidycmprsk)
library(ggpubr)

time_cumulative_incidence = 80

## Load data and Prepare for Analysis -----------------------------------

# Load the data set. 
df = read.csv("data/processed_data.csv") %>%
  select(-X) 

# Split Sanofi trials
df = df %>%
  mutate(trial = ifelse(
    trial == "Sanofi-1",
    ifelse(Bserostatus == 1, "Sanofi-1 non-naive", "Sanofi-1 naive"),
    ifelse(
      trial == "Sanofi-2",
      ifelse(Bserostatus == 1, "Sanofi-2 non-naive", "Sanofi-2 naive"),
      trial
    )
  )) 


# Rename the trials and put them in chronological order.
df = df %>%
  # Reorder trial factor.
  mutate(
    trial = forcats::fct_recode(
      trial,
      "Moderna (naive)" = "Moderna",
      "AstraZeneca (naive)" = "AstraZeneca",
      "Janssen (naive)" = "Janssen",
      "Novavax (naive)" = "Novavax",
      "Sanofi 1 (naive)" = "Sanofi-1 naive",
      "Sanofi 1 (non-naive)" = "Sanofi-1 non-naive",
      "Sanofi 2 (naive)" = "Sanofi-2 naive",
      "Sanofi 2 (non-naive)" = "Sanofi-2 non-naive"
    ),
    trial = fct_relevel(
      trial,
      "Moderna (naive)",
      "AstraZeneca (naive)",
      "Janssen (naive)",
      "Novavax (naive)",
      "Sanofi 1 (naive)",
      "Sanofi 1 (non-naive)",
      "Sanofi 2 (naive)",
      "Sanofi 2 (non-naive)"
    )
  )

# Add `Treatment` character variables indicating placebo versus vaccine.
df = df %>%
  mutate(Treatment = ifelse(treatment == 1, "Vaccine", "Placebo"))

# Compute pseudo-values. These are computed for each trial separately. We first
# fit a separate KM curve by trial and treatment group.
surv_fit_all = survfit(Surv(time_to_event, event)~strata(treatment, trial), data = df)
df$pseudo_value = 1 - pseudo(surv_fit_all, times = time_cumulative_incidence)

# Compute placebo cumulative incidence rate per trial. This will be used later
# on on some plots to standardize some things across trials.
cumulative_incidence_control_tbl = df %>%
  group_by(trial, treatment) %>%
  summarise(prob_infection_free = mean(pseudo_value)) %>%
  ungroup() %>%
  filter(treatment == 0) %>%
  select(-treatment)

# Define new variable indicating whether infection was observed by 80 days.
df = df %>%
  mutate(infection_120d = ifelse(time_to_event <= time_cumulative_incidence & event == 1, 1, 0))



# The case-cohort sampling weights are already present in this data set. We add
# a second set of weights where we ensure that each trial-treatment group
# combination gets the same total weight. We don't want one trial to have an
# undue influence on the predictions.
df = df %>%
  mutate(
    weight_prediction = case_cohort_weight_nAb
  )
df = df %>%
  left_join(df %>%
              group_by(treatment, trial) %>%
              summarize(total_weight = sum(weight_prediction)) %>%
              ungroup()) %>%
  mutate(weight_prediction = weight_prediction / total_weight)


# Convert the data set to long format. This simplifies plotting.
df_long = bind_rows(
  df %>%
    mutate(
      surrogate = "bAb",
      titer = bindSpike,
      Delta = Delta_bAb,
      sample_weight = case_cohort_weight_bAb
    ),
  df %>%
    mutate(
      surrogate = "nAb",
      titer = pseudoneutid50,
      Delta = Delta_nAb,
      sample_weight = case_cohort_weight_nAb
    ),
  df %>%
    mutate(
      surrogate = "Adjusted nAb",
      titer = pseudoneutid50_adjusted,
      Delta = Delta_nAb,
      sample_weight = case_cohort_weight_nAb
    )
)

# Figures -------------------------------------------------------------

## Distribution of Baseline Covariates --------------------------------

df %>%
  ggplot(aes(x = risk_score, y = trial, fill = Treatment)) +
  geom_boxplot(color = "black") +
  ylab("Trial") +
  xlab("Risk Score") +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(
  filename = "distribution-risk-score.pdf",
  device = "pdf",
  path = figures_dir,
  width = double_width,
  height = double_height, 
  units = unit
)

df %>%
  ggplot(aes(x = Age, y = trial, fill = Treatment)) +
  geom_boxplot(color = "black") +
  ylab("Trial") +
  xlab("Age (years)") +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(
  filename = "distribution-age.pdf",
  device = "pdf",
  path = figures_dir,
  width = double_width,
  height = double_height, 
  units = unit
)

## Distribution of Surrogates -----------------------------------------

# Distribution of the neutralization and binding Ab titers by trial. Note that
# we use the case-cohort sampling weights in these boxplots.
df_long %>%
  filter(!is.na(titer)) %>%
  ggplot(aes(x = titer, y = trial)) +
  geom_boxplot(aes(weight = sample_weight)) +
  facet_grid(Treatment~surrogate) +
  theme(legend.position = "none") +
  ylab("Trial")

ggsave(
  filename = "distribution-titers.pdf",
  device = "pdf",
  path = figures_dir,
  width = double_width,
  height = double_height, 
  units = unit
)

# Pairwise scatterplots of the distribution of the neutralization titer and the
# corresponding titer adjusted to circulating variants.
df %>%
  filter(Delta_nAb)  %>% 
  ggplot(aes(
    x = pseudoneutid50,
    y = pseudoneutid50_adjusted,
    color = Treatment
  )) +
  geom_point(alpha = 0.25) +
  facet_wrap( ~ trial) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Neutralizing antibody titer against the vaccine strain") +
  ylab("Neutralizing antibody titer against circulating strains") +
  theme(legend.position = "bottom")

ggsave(
  filename = "scatterplot-titer-adjustment.pdf",
  device = "pdf",
  path = figures_dir,
  width = double_width,
  height = double_height, 
  units = unit
)

## Cumulative Incidence Curves --------------------------------------------

df %>%
  group_by(trial) %>%
  summarise(
    cum_inc = cuminc(
      formula = Surv(time_to_event, factor(event)) ~ Treatment,
      data = pick(everything())
    ) %>%
      list()
  )


figure = ggarrange(
  plotlist = df %>%
    # Cut the data off at 120 days.
    mutate(
      event = ifelse(time_to_event > 120, 0, event),
      time_to_event = ifelse(time_to_event > 120, 120, time_to_event)
    ) %>%
    group_by(trial) %>%
    summarise(
      cum_inc_plot = list(
        cuminc(
          formula = Surv(time_to_event, factor(event)) ~ Treatment,
          data = pick(everything())
        ) %>%
          ggcuminc(outcome = "1") +
          add_confidence_interval() +
          scale_ggsurvfit() +
          ggtitle(trial) +
          coord_cartesian(xlim = c(0, time_cumulative_incidence + 40)) +
          scale_x_continuous(breaks = c(0, 40, 80, 120)) +
          ylab("") +
          xlab("") +
          theme(
            plot.margin = unit(c(0.2, 0.5, 0, -0.5), 'lines'),
            plot.title = element_text(hjust = 1)
          )
      )
    ) %>%
    pull(cum_inc_plot),
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

annotate_figure(
  figure,
  left = text_grob("Cumulative Incidence", rot = 90),
  bottom = text_grob("Time (days)", vjust = -4.5)
)


ggsave(
  filename = "cumulative-incidence.pdf",
  device = "pdf",
  path = figures_dir,
  width = double_width,
  height = double_height, 
  units = unit
)



# Tables ------------------------------------------------------------------

## Case-Cohort Sampling and Events ----------------------------------------

# Compute the number of infections by trial and treatment arm.
table_infections = df %>%
  mutate(
    treatment = ifelse(
      treatment == 1,
      "Vaccine",
      "Placebo"
    )
  ) %>%
  group_by(trial, treatment) %>%
  summarise(
    n = n(),
    n_event = sum(event),
    n_not_event = n - n_event
  ) %>%
  pivot_wider(names_from = "treatment", values_from = c("n", "n_event", "n_not_event"), names_sep = " - ")

# Compute the number of subjects sampled in the case-cohort sampling by trial
# and infection status. We only select the vaccine recipients here because the
# titer is not measured for placebo recipients.
table_case_cohort_bAb = df %>%
  mutate(
    infected = ifelse(
      event,
      "Cases",
      "Non-Cases"
    )
  ) %>%
  group_by(trial, infected, Treatment) %>%
  summarise(
    n_Delta = sum(Delta_bAb)
  ) %>%
  pivot_wider(names_from = c("infected", "Treatment"), values_from = "n_Delta", names_sep = " - ")

table_case_cohort_nAb = df %>%
  mutate(
    infected = ifelse(
      event,
      "Cases",
      "Non-Cases"
    )
  ) %>%
  group_by(trial, infected, Treatment) %>%
  summarise(
    n_Delta = sum(Delta_nAb)
  ) %>%
  pivot_wider(names_from = c("infected", "Treatment"), values_from = "n_Delta", names_sep = " - ")


# Save the tables to a csv file.
write.csv(
  table_infections,
  file = paste0(tables_dir, "/infections-by-treatment.csv")
)

write.csv(
  table_case_cohort_bAb,
  file = paste0(tables_dir, "/case-cohort-bAb.csv")
)

write.csv(
  table_case_cohort_nAb,
  file = paste0(tables_dir, "/case-cohort-nAb.csv")
)
         