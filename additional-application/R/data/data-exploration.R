# Load packages.
library(tidyverse)

# Data Preparation --------------------------------------------------------

full_data_tbl_wide = readRDS("R/data/intermediate-objects/full_data_tbl_wide.rds")
full_data_tbl_long = readRDS("R/data/intermediate-objects/full_data_tbl_long.rds")

# Exploration -------------------------------------------------------------

## Description of the Clusters --------------------------------------------

# Distribution of the cluster sizes. 
full_data_tbl_wide %>%
  group_by(CTR1N) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = n)) +
  geom_histogram(fill = "gray", color = "black") 
ggsave(
  filename = paste0("figures/data-exploration/", "histogram-cluster-sizes", ".pdf"),
  device = "pdf",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)

# Distribution of the number of events in each cluster, stratified by the event
# type.
full_data_tbl_wide %>%
  pivot_longer(
    cols = c(C4DTH, C4CVDTH, C4MI, C4STRK, C4USA, C4CORON, C4CVMB, C4CVMM, C4CVMB_DTH),
    names_to = "Endpoint",
    values_to = "censored"
  ) %>%
  group_by(Endpoint, CTR1N) %>%
  summarize(no_of_events = sum(censored == 0, na.rm = TRUE)) %>%
  ggplot(aes(x = no_of_events)) + 
  geom_bar(color = "black", fill = "gray") +
  facet_wrap(~Endpoint, ncol = 3) +
  xlab("Number of Events")
ggsave(
  filename = paste0("figures/data-exploration/", "histogram-number-of-events", ".pdf"),
  device = "pdf",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)

# Plot the number of of events in each country. 
full_data_tbl_wide %>%
  pivot_longer(
    cols = c(C4DTH, C4CVDTH, C4MI, C4STRK, C4USA, C4CORON, C4CVMB, C4CVMM),
    names_to = "Endpoint",
    values_to = "censored"
  ) %>%
  group_by(Endpoint, COU1A) %>%
  summarize(no_of_events = sum(censored == 0, na.rm = TRUE)) %>%
  na.omit() %>%
  ggplot(aes(x = COU1A, y = no_of_events, colour = COU1A)) + 
  geom_point() +
  scale_y_continuous(name = "Number of Events", transform = "log10") +
  facet_wrap(~Endpoint, ncol = 3) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Number of Events")
ggsave(
  filename = paste0("figures/data-exploration/", "number-of-events-COU1A", ".pdf"),
  device = "pdf",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)

# Also save the above numbers in a table. 
sink(file = "tables/data-exploration/number-of-events-COU1A.txt")
cat("Number of events by endpoint and country:\n\n")
full_data_tbl_wide %>%
  pivot_longer(
    cols = c(C4DTH, C4CVDTH, C4MI, C4STRK, C4USA, C4CORON, C4CVMB, C4CVMM, C4CVMB_DTH),
    names_to = "Endpoint",
    values_to = "censored"
  ) %>%
  group_by(Endpoint, COU1A) %>%
  summarize(no_of_events = sum(censored == 0, na.rm = TRUE)) %>%
  na.omit() %>%
  pivot_wider(
    names_from = "Endpoint",
    values_from = "no_of_events"
  )
sink()

# Number of clusters with five or more events for each event type.
sink(file = "tables/data-exploration/no-clusters-five-or-more-events.txt")
cat("Number of cluster with five or more events of a specific type:\n\n")
full_data_tbl_wide %>%
  pivot_longer(
    cols = c(C4DTH, C4CVDTH, C4MI, C4STRK, C4USA, C4CORON, C4CVMB, C4CVMM),
    names_to = "Endpoint",
    values_to = "censored"
  ) %>%
  group_by(Endpoint, CTR1N) %>%
  summarize(no_of_events = sum(censored == 0, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Endpoint) %>%
  summarize(no_clusters_5_events = sum(no_of_events >= 5)) %>%
  pivot_wider(names_from = c(Endpoint), values_from = c(no_clusters_5_events))
sink()

## Time-to-Event Outcomes ---------------------------------------------------

event_time_names = c("T2DTH", "T2CVDTH", "T2MI", "T2STRK", "T2USA", "T2CORON", "T2CVMB", "T2CVMM", "T2CVMB_DTH")
censoring_indicator_names = c("C4DTH", "C4CVDTH", "C4MI", "C4STRK", "C4USA", "C4CORON", "C4CVMB", "C4CVMM", "C4CVMB_DTH")

library(GGally)
library(survival)
library(gridExtra)
# Function for the KM estimate.
km_est = function(censoring_indicator, event_time) {
  survfit(Surv(event_time, event_indicator) ~ 1, data = full_data_tbl_wide %>%
            mutate(event_indicator = 1 - .data[[censoring_indicator]],
                   event_time = .data[[event_time]]))
}

km_plot_list = purrr::map2(
  .x = censoring_indicator_names,
  .y = event_time_names,
  .f = function(.x, .y) {
    km_fit = km_est(.x, .y)
    ggsurv(km_fit, plot.cens = FALSE, CI = FALSE) +
      xlab("Time (days)") +
      coord_cartesian(xlim = c(0, 1600), ylim = c(0.5, 1)) +
      ggtitle(.y)
  }
)
g = do.call(arrangeGrob, km_plot_list)
ggsave(
  filename = paste0("figures/data-exploration/", "km-estimates", ".pdf"),
  plot = g,
  device = "pdf",
  width = double_width,
  height = double_height,
  units = "cm",
  dpi = res
)

## Surrogate Endpoints ------------------------------------------------------

# Distribution of blood pressure measurements and sitting pulse as a function of
# time.
full_data_tbl_long %>%
  select(c(SBPAVG3N, DBPAVG3N, months)) %>%
  pivot_longer(
    cols = c(SBPAVG3N, DBPAVG3N),
    names_to = "Type",
    values_to = "Blood Pressure"
  ) %>%
  mutate(Type = case_match(Type, "DBPAVG3N" ~ "Diastolic", "SBPAVG3N" ~ "Systolic")) %>%
  ggplot(aes(x = months, y = `Blood Pressure`, color = Type)) +
  geom_jitter(height = 0,
              width = 0.1,
              alpha = 0.01) +
  geom_smooth(se = FALSE) +
  xlab("Time (months)") +
  scale_color_brewer(
    type = "div",
    palette = "Dark2",
    breaks = c("Systolic", "Diastolic"),
    name = NULL
  ) +
  theme(legend.position = "bottom")
ggsave(
  filename = paste0("figures/data-exploration/", "evolution-bp", ".png"),
  device = "png",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)
full_data_tbl_long %>%
  ggplot(aes(x = months, y = STNPLS1N)) +
  geom_jitter(height = 0, width = 0.1, alpha = 0.01) +
  geom_smooth(se = FALSE) +
  ylab("Sitting Pulse") +
  xlab("Time (months)")
ggsave(
  filename = paste0("figures/data-exploration/", "evolution-stnpls", ".png"),
  device = "png",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)

# Evolution of the blood pressure in a random selection of patients. 
set.seed(1)
# Select 50 patients randomly. A black dot is added for the last measurement for
# each patient.
subset_tbl = full_data_tbl_long %>%
  filter(SID1A %in% sample(unique(SID1A), size = 50))
subset_tbl %>%
  pivot_longer(
    cols = c(SBPAVG3N, DBPAVG3N),
    names_to = "Type",
    values_to = "Blood Pressure"
  ) %>%
  mutate(Type = case_match(Type, "DBPAVG3N" ~ "Diastolic", "SBPAVG3N" ~ "Systolic")) %>%
  ggplot(aes(
    x = months,
    y = `Blood Pressure`,
    color = Type,
    group = SID1A
  )) +
  geom_line(aes(group = paste(SID1A, Type)), alpha = 0.25) +
  geom_point(
    data = subset_tbl %>% pivot_longer(
      cols = c(SBPAVG3N, DBPAVG3N),
      names_to = "Type",
      values_to = "Blood Pressure"
    ) %>%
      mutate(Type = case_match(
        Type, "DBPAVG3N" ~ "Diastolic", "SBPAVG3N" ~ "Systolic"
      )) %>%
      group_by(SID1A, Type) %>%
      slice_max(order_by = months)
  ) +
  xlab("Time (months)") +
  scale_color_brewer(
    type = "div",
    palette = "Dark2",
    breaks = c("Systolic", "Diastolic"),
    name = NULL
  ) +
  theme(legend.position = "bottom")
ggsave(
  filename = paste0("figures/data-exploration/", "individual-evolution-bp", ".pdf"),
  device = "pdf",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)
subset_tbl %>%
  ggplot(aes(x = months, y = STNPLS1N, group = SID1A)) +
  geom_line(alpha = 0.25) +
  geom_point(data = subset_tbl %>%
               group_by(SID1A) %>%
               slice_max(order_by = months)) +
  xlab("Time (months)") +
  scale_color_brewer(type = "div", palette = "Dark2", breaks = c("Systolic", "Diastolic"))
ggsave(
  filename = paste0("figures/data-exploration/", "individual-evolution-stnpls", ".pdf"),
  device = "pdf",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)

# Proportion of missing values for blood pressure and sitting pulse as a
# function of time.

# First, we determine the possible measurement occasions.
possible_timings = full_data_tbl_long %>%
  filter(months >= 0) %>%
  pull(months) %>%
  as.integer() %>%
  unique()
completed_full_data_tbl_long = full_data_tbl_long %>%
  pivot_longer(
    cols = c("SBPAVG3N", "DBPAVG3N", "STNPLS1N"),
    names_to = "surrogate",
    values_to = "value"
  ) %>%
  select(SID1A, surrogate, value, months) %>%
  mutate(months = as.integer(months)) %>%
  group_by(SID1A, surrogate) %>%
  reframe(pick(everything()) %>% right_join(
    tibble(months = possible_timings),
    relationship = "many-to-one",
    by = "months"
  ))
  
completed_full_data_tbl_long %>%
  group_by(months, surrogate) %>%
  summarise(prop_missing = mean(!is.na(value))) %>%
  ggplot(aes(x = months, y = prop_missing)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(name = "Time (Months)", breaks = c(0, 12, 24, 36, 48)) +
  scale_y_continuous(name = "Proportion Missing") +
  facet_grid(. ~ surrogate)
ggsave(
  filename = paste0("figures/data-exploration/", "proportion-missing", ".pdf"),
  device = "pdf",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)
