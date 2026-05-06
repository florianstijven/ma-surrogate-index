library(tidyverse)

# Data Preparation --------------------------------------------------------

# Read in SAS-data sets.
baseline_tbl = haven::read_sas("data/baseline.sas7bdat")
death_tbl = haven::read_sas("data/death.sas7bdat")
demog_tbl = haven::read_sas("data/demog-w-countries.sas7bdat") %>%
  select(-TRTREG1C) # Drop treatment indicator because this variables is already in baseline_tbl
endpoint_tbl = haven::read_sas("data/endpoint.sas7bdat")
bp_tbl = haven::read_sas("data/bp.sas7bdat")

# The bp data set also contains unscheduled measurements. We will ignore these
# for now. We also drop the date of the measurements and the report number; we
# don't need these variables.
bp_tbl = bp_tbl %>%
  filter(VISNAM1A != "Unscheduled Visit") %>%
  select(-VIS1N, -VIS1D)

# The bp data set contains the measurement at day one twice. We first check
# whether the values are identical. If so, we can delete one of the measurements
# at day 1.

# The patient with id 0902_00011 contains has two rows for Day 1: one with
# missing values and one with observed values. We delete the former row.
bp_tbl = bp_tbl %>%
  filter(!(SID1A == "0902_00011" & is.na(SBPAVG3N)))
# The within-patient variance of the duplicated measurements at Day 1 is zero.
# We can the safely drop one row per patient.
bp_tbl %>%
  filter(VISNAM1A == "Day 1") %>%
  group_by(SID1A) %>%
  summarize(var_sbp = ifelse(n() > 1L, var(SBPAVG3N, na.rm = TRUE), 0),
            var_dbp = ifelse(n() > 1L, var(SBPAVG3N, na.rm = TRUE), 0),
            var_pls = ifelse(n() > 1L, var(SBPAVG3N, na.rm = TRUE), 0)) %>%
  ungroup() %>%
  summarize(identical_sbp = all(var_sbp == 0),
            identical_dbp = all(var_dbp == 0),
            identical_pls = all(var_pls == 0))

bp_tbl = bp_tbl %>%
  group_by(SID1A, VISNAM1A) %>%
  slice_head(n = 1) %>%
  ungroup()

# There are now no duplicates left in the bp data.
bp_tbl %>%
  group_by(SID1A, VISNAM1A) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  summarize(no_duplicates = all(n == 1L))

# There are 111 more patients in the baseline data set than in the demog data
# set. These are all patients with no value for the treatment indicator.
baseline_tbl %>%
  filter(!(SID1A %in% demog_tbl$SID1A)) 
baseline_tbl %>%
  filter(!(SID1A %in% demog_tbl$SID1A)) %>%
  summarize(all(TRTREG1C == ""))
# The above 111 patients are dropped from the data set.
baseline_tbl = baseline_tbl %>%
  filter(TRTREG1C != "")

# The blood pressure variables in the demog data set correspond to the endpoints
# at Day 1 in the endpoints data set. However, these measurements differ for 14
# patients.
full_join(
  x = demog_tbl,
  y = bp_tbl %>%
    filter(VISNAM1A == "Day 1"),
  by = "SID1A" 
) %>%
  summarize(sum(SBPAVG3N.x != SBPAVG3N.y, na.rm = TRUE),
            sum(DBPAVG3N.x != DBPAVG3N.y, na.rm = TRUE),
            sum(STNPLS1N.x != STNPLS1N.y, na.rm = TRUE))
# We only keep the variables from the bp data set.
demog_tbl = demog_tbl %>%
  select(-SBPAVG3N, -DBPAVG3N, -STNPLS1N)


# Join all data sets into a single wide data set.
full_data_tbl_wide = full_join(x = baseline_tbl, y = demog_tbl, by = "SID1A") %>%
  full_join(y = death_tbl, by = "SID1A") %>%
  full_join(y = endpoint_tbl, by = "SID1A") %>%
  full_join(
    y = bp_tbl %>% pivot_wider(
      id_cols = "SID1A",
      names_from = "VISNAM1A",
      values_from = c("SBPAVG3N", "DBPAVG3N", "STNPLS1N"),
      names_sep = "_",
      values_fill = NA
    ),
    by = "SID1A"
  )

# Vector of the censoring indicator names for all time-to-event variables
censoring_indicator_names = c("C4DTH",
                              "C4CVDTH",
                              "C4MI",
                              "C4STRK",
                              "C4USA",
                              "C4CORON",
                              "C4CVMB",
                              "C4CVMM")
# Vector of time-to-event variable names for all time-to-event variables
event_time_names = c("T2DTH",
                     "T2CVDTH",
                     "T2MI",
                     "T2STRK",
                     "T2USA",
                     "T2CORON",
                     "T2CVMB",
                     "T2CVMM")

# Negative survival times are set to NA.
full_data_tbl_wide = full_data_tbl_wide %>%
  mutate(across(all_of(event_time_names), ~ ifelse(.x <= 0, NA, .x)))

# Join all data sets into a single long data set.
full_data_tbl_long = full_join(x = baseline_tbl, y = demog_tbl, by = "SID1A") %>%
  full_join(y = death_tbl, by = "SID1A") %>%
  full_join(y = endpoint_tbl, by = "SID1A") %>%
  full_join(
    y = bp_tbl,
    by = "SID1A"
  )
# Add number of weeks to each row of the long data set. 
full_data_tbl_long = full_data_tbl_long %>%
  mutate(
    months = case_match(
      VISNAM1A,
      "Day 1" ~ 1 / (365.25 / 12),
      "Month 1" ~ 1,
      "Month 2" ~ 2,
      "Month 3" ~ 3,
      "Month 4" ~ 4,
      "Month 6" ~ 6,
      "Year 1" ~ 12,
      "1 Year 6 Months" ~ 18,
      "2 Years" ~ 24,
      "2 Years 6 Months" ~ 30,
      "3 Years" ~ 36,
      "3 Years 6 Months" ~ 42,
      "4 Years" ~ 48,
      "5 Years" ~ 60,
      "Week -2" ~ -14 / (365.25 / 12), 
      .default = NA
    ),
    days = months * (365.25 / 12)
  )

# The original time-to-event data are competing risks data. This introduces
# additional difficulties. We will therefore also consider alternative composite
# endpoints that are not subject to competing risks. Specifically, we will
# consider time to death, and time to cardiovascular event or death.
full_data_tbl_wide = full_data_tbl_wide %>%
  mutate(T2CVMB_DTH = pmin(T2CVMM, T2DTH),
         C4CVMB_DTH = pmin(C4CVMM, C4DTH)) %>%
  mutate(RGN1N = factor(RGN1N, levels = 1:2, labels = c("US", "Nordic")),
         COU1A = as.factor(COU1A))

full_data_tbl_long = full_data_tbl_long %>%
  mutate(T2CVMB_DTH = pmin(T2CVMM, T2DTH),
         C4CVMB_DTH = pmin(C4CVMM, C4DTH)) %>%
  mutate(RGN1N = factor(RGN1N, levels = 1:2, labels = c("US", "Nordic")),
         COU1A = as.factor(COU1A))

# Save Results -----------------------------------------------------------------

saveRDS(
  full_data_tbl_wide,
  file = "R/data/intermediate-objects/full_data_tbl_wide.rds"
)

saveRDS(
  full_data_tbl_long,
  file = "R/data/intermediate-objects/full_data_tbl_long.rds"
)
