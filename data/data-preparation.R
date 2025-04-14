# Load required packages
library(tidyr)
library(dplyr)

# Cleaning up the CrossProtocolData.csv ---------------------------------------

# Set the universal lower limits of detection for binding and neutralizing
# antibody titers.
llod_spike = log(10.84, base = 10)
llod_neut = log(2.61, base = 10)
# Observed values below the LLOD will be truncated the the LLOD divided by 2.
llod_spike_truncated = log(10.84 / 2, base = 10)
llod_neut_truncated = log(2.61 / 2, base = 10)

# Load the data set. Wstratum refers to the strata that determine the
# probability of being sampled for measuring S. 
df = read.csv("data/CrossProtocolData.csv") %>% 
  mutate(Wstratum = as.factor(Wstratum)) %>%
  droplevels() 

# Set all values of placebo to the LLOD divided by 2.
df <- df %>% mutate(bindSpike = ifelse(treatment == 0, llod_spike_truncated, bindSpike))
df <- df %>% mutate(pseudoneutid50 = ifelse(treatment == 0, llod_neut_truncated, pseudoneutid50))

# Patients with measured titers, who have a titer below the LLOD, will have
# their measurement truncated to the LLOD divided by 2.
df <- df %>% mutate(bindSpike = ifelse(Delta_bAb == T &
                                         bindSpike < llod_spike, llod_spike_truncated, bindSpike))
df <- df %>% mutate(
  pseudoneutid50 = ifelse(
    Delta_nAb == T & pseudoneutid50 < llod_neut,
    llod_neut_truncated,
    pseudoneutid50
  )
)

# Drop variables which are currently not used in any analyses.
df = df %>%
  select(
    -X,
    -CalendarDateEnrollment,
    -phase1_ind,
    -WhiteNonHispanic,
    -Country
  )

# For the non-JnJ trials, BMI should be categorized into <18.5, [18.5, 25), [25,
# 30), >= 30.
df = df %>%
  mutate(
    BMI_underweight = as.integer(ifelse(
      stringr::str_detect(protocol, "p3005"), BMI == 1, BMI < 18.5
    )),
    BMI_normal = as.integer(ifelse(
      stringr::str_detect(protocol, "p3005"),
      BMI == 2,
      BMI >=  18.5 &
        BMI < 25
    )),
    BMI_overweight = as.integer(ifelse(
      stringr::str_detect(protocol, "p3005"), BMI == 3, BMI >= 25 &
        BMI < 30
    )),
    BMI_obese = as.integer(ifelse(
      stringr::str_detect(protocol, "p3005"), BMI == 4, BMI >= 30
    )),
    BMI_underweight_normal = as.integer(ifelse(
      stringr::str_detect(protocol, "p3005"), BMI <= 2, BMI < 25
    ))
  )


# Adding the Circulating Variant Information -------------------------------

# Load the GISAID data set. 
gisaid_location = "S:/p300x_CrossProtocol/analysis/immunocorrelates/PartA_Blinded_Phase/data/covpn_meta_with_gisaid_usg_output_v6.csv"
gisaid_tbl = read.csv(gisaid_location)

# Load the table with abrogation coefficients.
abrogation_tbl = read.csv("data/GMT_ratios.csv")
# We set the missing values to -1000. This is ok because if the abrogation
# coefficient is missing for a trial, the corresponding variant should not be
# circulating. Setting it to a negative number ensures that we can just multiply
# with a zero proportion. If somehow, it would be multiplied with a non-zero
# proportion, we would notice it because the result may not be negative. 
abrogation_tbl = abrogation_tbl %>%
  mutate(across(reference:BA.4.5, ~ ifelse(is.na(.x), -1e10, .x)))

# Add a key to the joined data set to link it with the abrogation data.
data_joined = df %>%
  mutate(
    abrogation_key = trial,
    abrogation_key = ifelse(
      abrogation_key == "Sanofi",
      ifelse(
        treatment == 1,
        "Sanofi 1 Non-Naive Vaccine",
        "Sanofi 1 Non-Naive Placebo"
      ),
      abrogation_key
    )
  )

# Add the abrogation coefficients to data_joined.
data_joined = data_joined %>%
  left_join(abrogation_tbl, by = join_by(abrogation_key == trial))

# Join the clinical data with the GISAID data.
data_joined = data_joined %>%
  left_join(gisaid_tbl %>% select(-protocol, -SUBJID),
            by = "USUBJID",
            relationship = "one-to-one") %>%
  # We drop the information for omicron because we have the same information for
  # the omicron subtypes BA.1, BA.2, and BA.4.5.
  select(-prop.omicron, -omicron)


# For each trial, we set the proportion circulating to zero for every variant
# which was not known to be relevant for the given study.
data_joined = data_joined %>%
  mutate(total_proportion_row = rowSums(select(., prop.reference:prop.ba.4.5))) %>%
  mutate(
    prop.epsilon  = ifelse(epsilon < 0, 0, prop.epsilon),
    prop.gamma  = ifelse(gamma < 0, 0, prop.gamma),
    prop.zeta  = ifelse(zeta < 0, 0, prop.zeta),
    prop.alpha  = ifelse(alpha < 0, 0, prop.alpha),
    prop.delta  = ifelse(delta < 0, 0, prop.delta),
    prop.lambda  = ifelse(lambda < 0, 0, prop.lambda),
    prop.iota  = ifelse(iota < 0, 0, prop.iota),
    prop.beta  = ifelse(beta < 0, 0, prop.beta),
    prop.mu  = ifelse(mu < 0, 0, prop.mu),
    prop.ba.1  = ifelse(BA.1 < 0, 0, prop.ba.1),
    prop.ba.2  = ifelse(BA.2 < 0, 0, prop.ba.2),
    prop.ba.4.5  = ifelse(BA.4.5 < 0, 0, prop.ba.4.5)
  ) %>%
  # Make sure that the adjustment proportions sum to one.
  mutate(total_proportion_row_new = rowSums(select(., prop.reference:prop.ba.4.5))) %>%
  mutate(across(prop.reference:prop.ba.4.5, ~ .x / total_proportion_row_new))




# Compute the subject-specific adjustment factor.
data_joined = data_joined %>%
  mutate(
    abrogation_coefficient = prop.reference * reference +
      prop.epsilon * epsilon +
      prop.gamma * gamma +
      prop.zeta * zeta +
      prop.alpha * alpha +
      prop.delta * delta +
      prop.lambda * lambda +
      prop.iota * iota +
      prop.beta * beta +
      prop.mu * mu +
      prop.ba.1 * BA.1 +
      prop.ba.2 * BA.2 +
      prop.ba.4.5 * BA.4.5
  )


# Save processed data set.
write.csv(df, "data/processed_data.csv")
