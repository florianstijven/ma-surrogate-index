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

# Join the clinical data with the GISAID data.
data_joined = df %>%
  rename(SUBJID = Ptid) %>%
  left_join(gisaid_tbl %>% select(-protocol), by = "SUBJID")

# Add a key to the joined data set to link it with the abrogation data.
data_joined = data_joined %>%
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

# Compute subject-specific adjustment factor.
adjustment_factor_f = function(prop.epsilon,
                               prop.gamma ,
                               prop.zeta,
                               prop.alpha,
                               prop.delta,
                               prop.lambda,
                               prop.iota,
                               prop.beta,
                               prop.mu,
                               prop.omicron ,
                               prop.BA.1,
                               prop.BA.2,
                               prop.BA.4.5,
                               abrogation_key) {
  # Extract the abrogation factors for the correct trial. 
  abrogation_coefs_vec = abrogation_tbl[abrogation_tbl$trial == abrogation_key, ] %>% as.double()
  
  # Put the circulating variants into a single vector.
  prop_vec = c(
    prop.epsilon,
    prop.gamma ,
    prop.zeta,
    prop.alpha,
    prop.delta,
    prop.lambda,
    prop.iota,
    prop.beta,
    prop.mu,
    prop.omicron,
    prop.BA.1,
    prop.BA.2,
    prop.BA.4.5
  )
  
  # Compute the general abrogation coefficient.
  sum(prop_vec * (1 / abrogation_coefs_vec), na.rm = TRUE)
}
data_joined = data_joined %>%
  mutate(
    adjustment_factor = purrr::pmap_dbl(
      .l = list(
        prop.epsilon = prop.epsilon,
        prop.gamma = prop.gamma,
        prop.zeta = prop.zeta,
        prop.alpha = prop.alpha,
        prop.delta = prop.delta,
        prop.lambda = prop.lambda,
        prop.iota = prop.iota,
        prop.beta = prop.beta,
        prop.mu = prop.mu,
        prop.omicron = prop.omicron,
        prop.BA.1 = prop.ba.1,
        prop.BA.2 = prop.ba.2,
        prop.BA.4.5 = prop.ba.4.5,
        abrogation_key = abrogation_key
      ),
      .f = adjustment_factor_f
    )
  )

# It would be way more efficient to just join the data sets together and have
# the abrogation coefficients in columns for each subject. 






# Save processed data set.
write.csv(df, "data/processed_data.csv")
