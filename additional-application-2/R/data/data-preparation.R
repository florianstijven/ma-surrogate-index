library(tidyverse)
set.seed(1)
# Data Preparation --------------------------------------------------------

# Read in synthetic data from COVID-19 vaccine application.
synthetic_data_original = read.csv("data/processed_data_synthetic.csv")

# Add high-dimensional surrogate marker --------------------------------------

# Number of surrogate markers to add.
n_surrogate_markers = 250
# Number of surrogate markers with non-zero association with the clinical outcome.
n_non_zero_surrogate_markers = 10
# Effect size for active surrogate markers.
delta = 0.2

# Generate surrogate markers from standard normal distributions and add them to
# the data set. For the active markers, we add delta if the patient did not have
# the event.

# Generate active surrogate markers.
for (k in 1:n_non_zero_surrogate_markers) {
  surrogate_marker_name = paste0("surrogate_marker_", k)
  synthetic_data_original[[surrogate_marker_name]] = rnorm(nrow(synthetic_data_original)) - ifelse(synthetic_data_original$event == 0, delta, 0)
}
# Generate inactive surrogate markers.
for (k in (n_non_zero_surrogate_markers + 1):n_surrogate_markers) {
  surrogate_marker_name = paste0("surrogate_marker_", k)
  synthetic_data_original[[surrogate_marker_name]] = rnorm(nrow(synthetic_data_original))
}


# Save new data set -----------------------------------------------------------------

saveRDS(synthetic_data_original, "additional-application-2/data/processed_data_synthetic_high_d.rds")
