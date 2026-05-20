# Setup ------------------------------------------------------------------
# Load libraries
library(tidyverse)
library(scales)
library(splines)
library(WeightedROC)
library(mgcv)
library(future)
library(furrr)
library(survival)

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore")
} else {
  plan(multisession)
}

# Extract arguments for analysis.
args = commandArgs(trailingOnly = TRUE)

surr_indices_tbl_location = "additional-application/R/meta-analytic/intermediate-objects/landmark_predictions_tbl.rds"
# Load data set with IPD and estimated surrogate index for every subject.
ipd_surr_indices_tbl = readRDS(file = surr_indices_tbl_location)


## Prediction Model Performance -------------------------------------------
roc_tbl = ipd_surr_indices_tbl %>%
  filter(!is.na(surrogate_index) & !(is.na(event_status))) %>%
  group_by(method, surrogate, COU1A) %>%
  reframe(
    WeightedROC(
      guess = surrogate_index,
      label = event_status,
      weight = ipcw
    )
  )

# Make and save the plots for all combinations of `surrogate` and `weighting`.
roc_ggplots = roc_tbl %>%
  group_by(analysis_set, surrogate) %>%
  summarise(data = list(pick(everything()))) %>%
  ungroup() %>%
  mutate(ggplot_object = purrr::pmap(
    .l = list(
      data = data,
      analysis_set = analysis_set,
      surrogate = surrogate
    ),
    .f = function(data, analysis_set, surrogate) {
      surrogate_chr = switch(
        surrogate,
        bindSpike = "IgG Spike",
        pseudoneutid50 = "nAb ID50",
        pseudoneutid50_adjusted = "adjusted nAb ID50"
      )
      analysis_set_chr = switch(
        analysis_set,
        first_four = "first four trials",
        naive_only = "all trials, naive subjects only",
        mixed = "all trials"
      )
      subtitle = paste0("Surrogate index for ",
                        surrogate_chr,
                        " using ",
                        analysis_set_chr)
      data %>%
        mutate(
          weighting = ifelse(weighting == "normalized", "Normalized", "Unnormalized"),
          method = ifelse(method == "gam", "GAM logistic regression", ifelse(method == "cox", "Cox Model", "SuperLearner"))
        ) %>%
        ggplot(aes(color = method)) +
        geom_path(aes(FPR, TPR)) +
        facet_wrap(~ trial) +
        theme(legend.position = "bottom", legend.box = "vertical") +
        scale_color_discrete(name = "Method") +
        geom_abline(intercept = 0, slope = 1) +
        scale_color_discrete(name = "Estimated Surrogate Index") +
        ggtitle("ROC for estimated surrogate index") +
        labs(subtitle = subtitle) +
        theme(legend.spacing.y = unit(0, "cm"),
              legend.box.spacing = unit(0, "cm"))
    }
  ))

roc_ggplots %>%
  rowwise(analysis_set, surrogate) %>%
  summarise(
    ggsave(
      plot = ggplot_object,
      paste0("roc-", surrogate, "-", analysis_set, ".pdf"),
      path = figures_dir,
      height = double_height,
      width = double_width,
      device = "pdf",
      units = "cm"
    )
  )


roc_tbl %>%
  group_by(method, surrogate, weighting, trial, analysis_set, include_risk_score) %>%
  filter(method != "untransformed surrogate") %>%
  summarise(AUC = WeightedAUC(pick(c(
    TPR, FPR, FP, FN, threshold
  )))) %>%
  pivot_wider(names_from = "method", values_from = "AUC") %>%
  write.csv(paste0(tables_dir, "/auc-surrogate-indices.csv"))

rm("roc_tbl")