args = commandArgs(trailingOnly=TRUE)
# This R script takes only one argument: the variable used as clustering unit.
clustering_variable_chr = args[1]

# Load packages.
library(tidyverse)

# Data Preparation --------------------------------------------------------
landmark_times = 365.25 * c(1 / 12, 2 / 12, 3 / 12, 4 / 12, 5 / 12, 0.5, 1:2)

# Load intermediate results.
intermediate_object_location = paste0("R/data-fusion/intermediate-objects/",
                                      clustering_variable_chr,
                                      "/")
bias_inference_tbl = readRDS(file = paste0(intermediate_object_location, "bias_inference_tbl.rds")) %>%
  select(-np_boot_object)
bias_inference_tbl = bias_inference_tbl %>%
  mutate(
    lower_perc_ci = purrr::map_dbl(.x = perc_ci, .f = 1),
    upper_perc_ci = purrr::map_dbl(.x = perc_ci, .f = 2),
    est_surv_diff_se = purrr::map_dbl(
      .x = vcov_matrix,
      .f = function(x)
        sqrt(x[1, 1])
    )
  ) %>%
  mutate(
    landmark_month = 12 * landmark_time / 365.25,
    landmark_month_chr = factor(
      x = landmark_month,
      levels = 12 * landmark_times / 365.25,
      labels = paste0("Month ", 12 * landmark_times / 365.25)
    )
  )






# Bias Plots --------------------------------------------------------------

bias_plot = function(endpoint_plot,
                     include_year_four,
                     landmark_months) {
  landmark_months = unlist(landmark_months)
  if (include_year_four) {
    landmark_upper = 37
    times_upper = 365.25 * 4 + 1
  }  else {
    landmark_upper = 36
    times_upper = 365.25 * 4
  }
  p = bias_inference_tbl %>%
    filter(
      endpoint == endpoint_plot,
      landmark_month < landmark_upper,
      pred_time < times_upper,
      landmark_month %in% landmark_months
    )  %>%
    mutate(year = paste0("Year ", as.integer(pred_time/ 365.25))) %>%
    ggplot(
      aes(
        x = est_surv_diff,
        y = pred_surv_diff - est_surv_diff,
        color = clustering_variable,
        group = clustering_variable
      )
    ) +
    geom_point() +
    geom_errorbar(
      aes(
        ymin = pred_surv_diff - est_surv_diff - 1.96 * est_diff_se,
        ymax = pred_surv_diff - est_surv_diff + 1.96 * est_diff_se
      ),
      alpha = 0.50,
      linewidth = 0.35
    ) +
    geom_errorbarh(
      aes(
        xmin = est_surv_diff - 1.96 * est_surv_diff_se,
        xmax = est_surv_diff + 1.96 * est_surv_diff_se
      ),
      alpha = 0.50,
      linewidth = 0.35
    ) +
    geom_hline(aes(yintercept = 0), alpha = 0.5) +
    geom_vline(aes(xintercept = 0), alpha = 0.5) +
    scale_x_continuous(name = "Estimated Difference in Survival Probability", breaks = c(-0.05, 0, 0.05)) +
    scale_y_continuous(name = "Estimated Bias") +
    facet_grid(year ~ landmark_month_chr) +
    coord_cartesian(xlim = c(-0.05, 0.075)) +
    scale_color_discrete(name = "Clustering Unit") +
    theme(legend.position = "bottom", legend.box = "vertical")
  p = p + ggtitle("CVMB_DTH")
  p
  
  include_year_four_name = ifelse(include_year_four, "-4-years", "")
  filename = paste0(
    "figures/data-fusion/",
    clustering_variable_chr,
    "/bias-plot-",
    endpoint_plot,
    include_year_four_name,
    "-",
    landmark_months[1],
    "-",
    landmark_months[length(landmark_months)],
    ".pdf"
  )
  ggsave(
    filename = filename,
    device = "pdf",
    width = double_width,
    height = double_height,
    units = "cm",
    dpi = res
  )
  return(1)
}

# Enumerate all plotting options in a tibble. Call plotting function for each
# combinations of plotting options.
expand_grid(
  endpoint_plot = c("CVMB_DTH", "DTH"),
  include_year_four = c(FALSE),
  landmark_months = list(c(1, 2, 3, 4), c(5, 6, 12, 24))
) %>%
  group_by(endpoint_plot, include_year_four, landmark_months) %>%
  summarise(bias_plot(endpoint_plot, include_year_four, landmark_months))



# Scatter Plots with Confidence Ellipses ----------------------------------

scatter_plot_ellipse = function(endpoint_plot,
                                include_year_four,
                                landmark_months) {
  landmark_months = unlist(landmark_months)
  if (include_year_four) {
    landmark_upper = 37
    times_upper = 365.25 * 4 + 1
  }  else {
    landmark_upper = 36
    times_upper = 365.25 * 4
  }
  joint_confidence_regions_tbl = bias_inference_tbl %>%
    filter(
      endpoint == endpoint_plot,
      landmark_month < landmark_upper,
      pred_time < times_upper,
      landmark_month %in% landmark_months
    )  %>%
    mutate(year = paste0("Year ", as.integer(pred_time / 365.25))) %>%
    rowwise(clustering_variable, year, landmark_month_chr) %>%
    reframe(
      car::ellipse(
        center = c(est_surv_diff[1], pred_surv_diff),
        shape = vcov_matrix,
        radius = qchisq(p = 0.95, df = 2),
        draw = FALSE
      ) %>%
        as_tibble
    )
  colnames(joint_confidence_regions_tbl)[4:5] = c("x", "y")
  
  # Limits should depend on clustering variable and endpoint.
  # if (endpoint_plot == "CVMB_DTH" & clustering_variable_chr == "COU1A") {
  #   lims = c(-0.10, 0.15)
  # } else if (endpoint_plot == "CVMB_DTH" & clustering_variable_chr == "RGN1N") {
  #   lims = c(-0.05, 0.10)
  # } else if (endpoint_plot == "DTH" & clustering_variable_chr == "RGN1N") {
  #   lims = c(-0.05, 0.75)
  # }
  lims = c(-0.075, 0.15)
  
  p = bias_inference_tbl %>%
    filter(
      endpoint == endpoint_plot,
      landmark_month < landmark_upper,
      pred_time < times_upper,
      landmark_month %in% landmark_months
    )  %>%
    mutate(year = paste0("Year ", as.integer(pred_time / 365.25))) %>%
    ggplot(
      aes(
        x = est_surv_diff,
        y = pred_surv_diff,
        color = clustering_variable,
        group = clustering_variable
      )
    ) +
    geom_hline(aes(yintercept = 0), alpha = 0.5) +
    geom_vline(aes(xintercept = 0), alpha = 0.5) +
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dashed") +
    geom_path(aes(
      x = x,
      y = y,
      color = clustering_variable,
      group = clustering_variable
    ),
    data = joint_confidence_regions_tbl,
    linetype = "dashed") +
    geom_point() +
    scale_x_continuous(name = "Estimated Difference in Survival Probability") +
    scale_y_continuous(name = "Predicted Difference in Survival Probability") +
    facet_grid(year ~ landmark_month_chr) +
    coord_cartesian(xlim = lims, ylim = lims) +
    scale_color_discrete(name = "Clustering Unit") +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      axis.text.x = element_text(angle = 90)
    )
  p = p + ggtitle(endpoint_plot)
  p
  
  include_year_four_name = ifelse(include_year_four, "-4-years", "")
  filename = paste0(
    "figures/data-fusion/",
    clustering_variable_chr,
    "/scatter-plot-",
    endpoint_plot,
    include_year_four_name,
    "-",
    landmark_months[1],
    "-",
    landmark_months[length(landmark_months)],
    ".pdf"
  )
  ggsave(
    filename = filename,
    device = "pdf",
    width = double_width,
    height = double_height,
    units = "cm", 
    dpi = res
  )
  return(1)
}

# Enumerate all plotting options in a tibble. Call plotting function for each
# combinations of plotting options.
expand_grid(
  endpoint_plot = c("CVMB_DTH", "DTH"),
  include_year_four = c(FALSE),
  landmark_months = list(c(1, 2, 3, 4), c(5, 6, 12, 24))
) %>%
  group_by(endpoint_plot, include_year_four, landmark_months) %>%
  summarise(scatter_plot_ellipse(endpoint_plot, include_year_four, landmark_months))

