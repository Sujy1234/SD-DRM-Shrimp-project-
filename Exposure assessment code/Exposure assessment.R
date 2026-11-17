# Exposure assessment

library(mc2d)
library(ggplot2)
library(patchwork)

set.seed(123)

# Simulation settings
n_inner <- 1000
n_outer <- 500

ndvar(n_inner)
ndunc(n_outer)

# Input values
B_min <- 75
B_max <- 1100

M_shape <- 0.8760928
M_scale <- 43.276368

F_min <- 0.0018
F_max <- 0.0233

# Input variability
B_retail_v <- mcstoc(runif, type = "V", min = B_min, max = B_max)
M_v <- mcstoc(rgamma, type = "V", shape = M_shape, scale = M_scale)
F_pathogenic_v <- mcstoc(runif, type = "V", min = F_min, max = F_max)

# Daily pathogenic dose
D_raw_v <- B_retail_v * M_v * F_pathogenic_v

# Model summary
exposure_model <- mc(B_retail_v, M_v, F_pathogenic_v, D_raw_v)

print(exposure_model)
print(summary(exposure_model, probs = c(0.025, 0.50, 0.975)))

# Extract variability simulations
simulation_data <- data.frame(
  inner_id = seq_len(n_inner),
  B_retail = as.vector(unmc(B_retail_v, drop = TRUE)),
  M = as.vector(unmc(M_v, drop = TRUE)),
  F_pathogenic = as.vector(unmc(F_pathogenic_v, drop = TRUE)),
  D_raw = as.vector(unmc(D_raw_v, drop = TRUE))
)

# Add dose on the log10 scale
simulation_data$Log10_D_raw <- log10(simulation_data$D_raw)

# Save simulation values
write.csv(
  simulation_data,
  "raw_shrimp_mc2d_simulation_data.csv",
  row.names = FALSE
)

# Descriptive summary
summarise_parameter <- function(x) {
  percentiles <- quantile(
    x,
    probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  data.frame(
    Mean = mean(x, na.rm = TRUE),
    P2.5 = percentiles[1],
    P25 = percentiles[2],
    Median = percentiles[3],
    P75 = percentiles[4],
    P97.5 = percentiles[5],
    Minimum = min(x, na.rm = TRUE),
    Maximum = max(x, na.rm = TRUE)
  )
}

parameter_labels <- c(
  B_retail = "Retail concentration (MPN/g)",
  M = "Serving mass (g/shrimp-consumption day)",
  F_pathogenic = "Pathogenic fraction (proportion)",
  D_raw = "Pathogenic dose (MPN/shrimp-consumption day)",
  Log10_D_raw = "Log10 pathogenic dose [log10(MPN/shrimp-consumption day)]"
)

parameter_summary <- do.call(
  rbind,
  lapply(names(parameter_labels), function(parameter) {
    result <- summarise_parameter(simulation_data[[parameter]])
    result$Parameter <- parameter_labels[[parameter]]
    result
  })
)

parameter_summary <- parameter_summary[
  ,
  c("Parameter", "Mean", "P2.5", "P25", "Median",
    "P75", "P97.5", "Minimum", "Maximum")
]

rownames(parameter_summary) <- NULL

print(parameter_summary)

write.csv(
  parameter_summary,
  "raw_shrimp_input_parameter_summary.csv",
  row.names = FALSE
)

# Plotting data
plot_data <- data.frame(Log10_dose = simulation_data$Log10_D_raw)

q_dose <- quantile(
  plot_data$Log10_dose,
  probs = c(0.025, 0.50, 0.975),
  na.rm = TRUE
)

p2.5_log10_dose <- unname(q_dose[1])
median_log10_dose <- unname(q_dose[2])
p97.5_log10_dose <- unname(q_dose[3])

plot_style <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.tag = element_text(face = "bold", size = 12)
  )

# Histogram
dose_plot <- ggplot(plot_data, aes(x = Log10_dose)) +
  annotate(
    "rect",
    xmin = p2.5_log10_dose,
    xmax = p97.5_log10_dose,
    ymin = -Inf, ymax = Inf,
    fill = "#BBD1E3",
    alpha = 0.9
  ) +
  geom_histogram(
    binwidth = 0.10,
    boundary = 0,
    colour = "black",
    fill = "#7898B0",
    alpha = 0.65
  ) +
  geom_vline(
    xintercept = c(p2.5_log10_dose, p97.5_log10_dose),
    linetype = "dashed",
    colour = "#2C5D8A",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = median_log10_dose,
    colour = "black",
    linewidth = 0.9
  ) +
  labs(
    x = expression(Log[10] * " pathogenic dose (MPN per shrimp-consumption day)"),
    y = "Frequency"
  ) +
  plot_style

# Cumulative distribution function
cdf_plot <- ggplot(plot_data, aes(x = Log10_dose)) +
  stat_ecdf(
    geom = "step",
    linewidth = 0.9,
    colour = "black"
  ) +
  annotate(
    "segment",
    x = -Inf, xend = median_log10_dose,
    y = 0.5, yend = 0.5,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = median_log10_dose, xend = median_log10_dose,
    y = 0, yend = 0.5,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  annotate(
    "label",
    x = median_log10_dose + 0.10,
    y = 0.47,
    label = paste0("Median = ", round(median_log10_dose, 3)),
    hjust = 0,
    size = 3.5
  ) +
  scale_y_continuous(
    labels = scales::percent,
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = expression(Log[10] * " pathogenic dose (MPN per shrimp-consumption day)"),
    y = "Cumulative probability"
  ) +
  plot_style
# Combine figures
combined_dose_plot <- dose_plot + cdf_plot +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

print(combined_dose_plot)

ggsave(
  "dose_variability_histogram_cdf_mc2d.png",
  combined_dose_plot,
  width = 12,
  height = 5,
  dpi = 600,
  bg = "white"
)