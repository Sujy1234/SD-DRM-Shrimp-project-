# Risk assessment using Monte Carlo simulation

# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(triangle)

# Set seed
set.seed(123)

# Number of simulations
n_sim <- 100000


# Common inputs

# Uncertain prevalence of V. parahaemolyticus-positive shrimp
# Based on 10 positive and 75 negative samples
# Beta(1,1) prior gives Beta(11,76) posterior

P_pos <- rbeta(
  n = n_sim,
  shape1 = 11,
  shape2 = 76
)

# Indicator of whether each simulated shrimp sample is positive

I_pos <- rbinom(
  n = n_sim,
  size = 1,
  prob = P_pos
)

# Concentration among positive shrimp samples
# LogUniform(75,1100) MPN/g

C_pos <- exp(
  runif(
    n = n_sim,
    min = log(75),
    max = log(1100)
  )
)

# Unconditional retail concentration
# Zero for negative samples

C_retail <- I_pos * C_pos

# Daily shrimp consumption among shrimp consumers
# Arithmetic-scale mean = 36.38 g/day
# Estimated arithmetic-scale SD = 87.2 g/day
# Corresponding natural-log parameters:
# meanlog = 2.6396
# sdlog = 1.3816

M <- rlnorm(
  n = n_sim,
  meanlog = 2.6396,
  sdlog = 1.3816
)

# Raw shrimp mass handled
# Assumes that 13% of the initial mass is inedible or discarded

M_h <- M / (1 - 0.13)

# Pathogenic fraction

F_pathogenic <- runif(
  n = n_sim,
  min = 0.002,
  max = 0.02
)


# Scenario 1: Raw shrimp consumption

D_raw <- C_retail *
  M *
  F_pathogenic


# Scenario 2: Undercooked shrimp consumption

# Log10 reduction due to partial or insufficient cooking

R_under <- runif(
  n = n_sim,
  min = 0.7,
  max = 2.5
)

# Proportion remaining after cooking

Remaining_after_cooking <- 10^(-R_under)

# Pathogenic dose from undercooked shrimp

D_under <- C_retail *
  M *
  Remaining_after_cooking *
  F_pathogenic


# Scenario 3: Hand-mediated cross-contamination

# Shrimp-to-hand transfer
# X1 is modelled on the log10 percentage scale

X1 <- rlogis(
  n = n_sim,
  location = 1.59,
  scale = 0.14
)

T_shrimp_hand <- 10^X1 / 100

# Transfer fraction cannot exceed 1

T_shrimp_hand <- pmin(T_shrimp_hand, 1)

# Hand-to-food transfer
# X2 is modelled on the log10 percentage scale

X2 <- rtriangle(
  n = n_sim,
  a = -1.61,
  b = 1.32,
  c = 0.12
)

T_hand_food <- 10^X2 / 100

# Transfer fraction cannot exceed 1

T_hand_food <- pmin(T_hand_food, 1)

# Pathogenic dose transferred through contaminated hands

D_hand <- C_retail *
  M_h *
  T_shrimp_hand *
  T_hand_food *
  F_pathogenic


# Save full Monte Carlo simulation data

full_simulation_data <- data.frame(
  sim_id = seq_len(n_sim),
  P_pos = P_pos,
  I_pos = I_pos,
  C_pos = C_pos,
  C_retail = C_retail,
  M = M,
  M_h = M_h,
  F_pathogenic = F_pathogenic,
  R_under = R_under,
  Remaining_after_cooking = Remaining_after_cooking,
  X1_shrimp_hand = X1,
  T_shrimp_hand = T_shrimp_hand,
  X2_hand_food = X2,
  T_hand_food = T_hand_food,
  D_raw = D_raw,
  D_under = D_under,
  D_hand = D_hand
)

write.csv(
  full_simulation_data,
  file = "full_simulation_data_three_scenarios.csv",
  row.names = FALSE
)


# Prepare dose data

dose_data <- data.frame(
  Raw_shrimp = D_raw,
  Undercooked_shrimp = D_under,
  Hand_cross_contamination = D_hand
)

dose_long <- dose_data %>%
  pivot_longer(
    cols = everything(),
    names_to = "Scenario",
    values_to = "Dose"
  ) %>%
  mutate(
    Scenario = recode(
      Scenario,
      Raw_shrimp = "Raw shrimp consumption",
      Undercooked_shrimp = "Undercooked shrimp consumption",
      Hand_cross_contamination =
        "Hand-mediated cross-contamination"
    ),
    Scenario = factor(
      Scenario,
      levels = c(
        "Raw shrimp consumption",
        "Undercooked shrimp consumption",
        "Hand-mediated cross-contamination"
      )
    )
  )


# Summarise unconditional dose distributions

dose_summary <- dose_long %>%
  group_by(Scenario) %>%
  summarise(
    Probability_zero = mean(Dose == 0),
    Mean_dose = mean(Dose),
    Median_dose = median(Dose),
    P2.5 = quantile(Dose, 0.025),
    P25 = quantile(Dose, 0.25),
    P75 = quantile(Dose, 0.75),
    P97.5 = quantile(Dose, 0.975),
    Maximum = max(Dose),
    .groups = "drop"
  )

print(dose_summary)


# Summarise positive-dose simulations only

positive_dose_summary <- dose_long %>%
  filter(Dose > 0) %>%
  group_by(Scenario) %>%
  summarise(
    Mean_positive_dose = mean(Dose),
    Median_positive_dose = median(Dose),
    P2.5_positive = quantile(Dose, 0.025),
    P25_positive = quantile(Dose, 0.25),
    P75_positive = quantile(Dose, 0.75),
    P97.5_positive = quantile(Dose, 0.975),
    .groups = "drop"
  )

print(positive_dose_summary)


# Plot positive dose distributions

# Zero doses cannot be displayed on a log10 scale.
# The histogram therefore displays positive doses only.

dose_plot_data <- dose_long %>%
  filter(Dose > 0) %>%
  mutate(
    Log10_dose = log10(Dose)
  )

dose_histogram <- ggplot(
  data = dose_plot_data,
  aes(x = Log10_dose)
) +
  geom_histogram(
    binwidth = 0.10,
    boundary = 0,
    colour = "black",
    fill = "grey70"
  ) +
  facet_wrap(
    ~ Scenario,
    ncol = 1,
    scales = "free_y"
  ) +
  labs(
    x = expression(
      Log[10] *
        " pathogenic dose among positive exposures (MPN/day)"
    ),
    y = "Frequency"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(
      fill = "white",
      colour = "black"
    ),
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

print(dose_histogram)


# Save histogram

ggsave(
  filename = "dose_distribution_all_scenarios_updated.png",
  plot = dose_histogram,
  width = 7,
  height = 7,
  units = "in",
  dpi = 600
)
