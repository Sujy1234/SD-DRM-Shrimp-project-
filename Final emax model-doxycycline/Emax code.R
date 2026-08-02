# Doxycycline Emax and EC50 estimation
# Hill coefficient fixed at 1 with residual-bootstrap confidence intervals

# Packages

library(dplyr)

# Data and model settings

dat <- read.csv("Antibiotic data with control.csv", check.names = FALSE) %>%
  filter(Antibiotic == "Doxycycline", Time_h <= 6)

Hill <- 1
n_boot <- 1000
set.seed(123)

# Functions

estimate_rates <- function(df) {
  df %>%
    group_by(Concentration_mg_L) %>%
    group_modify(~ {
      fit <- lm(logCFU ~ Time_h, data = .x)
      tibble(E_C = unname(coef(fit)["Time_h"]),
             R_squared = summary(fit)$r.squared)
    }) %>%
    ungroup()
}

fit_emax <- function(rates, E0) {
  d <- filter(rates, Concentration_mg_L > 0)
  concentration <- d$Concentration_mg_L
  effect <- E0 - d$E_C
  
  profile_rss <- function(log_EC50) {
    EC50 <- exp(log_EC50)
    fraction <- concentration^Hill / (EC50^Hill + concentration^Hill)
    Emax <- max(sum(fraction * effect) / sum(fraction^2), 0)
    sum((effect - Emax * fraction)^2)
  }
  
  fit <- optimize(
    profile_rss,
    interval = log(c(min(concentration) * 1e-6,
                     max(concentration) * 1e6)),
    tol = 1e-10
  )
  
  EC50 <- exp(fit$minimum)
  fraction <- concentration^Hill / (EC50^Hill + concentration^Hill)
  Emax <- max(sum(fraction * effect) / sum(fraction^2), 0)
  
  c(Emax = Emax, EC50 = EC50)
}

bootstrap_data <- function(df) {
  bind_rows(lapply(split(df, df$Concentration_mg_L), function(g) {
    fit <- lm(logCFU ~ Time_h, data = g)
    residuals_centred <- residuals(fit) - mean(residuals(fit))
    g$logCFU <- fitted(fit) +
      sample(residuals_centred, nrow(g), replace = TRUE)
    g
  }))
}

# Concentration-specific rates

rates <- estimate_rates(dat)

cat("\nConcentration-specific net growth or killing rates\n")
print(rates)

E0 <- rates %>%
  filter(Concentration_mg_L == 0) %>%
  pull(E_C)

# Emax model

parameters <- fit_emax(rates, E0)

Emax_per_day <- unname(parameters["Emax"]) * 24 * log(10)
EC50_ug_mL <- unname(parameters["EC50"])

fitted_rates <- rates %>%
  mutate(
    Predicted_E_C = E0 -
      unname(parameters["Emax"]) * Concentration_mg_L^Hill /
      (EC50_ug_mL^Hill + Concentration_mg_L^Hill)
  )

cat("\nObserved and fitted concentration-specific rates\n")
print(fitted_rates)

cat("\nPoint estimates for the SD-DRM\n")
cat("Emax =", round(Emax_per_day, 4), "per day\n")
cat("EC50 =", round(EC50_ug_mL, 4), "µg/mL\n")
cat("Hill coefficient =", Hill, "(fixed)\n")

# Residual bootstrap

boot <- replicate(n_boot, {
  tryCatch({
    boot_rates <- suppressWarnings(estimate_rates(bootstrap_data(dat)))
    boot_E0 <- boot_rates %>%
      filter(Concentration_mg_L == 0) %>%
      pull(E_C)
    
    fit_emax(boot_rates, boot_E0)
  }, error = function(e) c(Emax = NA, EC50 = NA))
})

boot <- as.data.frame(t(boot))
boot <- boot[complete.cases(boot), , drop = FALSE]

pd_boot <- data.frame(
  Bootstrap_run = seq_len(nrow(boot)),
  Emax_per_day = boot$Emax * 24 * log(10),
  EC50_ug_mL = boot$EC50
)

Emax_ci <- quantile(pd_boot$Emax_per_day,
                    c(0.025, 0.975),
                    names = FALSE)

EC50_ci <- quantile(pd_boot$EC50_ug_mL,
                    c(0.025, 0.975),
                    names = FALSE)

cat("\nBootstrap results\n")
cat("Successful runs =", nrow(pd_boot), "of", n_boot, "\n")
cat("Emax 95% CI =", round(Emax_ci[1], 4), "to",
    round(Emax_ci[2], 4), "per day\n")
cat("EC50 95% CI =", round(EC50_ci[1], 4), "to",
    round(EC50_ci[2], 4), "µg/mL\n")

# Final results

results <- data.frame(
  Antibiotic = "Doxycycline",
  Emax_per_day = Emax_per_day,
  Emax_LCL_per_day = Emax_ci[1],
  Emax_UCL_per_day = Emax_ci[2],
  EC50_ug_mL = EC50_ug_mL,
  EC50_LCL_ug_mL = EC50_ci[1],
  EC50_UCL_ug_mL = EC50_ci[2],
  Successful_bootstrap_runs = nrow(pd_boot)
)

cat("\nFinal doxycycline parameter estimates\n")
print(results)

# Save results

write.csv(
  results,
  "Doxycycline_SD_DRM_parameter_summary.csv",
  row.names = FALSE
)

write.csv(
  pd_boot,
  "Doxycycline_SD_DRM_empirical_bootstrap.csv",
  row.names = FALSE
)