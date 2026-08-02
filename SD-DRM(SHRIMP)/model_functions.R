# SD-DRM model functions with annual infection-risk simulation

# Stable calculation of log(1 - exp(-a))

log1mexpm <- function(a) {
  ans <- numeric(length(a))
  ans[a <= log(2)] <- log(-expm1(-a[a <= log(2)]))
  ans[a > log(2)] <- log1p(-exp(-a[a > log(2)]))
  ans
}


# Recycle scalar inputs to the required length

recycle_input <- function(x, n, name) {
  if (!(length(x) %in% c(1, n))) {
    stop(paste0(name, " must have length 1 or ", n, "."))
  }
  rep_len(x, n)
}


# Fit susceptible dose-response parameters for one PD effect

fit_ab_effect <- function(mu, effect, tfs) {
  r_s <- -log1mexpm((mu + effect) * tfs)
  r_s <- pmin(pmax(r_s[is.finite(r_s)], 1e-12), 1 - 1e-12)
  
  fit <- fitdistrplus::fitdist(r_s, "beta", method = "mle")
  
  c(
    alpha_s = unname(fit$estimate["shape1"]),
    beta_s = unname(fit$estimate["shape2"])
  )
}


# Estimate susceptible parameters for one Emax-EC50 combination

get_ab_beta <- function(alpha, beta, C, Emax, EC50, tfs,
                        seed = 0, nsim = 10000) {
  
  if (length(C) != 1 || length(Emax) != 1 || length(EC50) != 1) {
    stop("C, Emax, and EC50 must each contain one value.")
  }
  
  if (C == 0) return(c(alpha_s = alpha, beta_s = beta))
  
  set.seed(seed)
  
  r <- stats::rbeta(nsim, alpha, beta)
  mu <- -log1mexpm(r) / tfs
  effect <- Emax * C / (EC50 + C)
  
  fit_ab_effect(mu, effect, tfs)
}


# Create lookup table for paired pharmacodynamic parameters

make_ab_lookup <- function(alpha, beta, Emax, EC50, C_max, tfs,
                           n_grid = 250, seed = 0, nsim = 10000) {
  
  if (length(Emax) != length(EC50)) {
    stop("Emax and EC50 must be paired vectors of equal length.")
  }
  
  if (any(!is.finite(Emax)) || any(!is.finite(EC50)) ||
      any(Emax < 0) || any(EC50 <= 0)) {
    stop("Invalid Emax or EC50 values.")
  }
  
  set.seed(seed)
  
  r <- stats::rbeta(nsim, alpha, beta)
  mu <- -log1mexpm(r) / tfs
  
  effect_max <- 1.001 * max(Emax * C_max / (EC50 + C_max))
  
  if (effect_max <= 0) {
    return(data.frame(effect = 0, alpha_s = alpha, beta_s = beta))
  }
  
  effect_grid <- c(
    0,
    exp(seq(
      log(effect_max * 1e-6),
      log(effect_max),
      length.out = n_grid - 1
    ))
  )
  
  pars <- t(vapply(
    effect_grid[-1],
    function(x) fit_ab_effect(mu, x, tfs),
    numeric(2)
  ))
  
  data.frame(
    effect = effect_grid,
    alpha_s = c(alpha, pars[, "alpha_s"]),
    beta_s = c(beta, pars[, "beta_s"])
  )
}


# Interpolate susceptible dose-response parameters

get_ab_from_lookup <- function(effect, ab_lookup) {
  
  upper_limit <- max(ab_lookup$effect)
  
  if (max(effect) > upper_limit * (1 + 1e-10)) {
    stop("The pharmacodynamic effect exceeds the lookup-table range.")
  }
  
  effect <- pmin(pmax(effect, 0), upper_limit)
  
  data.frame(
    alpha_s = exp(approx(
      ab_lookup$effect,
      log(ab_lookup$alpha_s),
      xout = effect,
      rule = 2
    )$y),
    
    beta_s = exp(approx(
      ab_lookup$effect,
      log(ab_lookup$beta_s),
      xout = effect,
      rule = 2
    )$y)
  )
}


# Calculate SD-DRM infection risk per consumption event

calc_sd_drm <- function(dose, fr, C, alpha, beta,
                        Emax, EC50, tfs, ab_lookup = NULL) {
  
  n <- max(
    length(dose), length(fr), length(C),
    length(Emax), length(EC50)
  )
  
  dose <- recycle_input(dose, n, "dose")
  fr <- recycle_input(fr, n, "fr")
  C <- recycle_input(C, n, "C")
  Emax <- recycle_input(Emax, n, "Emax")
  EC50 <- recycle_input(EC50, n, "EC50")
  
  if (any(!is.finite(dose)) || any(dose < 0)) {
    stop("dose must contain finite values greater than or equal to zero.")
  }
  
  if (any(!is.finite(fr)) || any(fr < 0 | fr > 1)) {
    stop("fr must contain values between 0 and 1.")
  }
  
  if (any(!is.finite(C)) || any(C < 0) ||
      any(!is.finite(Emax)) || any(Emax < 0) ||
      any(!is.finite(EC50)) || any(EC50 <= 0)) {
    stop("Invalid C, Emax, or EC50 values.")
  }
  
  effect <- ifelse(C == 0, 0, Emax * C / (EC50 + C))
  
  if (is.null(ab_lookup)) {
    
    if (any(effect != effect[1])) {
      stop("ab_lookup is required when the PD effect varies between rows.")
    }
    
    pars <- get_ab_beta(
      alpha, beta, C[1], Emax[1], EC50[1], tfs
    )
    
    alpha_s <- rep(unname(pars["alpha_s"]), n)
    beta_s <- rep(unname(pars["beta_s"]), n)
    
  } else {
    
    pars <- get_ab_from_lookup(effect, ab_lookup)
    alpha_s <- pars$alpha_s
    beta_s <- pars$beta_s
  }
  
  Ns <- dose * (1 - fr)
  Nr <- dose * fr
  
  p_ext_s <- (1 + Ns / beta_s)^(-alpha_s)
  p_ext_r <- (1 + Nr / beta)^(-alpha)
  
  risk_total <- 1 - p_ext_s * p_ext_r
  risk_less_treatable <- 1 - p_ext_r
  risk_more_treatable <- p_ext_r * (1 - p_ext_s)
  
  data.frame(
    risk_total = pmin(pmax(risk_total, 0), 1),
    risk_less_treatable = pmin(pmax(risk_less_treatable, 0), 1),
    risk_more_treatable = pmin(pmax(risk_more_treatable, 0), 1),
    PD_effect = effect,
    
    status = ifelse(
      risk_less_treatable >= risk_more_treatable,
      "Less likely treatable",
      "More likely treatable"
    ),
    
    stringsAsFactors = FALSE
  )
}


# Simulate annual infection risk

simulate_annual_risk <- function(
    daily_risk,
    n_person_years = 100000,
    consumption_mean = 0.0369,
    consumption_se = 0.0031,
    days_per_year = 365,
    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(daily_risk)) stop("daily_risk must be numeric.")
  
  daily_risk <- daily_risk[
    is.finite(daily_risk) & daily_risk >= 0 & daily_risk <= 1
  ]
  
  if (length(daily_risk) == 0) {
    stop("daily_risk contains no valid infection probabilities.")
  }
  
  if (length(n_person_years) != 1 || n_person_years < 1 ||
      n_person_years != as.integer(n_person_years)) {
    stop("n_person_years must be a positive integer.")
  }
  
  if (length(days_per_year) != 1 || days_per_year < 1 ||
      days_per_year != as.integer(days_per_year)) {
    stop("days_per_year must be a positive integer.")
  }
  
  if (!is.finite(consumption_mean) || !is.finite(consumption_se) ||
      consumption_se < 0) {
    stop("Invalid consumption_mean or consumption_se.")
  }
  
  F_j <- stats::rnorm(n_person_years, consumption_mean, consumption_se)
  F_j <- pmin(pmax(F_j, 0), 1)
  
  K_j <- stats::rbinom(n_person_years, days_per_year, F_j)
  annual_risk <- numeric(n_person_years)
  
  for (j in which(K_j > 0)) {
    event_risks <- sample(daily_risk, K_j[j], replace = TRUE)
    annual_risk[j] <- -expm1(sum(log1p(-event_risks)))
  }
  
  data.frame(
    person_year = seq_len(n_person_years),
    F_j = F_j,
    K_j = K_j,
    annual_risk = annual_risk
  )
}


# Apply annual-risk simulation to an SD-DRM results data frame

annualize_sd_drm <- function(
    risk_data,
    risk_col = "risk_total",
    scenario = NA_character_,
    n_person_years = 100000,
    consumption_mean = 0.0369,
    consumption_se = 0.0031,
    days_per_year = 365,
    seed = NULL) {
  
  if (!is.data.frame(risk_data)) stop("risk_data must be a data frame.")
  
  if (!risk_col %in% names(risk_data)) {
    stop(paste0("'", risk_col, "' was not found in risk_data."))
  }
  
  out <- simulate_annual_risk(
    daily_risk = risk_data[[risk_col]],
    n_person_years = n_person_years,
    consumption_mean = consumption_mean,
    consumption_se = consumption_se,
    days_per_year = days_per_year,
    seed = seed
  )
  
  out$Scenario <- scenario
  out
}


# Summarize an annual infection-risk distribution

summarize_annual_risk <- function(annual_data,
                                  risk_col = "annual_risk") {
  
  if (!is.data.frame(annual_data)) stop("annual_data must be a data frame.")
  
  if (!risk_col %in% names(annual_data)) {
    stop(paste0("'", risk_col, "' was not found in annual_data."))
  }
  
  x <- annual_data[[risk_col]]
  x <- x[is.finite(x) & x >= 0 & x <= 1]
  
  if (length(x) == 0) stop("No valid annual-risk values were found.")
  
  scenario_value <- if ("Scenario" %in% names(annual_data)) {
    unique(as.character(annual_data$Scenario))
  } else {
    NA_character_
  }
  
  if (length(scenario_value) != 1) scenario_value <- NA_character_
  
  data.frame(
    Scenario = scenario_value,
    n_person_years = length(x),
    probability_zero_annual_risk = mean(x == 0),
    mean_annual_risk = mean(x),
    median_annual_risk = stats::median(x),
    P2.5_annual_risk = unname(stats::quantile(x, 0.025)),
    P25_annual_risk = unname(stats::quantile(x, 0.25)),
    P75_annual_risk = unname(stats::quantile(x, 0.75)),
    P97.5_annual_risk = unname(stats::quantile(x, 0.975)),
    
    mean_consumption_days = if ("K_j" %in% names(annual_data)) {
      mean(annual_data$K_j)
    } else {
      NA_real_
    },
    
    stringsAsFactors = FALSE
  )
}


# Add log10-transformed columns

add_log10_cols <- function(data, cols) {
  
  missing_cols <- setdiff(cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  for (v in cols) {
    data[[paste0("log10_", v)]] <- ifelse(
      data[[v]] > 0,
      log10(data[[v]]),
      NA_real_
    )
  }
  
  data
}


# Safely transform values to log10 scale

safe_log10 <- function(x) {
  ifelse(x > 0, log10(x), NA_real_)
}


# Summarize selected unconditional event-level risk estimates

make_selected_summary <- function(
    data,
    selected_value,
    scenario_exploration,
    selected_condition) {
  
  x <- subset(data, value == selected_value)
  
  if (nrow(x) == 0) stop("No rows matched selected_value.")
  
  data.frame(
    Scenario_exploration = scenario_exploration,
    Selected_condition = selected_condition,
    Exposure_scenario = as.character(x$Scenario),
    Probability_zero_risk = round(x$probability_zero, 4),
    log10_unconditional_mean_risk = round(safe_log10(x$mean), 3),
    log10_unconditional_median_risk = round(safe_log10(x$median), 3),
    Dominant_outcome = as.character(x$outcome),
    stringsAsFactors = FALSE
  )
}