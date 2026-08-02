# SD-DRM plotting functions


# Colours for dose-response curves

make_group_colours <- function(groups) {
  
  if (length(groups) == 3) {
    colours <- c("#F8766D", "#7CAE00", "#C77CFF")
  } else if (length(groups) == 4) {
    colours <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
  } else {
    colours <- scales::hue_pal()(length(groups))
  }
  
  names(colours) <- groups
  colours
}


# Common log10 y-axis for mean annual-risk plots

make_common_mean_axis <- function(values, minimum_decades = 3) {
  
  values <- values[is.finite(values) & values > 0]
  
  lower_exponent <- floor(log10(min(values)))
  upper_exponent <- ceiling(log10(max(values)))
  
  if (upper_exponent - lower_exponent < minimum_decades) {
    lower_exponent <- upper_exponent - minimum_decades
  }
  
  list(
    limits = c(
      10^lower_exponent,
      10^upper_exponent
    ),
    breaks = 10^seq(
      lower_exponent,
      upper_exponent,
      by = 1
    )
  )
}


# Annual dose-response plot

plot_dose_response_function <- function(
    data, scenario_levels, colour_title,
    subtitle_text, file_name) {
  
  data$Scenario <- factor(
    data$Scenario,
    levels = scenario_levels
  )
  
  data$group <- factor(
    data$group,
    levels = unique(data$group)
  )
  
  data$status <- factor(
    data$status,
    levels = c(
      "Less likely treatable",
      "More likely treatable"
    )
  )
  
  marker_points <- do.call(
    rbind,
    lapply(
      split(
        data,
        interaction(
          data$Scenario,
          data$group,
          drop = TRUE
        )
      ),
      function(x) {
        x[seq(1, nrow(x), by = 30), ]
      }
    )
  )
  
  colours <- make_group_colours(
    levels(data$group)
  )
  
  plot_title <- if (
    grepl("MIC", colour_title, ignore.case = TRUE)
  ) {
    "Annual dose-response curves by doxycycline concentration"
  } else {
    "Annual dose-response curves by resistant fraction"
  }
  
  p <- ggplot(
    data,
    aes(
      x = dose,
      y = mean_annual_risk,
      colour = group
    )
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(
      data = marker_points,
      aes(shape = status),
      size = 2.2,
      stroke = 0.8
    ) +
    facet_wrap(
      ~Scenario,
      ncol = 2,
      scales = "free"
    ) +
    scale_x_log10(
      breaks = scales::log_breaks(n = 6),
      labels = scales::trans_format(
        "log10",
        scales::math_format(10^.x)
      )
    ) +
    scale_y_log10(
      breaks = scales::log_breaks(n = 6),
      labels = scales::trans_format(
        "log10",
        scales::math_format(10^.x)
      )
    ) +
    scale_colour_manual(
      values = colours
    ) +
    scale_shape_manual(
      values = c(
        "Less likely treatable" = 17,
        "More likely treatable" = 16
      )
    ) +
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = "Pathogenic dose (MPN/consumption event)",
      y = "Mean annual infection risk",
      colour = colour_title,
      shape = "Outcome"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        size = 13,
        face = "bold"
      ),
      legend.position = "top",
      legend.box = "vertical",
      strip.text = element_text(
        size = 11,
        face = "bold"
      ),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = file.path("images", file_name),
    plot = p,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  invisible(p)
}


# Mean annual infection-risk plot

plot_mean_annual_risk_function <- function(
    data, scenario_levels, x_title,
    subtitle_text, file_name,
    x_breaks = NULL) {
  
  data <- data[
    is.finite(data$value) &
      is.finite(data$mean_annual_risk) &
      data$mean_annual_risk > 0,
  ]
  
  data$Scenario <- factor(
    data$Scenario,
    levels = scenario_levels
  )
  
  data$outcome <- factor(
    data$outcome,
    levels = c(
      "More likely treatable",
      "Less likely treatable"
    )
  )
  
  if (is.null(x_breaks)) {
    x_breaks <- sort(unique(data$value))
  }
  
  common_axis <- make_common_mean_axis(
    data$mean_annual_risk
  )
  
  segments <- do.call(
    rbind,
    lapply(
      split(
        data,
        data$Scenario,
        drop = TRUE
      ),
      function(x) {
        
        x <- x[order(x$value), ]
        
        segment_data <- x[-nrow(x), ]
        segment_data$xend <- x$value[-1]
        segment_data$yend <- x$mean_annual_risk[-1]
        
        segment_data
      }
    )
  )
  
  plot_title <- if (
    grepl("MIC", x_title, ignore.case = TRUE)
  ) {
    "Mean annual infection risk by doxycycline concentration"
  } else {
    "Mean annual infection risk by resistant fraction"
  }
  
  p <- ggplot(
    segments,
    aes(
      x = value,
      xend = xend,
      y = mean_annual_risk,
      yend = yend,
      colour = outcome
    )
  ) +
    geom_segment(
      linewidth = 1.1,
      lineend = "butt"
    ) +
    facet_wrap(
      ~Scenario,
      ncol = 2,
      scales = "fixed"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      minor_breaks = NULL,
      expand = expansion(
        mult = c(0.02, 0.02)
      )
    ) +
    scale_y_log10(
      limits = common_axis$limits,
      breaks = common_axis$breaks,
      labels = scales::trans_format(
        "log10",
        scales::math_format(10^.x)
      ),
      expand = expansion(
        mult = c(0.03, 0.06)
      )
    ) +
    scale_colour_manual(
      values = c(
        "More likely treatable" = "black",
        "Less likely treatable" = "red"
      )
    ) +
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = x_title,
      y = "Mean annual infection risk",
      colour = "Outcome"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        size = 13,
        face = "bold"
      ),
      legend.position = "top",
      strip.text = element_text(
        size = 11,
        face = "bold"
      ),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8)
    )
  
  print(p)
  
  ggsave(
    filename = file.path("images", file_name),
    plot = p,
    width = 11,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  invisible(p)
}


# Conditional median annual infection-risk plot

plot_conditional_median_annual_risk_function <- function(
    data, scenario_levels, x_title,
    subtitle_text, file_name,
    x_breaks = NULL) {
  
  data <- data[
    is.finite(data$value) &
      is.finite(data$median_positive_annual_risk) &
      data$median_positive_annual_risk > 0,
  ]
  
  data$Scenario <- factor(
    data$Scenario,
    levels = scenario_levels
  )
  
  data$outcome <- factor(
    data$outcome,
    levels = c(
      "More likely treatable",
      "Less likely treatable"
    )
  )
  
  if (is.null(x_breaks)) {
    x_breaks <- sort(unique(data$value))
  }
  
  segments <- do.call(
    rbind,
    lapply(
      split(
        data,
        data$Scenario,
        drop = TRUE
      ),
      function(x) {
        
        x <- x[order(x$value), ]
        
        segment_data <- x[-nrow(x), ]
        segment_data$xend <- x$value[-1]
        segment_data$yend <-
          x$median_positive_annual_risk[-1]
        
        segment_data
      }
    )
  )
  
  plot_title <- if (
    grepl("MIC", x_title, ignore.case = TRUE)
  ) {
    "Conditional median annual infection risk by residual doxycycline concentration"
  } else {
    "Conditional median annual infection risk by resistant fraction"
  }
  
  p <- ggplot(
    data,
    aes(x = value)
  ) +
    geom_ribbon(
      aes(
        ymin = p025_positive_annual_risk,
        ymax = p975_positive_annual_risk
      ),
      fill = "grey80"
    ) +
    geom_ribbon(
      aes(
        ymin = p125_positive_annual_risk,
        ymax = p875_positive_annual_risk
      ),
      fill = "grey65"
    ) +
    geom_ribbon(
      aes(
        ymin = p25_positive_annual_risk,
        ymax = p75_positive_annual_risk
      ),
      fill = "grey50"
    ) +
    geom_segment(
      data = segments,
      aes(
        x = value,
        xend = xend,
        y = median_positive_annual_risk,
        yend = yend,
        colour = outcome
      ),
      linewidth = 1.2,
      lineend = "butt"
    ) +
    facet_wrap(
      ~Scenario,
      ncol = 2,
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      minor_breaks = NULL,
      expand = expansion(
        mult = c(0.02, 0.02)
      )
    ) +
    scale_y_log10(
      breaks = scales::log_breaks(n = 5),
      labels = scales::trans_format(
        "log10",
        scales::math_format(10^.x)
      )
    ) +
    scale_colour_manual(
      values = c(
        "More likely treatable" = "black",
        "Less likely treatable" = "red"
      )
    ) +
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = x_title,
      y = "Annual infection risk",
      colour = "Outcome"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        size = 13,
        face = "bold"
      ),
      legend.position = "top",
      strip.text = element_text(
        size = 11,
        face = "bold"
      ),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8)
    )
  
  print(p)
  
  ggsave(
    filename = file.path("images", file_name),
    plot = p,
    width = 11,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  invisible(p)
}