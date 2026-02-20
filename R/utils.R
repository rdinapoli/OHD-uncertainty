# utils.R - Shared helper functions for phase modeling and visualization
#
# Purpose: Provides data preparation, posterior validation, and OxCal-style plotting
#
# Used by: 15_site_15_233.R, 16_ahu_nau_nau.R, 17_ordering_test.R, 18_phase_boundary_figures.R
#
# Contents:
#   Section 1 - Stan data preparation and validation
#     - validate_calibration_for_stan()
#     - validate_stan_data_phases()
#     - prepare_stan_data_ohd_ordered()
#     - prepare_stan_data_c14_ordered()
#     - prepare_stan_data_combined_ordered()
#     - print_stan_data_summary()
#
#   Section 2 - Posterior validation
#     - validate_posterior_ages()
#     - validate_phase_model_recovery()
#     - plot_posterior_validation()
#
#   Section 3 - OxCal-style plotting
#     - calibrate_single_date()
#     - ohd_unmodeled_distribution()
#     - compute_hpd()
#     - plot_oxcal_c14()
#     - plot_oxcal_ohd()
#     - plot_phase_model_oxcal()
#     - plot_all_dates_summary()
#     - generate_phase_model_diagnostic()
#     - plot_oxcal_multipanel()


# --- Section 1: Stan data preparation and validation ---
#
# Centralized Stan data preparation with validation for phase models.
# Reference: Adapted from marine_date_sim/R/utils/prepare_stan_data.R

#' Validate calibration data for Stan models
#'
#' @param cal_ages Matrix or list of calibrated age grids (BP)
#' @param cal_probs Matrix or list of probability densities
#' @param N_grid Expected grid length
#' @return TRUE if valid, stops with error otherwise
validate_calibration_for_stan <- function(cal_ages, cal_probs, N_grid) {
  # Check dimensions match
  if (is.list(cal_ages)) {
    n_samples <- length(cal_ages)
    if (length(cal_probs) != n_samples) {
      stop("FATAL: cal_ages and cal_probs have different numbers of samples")
    }

    for (i in seq_len(n_samples)) {
      if (length(cal_ages[[i]]) != N_grid) {
        stop(sprintf("FATAL: cal_ages[[%d]] has length %d, expected %d",
                     i, length(cal_ages[[i]]), N_grid))
      }
      if (length(cal_probs[[i]]) != N_grid) {
        stop(sprintf("FATAL: cal_probs[[%d]] has length %d, expected %d",
                     i, length(cal_probs[[i]]), N_grid))
      }

      # Check probabilities sum to approximately 1
      prob_sum <- sum(cal_probs[[i]])
      if (abs(prob_sum - 1.0) > 0.01) {
        warning(sprintf("cal_probs[[%d]] sums to %.4f (expected ~1.0)", i, prob_sum))
      }
    }
  }

  TRUE
}


#' Validate Stan data for phase models
#'
#' @param stan_data List of Stan data
#' @return TRUE if valid, stops with error otherwise
validate_stan_data_phases <- function(stan_data) {
  # Check N and N_phases exist
  if (is.null(stan_data$N) && is.null(stan_data$N_c14) && is.null(stan_data$N_ohd)) {
    stop("FATAL: stan_data must contain N, N_c14, or N_ohd")
  }

  if (is.null(stan_data$N_phases)) {
    stop("FATAL: stan_data must contain N_phases")
  }

  if (stan_data$N_phases < 2) {
    stop("FATAL: N_phases must be >= 2")
  }

  # Check phase_id values are valid
  if (!is.null(stan_data$phase_id)) {
    if (any(stan_data$phase_id < 1 | stan_data$phase_id > stan_data$N_phases)) {
      stop(sprintf("FATAL: phase_id values must be in [1, %d]", stan_data$N_phases))
    }
  }

  if (!is.null(stan_data$phase_id_c14)) {
    if (any(stan_data$phase_id_c14 < 1 | stan_data$phase_id_c14 > stan_data$N_phases)) {
      stop(sprintf("FATAL: phase_id_c14 values must be in [1, %d]", stan_data$N_phases))
    }
  }

  if (!is.null(stan_data$phase_id_ohd)) {
    if (any(stan_data$phase_id_ohd < 1 | stan_data$phase_id_ohd > stan_data$N_phases)) {
      stop(sprintf("FATAL: phase_id_ohd values must be in [1, %d]", stan_data$N_phases))
    }
  }

  # Check age bounds
  if (!is.null(stan_data$age_min) && !is.null(stan_data$age_max)) {
    if (stan_data$age_min >= stan_data$age_max) {
      stop(sprintf("FATAL: age_min (%.1f) must be less than age_max (%.1f)",
                   stan_data$age_min, stan_data$age_max))
    }
  }

  # Check OHD data if present
  if (!is.null(stan_data$OHD_age_BP)) {
    if (any(is.na(stan_data$OHD_age_BP))) {
      stop("FATAL: OHD_age_BP contains NA values")
    }
    if (any(stan_data$OHD_SD <= 0)) {
      stop("FATAL: OHD_SD must be positive")
    }
  }

  TRUE
}


#' Prepare Stan data for OHD ordered phase model
#'
#' @param ohd_data Data frame with OHD dates
#' @param age_min Minimum age bound (BP)
#' @param age_max Maximum age bound (BP)
#' @param phase_col Column name for phase assignment
#' @param age_col Column name for OHD age (BP)
#' @param sd_col Column name for OHD SD
#' @return List of Stan data
prepare_stan_data_ohd_ordered <- function(ohd_data, age_min, age_max,
                                          phase_col = "Phase",
                                          age_col = "OHD_BP",
                                          sd_col = "OHD_SD") {
  N <- nrow(ohd_data)
  N_phases <- length(unique(ohd_data[[phase_col]]))

  stan_data <- list(
    N = N,
    N_phases = N_phases,
    phase_id = as.integer(ohd_data[[phase_col]]),
    OHD_age_BP = ohd_data[[age_col]],
    OHD_SD = ohd_data[[sd_col]],
    age_min = age_min,
    age_max = age_max
  )

  validate_stan_data_phases(stan_data)

  stan_data
}


#' Prepare Stan data for radiocarbon ordered phase model
#'
#' @param c14_data Data frame with radiocarbon dates
#' @param cal_ages List of calibrated age grids (one per sample)
#' @param cal_probs List of probability densities (one per sample)
#' @param age_min Minimum age bound (BP)
#' @param age_max Maximum age bound (BP)
#' @param phase_col Column name for phase assignment
#' @return List of Stan data
prepare_stan_data_c14_ordered <- function(c14_data, cal_ages, cal_probs,
                                          age_min, age_max,
                                          phase_col = "Phase") {
  N <- nrow(c14_data)
  N_phases <- length(unique(c14_data[[phase_col]]))
  N_grid <- length(cal_ages[[1]])

  # Validate calibration data
  validate_calibration_for_stan(cal_ages, cal_probs, N_grid)

  stan_data <- list(
    N = N,
    N_phases = N_phases,
    phase_id = as.integer(c14_data[[phase_col]]),
    N_grid = N_grid,
    cal_ages = cal_ages,
    cal_probs = cal_probs,
    age_min = age_min,
    age_max = age_max
  )

  validate_stan_data_phases(stan_data)

  stan_data
}


#' Prepare Stan data for combined radiocarbon + OHD phase model
#'
#' @param c14_data Data frame with radiocarbon dates
#' @param ohd_data Data frame with OHD dates
#' @param cal_ages List of calibrated age grids (one per C14 sample)
#' @param cal_probs List of probability densities (one per C14 sample)
#' @param age_min Minimum age bound (BP)
#' @param age_max Maximum age bound (BP)
#' @param phase_col Column name for phase assignment
#' @param ohd_age_col Column name for OHD age (BP)
#' @param ohd_sd_col Column name for OHD SD
#' @return List of Stan data
prepare_stan_data_combined_ordered <- function(c14_data, ohd_data,
                                                cal_ages, cal_probs,
                                                age_min, age_max,
                                                phase_col = "Phase",
                                                ohd_age_col = "OHD_age_BP",
                                                ohd_sd_col = "OHD_SD") {
  N_c14 <- nrow(c14_data)
  N_ohd <- nrow(ohd_data)
  N_grid <- length(cal_ages[[1]])

  # Determine N_phases from combined data
  all_phases <- c(c14_data[[phase_col]], ohd_data[[phase_col]])
  N_phases <- length(unique(all_phases))

  # Validate calibration data
  validate_calibration_for_stan(cal_ages, cal_probs, N_grid)

  stan_data <- list(
    N_c14 = N_c14,
    N_ohd = N_ohd,
    N_phases = N_phases,
    phase_id_c14 = as.integer(c14_data[[phase_col]]),
    phase_id_ohd = as.integer(ohd_data[[phase_col]]),
    N_grid = N_grid,
    cal_ages = cal_ages,
    cal_probs = cal_probs,
    OHD_age_BP = ohd_data[[ohd_age_col]],
    OHD_SD = ohd_data[[ohd_sd_col]],
    age_min = age_min,
    age_max = age_max
  )

  validate_stan_data_phases(stan_data)

  stan_data
}


#' Print summary of Stan data
#'
#' @param stan_data List of Stan data
print_stan_data_summary <- function(stan_data) {
  cat("Stan Data Summary:\n")

  if (!is.null(stan_data$N)) {
    cat(sprintf("  N samples: %d\n", stan_data$N))
  }
  if (!is.null(stan_data$N_c14)) {
    cat(sprintf("  N C14 samples: %d\n", stan_data$N_c14))
  }
  if (!is.null(stan_data$N_ohd)) {
    cat(sprintf("  N OHD samples: %d\n", stan_data$N_ohd))
  }

  cat(sprintf("  N phases: %d\n", stan_data$N_phases))
  cat(sprintf("  Age range: %.1f - %.1f BP\n", stan_data$age_min, stan_data$age_max))

  if (!is.null(stan_data$N_grid)) {
    cat(sprintf("  Calibration grid size: %d\n", stan_data$N_grid))
  }

  cat("\n")
}


# --- Section 2: Posterior validation ---
#
# Posterior validation functions for phase models.
# Reference: Adapted from marine_date_sim/R/utils/validate_posterior.R

#' Validate posterior ages from a phase model
#'
#' Checks for boundary behavior, reasonable CI widths, and (optionally)
#' agreement with true values for synthetic data validation.
#'
#' @param fit A stanfit object
#' @param true_ages Optional vector of true ages (for synthetic data)
#' @param age_min Minimum age bound
#' @param age_max Maximum age bound
#' @param model_type Type of model ("radiocarbon", "ohd", or "combined")
#' @param age_param Name of the age parameter (default "true_age")
#' @param stop_on_fail Whether to stop on validation failure
#' @return List with validation results
validate_posterior_ages <- function(fit, true_ages = NULL,
                                    age_min, age_max,
                                    model_type = "ohd",
                                    age_param = "true_age",
                                    stop_on_fail = FALSE) {
  # Extract posterior
  posterior <- as.data.frame(fit)

  # Find age columns
  age_cols <- grep(sprintf("^%s\\[", age_param), names(posterior), value = TRUE)

  if (length(age_cols) == 0) {
    msg <- sprintf("No parameters matching '%s[' found in posterior", age_param)
    if (stop_on_fail) stop(msg)
    return(list(passed = FALSE, errors = msg))
  }

  ages <- as.matrix(posterior[, age_cols])
  n_events <- ncol(ages)

  # Initialize results
  warnings <- character(0)
  errors <- character(0)

  # Check for ages at boundaries
  boundary_tolerance <- (age_max - age_min) * 0.01  # 1% of range

  n_at_min <- sum(apply(ages, 2, function(x) {
    quantile(x, 0.5) < age_min + boundary_tolerance
  }))

  n_at_max <- sum(apply(ages, 2, function(x) {
    quantile(x, 0.5) > age_max - boundary_tolerance
  }))

  if (n_at_min > 0) {
    warnings <- c(warnings,
                  sprintf("%d events have posterior median near age_min", n_at_min))
  }

  if (n_at_max > 0) {
    warnings <- c(warnings,
                  sprintf("%d events have posterior median near age_max", n_at_max))
  }

  # Calculate CI widths
  ci_widths <- apply(ages, 2, function(x) {
    diff(quantile(x, c(0.025, 0.975)))
  })

  mean_ci_width <- mean(ci_widths)
  max_ci_width <- max(ci_widths)

  # Check for unreasonably narrow CIs (possible convergence issue)
  if (any(ci_widths < 1)) {
    warnings <- c(warnings,
                  sprintf("%d events have CI width < 1 year (suspiciously narrow)",
                          sum(ci_widths < 1)))
  }

  # Check for unreasonably wide CIs
  expected_max_width <- age_max - age_min
  if (any(ci_widths > expected_max_width * 0.9)) {
    warnings <- c(warnings,
                  sprintf("%d events have CI width > 90%% of age range (uninformative)",
                          sum(ci_widths > expected_max_width * 0.9)))
  }

  # If true ages provided, calculate coverage and bias
  bias <- NA_real_
  coverage <- NA_real_

  if (!is.null(true_ages)) {
    if (length(true_ages) != n_events) {
      errors <- c(errors,
                  sprintf("true_ages length (%d) != n_events (%d)",
                          length(true_ages), n_events))
    } else {
      # Calculate bias
      posterior_means <- colMeans(ages)
      bias <- mean(posterior_means - true_ages)

      # Calculate coverage (% of true values within 95% CI)
      in_ci <- sapply(seq_len(n_events), function(i) {
        ci <- quantile(ages[, i], c(0.025, 0.975))
        true_ages[i] >= ci[1] & true_ages[i] <= ci[2]
      })
      coverage <- mean(in_ci)

      # Check coverage
      if (coverage < 0.85) {
        warnings <- c(warnings,
                      sprintf("Coverage = %.1f%% (expected ~95%%)", coverage * 100))
      }

      # Check bias
      if (abs(bias) > mean_ci_width / 4) {
        warnings <- c(warnings,
                      sprintf("Mean bias = %.1f years (>25%% of mean CI width)", bias))
      }
    }
  }

  # Determine pass/fail
  passed <- length(errors) == 0

  list(
    passed = passed,
    n_events = n_events,
    mean_ci_width = mean_ci_width,
    max_ci_width = max_ci_width,
    n_at_boundary = n_at_min + n_at_max,
    bias = bias,
    coverage = coverage,
    warnings = warnings,
    errors = errors
  )
}


#' Validate phase model parameter recovery
#'
#' For synthetic data with known true values.
#'
#' @param fit A stanfit object
#' @param true_boundaries True phase boundaries
#' @param true_ages True event ages
#' @param stop_on_fail Whether to stop on validation failure
#' @return List with validation results
validate_phase_model_recovery <- function(fit, true_boundaries, true_ages,
                                          stop_on_fail = FALSE) {
  posterior <- as.data.frame(fit)

  # Validate boundaries
  boundary_cols <- grep("^phase_boundary\\[", names(posterior), value = TRUE)

  if (length(boundary_cols) != length(true_boundaries)) {
    msg <- sprintf("Number of boundary parameters (%d) != true_boundaries (%d)",
                   length(boundary_cols), length(true_boundaries))
    if (stop_on_fail) stop(msg)
    return(list(passed = FALSE, errors = msg))
  }

  boundaries <- as.matrix(posterior[, boundary_cols])

  # Check boundary recovery
  boundary_means <- colMeans(boundaries)
  boundary_errors <- boundary_means - true_boundaries

  boundary_in_ci <- sapply(seq_along(true_boundaries), function(i) {
    ci <- quantile(boundaries[, i], c(0.025, 0.975))
    true_boundaries[i] >= ci[1] & true_boundaries[i] <= ci[2]
  })

  boundary_coverage <- mean(boundary_in_ci)

  # Validate ages
  age_validation <- validate_posterior_ages(
    fit, true_ages = true_ages,
    age_min = min(true_ages) - 100,
    age_max = max(true_ages) + 100,
    stop_on_fail = FALSE
  )

  list(
    passed = boundary_coverage >= 0.8 & age_validation$passed,
    boundary_coverage = boundary_coverage,
    boundary_bias = mean(boundary_errors),
    boundary_rmse = sqrt(mean(boundary_errors^2)),
    age_coverage = age_validation$coverage,
    age_bias = age_validation$bias,
    age_mean_ci_width = age_validation$mean_ci_width
  )
}


#' Plot posterior validation results
#'
#' @param fit A stanfit object
#' @param true_ages Optional vector of true ages
#' @param age_min Minimum age bound
#' @param age_max Maximum age bound
#' @param title Plot title
#' @return A ggplot object
plot_posterior_validation <- function(fit, true_ages = NULL,
                                      age_min, age_max,
                                      title = "Posterior Validation") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for this function")
  }

  posterior <- as.data.frame(fit)
  age_cols <- grep("^true_age\\[", names(posterior), value = TRUE)
  ages <- as.matrix(posterior[, age_cols])

  n_events <- ncol(ages)

  # Create summary data frame
  df <- tibble(
    event = seq_len(n_events),
    mean = colMeans(ages),
    ci_lower = apply(ages, 2, quantile, 0.025),
    ci_upper = apply(ages, 2, quantile, 0.975)
  )

  if (!is.null(true_ages)) {
    df$true_age <- true_ages
    df$in_ci <- df$true_age >= df$ci_lower & df$true_age <= df$ci_upper
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = event)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                           width = 0.3, color = "gray50") +
    ggplot2::geom_point(ggplot2::aes(y = mean), size = 2)

  if (!is.null(true_ages)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(y = true_age, color = in_ci),
                                  shape = 4, size = 3)
  }

  p <- p +
    ggplot2::labs(title = title, x = "Event", y = "Age (BP)") +
    ggplot2::theme_minimal()

  p
}


# --- Section 3: OxCal-style plotting ---
#
# OxCal-style diagnostic plots for phase model validation.
#
# OxCal Conventions:
# - Light distribution: Unmodeled (raw calibration or prior)
# - Dark distribution: Modeled (posterior from Bayesian model)
# - Vertical line: True age (in simulation context)
# - Bracket: 95% HPD interval
#
# Key differences from marine_date_sim:
# - OHD dates: Use normal distribution as "unmodeled" (measurement +/- SD)
# - C14 dates: Calibrate against SHCal20 (Southern Hemisphere for Easter Island)
# - Combined visualization of both date types

#' Calibrate a single C14 date against a calibration curve
#'
#' Computes the unmodeled (prior) probability distribution over calendar ages
#' for a single radiocarbon measurement.
#'
#' @param c14_age Measured C14 age (BP)
#' @param c14_error Measurement uncertainty (1 sigma)
#' @param curve_cal Vector of calendar ages from calibration curve
#' @param curve_c14 Vector of C14 ages from calibration curve
#' @param curve_sigma Vector of curve uncertainties
#' @param cal_range Optional vector of calendar ages to evaluate
#' @param normalize If TRUE, normalize to sum to 1
#'
#' @return Data frame with columns: cal_age, probability
calibrate_single_date <- function(c14_age, c14_error,
                                   curve_cal, curve_c14, curve_sigma,
                                   cal_range = NULL,
                                   normalize = TRUE) {

  # Use curve range if not specified
  if (is.null(cal_range)) {
    cal_range <- curve_cal
  }

  # Interpolate curve to cal_range
  expected_c14 <- approx(curve_cal, curve_c14, xout = cal_range, rule = 2)$y
  interp_sigma <- approx(curve_cal, curve_sigma, xout = cal_range, rule = 2)$y

  # Total uncertainty at each calendar age
  total_sigma <- sqrt(c14_error^2 + interp_sigma^2)

  # Calibrated probability (Gaussian likelihood)
  prob <- dnorm(c14_age, mean = expected_c14, sd = total_sigma)

  # Normalize to probability distribution
  if (normalize && sum(prob) > 0) {
    prob <- prob / sum(prob)
  }

  data.frame(
    cal_age = cal_range,
    probability = prob
  )
}


#' Create unmodeled distribution for OHD date
#'
#' OHD dates are already in calendar years, so the "unmodeled" distribution
#' is simply a normal distribution centered on the OHD measurement.
#'
#' @param ohd_age OHD age (BP)
#' @param ohd_sd OHD uncertainty (1 sigma)
#' @param cal_range Vector of calendar ages to evaluate
#' @param normalize If TRUE, normalize to sum to 1
#'
#' @return Data frame with columns: cal_age, probability
ohd_unmodeled_distribution <- function(ohd_age, ohd_sd,
                                        cal_range,
                                        normalize = TRUE) {

  # OHD measurement as normal distribution
  prob <- dnorm(cal_range, mean = ohd_age, sd = ohd_sd)

  # Normalize
  if (normalize && sum(prob) > 0) {
    prob <- prob / sum(prob)
  }

  data.frame(
    cal_age = cal_range,
    probability = prob
  )
}


#' Compute highest posterior density (HPD) interval
#'
#' @param cal_ages Vector of calendar ages
#' @param probs Vector of probabilities (must sum to ~1)
#' @param level Credible level (default 0.95)
#'
#' @return List with lower, upper bounds and the level
compute_hpd <- function(cal_ages, probs, level = 0.95) {
  # Sort by probability (descending)
  ord <- order(probs, decreasing = TRUE)
  sorted_probs <- probs[ord]
  sorted_ages <- cal_ages[ord]

  # Accumulate until we reach desired level
  cumprob <- cumsum(sorted_probs)
  in_hpd <- cumprob <= level

  # The HPD region includes all points up to the level
  hpd_ages <- sorted_ages[in_hpd]

  list(
    lower = min(hpd_ages),
    upper = max(hpd_ages),
    level = level
  )
}


#' Create OxCal-style plot for a single C14 event
#'
#' @param event_id Event identifier
#' @param c14_age Measured C14 age
#' @param c14_error Measurement uncertainty
#' @param curve_cal Calibration curve calendar ages
#' @param curve_c14 Calibration curve C14 ages
#' @param curve_sigma Calibration curve uncertainties
#' @param posterior_samples Vector of posterior samples for this event
#' @param true_age True calendar age (if known)
#' @param phase_boundaries Optional vector of phase boundary ages
#' @param cal_min Minimum calendar age for plot
#' @param cal_max Maximum calendar age for plot
#'
#' @return ggplot object
plot_oxcal_c14 <- function(event_id,
                            c14_age, c14_error,
                            curve_cal, curve_c14, curve_sigma,
                            posterior_samples,
                            true_age = NULL,
                            phase_boundaries = NULL,
                            cal_min = NULL, cal_max = NULL) {

  # Determine plot range from posterior if not specified
  if (is.null(cal_min)) cal_min <- min(posterior_samples) - 100
  if (is.null(cal_max)) cal_max <- max(posterior_samples) + 100

  # Create fine grid for calibration
  cal_grid <- seq(cal_min, cal_max, length.out = 500)

  # Get unmodeled (raw calibrated) distribution
  unmodeled <- calibrate_single_date(
    c14_age, c14_error,
    curve_cal, curve_c14, curve_sigma,
    cal_range = cal_grid
  )

  # Get modeled (posterior) distribution via KDE
  posterior_dens <- density(posterior_samples, from = cal_min, to = cal_max, n = 500)
  modeled <- data.frame(
    cal_age = posterior_dens$x,
    probability = posterior_dens$y / sum(posterior_dens$y)
  )

  # Scale distributions to same height for visual comparison
  max_prob <- max(c(unmodeled$probability, modeled$probability))
  unmodeled$prob_scaled <- unmodeled$probability / max_prob
  modeled$prob_scaled <- modeled$probability / max_prob

  # Compute HPD intervals
  unmodeled_hpd <- compute_hpd(unmodeled$cal_age, unmodeled$probability)
  modeled_hpd <- compute_hpd(modeled$cal_age, modeled$probability)

  # Build plot
  p <- ggplot() +
    # Unmodeled distribution (light blue)
    geom_area(data = unmodeled,
              aes(x = cal_age, y = prob_scaled),
              fill = "lightblue", alpha = 0.5, color = "steelblue") +
    # Modeled distribution (dark gray)
    geom_area(data = modeled,
              aes(x = cal_age, y = prob_scaled),
              fill = "gray30", alpha = 0.6, color = "black") +
    # Phase boundaries
    {if (!is.null(phase_boundaries)) {
      geom_vline(xintercept = phase_boundaries,
                 linetype = "dashed", color = "darkgreen", linewidth = 0.5)
    }} +
    # True age (if known)
    {if (!is.null(true_age)) {
      geom_vline(xintercept = true_age,
                 linetype = "solid", color = "red", linewidth = 1)
    }} +
    # HPD brackets
    annotate("segment",
             x = unmodeled_hpd$lower, xend = unmodeled_hpd$upper,
             y = -0.05, yend = -0.05,
             color = "steelblue", linewidth = 2) +
    annotate("segment",
             x = modeled_hpd$lower, xend = modeled_hpd$upper,
             y = -0.12, yend = -0.12,
             color = "black", linewidth = 2) +
    # Labels
    scale_x_reverse(name = "Calendar Age (BP)") +
    scale_y_continuous(name = "Relative Probability", limits = c(-0.2, 1.1)) +
    labs(
      title = sprintf("C14 Event %s", event_id),
      subtitle = sprintf("C14: %d ± %d BP", round(c14_age), round(c14_error))
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )

  return(p)
}


#' Create OxCal-style plot for a single OHD event
#'
#' @param event_id Event identifier
#' @param ohd_age OHD age (BP)
#' @param ohd_sd OHD uncertainty
#' @param posterior_samples Vector of posterior samples for this event
#' @param true_age True calendar age (if known)
#' @param phase_boundaries Optional vector of phase boundary ages
#' @param cal_min Minimum calendar age for plot
#' @param cal_max Maximum calendar age for plot
#'
#' @return ggplot object
plot_oxcal_ohd <- function(event_id,
                            ohd_age, ohd_sd,
                            posterior_samples,
                            true_age = NULL,
                            phase_boundaries = NULL,
                            cal_min = NULL, cal_max = NULL) {

  # Determine plot range
  if (is.null(cal_min)) cal_min <- min(c(posterior_samples, ohd_age - 3*ohd_sd)) - 50
  if (is.null(cal_max)) cal_max <- max(c(posterior_samples, ohd_age + 3*ohd_sd)) + 50

  cal_grid <- seq(cal_min, cal_max, length.out = 500)

  # Get unmodeled (OHD measurement as normal)
  unmodeled <- ohd_unmodeled_distribution(ohd_age, ohd_sd, cal_grid)

  # Get modeled (posterior) distribution via KDE
  posterior_dens <- density(posterior_samples, from = cal_min, to = cal_max, n = 500)
  modeled <- data.frame(
    cal_age = posterior_dens$x,
    probability = posterior_dens$y / sum(posterior_dens$y)
  )

  # Scale distributions
  max_prob <- max(c(unmodeled$probability, modeled$probability))
  unmodeled$prob_scaled <- unmodeled$probability / max_prob
  modeled$prob_scaled <- modeled$probability / max_prob

  # Compute HPD intervals
  unmodeled_hpd <- compute_hpd(unmodeled$cal_age, unmodeled$probability)
  modeled_hpd <- compute_hpd(modeled$cal_age, modeled$probability)

  # Build plot (coral color scheme for OHD to distinguish from C14)
  p <- ggplot() +
    # Unmodeled distribution (light coral)
    geom_area(data = unmodeled,
              aes(x = cal_age, y = prob_scaled),
              fill = "mistyrose", alpha = 0.5, color = "coral") +
    # Modeled distribution (dark gray)
    geom_area(data = modeled,
              aes(x = cal_age, y = prob_scaled),
              fill = "gray30", alpha = 0.6, color = "black") +
    # Phase boundaries
    {if (!is.null(phase_boundaries)) {
      geom_vline(xintercept = phase_boundaries,
                 linetype = "dashed", color = "darkgreen", linewidth = 0.5)
    }} +
    # True age (if known)
    {if (!is.null(true_age)) {
      geom_vline(xintercept = true_age,
                 linetype = "solid", color = "red", linewidth = 1)
    }} +
    # HPD brackets
    annotate("segment",
             x = unmodeled_hpd$lower, xend = unmodeled_hpd$upper,
             y = -0.05, yend = -0.05,
             color = "coral", linewidth = 2) +
    annotate("segment",
             x = modeled_hpd$lower, xend = modeled_hpd$upper,
             y = -0.12, yend = -0.12,
             color = "black", linewidth = 2) +
    # Labels
    scale_x_reverse(name = "Calendar Age (BP)") +
    scale_y_continuous(name = "Relative Probability", limits = c(-0.2, 1.1)) +
    labs(
      title = sprintf("OHD Event %s", event_id),
      subtitle = sprintf("OHD: %d ± %d BP", round(ohd_age), round(ohd_sd))
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )

  return(p)
}


#' Create OxCal-style vertical strip plot for phase model
#'
#' Classic OxCal layout with events arranged vertically by stratigraphic order,
#' showing phase boundaries as colored bands.
#'
#' @param event_data Data frame with columns: event_id, phase, age_type ("C14" or "OHD"),
#'                   measurement, error, posterior_mean, ci_lower, ci_upper
#' @param phase_boundaries Vector of phase boundary posterior means
#' @param phase_boundary_ci Matrix with columns lower/upper for each boundary
#' @param title Plot title
#'
#' @return ggplot object
plot_phase_model_oxcal <- function(event_data,
                                    phase_boundaries,
                                    phase_boundary_ci = NULL,
                                    title = "Phase Model Results") {

  n_events <- nrow(event_data)
  n_phases <- length(phase_boundaries) + 1

  # Create phase background bands
  phase_bands <- data.frame(
    phase = 1:n_phases,
    xmin = c(max(phase_boundaries) + 50, phase_boundaries),
    xmax = c(phase_boundaries, min(phase_boundaries) - 50)
  )
  phase_bands$fill <- ifelse(phase_bands$phase %% 2 == 1, "gray95", "gray85")

  # Order events by phase (oldest=highest phase number at bottom, youngest=Phase 1 at top)
  # Within each phase, order by posterior mean (older first = higher y position)
  # This creates: Phase N (oldest) at y=1 (bottom), Phase 1 (youngest) at y=n_events (top)
  event_data <- event_data %>%
    arrange(desc(phase), posterior_mean) %>%
    mutate(y_pos = row_number())

  # Build plot
  p <- ggplot() +
    # Phase background bands
    geom_rect(data = phase_bands,
              aes(xmin = xmin, xmax = xmax, ymin = 0.5, ymax = n_events + 0.5),
              fill = phase_bands$fill, alpha = 0.5) +
    # Phase boundary lines
    geom_vline(xintercept = phase_boundaries,
               linetype = "dashed", color = "darkgreen", linewidth = 1) +
    # Phase boundary CIs (if provided)
    {if (!is.null(phase_boundary_ci)) {
      geom_rect(data = data.frame(
        xmin = phase_boundary_ci[, "lower"],
        xmax = phase_boundary_ci[, "upper"],
        ymin = 0.5,
        ymax = n_events + 0.5
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "darkgreen", alpha = 0.1)
    }} +
    # 95% CI bars
    geom_segment(data = event_data,
                 aes(x = ci_lower, xend = ci_upper, y = y_pos, yend = y_pos,
                     color = age_type),
                 linewidth = 1, alpha = 0.7) +
    # Posterior median points
    geom_point(data = event_data,
               aes(x = posterior_mean, y = y_pos, color = age_type),
               size = 2) +
    # Measurement points (hollow)
    geom_point(data = event_data,
               aes(x = measurement, y = y_pos, color = age_type),
               shape = 1, size = 2.5) +
    # Color scheme
    scale_color_manual(values = c("C14" = "steelblue", "OHD" = "coral"),
                       name = "Date Type") +
    # Axis formatting
    scale_x_reverse(name = "Calendar Age (BP)") +
    scale_y_continuous(name = "Event",
                       breaks = event_data$y_pos,
                       labels = event_data$event_id) +
    labs(
      title = title,
      subtitle = paste0("Filled points = Posterior median | ",
                        "Open circles = Measurement | ",
                        "Bars = 95% CI")
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  # Add phase labels on right side
  phase_label_data <- event_data %>%
    group_by(phase) %>%
    summarise(y_mid = mean(y_pos), .groups = "drop")

  # Get x-axis range for label placement
  x_range <- range(c(event_data$ci_lower, event_data$ci_upper), na.rm = TRUE)

  p <- p +
    geom_text(data = phase_label_data,
              aes(x = x_range[1] - (x_range[2] - x_range[1]) * 0.02,
                  y = y_mid,
                  label = paste0("Phase ", phase)),
              hjust = 1, size = 3, fontface = "bold")

  return(p)
}


#' Plot all dates summary panel (OxCal-style)
#'
#' Creates a single-panel figure showing all dates organized by phase.
#' Youngest phase (Phase 1) at top, oldest phase at bottom.
#' Older dates on left (high BP), younger dates on right (low BP).
#'
#' @param event_data Data frame with columns: event_id, phase, age_type,
#'        measurement, error, posterior_mean, ci_lower, ci_upper
#' @param phase_boundaries Vector of phase boundary means (in BP)
#' @param phase_boundary_ci Optional matrix of boundary CIs (columns: lower, upper)
#' @param title Plot title
#' @param show_measurements Whether to show original measurements (default TRUE)
#' @param show_phase_bands Whether to show alternating phase backgrounds (default TRUE)
#'
#' @return ggplot object
plot_all_dates_summary <- function(event_data,
                                    phase_boundaries = NULL,
                                    phase_boundary_ci = NULL,
                                    title = "All Dates by Phase",
                                    show_measurements = TRUE,
                                    show_phase_bands = TRUE) {

  n_events <- nrow(event_data)
  n_phases <- max(event_data$phase)

  # Order events: oldest phase (highest number) at bottom (y=1), youngest at top
  # Within phase: older ages at lower y position within the phase group
  event_data <- event_data %>%
    arrange(desc(phase), posterior_mean) %>%
    mutate(y_pos = row_number())

  # Determine x-axis range
  x_min <- min(event_data$ci_lower, na.rm = TRUE) - 50
  x_max <- max(event_data$ci_upper, na.rm = TRUE) + 50

  # Build plot
  p <- ggplot()

  # Add phase background bands if requested
  if (show_phase_bands && !is.null(phase_boundaries)) {
    phase_bands <- data.frame(
      phase = 1:n_phases,
      xmin = c(x_max, phase_boundaries),
      xmax = c(phase_boundaries, x_min)
    )
    phase_bands$fill <- ifelse(phase_bands$phase %% 2 == 1, "gray95", "gray85")

    p <- p +
      geom_rect(data = phase_bands,
                aes(xmin = xmin, xmax = xmax, ymin = 0.5, ymax = n_events + 0.5),
                fill = phase_bands$fill, alpha = 0.5)
  }

  # Add phase boundary lines
  if (!is.null(phase_boundaries)) {
    p <- p +
      geom_vline(xintercept = phase_boundaries,
                 linetype = "dashed", color = "darkgreen", linewidth = 0.8)

    # Add boundary CIs if provided
    if (!is.null(phase_boundary_ci)) {
      p <- p +
        geom_rect(data = data.frame(
          xmin = phase_boundary_ci[, "lower"],
          xmax = phase_boundary_ci[, "upper"],
          ymin = 0.5,
          ymax = n_events + 0.5
        ),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "darkgreen", alpha = 0.1)
    }
  }

  # Add 95% CI bars
  p <- p +
    geom_segment(data = event_data,
                 aes(x = ci_lower, xend = ci_upper, y = y_pos, yend = y_pos,
                     color = age_type),
                 linewidth = 1.5, alpha = 0.7)

  # Add posterior median points
  p <- p +
    geom_point(data = event_data,
               aes(x = posterior_mean, y = y_pos, color = age_type),
               size = 2.5)

  # Add measurement points (hollow) if requested
  if (show_measurements) {
    p <- p +
      geom_point(data = event_data,
                 aes(x = measurement, y = y_pos, color = age_type),
                 shape = 1, size = 3)
  }

  # Color scheme
  p <- p +
    scale_color_manual(values = c("C14" = "steelblue", "OHD" = "coral"),
                       name = "Date Type")

  # Axis formatting (reversed x for BP convention: older left, younger right)
  p <- p +
    scale_x_reverse(name = "Calendar Age (BP)",
                    expand = expansion(mult = c(0.05, 0.15))) +  # Extra space on right for labels
    scale_y_continuous(name = "",
                       breaks = event_data$y_pos,
                       labels = event_data$event_id,
                       expand = expansion(mult = 0.02))

  # Add phase labels on right margin
  phase_label_data <- event_data %>%
    group_by(phase) %>%
    summarise(
      y_min = min(y_pos),
      y_max = max(y_pos),
      y_mid = mean(y_pos),
      .groups = "drop"
    )

  p <- p +
    geom_text(data = phase_label_data,
              aes(x = x_min - (x_max - x_min) * 0.08,
                  y = y_mid,
                  label = paste0("Phase ", phase)),
              hjust = 0.5, size = 3.5, fontface = "bold", color = "gray30")

  # Add horizontal separators between phases
  phase_breaks <- event_data %>%
    group_by(phase) %>%
    summarise(y_max = max(y_pos), .groups = "drop") %>%
    filter(phase < n_phases) %>%
    pull(y_max)

  if (length(phase_breaks) > 0) {
    p <- p +
      geom_hline(yintercept = phase_breaks + 0.5,
                 linetype = "dotted", color = "gray50", linewidth = 0.5)
  }

  # Theme and labels
  p <- p +
    labs(
      title = title,
      subtitle = paste0("Filled = Posterior median | Open = Measurement | Bars = 95% CI\n",
                        "Top = Youngest (Phase 1) | Bottom = Oldest (Phase ", n_phases, ")")
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 8),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )

  return(p)
}


#' Generate comprehensive phase model diagnostic
#'
#' Creates a multi-panel figure with:
#' 1. OxCal-style vertical strip plot
#' 2. Phase boundary posteriors
#' 3. CI width comparison by date type
#'
#' @param fit Stan fit object
#' @param c14_data Data frame with C14 measurements
#' @param ohd_data Data frame with OHD measurements
#' @param model_name Model identifier
#' @param output_file Path for PNG output (optional)
#'
#' @return List of ggplot objects
generate_phase_model_diagnostic <- function(fit,
                                             c14_data = NULL,
                                             ohd_data = NULL,
                                             model_name = "Phase Model",
                                             output_file = NULL) {

  posterior <- as.data.frame(fit)

  # Extract phase boundaries
  boundary_cols <- grep("^phase_boundary\\[", names(posterior), value = TRUE)
  n_boundaries <- length(boundary_cols)

  if (n_boundaries > 0) {
    boundaries <- as.matrix(posterior[, boundary_cols])
    boundary_means <- colMeans(boundaries)
    boundary_ci <- apply(boundaries, 2, quantile, probs = c(0.025, 0.975))
  } else {
    boundary_means <- NULL
    boundary_ci <- NULL
  }

  # Build event data frame
  event_data <- data.frame()

  # Add C14 events
  if (!is.null(c14_data)) {
    c14_age_cols <- grep("^true_age\\[|^true_age_c14\\[", names(posterior), value = TRUE)
    if (length(c14_age_cols) > 0) {
      c14_ages <- as.matrix(posterior[, c14_age_cols])
      for (i in 1:nrow(c14_data)) {
        event_data <- rbind(event_data, data.frame(
          event_id = paste0("C14_", i),
          phase = c14_data$Phase[i],
          age_type = "C14",
          measurement = c14_data$C14_Age[i],
          error = c14_data$C14_SD[i],
          posterior_mean = mean(c14_ages[, i]),
          ci_lower = quantile(c14_ages[, i], 0.025),
          ci_upper = quantile(c14_ages[, i], 0.975)
        ))
      }
    }
  }

  # Add OHD events
  if (!is.null(ohd_data)) {
    ohd_age_cols <- grep("^true_age\\[|^true_age_ohd\\[", names(posterior), value = TRUE)
    if (length(ohd_age_cols) > 0) {
      ohd_ages <- as.matrix(posterior[, ohd_age_cols])
      for (i in 1:nrow(ohd_data)) {
        event_data <- rbind(event_data, data.frame(
          event_id = paste0("OHD_", i),
          phase = ohd_data$Phase[i],
          age_type = "OHD",
          measurement = ohd_data$OHD_BP[i],
          error = ohd_data$OHD_SD[i],
          posterior_mean = mean(ohd_ages[, i]),
          ci_lower = quantile(ohd_ages[, i], 0.025),
          ci_upper = quantile(ohd_ages[, i], 0.975)
        ))
      }
    }
  }

  # Compute CI widths
  event_data$ci_width <- event_data$ci_upper - event_data$ci_lower

  # --- Plot 1: OxCal-style strip plot ---
  p1 <- plot_phase_model_oxcal(
    event_data = event_data,
    phase_boundaries = boundary_means,
    phase_boundary_ci = if (!is.null(boundary_ci)) t(boundary_ci) else NULL,
    title = sprintf("%s: Event Ages", model_name)
  )

  # --- Plot 2: Phase boundary posteriors ---
  if (n_boundaries > 0) {
    boundary_df <- as.data.frame(boundaries) %>%
      tidyr::pivot_longer(everything(), names_to = "boundary", values_to = "age") %>%
      mutate(boundary_num = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", boundary)))

    p2 <- ggplot(boundary_df, aes(x = age, y = factor(boundary_num))) +
      ggridges::geom_density_ridges(fill = "darkgreen", alpha = 0.5, color = "darkgreen") +
      scale_x_reverse(name = "Calendar Age (BP)") +
      labs(
        title = sprintf("%s: Phase Boundaries", model_name),
        y = "Boundary"
      ) +
      theme_minimal()
  } else {
    p2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No phase boundaries in model") +
      theme_void()
  }

  # --- Plot 3: CI width comparison ---
  p3 <- ggplot(event_data, aes(x = age_type, y = ci_width, fill = age_type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = c("C14" = "steelblue", "OHD" = "coral")) +
    labs(
      title = sprintf("%s: Precision Comparison", model_name),
      x = "Date Type",
      y = "95% CI Width (years)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  # --- Plot 4: All dates summary (single panel) ---
  p4 <- plot_all_dates_summary(
    event_data = event_data,
    phase_boundaries = boundary_means,
    phase_boundary_ci = if (!is.null(boundary_ci)) t(boundary_ci) else NULL,
    title = sprintf("%s: All Dates Summary", model_name)
  )

  # Combine and save
  if (!is.null(output_file)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined <- p1 / (p2 | p3) +
        patchwork::plot_layout(heights = c(2, 1))
      ggsave(output_file, combined, width = 12, height = 14, dpi = 300)
      cat(sprintf("Saved diagnostic to %s\n", output_file))

      # Also save the summary panel separately
      summary_file <- gsub("\\.png$", "_summary.png", output_file)
      ggsave(summary_file, p4, width = 10, height = max(8, nrow(event_data) * 0.3), dpi = 300)
      cat(sprintf("Saved summary to %s\n", summary_file))
    }
  }

  list(
    strip_plot = p1,
    boundary_plot = p2,
    precision_plot = p3,
    summary_plot = p4,
    event_data = event_data
  )
}


#' Create multi-panel OxCal-style plot showing unmodeled vs modeled
#'
#' @param c14_data Data frame with C14 measurements (optional)
#' @param ohd_data Data frame with OHD measurements (optional)
#' @param posterior_c14 Matrix of posterior samples for C14 ages (optional)
#' @param posterior_ohd Matrix of posterior samples for OHD ages (optional)
#' @param cal_curve List with cal_BP, c14_BP, sigma (for C14 calibration)
#' @param phase_boundaries Phase boundary ages (optional)
#' @param n_events Maximum events to show (default 8)
#' @param title Main title
#'
#' @return ggplot object (combined panels)
plot_oxcal_multipanel <- function(c14_data = NULL,
                                   ohd_data = NULL,
                                   posterior_c14 = NULL,
                                   posterior_ohd = NULL,
                                   cal_curve = NULL,
                                   phase_boundaries = NULL,
                                   n_events = 8,
                                   title = "OxCal-Style Diagnostic") {

  plots <- list()

  # Determine global age range
  all_ages <- c()
  if (!is.null(posterior_c14)) all_ages <- c(all_ages, as.vector(posterior_c14))
  if (!is.null(posterior_ohd)) all_ages <- c(all_ages, as.vector(posterior_ohd))
  if (!is.null(c14_data)) all_ages <- c(all_ages, c14_data$C14_Age)
  if (!is.null(ohd_data)) all_ages <- c(all_ages, ohd_data$OHD_BP)

  cal_min <- min(all_ages) - 100
  cal_max <- max(all_ages) + 100

  # Add C14 plots
  if (!is.null(c14_data) && !is.null(posterior_c14) && !is.null(cal_curve)) {
    n_c14 <- min(nrow(c14_data), ceiling(n_events / 2))
    indices <- round(seq(1, nrow(c14_data), length.out = n_c14))

    for (i in indices) {
      p <- plot_oxcal_c14(
        event_id = i,
        c14_age = c14_data$C14_Age[i],
        c14_error = c14_data$C14_SD[i],
        curve_cal = cal_curve$cal_BP,
        curve_c14 = cal_curve$c14_BP,
        curve_sigma = cal_curve$sigma,
        posterior_samples = posterior_c14[, i],
        phase_boundaries = phase_boundaries,
        cal_min = cal_min,
        cal_max = cal_max
      )
      plots[[length(plots) + 1]] <- p
    }
  }

  # Add OHD plots
  if (!is.null(ohd_data) && !is.null(posterior_ohd)) {
    n_ohd <- min(nrow(ohd_data), ceiling(n_events / 2))
    indices <- round(seq(1, nrow(ohd_data), length.out = n_ohd))

    for (i in indices) {
      p <- plot_oxcal_ohd(
        event_id = i,
        ohd_age = ohd_data$OHD_BP[i],
        ohd_sd = ohd_data$OHD_SD[i],
        posterior_samples = posterior_ohd[, i],
        phase_boundaries = phase_boundaries,
        cal_min = cal_min,
        cal_max = cal_max
      )
      plots[[length(plots) + 1]] <- p
    }
  }

  # Combine using patchwork
  if (length(plots) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = 2) +
      patchwork::plot_annotation(
        title = title,
        subtitle = paste0("Light fill = Unmodeled | Dark fill = Modeled | ",
                          "Blue = C14 | Coral = OHD"),
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
    return(combined)
  } else if (length(plots) > 0) {
    return(plots)
  } else {
    return(NULL)
  }
}
