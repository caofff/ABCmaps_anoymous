#' Estimate the low-fidelity false-positive rate from pilot simulations
#'
#' Compares low- and high-fidelity ABC decisions on a shared pilot particle
#' population. The function can select pilot thresholds automatically or use
#' thresholds supplied by the caller.
#'
#' @param Theta Pilot particles. A numeric vector is accepted for a
#'   one-dimensional parameter; otherwise supply an `N` by `d` matrix.
#' @param Dis_fun A function with arguments `theta`, `n_rep`, and `type`.
#'   It must return an `nrow(theta)` by `n_rep` discrepancy matrix.
#' @param n_L Number of low-fidelity replicates per pilot particle.
#' @param n_H Number of high-fidelity replicates per pilot particle.
#' @param weights Optional non-negative pilot weights. Uniform weights are used
#'   by default.
#' @param kappa Quantile used to select `epsilon` when `epsilon` is `NULL`.
#' @param epsilon_target Optional lower bound for the selected high-fidelity
#'   threshold.
#' @param epsilon Optional fixed high-fidelity threshold.
#' @param epsilon_L Optional fixed low-fidelity threshold. When it is `NULL`,
#'   the threshold is selected as the `1 - a_L` weighted quantile of minimum
#'   low-fidelity discrepancies under the pilot HF-ABC posterior weights.
#' @param a_L Target posterior-weighted low-fidelity false-negative rate used
#'   when selecting `epsilon_L` automatically.
#' @param cost_mode Cost source for `fpr_max`. Use `"fixed"` for supplied
#'   costs or `"measured"` for pilot wall-clock costs.
#' @param c_L Average cost of one low-fidelity simulation in fixed-cost mode.
#' @param c_H Average cost of one high-fidelity simulation in fixed-cost mode.
#' @param verbose Whether to report the selected thresholds and estimated rates.
#'
#' @return A list with a one-row `summary` data frame, a particle-level
#'   `particles` data frame, and the low- and high-fidelity discrepancy
#'   matrices.
#' @export
estimate_pilot_fpr <- function(Theta, Dis_fun, n_L, n_H,
                               weights = NULL,
                               kappa = 0.1,
                               epsilon_target = NULL,
                               epsilon = NULL,
                               epsilon_L = NULL,
                               a_L = 0.05,
                               cost_mode = c("fixed", "measured"),
                               c_L = 1,
                               c_H = 10,
                               verbose = FALSE) {
  if (!is.function(Dis_fun)) stop("Dis_fun must be a function.", call. = FALSE)
  n_L <- .scalar_integer(n_L, "n_L")
  n_H <- .scalar_integer(n_H, "n_H")
  kappa <- .scalar_probability(kappa, "kappa")
  a_L <- .scalar_probability(a_L, "a_L", include_zero = TRUE)
  cost_mode <- match.arg(cost_mode)

  if (is.null(dim(Theta))) {
    Theta <- matrix(Theta, ncol = 1L)
  } else {
    Theta <- as.matrix(Theta)
  }
  if (!is.numeric(Theta) || nrow(Theta) < 1L || ncol(Theta) < 1L ||
      any(!is.finite(Theta))) {
    stop("Theta must contain at least one finite numeric pilot particle.",
         call. = FALSE)
  }
  pilot_size <- nrow(Theta)

  if (is.null(weights)) {
    weights <- rep(1 / pilot_size, pilot_size)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != pilot_size || any(!is.finite(weights)) ||
        any(weights < 0) || sum(weights) <= 0) {
      stop("weights must be finite, non-negative, and match nrow(Theta).",
           call. = FALSE)
    }
    weights <- weights / sum(weights)
  }

  validate_optional_threshold <- function(value, name) {
    if (!is.null(value) &&
        (!is.numeric(value) || length(value) != 1L || !is.finite(value))) {
      stop(name, " must be NULL or a finite number.", call. = FALSE)
    }
    if (is.null(value)) NULL else as.numeric(value)
  }
  epsilon <- validate_optional_threshold(epsilon, "epsilon")
  epsilon_L <- validate_optional_threshold(epsilon_L, "epsilon_L")
  epsilon_target <- validate_optional_threshold(
    epsilon_target, "epsilon_target"
  )
  if (cost_mode == "fixed" &&
      (!is.numeric(c_L) || length(c_L) != 1L || !is.finite(c_L) ||
       c_L <= 0 || !is.numeric(c_H) || length(c_H) != 1L ||
       !is.finite(c_H) || c_H <= 0)) {
    stop("c_L and c_H must be finite positive numbers in fixed-cost mode.",
         call. = FALSE)
  }

  timed_discrepancy <- function(n_rep, type) {
    start <- Sys.time()
    value <- Dis_fun(Theta, n_rep, type = type)
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    value <- .discrepancy_matrix(
      value, pilot_size, n_rep, paste0("Dis_fun(type = '", type, "')")
    )
    list(value = value, elapsed = elapsed)
  }

  low_timing <- timed_discrepancy(n_L, "low")
  high_timing <- timed_discrepancy(n_H, "high")
  low_discrepancy <- low_timing$value
  high_discrepancy <- high_timing$value
  measured_c_L <- low_timing$elapsed / (pilot_size * n_L)
  measured_c_H <- high_timing$elapsed / (pilot_size * n_H)

  if (cost_mode == "measured") {
    if (!is.finite(measured_c_L) || measured_c_L <= 0 ||
        !is.finite(measured_c_H) || measured_c_H <= 0) {
      stop(
        paste0("Pilot timing did not resolve positive LF and HF costs; ",
               "use cost_mode = 'fixed' or a larger pilot."),
        call. = FALSE
      )
    }
    c_L <- measured_c_L
    c_H <- measured_c_H
  } else {
    c_L <- as.numeric(c_L)
    c_H <- as.numeric(c_H)
  }
  minimum_low <- .row_minimum(low_discrepancy)
  minimum_high <- .row_minimum(high_discrepancy)

  epsilon_selected <- is.null(epsilon)
  if (epsilon_selected) {
    epsilon <- .weighted_quantile(minimum_high, weights, kappa)
  }
  if (!is.null(epsilon_target)) {
    epsilon <- max(epsilon, epsilon_target)
  }

  high_match_count <- rowSums(high_discrepancy <= epsilon)
  high_pass <- high_match_count > 0L
  high_likelihood <- high_match_count / n_H
  high_posterior_weight <- weights * high_likelihood
  high_posterior_mass <- sum(high_posterior_weight)
  if (!is.finite(high_posterior_mass) || high_posterior_mass <= 0) {
    stop("The pilot produced no positive HF-ABC weight.", call. = FALSE)
  }
  high_posterior_weight <- high_posterior_weight / high_posterior_mass

  epsilon_L_selected <- is.null(epsilon_L)
  if (epsilon_L_selected) {
    epsilon_L <- .weighted_quantile(
      minimum_low, high_posterior_weight, 1 - a_L
    )
  }

  low_match_count <- rowSums(low_discrepancy <= epsilon_L)
  low_pass <- low_match_count > 0L
  false_positive <- !high_pass & low_pass
  false_negative <- high_pass & !low_pass
  true_positive <- high_pass & low_pass
  true_negative <- !high_pass & !low_pass

  high_reject_count <- sum(!high_pass)
  fpr <- if (high_reject_count > 0L) {
    sum(false_positive) / high_reject_count
  } else {
    NA_real_
  }
  high_reject_weight <- sum(weights[!high_pass])
  fpr_weighted <- if (high_reject_weight > 0) {
    sum(weights[false_positive]) / high_reject_weight
  } else {
    NA_real_
  }
  fnr_posterior <- sum(high_posterior_weight[false_negative])
  high_pass_count <- sum(high_pass)
  fnr_unweighted <- if (high_pass_count > 0L) {
    sum(false_negative) / high_pass_count
  } else {
    NA_real_
  }
  p_H <- mean(high_pass)
  fpr_max_denominator <- (1 - p_H) * n_H * c_H
  fpr_max <- if (fpr_max_denominator > 0) {
    1 - n_L * c_L / fpr_max_denominator
  } else {
    NA_real_
  }
  fpr_below_max <- is.finite(fpr) && is.finite(fpr_max) && fpr < fpr_max

  theta_columns <- as.data.frame(Theta, stringsAsFactors = FALSE)
  names(theta_columns) <- paste0("Theta", seq_len(ncol(Theta)))
  particles <- data.frame(
    particle = seq_len(pilot_size),
    theta_columns,
    weight = weights,
    min_low_discrepancy = minimum_low,
    min_high_discrepancy = minimum_high,
    low_match_count = low_match_count,
    high_match_count = high_match_count,
    high_likelihood = high_likelihood,
    high_posterior_weight = high_posterior_weight,
    low_pass = low_pass,
    high_pass = high_pass,
    true_positive = true_positive,
    false_positive = false_positive,
    true_negative = true_negative,
    false_negative = false_negative,
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    pilot_size = pilot_size,
    n_L = n_L,
    n_H = n_H,
    kappa = kappa,
    a_L = a_L,
    epsilon = epsilon,
    epsilon_L = epsilon_L,
    epsilon_selected = epsilon_selected,
    epsilon_L_selected = epsilon_L_selected,
    p_H = p_H,
    fpr = fpr,
    fpr_weighted = fpr_weighted,
    fpr_max = fpr_max,
    fpr_below_max = fpr_below_max,
    fnr_hf_posterior = fnr_posterior,
    fnr_unweighted = fnr_unweighted,
    high_acceptance_rate = sum(weights[high_pass]),
    low_acceptance_rate = sum(weights[low_pass]),
    high_posterior_ess = 1 / sum(high_posterior_weight^2),
    n_true_positive = sum(true_positive),
    n_false_positive = sum(false_positive),
    n_true_negative = sum(true_negative),
    n_false_negative = sum(false_negative),
    cost_mode = cost_mode,
    c_L = c_L,
    c_H = c_H,
    measured_c_L = measured_c_L,
    measured_c_H = measured_c_H,
    low_elapsed_seconds = low_timing$elapsed,
    high_elapsed_seconds = high_timing$elapsed,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message(sprintf(
      paste0("Pilot FPR: epsilon=%.6g, epsilon_L=%.6g, ",
             "FPR=%s, FPR_max=%s, posterior FNR=%.6g"),
      epsilon, epsilon_L,
      if (is.finite(fpr)) format(fpr, digits = 6L) else "NA",
      if (is.finite(fpr_max)) format(fpr_max, digits = 6L) else "NA",
      fnr_posterior
    ))
  }

  list(
    summary = summary,
    particles = particles,
    low_discrepancy = low_discrepancy,
    high_discrepancy = high_discrepancy
  )
}
