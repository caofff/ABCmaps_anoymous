.scalar_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum) {
    stop(name, " must be an integer greater than or equal to ", minimum, ".",
         call. = FALSE)
  }
  as.integer(x)
}

.scalar_probability <- function(x, name, include_zero = FALSE) {
  lower_ok <- if (include_zero) x >= 0 else x > 0
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      !lower_ok || x >= 1) {
    interval <- if (include_zero) "[0, 1)" else "(0, 1)"
    stop(name, " must lie in ", interval, ".", call. = FALSE)
  }
  as.numeric(x)
}

.validate_prior <- function(prior, require_density = TRUE) {
  if (!is.list(prior) || !is.function(prior$sample)) {
    stop("Prior must be a list containing a sample function.", call. = FALSE)
  }
  if (require_density && !is.function(prior$density)) {
    stop("Prior must contain a density function for MAPS.", call. = FALSE)
  }
  invisible(prior)
}

.sample_particles <- function(prior, n) {
  theta <- prior$sample(n)
  if (is.null(dim(theta))) {
    if (length(theta) != n) {
      stop("Prior$sample(N) returned an unexpected number of values.",
           call. = FALSE)
    }
    theta <- matrix(theta, nrow = n, ncol = 1L)
  } else {
    theta <- as.matrix(theta)
  }
  if (nrow(theta) != n || ncol(theta) < 1L || !is.numeric(theta) ||
      any(!is.finite(theta))) {
    stop("Prior$sample(N) must return a finite numeric N by d matrix.",
         call. = FALSE)
  }
  theta
}

.proposal_matrix <- function(theta, n, dimension) {
  if (is.null(dim(theta))) {
    if (length(theta) != n * dimension) {
      stop("proposal_fun returned an unexpected number of values.",
           call. = FALSE)
    }
    theta <- matrix(theta, nrow = n, ncol = dimension)
  } else {
    theta <- as.matrix(theta)
  }
  if (!identical(dim(theta), c(n, dimension)) || !is.numeric(theta) ||
      any(!is.finite(theta))) {
    stop("proposal_fun must return a finite numeric matrix with one row per input particle.",
         call. = FALSE)
  }
  theta
}

.discrepancy_matrix <- function(discrepancy, n, n_rep, label) {
  if (is.null(dim(discrepancy))) {
    if (length(discrepancy) != n * n_rep) {
      stop(label, " returned an unexpected number of discrepancies.",
           call. = FALSE)
    }
    discrepancy <- matrix(discrepancy, nrow = n, ncol = n_rep)
  } else {
    discrepancy <- as.matrix(discrepancy)
  }
  if (!identical(dim(discrepancy), c(n, n_rep)) ||
      !is.numeric(discrepancy) || any(!is.finite(discrepancy))) {
    stop(label, " must return a finite numeric nrow(theta) by n_rep matrix.",
         call. = FALSE)
  }
  discrepancy
}

.evaluate_discrepancy <- function(dis_fun, theta, n_rep, type) {
  n <- nrow(theta)
  if (n == 0L) return(matrix(numeric(), nrow = 0L, ncol = n_rep))
  value <- dis_fun(theta, n_rep, type = type)
  .discrepancy_matrix(value, n, n_rep, paste0("Dis_fun(type = '", type, "')"))
}

.prior_density <- function(prior, theta) {
  density <- as.numeric(prior$density(theta))
  if (length(density) != nrow(theta) || any(!is.finite(density)) ||
      any(density < 0)) {
    stop("Prior$density(theta) must return one finite, non-negative value per row.",
         call. = FALSE)
  }
  density
}

.normalize_weights <- function(weights, context) {
  weights <- as.numeric(weights)
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Invalid particle weights ", context, ".", call. = FALSE)
  }
  total <- sum(weights)
  if (!is.finite(total) || total <= 0) {
    stop("All particle weights are zero ", context, ".", call. = FALSE)
  }
  weights / total
}

.row_minimum <- function(x) {
  if (nrow(x) == 0L) return(numeric())
  apply(x, 1L, min)
}

.quantile_value <- function(x, probability) {
  unname(stats::quantile(x, probs = probability, names = FALSE))
}

.weighted_quantile <- function(x, weights, probability) {
  valid <- is.finite(x) & is.finite(weights) & weights > 0
  x <- as.numeric(x[valid])
  weights <- as.numeric(weights[valid])
  if (length(x) == 0L) {
    stop("The weighted quantile has no positive-weight values.", call. = FALSE)
  }

  ordering <- order(x)
  x <- x[ordering]
  weights <- weights[ordering] / sum(weights)
  cumulative_weight <- cumsum(weights)
  cumulative_weight[length(cumulative_weight)] <- 1
  x[which(cumulative_weight >= probability)[1L]]
}

.unique_particle_ess <- function(theta, weights) {
  keys <- apply(theta, 1L, function(row) {
    paste(format(row, digits = 17L, scientific = TRUE, trim = TRUE),
          collapse = "\r")
  })
  grouped <- rowsum(as.numeric(weights), keys, reorder = FALSE)[, 1L]
  grouped <- grouped / sum(grouped)
  1 / sum(grouped^2)
}

.systematic_resample <- function(weights, n) {
  positions <- stats::runif(1L, min = 0, max = 1 / n) +
    (seq_len(n) - 1L) / n
  pmin(findInterval(positions, cumsum(weights)) + 1L, n)
}

.history_frame <- function(theta, weights, iteration, fidelity) {
  result <- as.data.frame(theta, stringsAsFactors = FALSE)
  names(result) <- paste0("Theta", seq_len(ncol(theta)))
  result$weights <- as.numeric(weights)
  result$t <- as.integer(iteration)
  result$type <- fidelity
  result
}

.maps_proposals <- function(theta, weights, prior, proposal_fun) {
  active <- which(weights > 0)
  proposed <- proposal_fun(
    theta[active, , drop = FALSE], theta, weights
  )
  proposed <- .proposal_matrix(proposed, length(active), ncol(theta))
  proposed_density <- .prior_density(prior, proposed)
  valid <- proposed_density > 0

  list(
    theta = proposed[valid, , drop = FALSE],
    current_index = active[valid],
    density = proposed_density[valid]
  )
}
