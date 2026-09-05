#' Multifidelity adaptive particle sampling
#'
#' Runs MAPS with a low-fidelity screening stage and a high-fidelity ABC
#' correction stage.
#'
#' @param Prior A list with `sample(n)` and `density(theta)` functions.
#' @param Dis_fun A function with arguments `theta`, `n_rep`, and `type`.
#'   It must return an `nrow(theta)` by `n_rep` discrepancy matrix.
#' @param proposal_fun A proposal function accepting active particles, the
#'   complete particle matrix, and the current weights.
#' @param N Number of particles.
#' @param n_L Number of low-fidelity replicates per evaluated particle.
#' @param n_H Number of high-fidelity replicates per evaluated particle.
#' @param alpha High-fidelity threshold quantile.
#' @param alpha_L Low-fidelity threshold quantile.
#' @param epsilon_T Final high-fidelity ABC threshold.
#' @param epsilon_L_T Fixed lower bound for the low-fidelity threshold. Use
#'   `-Inf` to disable the lower bound.
#' @param gamma Resampling threshold as a fraction of `N`.
#' @param verbose Whether to report iteration diagnostics.
#' @param max_iter Maximum number of MAPS iterations.
#'
#' @return A list containing final particles and weights, threshold and cost
#'   traces, effective sample sizes, and particle history.
#' @export
MAPS <- function(Prior, Dis_fun,
                 proposal_fun = maps_gaussian_proposal,
                 N, n_L, n_H, alpha, alpha_L,
                 epsilon_T, epsilon_L_T = -Inf,
                 gamma = 0.5, verbose = FALSE, max_iter = 100L) {
  .validate_prior(Prior, require_density = TRUE)
  if (!is.function(Dis_fun)) stop("Dis_fun must be a function.", call. = FALSE)
  if (!is.function(proposal_fun)) {
    stop("proposal_fun must be a function.", call. = FALSE)
  }

  N <- .scalar_integer(N, "N")
  n_L <- .scalar_integer(n_L, "n_L")
  n_H <- .scalar_integer(n_H, "n_H")
  max_iter <- .scalar_integer(max_iter, "max_iter")
  alpha <- .scalar_probability(alpha, "alpha")
  alpha_L <- .scalar_probability(alpha_L, "alpha_L")
  gamma <- .scalar_probability(gamma, "gamma", include_zero = TRUE)
  if (!is.numeric(epsilon_T) || length(epsilon_T) != 1L ||
      !is.finite(epsilon_T)) {
    stop("epsilon_T must be a finite number.", call. = FALSE)
  }
  epsilon_T <- as.numeric(epsilon_T)
  if (!is.numeric(epsilon_L_T) || length(epsilon_L_T) != 1L ||
      is.na(epsilon_L_T) || epsilon_L_T == Inf) {
    stop("epsilon_L_T must be finite or -Inf.", call. = FALSE)
  }
  epsilon_L_T <- as.numeric(epsilon_L_T)

  run_start <- proc.time()[["elapsed"]]
  theta <- .sample_particles(Prior, N)
  initial_density <- .prior_density(Prior, theta)
  if (any(initial_density <= 0)) {
    stop("Prior$sample(N) generated particles with zero prior density.",
         call. = FALSE)
  }

  weights <- rep(1 / N, N)
  history <- list(.history_frame(theta, weights, 0L, "High"))
  epsilon_trace <- numeric()
  epsilon_low_trace <- numeric()
  low_count_trace <- integer()
  high_count_trace <- integer()
  ess_trace <- numeric()
  elapsed_trace <- numeric()
  low_count <- 0L
  high_count <- 0L

  low_discrepancy <- .evaluate_discrepancy(Dis_fun, theta, n_L, "low")
  low_count <- low_count + N
  minimum_low <- .row_minimum(low_discrepancy)
  epsilon_low <- max(
    .quantile_value(minimum_low, alpha_L), epsilon_L_T
  )
  weights <- .normalize_weights(
    weights * (minimum_low <= epsilon_low),
    "after initial low-fidelity screening"
  )

  if (.unique_particle_ess(theta, weights) < gamma * N) {
    index <- .systematic_resample(weights, N)
    theta <- theta[index, , drop = FALSE]
    minimum_low <- minimum_low[index]
    weights <- rep(1 / N, N)
  }

  proposals <- .maps_proposals(theta, weights, Prior, proposal_fun)
  if (nrow(proposals$theta) > 0L) {
    proposed_low <- .evaluate_discrepancy(
      Dis_fun, proposals$theta, n_L, "low"
    )
    low_count <- low_count + nrow(proposals$theta)
    proposed_minimum_low <- .row_minimum(proposed_low)
    current_density <- .prior_density(
      Prior, theta[proposals$current_index, , drop = FALSE]
    )
    probability <- pmin(1, proposals$density / current_density) *
      (proposed_minimum_low <= epsilon_low)
    accepted <- which(stats::runif(length(probability)) < probability)
    if (length(accepted) > 0L) {
      current <- proposals$current_index[accepted]
      theta[current, ] <- proposals$theta[accepted, , drop = FALSE]
      minimum_low[current] <- proposed_minimum_low[accepted]
    }
  }
  history[[length(history) + 1L]] <-
    .history_frame(theta, weights, 1L, "Low")

  active <- which(weights > 0)
  high_discrepancy <- matrix(Inf, nrow = N, ncol = n_H)
  high_discrepancy[active, ] <- .evaluate_discrepancy(
    Dis_fun, theta[active, , drop = FALSE], n_H, "high"
  )
  high_count <- high_count + length(active)
  minimum_high <- .row_minimum(high_discrepancy)
  epsilon <- max(
    .quantile_value(minimum_high[active], alpha), epsilon_T
  )
  high_matches <- rowSums(high_discrepancy[active, , drop = FALSE] <= epsilon)
  weights[active] <- weights[active] * high_matches / n_H
  weights <- .normalize_weights(
    weights, "after initial high-fidelity correction"
  )

  epsilon_trace[1L] <- epsilon
  epsilon_low_trace[1L] <- epsilon_low
  low_count_trace[1L] <- low_count
  high_count_trace[1L] <- high_count
  ess_trace[1L] <- .unique_particle_ess(theta, weights)
  elapsed_trace[1L] <- proc.time()[["elapsed"]] - run_start
  history[[length(history) + 1L]] <-
    .history_frame(theta, weights, 1L, "High")

  if (verbose) {
    message(sprintf(
      "MAPS iteration 1: epsilon=%.6g, epsilon_L=%.6g, ESS=%.2f",
      epsilon, epsilon_low, ess_trace[1L]
    ))
  }

  iteration <- 1L
  while (epsilon > epsilon_T && iteration < max_iter) {
    epsilon_old <- epsilon
    active <- which(weights > 0)
    iteration <- iteration + 1L

    epsilon_low <- max(
      .quantile_value(minimum_low[active], alpha_L), epsilon_L_T
    )
    weights <- .normalize_weights(
      weights * (minimum_low <= epsilon_low),
      paste0("after low-fidelity screening at iteration ", iteration)
    )

    if (.unique_particle_ess(theta, weights) < gamma * N) {
      index <- .systematic_resample(weights, N)
      theta <- theta[index, , drop = FALSE]
      high_discrepancy <- high_discrepancy[index, , drop = FALSE]
      minimum_high <- minimum_high[index]
      minimum_low <- minimum_low[index]
      weights <- rep(1 / N, N)
    }

    active <- which(weights > 0)
    proposals <- .maps_proposals(theta, weights, Prior, proposal_fun)
    low_pass_count <- 0L
    accepted_count <- 0L

    if (nrow(proposals$theta) > 0L) {
      proposed_low <- .evaluate_discrepancy(
        Dis_fun, proposals$theta, n_L, "low"
      )
      low_count <- low_count + nrow(proposals$theta)
      proposed_minimum_low <- .row_minimum(proposed_low)
      low_pass <- which(proposed_minimum_low <= epsilon_low)
      low_pass_count <- length(low_pass)

      if (low_pass_count > 0L) {
        proposed_high <- .evaluate_discrepancy(
          Dis_fun, proposals$theta[low_pass, , drop = FALSE], n_H, "high"
        )
        high_count <- high_count + low_pass_count
        proposed_old_matches <- rowSums(proposed_high <= epsilon_old)
        current <- proposals$current_index[low_pass]
        current_old_matches <- rowSums(
          high_discrepancy[current, , drop = FALSE] <= epsilon_old
        )
        current_density <- .prior_density(
          Prior, theta[current, , drop = FALSE]
        )
        probability <- pmin(
          1,
          proposals$density[low_pass] * proposed_old_matches /
            (current_density * current_old_matches)
        )
        probability[!is.finite(probability)] <- 0
        accepted_in_low <- which(
          stats::runif(length(probability)) < probability
        )

        if (length(accepted_in_low) > 0L) {
          accepted_proposal <- low_pass[accepted_in_low]
          accepted_current <- proposals$current_index[accepted_proposal]
          theta[accepted_current, ] <-
            proposals$theta[accepted_proposal, , drop = FALSE]
          high_discrepancy[accepted_current, ] <-
            proposed_high[accepted_in_low, , drop = FALSE]
          minimum_high[accepted_current] <-
            .row_minimum(proposed_high)[accepted_in_low]
          minimum_low[accepted_current] <-
            proposed_minimum_low[accepted_proposal]
          accepted_count <- length(accepted_current)
        }
      }
    }

    history[[length(history) + 1L]] <-
      .history_frame(theta, weights, iteration, "Low")

    active <- which(weights > 0)
    epsilon <- max(
      .quantile_value(minimum_high[active], alpha), epsilon_T
    )
    old_matches <- rowSums(
      high_discrepancy[active, , drop = FALSE] <= epsilon_old
    )
    new_matches <- rowSums(
      high_discrepancy[active, , drop = FALSE] <= epsilon
    )
    weights[active] <- weights[active] * new_matches / old_matches
    weights <- .normalize_weights(
      weights, paste0("after high-fidelity correction at iteration ", iteration)
    )

    epsilon_trace[iteration] <- epsilon
    epsilon_low_trace[iteration] <- epsilon_low
    low_count_trace[iteration] <- low_count
    high_count_trace[iteration] <- high_count
    ess_trace[iteration] <- .unique_particle_ess(theta, weights)
    elapsed_trace[iteration] <- proc.time()[["elapsed"]] - run_start
    history[[length(history) + 1L]] <-
      .history_frame(theta, weights, iteration, "High")

    if (verbose) {
      message(sprintf(
        paste0("MAPS iteration %d: epsilon=%.6g, epsilon_L=%.6g, ",
               "ESS=%.2f, LF-pass=%d, accepted=%d"),
        iteration, epsilon, epsilon_low, ess_trace[iteration],
        low_pass_count, accepted_count
      ))
    }
  }

  if (epsilon > epsilon_T) {
    stop("MAPS did not reach epsilon_T before max_iter.", call. = FALSE)
  }

  particle_history <- do.call(rbind, history)
  rownames(particle_history) <- NULL
  particle_history$type <- factor(
    particle_history$type, levels = c("Low", "High")
  )
  particle_history$group <- paste0(
    particle_history$type, particle_history$t
  )

  list(
    MF_Data_ite = particle_history,
    ESS = ess_trace,
    Eps = epsilon_trace,
    Eps_L = epsilon_low_trace,
    epsilon_L_T = epsilon_L_T,
    epsilon_L_floor_fixed = is.finite(epsilon_L_T),
    n_low = low_count_trace,
    n_high = high_count_trace,
    elapsed = elapsed_trace,
    Theta = theta,
    weights = weights
  )
}
