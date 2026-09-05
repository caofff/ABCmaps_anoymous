#' Adaptive sequential Monte Carlo ABC
#'
#' Runs ABC-ASMC with adaptive discrepancy thresholds, resampling, and
#' Metropolis-Hastings particle mutation.
#'
#' @param Prior A list with `sample(n)` and `density(theta)` functions.
#' @param Dis_fun A function with arguments `theta`, `n_rep`, and `type`.
#'   It must return an `nrow(theta)` by `n_rep` discrepancy matrix.
#' @param proposal_fun A proposal function accepting active particles, the
#'   complete particle matrix, and the current weights.
#' @param N Number of particles.
#' @param n_H Number of simulator replicates per evaluated particle.
#' @param alpha Adaptive threshold quantile.
#' @param epsilon_T Final ABC threshold.
#' @param gamma Resampling threshold as a fraction of `N`.
#' @param type Fidelity label passed to `Dis_fun`.
#' @param verbose Whether to report iteration diagnostics.
#' @param max_iter Maximum number of iterations.
#'
#' @return A list containing final particles and weights, threshold and
#'   simulation-count traces, effective sample sizes, elapsed times, and
#'   particle history.
#' @export
ABC_ASMC <- function(Prior, Dis_fun,
                     proposal_fun = maps_gaussian_proposal,
                     N, n_H, alpha, epsilon_T,
                     gamma = 0.5, type = "high",
                     verbose = FALSE, max_iter = 100L) {
  .validate_prior(Prior, require_density = TRUE)
  if (!is.function(Dis_fun)) stop("Dis_fun must be a function.", call. = FALSE)
  if (!is.function(proposal_fun)) {
    stop("proposal_fun must be a function.", call. = FALSE)
  }
  N <- .scalar_integer(N, "N")
  n_H <- .scalar_integer(n_H, "n_H")
  max_iter <- .scalar_integer(max_iter, "max_iter")
  alpha <- .scalar_probability(alpha, "alpha")
  gamma <- .scalar_probability(gamma, "gamma", include_zero = TRUE)
  if (!is.numeric(epsilon_T) || length(epsilon_T) != 1L ||
      !is.finite(epsilon_T)) {
    stop("epsilon_T must be a finite number.", call. = FALSE)
  }
  if (!is.character(type) || length(type) != 1L || !nzchar(type)) {
    stop("type must be a non-empty character string.", call. = FALSE)
  }

  run_start <- proc.time()[["elapsed"]]
  theta <- .sample_particles(Prior, N)
  initial_density <- .prior_density(Prior, theta)
  if (any(initial_density <= 0)) {
    stop("Prior$sample(N) generated particles with zero prior density.",
         call. = FALSE)
  }

  weights <- rep(1 / N, N)
  discrepancy <- .evaluate_discrepancy(Dis_fun, theta, n_H, type)
  evaluated_count <- N
  epsilon <- Inf
  iteration <- 0L
  epsilon_trace <- numeric()
  count_trace <- integer()
  ess_trace <- numeric()
  elapsed_trace <- numeric()

  make_history <- function(theta, weights, iteration) {
    result <- as.data.frame(theta, stringsAsFactors = FALSE)
    names(result) <- paste0("Theta", seq_len(ncol(theta)))
    result$weights <- as.numeric(weights)
    result$t <- as.integer(iteration)
    result
  }
  history <- list(make_history(theta, weights, 0L))

  while (epsilon > epsilon_T && iteration < max_iter) {
    iteration <- iteration + 1L
    epsilon_old <- epsilon
    active <- which(weights > 0)
    minimum_discrepancy <- .row_minimum(discrepancy)
    epsilon <- max(
      .quantile_value(minimum_discrepancy[active], alpha), epsilon_T
    )

    old_matches <- if (is.infinite(epsilon_old)) {
      rep(n_H, length(active))
    } else {
      rowSums(discrepancy[active, , drop = FALSE] <= epsilon_old)
    }
    new_matches <- rowSums(
      discrepancy[active, , drop = FALSE] <= epsilon
    )
    weights[active] <- weights[active] * new_matches / old_matches
    weights <- .normalize_weights(
      weights, paste0("after ABC-ASMC correction at iteration ", iteration)
    )

    if (.unique_particle_ess(theta, weights) < gamma * N) {
      index <- .systematic_resample(weights, N)
      theta <- theta[index, , drop = FALSE]
      discrepancy <- discrepancy[index, , drop = FALSE]
      weights <- rep(1 / N, N)
    }

    proposals <- .maps_proposals(theta, weights, Prior, proposal_fun)
    accepted_count <- 0L
    if (nrow(proposals$theta) > 0L) {
      proposed_discrepancy <- .evaluate_discrepancy(
        Dis_fun, proposals$theta, n_H, type
      )
      evaluated_count <- evaluated_count + nrow(proposals$theta)
      proposed_matches <- rowSums(proposed_discrepancy <= epsilon)
      current <- proposals$current_index
      current_matches <- rowSums(
        discrepancy[current, , drop = FALSE] <= epsilon
      )
      current_density <- .prior_density(
        Prior, theta[current, , drop = FALSE]
      )
      probability <- pmin(
        1,
        proposals$density * proposed_matches /
          (current_density * current_matches)
      )
      probability[!is.finite(probability)] <- 0
      accepted <- which(stats::runif(length(probability)) < probability)
      if (length(accepted) > 0L) {
        accepted_current <- current[accepted]
        theta[accepted_current, ] <-
          proposals$theta[accepted, , drop = FALSE]
        discrepancy[accepted_current, ] <-
          proposed_discrepancy[accepted, , drop = FALSE]
        accepted_count <- length(accepted)
      }
    }

    epsilon_trace[iteration] <- epsilon
    count_trace[iteration] <- evaluated_count
    ess_trace[iteration] <- .unique_particle_ess(theta, weights)
    elapsed_trace[iteration] <- proc.time()[["elapsed"]] - run_start
    history[[length(history) + 1L]] <-
      make_history(theta, weights, iteration)

    if (verbose) {
      message(sprintf(
        "ABC-ASMC iteration %d: epsilon=%.6g, ESS=%.2f, accepted=%d",
        iteration, epsilon, ess_trace[iteration], accepted_count
      ))
    }
  }

  if (epsilon > epsilon_T) {
    stop("ABC-ASMC did not reach epsilon_T before max_iter.", call. = FALSE)
  }

  particle_history <- do.call(rbind, history)
  rownames(particle_history) <- NULL

  list(
    smc_Data_ite = particle_history,
    ESS = ess_trace,
    Eps = epsilon_trace,
    n_high = count_trace,
    elapsed = elapsed_trace,
    Theta = theta,
    weights = weights
  )
}
