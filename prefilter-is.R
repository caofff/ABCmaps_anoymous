#' Prior-proposal prefiltered importance sampling
#'
#' Draws candidates from the prior, screens them with a low-fidelity
#' discrepancy, and evaluates the high-fidelity discrepancy only for survivors.
#'
#' @param Prior A list containing a `sample(n)` function.
#' @param Dis_fun A function with arguments `theta`, `n_rep`, and `type`.
#'   It must return an `nrow(theta)` by `n_rep` discrepancy matrix.
#' @param N Number of candidates sampled from the prior.
#' @param n_L Number of low-fidelity replicates per candidate.
#' @param n_H Number of high-fidelity replicates per surviving candidate.
#' @param epsilon High-fidelity ABC threshold.
#' @param epsilon_L Low-fidelity screening threshold.
#' @param chunk_size Number of candidates processed per chunk.
#' @param verbose Whether to report chunk progress.
#'
#' @return A list containing positive-weight particles, normalized weights,
#'   thresholds, and simulation counts.
#' @export
PrefilterIS <- function(Prior, Dis_fun, N, n_L, n_H,
                        epsilon, epsilon_L,
                        chunk_size = 5000L,
                        verbose = FALSE) {
  .validate_prior(Prior, require_density = FALSE)
  if (!is.function(Dis_fun)) stop("Dis_fun must be a function.", call. = FALSE)
  N <- .scalar_integer(N, "N")
  n_L <- .scalar_integer(n_L, "n_L")
  n_H <- .scalar_integer(n_H, "n_H")
  chunk_size <- .scalar_integer(chunk_size, "chunk_size")
  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || !is.numeric(epsilon_L) ||
      length(epsilon_L) != 1L || !is.finite(epsilon_L)) {
    stop("epsilon and epsilon_L must be finite numbers.", call. = FALSE)
  }

  theta <- .sample_particles(Prior, N)
  unnormalized_weights <- numeric(N)
  low_pass <- logical(N)
  high_evaluated <- 0L
  starts <- seq.int(1L, N, by = chunk_size)

  for (chunk_index in seq_along(starts)) {
    first <- starts[chunk_index]
    last <- min(first + chunk_size - 1L, N)
    index <- first:last
    theta_chunk <- theta[index, , drop = FALSE]
    low_discrepancy <- .evaluate_discrepancy(
      Dis_fun, theta_chunk, n_L, "low"
    )
    pass <- rowSums(low_discrepancy <= epsilon_L) > 0L
    low_pass[index] <- pass

    if (any(pass)) {
      survivors <- index[pass]
      high_discrepancy <- .evaluate_discrepancy(
        Dis_fun, theta[survivors, , drop = FALSE], n_H, "high"
      )
      high_evaluated <- high_evaluated + length(survivors)
      unnormalized_weights[survivors] <- rowSums(
        high_discrepancy <= epsilon
      )
    }

    if (verbose) {
      message(sprintf(
        "PF-IS chunk %d/%d: candidates=%d, LF survivors=%d",
        chunk_index, length(starts), length(index), sum(pass)
      ))
    }
  }

  positive <- is.finite(unnormalized_weights) & unnormalized_weights > 0
  if (!any(positive)) {
    stop("PF-IS produced no positive-weight particles.", call. = FALSE)
  }
  weights <- unnormalized_weights[positive]
  weights <- weights / sum(weights)

  list(
    Theta = theta[positive, , drop = FALSE],
    weights = weights,
    candidate_count = N,
    lf_survivor_count = sum(low_pass),
    lf_survival_rate = mean(low_pass),
    positive_weight_count = sum(positive),
    positive_weight_rate = mean(positive),
    low_simulations = N * n_L,
    high_simulations = high_evaluated * n_H,
    n_L = n_L,
    n_H = n_H,
    epsilon = as.numeric(epsilon),
    epsilon_L = as.numeric(epsilon_L)
  )
}
