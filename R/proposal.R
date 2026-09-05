#' Gaussian random-walk proposal for MAPS
#'
#' Constructs a Gaussian random-walk proposal with covariance twice the
#' weighted covariance of the current particle population.
#'
#' @param theta Active particles to perturb, as a matrix.
#' @param Theta Complete particle population, as a matrix.
#' @param weights Non-negative particle weights.
#'
#' @return A matrix with the same dimensions as `theta`.
#' @export
maps_gaussian_proposal <- function(theta, Theta, weights) {
  theta <- as.matrix(theta)
  Theta <- as.matrix(Theta)
  weights <- as.numeric(weights)

  if (ncol(theta) != ncol(Theta) || length(weights) != nrow(Theta) ||
      any(!is.finite(theta)) || any(!is.finite(Theta)) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("Invalid particles or weights supplied to maps_gaussian_proposal().",
         call. = FALSE)
  }

  weights <- weights / sum(weights)
  center <- colSums(Theta * weights)
  centered <- sweep(Theta, 2L, center, "-")
  denominator <- 1 - sum(weights^2)
  if (denominator > sqrt(.Machine$double.eps)) {
    covariance <- crossprod(centered * sqrt(weights)) / denominator
  } else {
    covariance <- matrix(0, nrow = ncol(Theta), ncol = ncol(Theta))
  }

  kernel_covariance <- 2 * covariance
  kernel_covariance <- (kernel_covariance + t(kernel_covariance)) / 2
  decomposition <- eigen(kernel_covariance, symmetric = TRUE)
  eigenvalues <- pmax(decomposition$values, 1e-10)
  transform <- diag(sqrt(eigenvalues), nrow = length(eigenvalues)) %*%
    t(decomposition$vectors)
  increments <- matrix(
    stats::rnorm(nrow(theta) * ncol(theta)),
    nrow = nrow(theta), ncol = ncol(theta)
  ) %*% transform

  theta + increments
}
