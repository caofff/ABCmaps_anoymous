library(ABCmaps)

set.seed(7)

prior <- list(
  sample = function(n) matrix(runif(n, -1, 1), ncol = 1L),
  density = function(theta) dunif(as.numeric(theta), -1, 1)
)

discrepancy <- function(theta, n_rep, type = "high") {
  theta <- as.numeric(theta)
  center <- if (type == "low") 0.20 else 0.25
  matrix(rep((theta - center)^2, n_rep), nrow = length(theta))
}

maps_fit <- MAPS(
  prior, discrepancy, N = 100, n_L = 1, n_H = 1,
  alpha = 0.5, alpha_L = 0.8, epsilon_T = 0.5,
  epsilon_L_T = 0.8, gamma = 0.5
)
stopifnot(
  is.matrix(maps_fit$Theta),
  nrow(maps_fit$Theta) == 100L,
  abs(sum(maps_fit$weights) - 1) < 1e-12,
  tail(maps_fit$Eps, 1L) == 0.5
)

asmc_fit <- ABC_ASMC(
  prior, discrepancy, N = 100, n_H = 1,
  alpha = 0.5, epsilon_T = 0.5, gamma = 0.5
)
stopifnot(
  is.matrix(asmc_fit$Theta),
  nrow(asmc_fit$Theta) == 100L,
  abs(sum(asmc_fit$weights) - 1) < 1e-12,
  tail(asmc_fit$Eps, 1L) == 0.5,
  tail(asmc_fit$n_high, 1L) >= 100L
)

pfis_fit <- PrefilterIS(
  prior, discrepancy, N = 200, n_L = 1, n_H = 1,
  epsilon = 0.5, epsilon_L = 0.8, chunk_size = 31
)
stopifnot(
  is.matrix(pfis_fit$Theta),
  nrow(pfis_fit$Theta) == length(pfis_fit$weights),
  abs(sum(pfis_fit$weights) - 1) < 1e-12,
  pfis_fit$low_simulations == 200L,
  pfis_fit$high_simulations <= 200L
)

pilot_theta <- matrix(seq_len(4L), ncol = 1L)
pilot_high <- rbind(
  c(0.1, 0.2, 0.3, 0.4),
  c(0.1, 1.0, 1.0, 1.0),
  rep(1.0, 4L),
  rep(1.0, 4L)
)
pilot_low <- matrix(c(1.0, 10.0, 0.5, 3.0), ncol = 1L)

pilot_discrepancy <- function(theta, n_rep, type = "high") {
  index <- as.integer(theta[, 1L])
  if (type == "high") {
    pilot_high[index, , drop = FALSE]
  } else {
    pilot_low[index, , drop = FALSE]
  }
}

pilot_fit <- estimate_pilot_fpr(
  Theta = pilot_theta,
  Dis_fun = pilot_discrepancy,
  n_L = 1,
  n_H = 4,
  weights = c(0.1, 0.2, 0.3, 0.4),
  kappa = 0.5,
  epsilon = 0.5,
  a_L = 0.4,
  cost_mode = "fixed",
  c_L = 1,
  c_H = 10
)

stopifnot(
  isTRUE(all.equal(pilot_fit$summary$epsilon_L, 1)),
  isTRUE(all.equal(pilot_fit$summary$p_H, 1 / 2)),
  isTRUE(all.equal(pilot_fit$summary$fpr, 1 / 2)),
  isTRUE(all.equal(pilot_fit$summary$fpr_weighted, 3 / 7)),
  isTRUE(all.equal(pilot_fit$summary$fpr_max, 0.95)),
  isTRUE(all.equal(pilot_fit$summary$fnr_hf_posterior, 1 / 3)),
  isTRUE(pilot_fit$summary$fpr_below_max),
  identical(pilot_fit$summary$cost_mode, "fixed"),
  isTRUE(all.equal(sum(pilot_fit$particles$high_posterior_weight), 1))
)

timed_pilot_discrepancy <- function(theta, n_rep, type = "high") {
  Sys.sleep(if (type == "low") 0.002 * n_rep else 0.004 * n_rep)
  pilot_discrepancy(theta, n_rep, type)
}

timed_pilot_fit <- estimate_pilot_fpr(
  Theta = pilot_theta,
  Dis_fun = timed_pilot_discrepancy,
  n_L = 1,
  n_H = 4,
  weights = c(0.1, 0.2, 0.3, 0.4),
  kappa = 0.5,
  epsilon = 0.5,
  a_L = 0.4,
  cost_mode = "measured"
)

stopifnot(
  identical(timed_pilot_fit$summary$cost_mode, "measured"),
  timed_pilot_fit$summary$c_L > 0,
  timed_pilot_fit$summary$c_H > 0,
  isTRUE(all.equal(
    timed_pilot_fit$summary$c_L,
    timed_pilot_fit$summary$measured_c_L
  )),
  isTRUE(all.equal(
    timed_pilot_fit$summary$c_H,
    timed_pilot_fit$summary$measured_c_H
  )),
  is.finite(timed_pilot_fit$summary$fpr_max)
)

toy_environment <- new.env(parent = globalenv())
example_directory <- system.file("examples", package = "ABCmaps")
toy_environment$N <- 200
toy_environment$pilot_N <- 200
toy_environment$num_rep <- 1
toy_environment$output_dir <- tempfile("maps-toy-output-")
source(
  file.path(example_directory, "toy-example.R"),
  local = toy_environment
)
stopifnot(
  nrow(toy_environment$comparison) == 2L,
  identical(
    toy_environment$comparison$method,
    c("ABC-ASMC", "MAPS")
  ),
  all(is.finite(toy_environment$comparison$ISE)),
  all(is.finite(toy_environment$comparison$time))
)
