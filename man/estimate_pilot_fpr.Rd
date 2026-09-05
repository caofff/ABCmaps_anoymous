\name{estimate_pilot_fpr}
\alias{estimate_pilot_fpr}
\title{Estimate a Low-Fidelity False-Positive Rate from Pilot Simulations}
\usage{
estimate_pilot_fpr(
  Theta,
  Dis_fun,
  n_L,
  n_H,
  weights = NULL,
  kappa = 0.1,
  epsilon_target = NULL,
  epsilon = NULL,
  epsilon_L = NULL,
  a_L = 0.05,
  cost_mode = c("fixed", "measured"),
  c_L = 1,
  c_H = 10,
  verbose = FALSE
)
}
\arguments{
  \item{Theta}{Pilot particles as a numeric vector for a one-dimensional
  parameter or an \code{N} by \code{d} numeric matrix.}
  \item{Dis_fun}{A discrepancy function accepting \code{theta},
  \code{n_rep}, and \code{type}.}
  \item{n_L}{Low-fidelity replicates per pilot particle.}
  \item{n_H}{High-fidelity replicates per pilot particle.}
  \item{weights}{Optional non-negative pilot weights. Uniform weights are used
  by default.}
  \item{kappa}{Quantile used to select \code{epsilon} when it is not supplied.}
  \item{epsilon_target}{Optional lower bound for the selected high-fidelity
  threshold.}
  \item{epsilon}{Optional fixed high-fidelity threshold.}
  \item{epsilon_L}{Optional fixed low-fidelity threshold.}
  \item{a_L}{Target posterior-weighted false-negative rate used when selecting
  \code{epsilon_L}.}
  \item{cost_mode}{Cost source for \code{fpr_max}. Use \code{"fixed"} for
  supplied costs or \code{"measured"} for pilot wall-clock costs.}
  \item{c_L}{Average cost of one low-fidelity simulation in fixed-cost mode.}
  \item{c_H}{Average cost of one high-fidelity simulation in fixed-cost mode.}
  \item{verbose}{Whether to report selected thresholds and estimated rates.}
}
\value{
A list containing a one-row \code{summary} data frame, a particle-level
\code{particles} data frame, and the two pilot discrepancy matrices.
}
\description{
Evaluates low- and high-fidelity ABC decisions on a shared pilot population.
The false-positive rate is the fraction of high-fidelity-rejected pilot
particles that pass the low-fidelity screen. The function also reports the
cost-based upper bound \code{fpr_max} and the high-fidelity-posterior-weighted
false-negative rate.
}
