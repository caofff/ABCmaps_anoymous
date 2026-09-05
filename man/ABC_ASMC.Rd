\name{ABC_ASMC}
\alias{ABC_ASMC}
\title{Adaptive Sequential Monte Carlo ABC}
\usage{
ABC_ASMC(
  Prior,
  Dis_fun,
  proposal_fun = maps_gaussian_proposal,
  N,
  n_H,
  alpha,
  epsilon_T,
  gamma = 0.5,
  type = "high",
  verbose = FALSE,
  max_iter = 100L
)
}
\arguments{
  \item{Prior}{A list with \code{sample(n)} and \code{density(theta)}
  functions.}
  \item{Dis_fun}{A discrepancy function accepting \code{theta},
  \code{n_rep}, and \code{type}.}
  \item{proposal_fun}{A function that perturbs active particles.}
  \item{N}{Number of particles.}
  \item{n_H}{Simulator replicates per evaluated particle.}
  \item{alpha}{Adaptive threshold quantile in \code{(0, 1)}.}
  \item{epsilon_T}{Final ABC threshold.}
  \item{gamma}{Resampling threshold as a fraction of \code{N}.}
  \item{type}{Fidelity label passed to \code{Dis_fun}.}
  \item{verbose}{Whether to report iteration diagnostics.}
  \item{max_iter}{Maximum number of iterations.}
}
\value{
A list with final particles in \code{Theta}, normalized \code{weights},
thresholds in \code{Eps}, cumulative evaluated-particle counts in
\code{n_high}, effective sample sizes, elapsed times, and particle history in
\code{smc_Data_ite}.
}
\description{
Runs adaptive sequential Monte Carlo ABC with quantile-based thresholds,
effective-sample-size resampling, and Metropolis-Hastings mutation.
}
