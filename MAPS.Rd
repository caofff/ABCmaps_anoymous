\name{MAPS}
\alias{MAPS}
\title{Multifidelity Adaptive Particle Sampling}
\usage{
MAPS(
  Prior,
  Dis_fun,
  proposal_fun = maps_gaussian_proposal,
  N,
  n_L,
  n_H,
  alpha,
  alpha_L,
  epsilon_T,
  epsilon_L_T = -Inf,
  gamma = 0.5,
  verbose = FALSE,
  max_iter = 100L
)
}
\arguments{
  \item{Prior}{A list with \code{sample(n)} and \code{density(theta)} functions.}
  \item{Dis_fun}{A discrepancy function accepting \code{theta}, \code{n_rep}, and \code{type}.}
  \item{proposal_fun}{A function that perturbs active particles.}
  \item{N}{Number of particles.}
  \item{n_L}{Low-fidelity replicates per evaluated particle.}
  \item{n_H}{High-fidelity replicates per evaluated particle.}
  \item{alpha}{High-fidelity threshold quantile.}
  \item{alpha_L}{Low-fidelity threshold quantile.}
  \item{epsilon_T}{Final high-fidelity ABC threshold.}
  \item{epsilon_L_T}{Fixed lower bound for the low-fidelity threshold.}
  \item{gamma}{Resampling threshold as a fraction of \code{N}.}
  \item{verbose}{Whether to report iteration diagnostics.}
  \item{max_iter}{Maximum number of iterations.}
}
\value{
A list with final particles in \code{Theta}, normalized \code{weights},
threshold and simulation-count traces, effective sample sizes, elapsed times,
and particle history in \code{MF_Data_ite}.
}
\description{
Runs MAPS using low-fidelity screening followed by high-fidelity ABC
correction.
}
