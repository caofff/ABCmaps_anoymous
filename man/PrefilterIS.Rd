\name{PrefilterIS}
\alias{PrefilterIS}
\title{Prior-Proposal Prefiltered Importance Sampling}
\usage{
PrefilterIS(
  Prior,
  Dis_fun,
  N,
  n_L,
  n_H,
  epsilon,
  epsilon_L,
  chunk_size = 5000L,
  verbose = FALSE
)
}
\arguments{
  \item{Prior}{A list containing a \code{sample(n)} function.}
  \item{Dis_fun}{A discrepancy function accepting \code{theta}, \code{n_rep}, and \code{type}.}
  \item{N}{Number of prior candidates.}
  \item{n_L}{Low-fidelity replicates per candidate.}
  \item{n_H}{High-fidelity replicates per surviving candidate.}
  \item{epsilon}{High-fidelity ABC threshold.}
  \item{epsilon_L}{Low-fidelity screening threshold.}
  \item{chunk_size}{Candidates processed per chunk.}
  \item{verbose}{Whether to report chunk progress.}
}
\value{
A list with positive-weight particles, normalized weights, thresholds,
survival rates, and simulation counts.
}
\description{
Samples candidates from the prior, screens them with the low-fidelity
discrepancy, and evaluates high fidelity only for survivors.
}
