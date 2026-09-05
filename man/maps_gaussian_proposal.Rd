\name{maps_gaussian_proposal}
\alias{maps_gaussian_proposal}
\title{Gaussian Random-Walk Proposal for MAPS}
\usage{
maps_gaussian_proposal(theta, Theta, weights)
}
\arguments{
  \item{theta}{Active particles to perturb, as a matrix.}
  \item{Theta}{Complete particle population, as a matrix.}
  \item{weights}{Non-negative particle weights.}
}
\value{A matrix with the same dimensions as \code{theta}.}
\description{
Uses a Gaussian random walk whose covariance is twice the weighted particle
covariance.
}
