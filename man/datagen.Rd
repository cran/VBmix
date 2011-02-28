\name{datagen}
\alias{datagen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
datagen
}
\description{
generates data from a random multivariate Gaussian, and adds redundant dimensions by random linear combinations with noise.
}
\usage{
datagen(dreal = 2, deff = 6, npts = 200, noise = 0.1, genmean = rep(0, dreal), genspan = 6, iso = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dreal}{
dimensionality of the multivariate Gaussian.
}
  \item{deff}{
dimensionality of the returned sample.
}
  \item{npts}{
number of elements to be sampled.
}
  \item{noise}{
noise magnitude for the linear combination.
}
  \item{genmean}{
mean of the multivariate Gaussian.
}
  \item{genspan}{
maximal magnitude of the diagonal elements in the covariance matrix. Non-diagonal elements are sampled under constraints of positive-definiteness.
}
  \item{iso}{
sample from an isotropic multivariate Gaussian (i.e. diagonal covariance matrix).
}
}
\value{
matrix of sampled row-elements.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- datagen()
}

