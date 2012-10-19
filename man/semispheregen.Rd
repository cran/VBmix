\name{semispheregen}
\alias{semispheregen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
semispheregen
}
\description{
sample data points along a semi-sphere.
}
\usage{
semispheregen(npts = 200, radius = 10, noise = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npts}{
number of elements to be sampled.
}
  \item{radius}{
radius of the sphere.
}
  \item{noise}{
additive gaussian white noise to the sampled points.
}
}
\value{
matrix of row-elements with the sampled elements.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- semispheregen()
}

