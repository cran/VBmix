\name{circlegen}
\alias{circlegen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
circlegen
}
\description{
generate data elements along a 2D circle with additional noise.
}
\usage{
circlegen(npts = 200, radius = 10, noise = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npts}{
number of elements to generate.
}
  \item{radius}{
radius of the circle.
}
  \item{noise}{
determines the width of the circle stroke.
}
}
\value{
matrix of sampled row-elements.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- circlegen()
}

