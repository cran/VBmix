\name{spiralgen}
\alias{spiralgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
spiralgen
}
\description{
generates data elements along a spiral with additional noise.
}
\usage{
spiralgen(radius = 10, n = 1000, laps = 2, noise = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{radius}{
determines the radius of a spiral revolution.
}
  \item{n}{
number of elements to generate.
}
  \item{laps}{
number of revolutions of the spiral.
}
  \item{noise}{
determines the width of the spiral stroke.
}
}
\value{
matrix of sampled row-elements.
}
\author{
Pierrick Bruneau
}
\seealso{
datagen circlegen
}
\examples{
temp <- spiralgen()
}

