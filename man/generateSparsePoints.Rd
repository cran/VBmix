\name{generateSparsePoints}
\alias{generateSparsePoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generateSparsePoints
}
\description{
generates a set of points pairwise-separated by a minimal distance. Is not guaranteed to converge: when maxit is reached, current points are returned.
}
\usage{
generateSparsePoints(npoints, dim = 2, span = 10, mindist = 2, maxit = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npoints}{
number of points to generate (i.e. in a matrix with elements as rows.
}
  \item{dim}{
number of variables to generate.
}
  \item{span}{
[-span, span] is used as bounds to uniform sampling for all variables.
}
  \item{mindist}{
minimal distance that each element should have with all others. the "control" C routine is used to perform this verification. All points that do not respect this constraint are resampled.
}
  \item{maxit}{
maximal number of iterations before current elements are returned.
}
}
\value{
matrix with well separated elements as its rows.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- generateSparsePoints(10)
}

