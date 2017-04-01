\name{covgen}
\alias{covgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
covgen
}
\description{
generates random definite positive matrices (i.e. valid covariance matrices).
}
\usage{
covgen(d = 2, bounds = c(1, 5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
rank of the square matrix to be returned.
}
  \item{bounds}{
minima and maximal values for diagonal values.
}
}
\value{
random definite positive matrix
}
\author{
Pierrick Bruneau
}
\note{
Matrix cells are sampled with an heuristic not guaranteed to lead to definite positiveness: this characteristic is only controlled before function return. If positive definite after control, the matrix is returned. If not, an error message is issued.
}
\seealso{
randomGmm
}
\examples{
temp <- covgen()
}

