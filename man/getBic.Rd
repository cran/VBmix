\name{getBic}
\alias{getBic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getBic
}
\description{
computes BIC criterion (see references) for a specific GMM and data set.
}
\usage{
getBic(gmm, dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmm}{
GMM object.
}
  \item{dat}{
matrix of row-elements.
}
}
\value{
BIC estimate.
}
\references{
Schwarz, G. (1978) _Estimating the dimension of a model_ The Annals of Statistics, Volume 6, Pages 461-464.
}
\author{
Pierrick Bruneau
}
\seealso{
getDataLikelihood varbayes
}
\examples{
temp <- getBic(gmmpen[[1]], pendat)
}

