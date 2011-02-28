\name{getDataLikelihood}
\alias{getDataLikelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getDataLikelihood
}
\description{
gets log-likelihoods associated to a matrix of row-elements.
}
\usage{
getDataLikelihood(gmm, dat)
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
numeric vector of log-likelihoods.
}
\author{
Pierrick Bruneau
}
\seealso{
getBic gmmgen
}
\examples{
temp <- getDataLikelihood(gmmpen[[3]], pendat)
}

