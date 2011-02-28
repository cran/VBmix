\name{binnedEntropy}
\alias{binnedEntropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
binnedEntropy
}
\description{
uses bins to approximate the empirical entropy of a variable.
}
\usage{
binnedEntropy(v, nbins = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
a numeric vector.
}
  \item{nbins}{
number of bins used to estimate the entropy.
}
}
\value{
entropy value.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- binnedEntropy(irisdata[,1])
}

