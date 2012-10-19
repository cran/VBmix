\name{gmmgen}
\alias{gmmgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gmmgen
}
\description{
sample elements from a GMM.
}
\usage{
gmmgen(mod, nitem)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
GMM sampled from.
}
  \item{nitem}{
number of elements to be sampled.
}
}
\value{
nitem x d matrix with elements as rows.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- gmmgen(gmmpen[[1]], 50)
}

