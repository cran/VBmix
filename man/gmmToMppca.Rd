\name{gmmToMppca}
\alias{gmmToMppca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gmmToMppca
}
\description{
uses eigen decompositions to convert a GMM to a MPPCA model.
}
\usage{
gmmToMppca(model, alpha = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
GMM to be converted.
}
  \item{alpha}{
GMM are associated to weights, and MPPCA models to population sizes. alpha is the chosen population size for the output MPPCA.
}
}
\value{
converted MPPCA model.
}
\author{
Pierrick Bruneau
}
\seealso{
mppcaToGmm
}
\examples{
temp <- gmmToMppca(gmmpen[[3]])
}

