\name{gaussianKL}
\alias{gaussianKL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gaussianKL
}
\description{
computes KL( N(0, Sigma_0) || N(0, Sigma_1) ).
}
\usage{
gaussianKL(N0, N1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N0}{
Sigma_0
}
  \item{N1}{
Sigma_1
}
}
\value{
KL value.
}
\author{
Pierrick Bruneau
}
\seealso{
klmc
}
\examples{
temp <- gaussianKL(gmmpen[[1]]$cov[[1]], gmmpen[[1]]$cov[[2]])
}

