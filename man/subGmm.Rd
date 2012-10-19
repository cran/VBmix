\name{subGmm}
\alias{subGmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
subGmm
}
\description{
select a subset of components and dimensions from an input GMM.
}
\usage{
subGmm(model, dims = c(1, 2), inds = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
GMM from which to extract subsets.
}
  \item{dims}{
numeric vector of the extracted dimensions.
}
  \item{inds}{
numeric vector of selected components indices. If NULL, all components are selected.
}
}
\value{
subset of input GMM.
}
\author{
Pierrick Bruneau
}
\seealso{
newGmm
}
\examples{
temp <- subGmm(gmmpen[[1]], inds=1:3)
}
