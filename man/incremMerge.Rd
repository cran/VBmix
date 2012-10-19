\name{incremMerge}
\alias{incremMerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
incremMerge
}
\description{
updates a reference MPPCA model with an input distribution.
}
\usage{
incremMerge(modref, newmod, k = 200, nit = 100, quick = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modref}{
reference MPPCA to update.
}
  \item{newmod}{
new MPPCA to incorporate.
}
  \item{k}{
number of components of the output variational posterior.
}
  \item{nit}{
number of iterations used in the mmppca call that performs the update.
}
  \item{quick}{
boolean parameter transmitted to the subMppca routine that shrinks the output variational posterior.
}
}
\value{
updated variational posterior.
}
\author{
Pierrick Bruneau
}
\seealso{
mppca mmppca
}
\examples{
# commented for packaging needs (requires approx. 5s)
#temp <- incremMerge(pcapen[[1]], pcapen[[2]], quick=T)
}

