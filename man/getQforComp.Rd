\name{getQforComp}
\alias{getQforComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getQforComp
}
\description{
gets the rank associated with a properly aligned factor matrix.
}
\usage{
getQforComp(loadings, tau = 1, verbose = FALSE, quick = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loadings}{
aligned factor matrix.
}
  \item{tau}{
diagonal noise used for KL computations.
}
  \item{verbose}{
if TRUE maximal info is displayed.
}
  \item{quick}{
if TRUE, column norm values are used instead of KL computations (less accurate but faster).
}
}
\value{
rank associated with loadings.
}
\references{
Bruneau, P., Gelgon, M. and Picarougne, F. (2011) _Component-level aggregation of probabilistic PCA mixtures using variational-Bayes_, Tech Report, http://hal.archives-ouvertes.fr/docs/00/56/72/99/PDF/techrep.pdf.
}
\author{
Pierrick Bruneau
}
\seealso{
newMppca mppca
}
\examples{
temp <- getQforComp(pcapen[[1]]$wmean[[2]], quick=TRUE)
}

