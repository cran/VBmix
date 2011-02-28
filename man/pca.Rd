\name{pca}
\alias{pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pca
}
\description{
transforms a data set, and returns coordinates in the principal basis.
}
\usage{
pca(dat, ncomp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
matrix of row-elements.
}
  \item{ncomp}{
number of retained variables in the output result. If NULL, all transformed variables are returned.
}
}
\value{
matrix of transformed row-elements.
}
\references{
Tipping, M. E. and Bishop, C. M. (1999) _Probabilistic principal component analysis_ ,Journal of the Royal
Statistical Society - B Series, Volume 61, Number 3, Pages 611-622.

Bruneau, P., Gelgon, M. and Picarougne, F. (2010) _Aggregation of probabilistic PCA mixtures with
a variational-Bayes technique over parameters_, ICPR'10.

Bruneau, P., Gelgon, M. and Picarougne, F. (2011) _Component-level aggregation of probabilistic PCA mixtures using variational-Bayes_, Tech Report, http://hal.archives-ouvertes.fr/docs/00/56/72/99/PDF/techrep.pdf.
}
\author{
Pierrick Bruneau
}
\seealso{
mppca
}
\examples{
temp <- pca(irisdata, 3)
}
