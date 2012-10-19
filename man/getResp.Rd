\name{getResp}
\alias{getResp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getResp
}
\description{
get posterior responsibilities of elements in a data set, according to a posterior MPPCA distribution.
}
\usage{
getResp(data, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
matrix of row-elements.
}
  \item{model}{
posterior MPPCA.
}
}
\value{
nxk matrix (with n the number of row-elements, and k the number of components in the MPPCA) of membership probabilities. (i.e. Z in references)
}
\references{
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
temp <- getResp(pendat, pcapen[[1]])
}

