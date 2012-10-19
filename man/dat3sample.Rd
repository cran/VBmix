\name{dat3sample}
\alias{dat3sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
dat3sample
}
\description{
generates data elements according to SYN3 process (sample along a 2D circle with additive noise, and linearly transform to higher dimensional space with further noise addition, see reference).
}
\usage{
dat3sample(nelts, radius, noise, transform, oldbounds = NULL, newbounds = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nelts}{
number of elements to generate.
}
  \item{radius}{
radius of the sphere to sample from.
}
  \item{noise}{
additive noise magnitude.
}
  \item{transform}{
matrix defining linear transform.
}
  \item{oldbounds}{
optional argument for sample rescaling. If not NULL, transmitted to setDomain as oldspan.
}
  \item{newbounds}{
optional argument for sample rescaling. If not NULL, transmitted to setDomain as newspan.
}
}
\value{
matrix of sampled row-elements.
}
\references{
runeau, P., Gelgon, M. and Picarougne, F. (2011) _Component-level aggregation of probabilistic PCA mixtures using variational-Bayes_, Tech Report, http://hal.archives-ouvertes.fr/docs/00/56/72/99/PDF/techrep.pdf.
}
\author{
Pierrick Bruneau
}
\seealso{
dat1sample dat2sample
}
\examples{
temp <- dat3sample(500, 10, 1, generate2Dtransform())
}

