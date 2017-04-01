\name{dat2sample}
\alias{dat2sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
dat2sample
}
\description{
generates data elements according to SYN2 process (sample along a semi-sphere with additive noise, see reference).
}
\usage{
dat2sample(nelts, radius, noise, oldbounds = NULL, newbounds = NULL)
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
Bruneau, P., Gelgon, M. and Picarougne, F. (2011) _Component-level aggregation of probabilistic PCA mixtures using variational-Bayes_, Tech Report, http://hal.archives-ouvertes.fr/docs/00/56/72/99/PDF/techrep.pdf.
}
\author{
Pierrick Bruneau
}
\seealso{
dat1sample dat3sample
}
\examples{
temp <- dat2sample(500, 10, 1)
}

