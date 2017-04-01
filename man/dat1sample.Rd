\name{dat1sample}
\alias{dat1sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
dat1sample
}
\description{
generates data elements according to SYN1 process (sample from a 2D GMM, linearly transformed with additive noise, see reference).
}
\usage{
dat1sample(nelts, gmm, noise, transform=generate2Dtransform(2), 
	oldbounds = NULL, newbounds = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nelts}{
number of elements to generate.
}
  \item{gmm}{
2D GMM to be sampled from.
}
  \item{noise}{
additive noise magnitude.
}
  \item{transform}{
matrix defining linear transform. Defaults to I.
}
  \item{oldbounds}{
optional argument for sample rescaling. If not NULL, transmitted to setDomain as oldspan.
}
  \item{newbounds}{
optional argument for sample rescaling. If not NULL, transmitted to setDomain as newspan.
}
}
\value{
matrix of sampled row-elements
}
\references{
Bruneau, P., Gelgon, M. and Picarougne, F. (2011) _Component-level aggregation of probabilistic PCA mixtures using variational-Bayes_, Tech Report, http://hal.archives-ouvertes.fr/docs/00/56/72/99/PDF/techrep.pdf.
}
\author{
Pierrick Bruneau
}
\seealso{
dat2sample dat3sample
}
\examples{

temp <- dat1sample(500, randomGmm(), 1, generate2Dtransform())

}




