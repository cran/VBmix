\name{newMppca}
\alias{newMppca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
newMppca
}
\description{
creates an empty posterior MPPCA data structure.
}
\usage{
newMppca()
}
%- maybe also 'usage' for other objects documented here.
\value{
list object with the following members:
\item{alpha}{numeric vector for bayesian alpha parameter.}
\item{numoment}{list of numeric vectors, containing E[nu_(kj)] parameters.}
\item{nub}{list of numeric vectors, containing b_(kj) parameters for nu.}
\item{taumoment}{numeric vector for tau parameter. NB: all set identically and statically to 1, as in [Bruneau 2011] a single static tau parameter is used.}
\item{taua}{numeric vector for a_k parameters for tau.}
\item{taub}{numeric vector for b_k parameters for tau.}
\item{wmean}{list of matrices containing E[Lambda_k] parameters.}
\item{wsigma}{list of matrices containing Cov(Lambda_k^(i.)).}
\item{xsigma}{list of matrices containing Cov(x_k).}
\item{mumean}{list of numeric vectors, containing means of the MPPCA model.}
\item{musigma}{list of matrices with covariances for the mean estimates.}
\item{mustar}{list of numeric vectors, containing prior means of the MPPCA model, used for initialisation.}
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
mppca mmppca
}
\examples{
temp <- newMppca()
}
