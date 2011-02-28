\name{mppcaToGmm}
\alias{mppcaToGmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mppcaToGmm
}
\description{
converts a MPPCA model to a GMM.
}
\usage{
mppcaToGmm(model, notau = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
MPPCA model to be converted.
}
  \item{notau}{
if TRUE, covariances are built with Lambda.Lambda^T without adding tau.
}
}
\value{
GMM after conversion.
}
\references{
Tipping, M. E. and Bishop, C. M. (1999) _Probabilistic principal component analysis_ ,Journal of the Royal
Statistical Society - B Series, Volume 61, Number 3, Pages 611-622.
}
\author{
Pierrick Bruneau
}

\seealso{
mppca varbayes
}
\examples{
temp <- mppcaToGmm(pcapen[[1]])
}
