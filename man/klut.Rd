\name{klut}
\alias{klut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
klut
}
\description{
compute Unscented Transform approximation to KL divergence between GMM.
}
\usage{
klut(mod1, mod2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod1}{
GMM parameter to KL(mod1 || mod2).
}
  \item{mod2}{
GMM parameter to KL(mod1 || mod2).
}
}
\value{
KL value.
}
\references{
Goldberger, J., Gordon, and Greenspan, H. (2003) _An efficient image similarity measure based on
approximations of KL-divergence between two Gaussian mixtures_ ICCV Proceedings, Volume 1, Pages 487-493.
}
\author{
Pierrick Bruneau
}
\seealso{
klmc
}
\examples{
temp <- klut(gmmpen[[1]], gmmpen[[2]])
}

