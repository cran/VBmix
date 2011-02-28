\name{klmc}
\alias{klmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
klmc
}
\description{
computes Monte Carlo estimate of KL divergence between GMM.
}
\usage{
klmc(mod1, mod2, nsamp = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod1}{
GMM parameter to KL(mod1 || mod2).
}
  \item{mod2}{
GMM parameter to KL(mod1 || mod2).
}
  \item{nsamp}{
number of samples used to build estimate.
}
}
\value{
KL value.
}
\author{
Pierrick Bruneau
}
\seealso{
jsmc gaussianKL
}
\examples{
temp <- klmc(gmmpen[[1]], gmmpen[[2]])
}

