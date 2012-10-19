\name{jsmc}
\alias{jsmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
jsmc
}
\description{
computes Monte Carlo estimate of Jensen-Shannon (JS) divergence between GMM.
}
\usage{
jsmc(mod1, mod2, nsamp = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod1}{
GMM parameter to JS(mod1 || mod2).
}
  \item{mod2}{
GMM parameter to JS(mod1 || mod2).
}
  \item{nsamp}{
number of samples used to build estimate.
}
}
\value{
JS divergence value.
}
\author{
Pierrick Bruneau
}
\seealso{
klmc gaussianKL
}
\examples{
temp <- jsmc(gmmpen[[1]], gmmpen[[2]])
}

