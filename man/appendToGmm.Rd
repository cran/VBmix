\name{appendToGmm}
\alias{appendToGmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
appendToGmm
}
\description{
concatenates mod2 to mod1.
}
\usage{
appendToGmm(mod1, mod2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod1}{
GMM to which mod2 is appended.
}
  \item{mod2}{
GMM appended to mod1.
}
}
\value{
 GMM with concatenated models, with a set accordingly.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- appendToGmm(gmmpen[[1]], gmmpen[[2]])
}
