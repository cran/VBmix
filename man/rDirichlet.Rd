\name{rDirichlet}
\alias{rDirichlet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rDirichlet
}
\description{
samples from the Dirichlet distribution.
}
\usage{
rDirichlet(K, alpha = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{
order of the sample.
}
  \item{alpha}{
alpha parameter of the distribution (i.e. alpha repeated K times).
}
}
\value{
numeric vector, which values are in [0,1] and sum to 1.
}
\author{
Pierrick Bruneau
}
\seealso{
dDirichlet
}
\examples{
temp <- rDirichlet(4)
}

