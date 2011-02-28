\name{gmmdensity}
\alias{gmmdensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gmmdensity
}
\description{
get densities of a set of elements w.r.t a GMM.
}
\usage{
gmmdensity(mod, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
reference GMM.
}
  \item{data}{
matrix of row-elements.
}
}
\value{
numeric vector containing densities.
}
\author{
Pierrick Bruneau
}
\seealso{
gmmgen
}
\examples{
temp <- gmmgen(gmmpen[[1]], 50)
temp2 <- gmmdensity(gmmpen[[1]], temp[[1]])
}

