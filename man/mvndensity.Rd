\name{mvndensity}
\alias{mvndensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mvndensity
}
\description{
get densities of a set of elements w.r.t a multivariate normal.
}
\usage{
mvndensity(mean, cov, data, rescaled=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mean}{
numeric vector, mean of the multivariate normal.
}
  \item{cov}{
covariance matrix of the multivariate normal.
}
  \item{data}{
matrix of row-elements.
}
\item{rescaled}{
	if TRUE, a variant accouting for data dimensionality is computed.
}
}
\value{
numeric vector containing densities.
}
\author{
Pierrick Bruneau
}
\seealso{
mvngen
}
\examples{
temp <- mvngen(c(0, 0), diag(2), 5)
mvndensity(c(0,0), diag(2), temp)
#[1] 0.137188286 0.032318242 0.005181099 0.047312602 0.033178600
}

