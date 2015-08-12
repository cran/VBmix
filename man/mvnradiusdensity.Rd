\name{mvnradiusdensity}
\alias{mvnradiusdensity}
\title{
mvnradiusdensity
}
\description{
get densities of a set of squared radii, i.e. obtained from a Mahalanobis distance computed externally wrt an inverse covariance matrix.
}
\usage{
mvnradiusdensity(cov, radii)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov}{
Covariance matrix from which we compute the determinant.
}
  \item{radii}{
Radii wrt which we directly take the density values.
}
}
\value{
numeric vector containing densities.
}
\author{
Pierrick Bruneau
}
\seealso{
mvngen mvndensity
}
\examples{
temp <- mvngen(c(0, 0), diag(2), 5)
R <- chol(solve(diag(2))) 
# trivial here, as Cholesky R of I(-1) is I
temp <- temp %*% t(R) %*% R %*% t(temp)
mvnradiusdensity(diag(2), diag(temp))
}

