\name{mvngen}
\alias{mvngen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mvngen
}
\description{
sample nitem elements from N(mean, cov).
}
\usage{
mvngen(mean, cov, nitem)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mean}{
numeric vector.
}
  \item{cov}{
covariance matrix.
}
  \item{nitem}{
number of items to generate.
}
}
\value{
nitem x d matrix with elements as rows (further denoted as a matrix of row-elements).
}
\author{
Pierrick Bruneau
}
\examples{
mvngen(c(0, 0), diag(2), 5)
#       [,1]        [,2]
#[1,] -0.09898211  1.4516438
#[2,]  0.20814926 -0.1233861
#[3,]  0.18410071  0.5995621
#[4,]  0.65994562  0.8328315
#[5,]  2.33098055 -0.5943117
}

