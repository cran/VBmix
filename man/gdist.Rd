\name{gdist}
\alias{gdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise distance between groups
}
\description{
dist computes distances between every pair of elements in a group: gdist takes two arguments, and returns the matrix of distances wrt every possible pair with one argument from each group. General Mahalanobis metrics are also allowed.
}
\usage{
gdist(g1, g2, metric = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g1}{
  n1 x d matrix with n1 data elements.
}
  \item{g2}{
  n2 x d matrix with n2 data elements.
}
  \item{metric}{
  If NULL, defaults to identity (i.e. Euclidean distance). A d x d matrix is assumed as a sample covariance matrix, and its inverse is used to compute distances (i.e. Mahalanobis distance). Likewise, a list of n2 d x d matrices can be provided, yielding distances specific to each row in g2. 
}
}
\details{
This function is especially useful to make algorithms such as k-means (or mkmeans in the package) more efficient - rows in g1 are then generally the data set, and in g2 respectively cluster centres.
}
\value{
n1 x n2 matrix of distances
}

\author{
P. Bruneau
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{dist}}, \code{\link{mkmeans}}
}
\examples{
dists <- gdist(irisdata, irisdata[c(1,11,21),])
}
