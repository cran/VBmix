\name{mkmeans}
\alias{mkmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mahalanobis K-means
}
\description{
K-means variant that uses a class-wise Mahalanobis metric. The implementation follows somewhat Lloyd's, with class-wise covariance computation step following that of centres.
}
\usage{
mkmeans(dat, k, maxiter = 100, seeds = NULL, prior = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Matrix with n rows and d columns of n d-dimensional data elements to cluster.
}
  \item{k}{
Number of clusters in the output.
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{seeds}{
Optional indexes of initial centres taken in the input data. If NULL, uniform sampling is used.
}
  \item{prior}{
  Prior population size used for regularizing components.
  }
}
\details{
K-means is characterized by the use of identity as the metric. To remain close to this in spirit, each class-wise covariance matrix is normalized after computation so that is trace equals d. This avoids excessively unbalanced classes, while facilitating the case where the support of a given cluster is less than 2 - covariance cannot be computed in this case. Covariance then defaults to identity. Also to prevent degeneracies when 2 < cluster size < d, a regularization term proportional to sample data features is added to the covariance diagonal.

The returned value follows the GMM data structure (i.e., as returned by e.g. varbayes() and newGmm())
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
\item{labels}{Cluster labels taking values in 1..k}
\item{w}{Numeric vector of cluster weights}
\item{mean}{List of mean vectors}
\item{cov}{List of covariance matrices}
}
\author{
P. Bruneau
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{newGmm}}, \code{\link{varbayes}}
}
\examples{
	mod <- mkmeans(irisdata, 3)
}
