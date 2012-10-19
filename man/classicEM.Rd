\name{classicEM}
\alias{classicEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
classicEM
}
\description{
estimates a GMM on data using EM algorithm. A lower bound is calculated and monitored at each iteration.
}
\usage{
classicEM(data, k, thres = 0.1, maxit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
matrix of row-elements.
}
  \item{k}{
maximal number of components in the GMM. In case of degeneracies, the final model size may be less than 0.
}
  \item{thres}{
threshold for lower bound variations between 2 iterations. Convergence is decided if this variation is below thres.
}
  \item{maxit}{
 if NULL, the stopping criterion is related to thres. If not NULL, maxit iterations are performed.
}
}
\value{
estimated GMM with at most k components, with labels containing associated labels for data in addition.
}
\references{
Bishop, C. M. (2006) _Pattern Recognition and Machine Learning_, Chapter 9, Springer.
}
\author{
Pierrick Bruneau
}
\seealso{
newGmm varbayes
}
\examples{
temp <- classicEM(irisdata, 4)
}

