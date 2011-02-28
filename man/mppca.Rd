\name{mppca}
\alias{mppca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mppca
}
\description{
estimates the variational posterior distribution of a MPPCA on a data set. A lower bound is calculated and monitored at each iteration. This posterior can be used for various purposes (e.g. MC proposal distribution). It can be transformed using mppcaToGmm and subMppca, outputing a GMM.
}
\usage{
mppca(data, ncomp, thres = 0.1, maxit = NULL, qmax = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
matrix of row-elements.
}
  \item{ncomp}{
number of components in the posterior.
}
  \item{thres}{
threshold for lower bound variations between 2 iterations. Convergence is decided if this variation is below thres.
}
  \item{maxit}{
if NULL, the stopping criterion is related to thres. If not NULL, maxit iterations are performed.
}
  \item{qmax}{
maximal rank of the posterior factor matrices. If NULL, is set to d-1.
}
}
\value{
estimated posterior MPPCA with ncomp components.
}
\references{
Beal, M. J. (2003) _Variational Algorithms for approximate inference_, PhD Thesis, University of London.
}
\author{
Pierrick Bruneau
}
\seealso{
newMppca mppcaToGmm subMppca
}
\examples{
temp <- mppca(pendat, 100, maxit=30, qmax=8)
}

