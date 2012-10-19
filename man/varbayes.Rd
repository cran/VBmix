\name{varbayes}
\alias{varbayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
varbayes
}
\description{
estimates the variational posterior distribution of a GMM on data using the variational EM algorithm (see references). A lower bound is calculated and monitored at each iteration. This posterior can be used for various purposes (e.g. MC proposal distribution). It can be transformed using extractSimpleModel, outputing a GMM.
}
\usage{
varbayes(data, ncomp, thres = 0.1, maxit = NULL)
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
}
\value{
estimated posterior with ncomp components. Structured in a list object as follows:
\item{alpha}{hyperparameters influencing the active components in the posterior.}
\item{beta}{hyperparameters regarding shaping of the Normal-Wishart posteriors.}
\item{nu}{hyperparameters regarding shaping of the Normal-Wishart posteriors.}
\item{mean}{hyperparameters regarding shaping of the Normal-Wishart posteriors.}
\item{wish}{hyperparameters regarding shaping of the Normal-Wishart posteriors.}
}
\references{
Bishop, C. M. (2006) _Pattern Recognition and Machine Learning_, Chapter 10, Springer.
}
\author{
Pierrick Bruneau
}
\seealso{
classicEM extractSimpleModel
}
\examples{
temp <- varbayes(irisdata, 20)
}

