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
A list object, with the following items:

\item{model}{posterior variational distribution.}
\item{data}{a copy of the data parameter.}
\item{nk}{counts, for each iteration, of the population modeled by each Gaussian component.}
\item{agitation}{agitation measures (see Beal 2003 for explanation) for each iteration and Gaussian component.}
\item{bound}{latest monitored bound value (convergence criterion maximized throughout the process).}

The model item is structured in a list as follows:

\item{alpha}{hyperparameters influencing the active components in the posterior.}
\item{beta}{hyperparameters regarding shaping of the Normal-Wishart posteriors.}
\item{nu}{hyperparameters regarding shaping of the Normal-Wishart posteriors.}
\item{mean}{hyperparameters regarding shaping of the Normal-Wishart posteriors.}
\item{wish}{hyperparameters regarding shaping of the Normal-Wishart posteriors.}
}
\references{
Bishop, C. M. (2006) \emph{Pattern Recognition and Machine Learning}, Chapter 10, Springer.

Beal, M. J. (2003) \emph{Variational Algorithms for approximate inference}, PhD thesis, University of London.
}
\author{
Pierrick Bruneau
}
\seealso{
EM extractSimpleModel
}
\examples{
temp <- varbayes(irisdata, 20)
}

