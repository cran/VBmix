\name{vbconstr}
\alias{vbconstr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
vbconstr
}
\description{
estimates the variational posterior distribution of a GMM that aggregates a constrained collection of GMM.  A lower bound is calculated and monitored at each iteration. This posterior can be used for various purposes (e.g. MC proposal distribution). It can be transformed using extractSimpleModel, outputing a GMM.
}
\usage{
vbconstr(models, ncomp, thres = 0.1, maxit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
GMM made with the weighted sum of the collection of GMM to aggregate. a is used to model constraints between components in this GMM.
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
estimated posterior with ncomp components.
}
\references{
Bruneau, P., Gelgon, M., and Picarougne, F. (2010) _Parsimonious reduction of Gaussian mixture
models with a variational-Bayes approach_, Pattern Recognition, Volume 43, Pages 850-858.
}
\author{
Pierrick Bruneau
}
\seealso{
vbcomp extractSimpleModel
}
\examples{
temp1 <- newGmm()
for(i in 1:10) temp1 <- appendToGmm(temp1, gmmpen[[i]])
temp2 <- vbconstr(temp1, 50)
}

