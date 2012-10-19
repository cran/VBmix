\name{mmppca}
\alias{mmppca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mmppca
}
\description{
estimates the variational posterior distribution of a MPPCA that aggregates a collection of input MPPCA models. A lower bound is calculated and monitored at each iteration. This posterior can be used for various purposes (e.g. MC proposal distribution). It can be transformed using mppcaToGmm and subMppca, outputing a GMM. The maximal rank of output factor matrices is determined by the inputs.
}
\usage{
mmppca(mods, ncomp, thres = 0.1, maxit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mods}{
input MPPCA that concatenates the set of components to aggregate.
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
estimated posterior MPPCA with ncomp components.
}
\references{
Bruneau, P., Gelgon, M. and Picarougne, F. (2010) _Aggregation of probabilistic PCA mixtures with
a variational-Bayes technique over parameters_, ICPR'10.

Bruneau, P., Gelgon, M. and Picarougne, F. (2011) _Component-level aggregation of probabilistic PCA mixtures using variational-Bayes_, Tech Report, http://hal.archives-ouvertes.fr/docs/00/56/72/99/PDF/techrep.pdf.
}
\author{
Pierrick Bruneau
}
\seealso{
newMppca mppca subMppca
}
\examples{
temp <- newMppca()
for(i in 1:3) temp <- appendToMppca(temp, pcapen[[i]])
temp2 <- mmppca(temp, 50, maxit=30)
}

