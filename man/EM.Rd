\name{EM}
\alias{EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM
}
\description{
estimates a GMM on data using EM algorithm.
}
\usage{
EM(data, ncomp, model=c("general", "diagonal", "spherical"), class=FALSE, 
	thres = 0.1, maxit = NULL, rbic=FALSE, debug=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
matrix of row-elements.
}
  \item{ncomp}{
maximal number of components in the GMM. In case of degeneracies, the final model size may be less than ncomp.
}

  \item{model}{
Hypothesis on the modfel to estimate: "general", "diagonal" or "spherical" covariance matrices.
}

\item{class}{
If TRUE, hard allocate elements in the E step (see CEM variant in Biernacki et al.). If FALSE, compute soft responsibilities as in usual EM algorithm.
}

  \item{thres}{
threshold for lower bound variations between 2 iterations. Convergence is decided if this variation is below thres.
}
  \item{maxit}{
 if NULL, the stopping criterion is related to thres. If not NULL, maxit iterations are performed.
}
\item{rbic}{
if FALSE, output BIC criterion associated to the obtained GMM. If TRUE, use a variant that accounts for the dimensionality of the model.
}
\item{debug}{
if TRUE, display debug markers.
}
}
\value{
estimated GMM with at most ncomp components, with labels containing associated labels for data in addition.
\item{labels}{Cluster labels taking values in 1..k}
\item{w}{Numeric vector of cluster weights}
\item{mean}{List of mean vectors}
\item{cov}{List of covariance matrices}
\item{likelihood}{Likelihood value of the model}
\item{bic}{BIC criterion of the model}
}
\references{
Bishop, C. M. (2006) _Pattern Recognition and Machine Learning_, Chapter 9, Springer.
Biernacki, C. et al. _Model-based cluster and discriminant analysis with the MIXMOD software_, Computational Statistics and Data Analysis 51.2 (2006): 587-600.
}
\author{
Pierrick Bruneau
}
\seealso{
newGmm varbayes
}
\examples{
temp <- EM(irisdata, 4)
}

