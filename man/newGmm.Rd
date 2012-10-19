\name{newGmm}
\alias{newGmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
newGmm
}
\description{
creates an empty GMM data structure.
}
\usage{
newGmm()
}
%- maybe also 'usage' for other objects documented here.
\value{
list object with the following members:
\item{w}{numeric vector containing the component weights of the mixture model.}
\item{mean}{list with respective means (numeric vectors) as elements.}
\item{cov}{list with respective covariance matrices as elements.}
\item{a}{constraints between components, encoded in a numeric vector. One value per component. 2 components associated to the same value are said to be from the same origin. Used in vbconstr.}
}
\author{
Pierrick Bruneau
}
\seealso{
varbayes
vbconstr
}
\examples{
temp <- newGmm()
}
