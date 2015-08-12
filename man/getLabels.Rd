\name{getLabels}
\alias{getLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getLabels
}
\description{
gets numeric labels that associates a data set and a GMM.
}
\usage{
getLabels(model, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
GMM.
}
  \item{data}{
matrix of row-elements.
}
}
\value{
vector of numeric labels, that take values of the respective component indexes in the GMM.
}
\author{
Pierrick Bruneau
}
\seealso{
newGmm
}
\examples{
temp <- EM(irisdata, 4)
temp2 <- getLabels(temp, irisdata)
}

