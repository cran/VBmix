\name{getVarbayesResp}
\alias{getVarbayesResp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getVarbayesResp
}
\description{
gets posterior responsibilities for a data set, according to the variational posterior of a GMM.
}
\usage{
getVarbayesResp(data, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
matrix of row-elements.
}
  \item{model}{
variational posterior of a GMM
}
}
\value{
responsibility matrix (Z in references) resulting from the parameters.
}
\references{
Bishop, C. M. (2006) _Pattern Recognition and Machine Learning_ Chap. 10, Springer.
}
\author{
Pierrick Bruneau
}
\seealso{
getResp ZtoLabels
}
\examples{
temp <- getVarbayesResp(pendat, vbpen[[2]])
}

