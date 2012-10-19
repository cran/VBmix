\name{ZtoLabels}
\alias{ZtoLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ZtoLabels
}
\description{
converts a responsibility matrix (Z in references) to a vector of numeric labels.
}
\usage{
ZtoLabels(resp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{
responsibility matrix to convert.
}
}
\value{
labels vector.
}
\references{
Bishop, C. M. (2006) _Pattern Recognition and Machine Learning_ Chap. 9, Springer.
}
\author{
Pierrick Bruneau
}
\seealso{
getResp getVarbayesResp
}
\examples{
temp <- getResp(pendat, pcapen[[2]])
temp2 <- ZtoLabels(temp)
}

