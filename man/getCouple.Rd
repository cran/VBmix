\name{getCouple}
\alias{getCouple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getCouple
}
\description{
computes classification error function described in references, a.k.a couple error. In brief, evaluates how elements are gathered similarly, irrespectively of exact label values (adapted to clustering).
}
\usage{
getCouple(vec1, vec2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec1}{
 vector of numeric labels.
}
  \item{vec2}{
 vector of numeric labels.
}
}
\value{
classification error in [0,1].
}
\references{
Fowlkes, E. B. and Mallows, C. L. (1983) _A method for comparing two hierarchical clusterings_, J. Am.
Stat. Assoc., Volume 78, Pages 553-569.

Picarougne, F. ,Azzag, H., Venturini, G. and Guinot, C. (2007) _A new approach of data clustering using a
flock of agents_, Evolutionary Computation, Volume 15, Number 3, Pages 345-367.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- EM(irisdata, 4)
getCouple(temp$labels, irislabels)
}

