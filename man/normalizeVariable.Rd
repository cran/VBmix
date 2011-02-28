\name{normalizeVariable}
\alias{normalizeVariable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
normalizeVariable
}
\description{
normalizes a variable (numeric vector) in [0,1].
}
\usage{
normalizeVariable(v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
a numeric vector.
}
}
\value{
normalized numeric vector.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- normalizeVariable(irisdata[,1])
}

