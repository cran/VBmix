\name{setDomain}
\alias{setDomain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
setDomain
}
\description{
performs linear rescaling of given data.
}
\usage{
setDomain(dat, span = 10, oldspan = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
data to rescale. matrix object, with elements as rows, and variables as columns (i.e. variables are rescaled).
}
  \item{span}{
new domain to which dat is rescaled. If type is numeric and length = 1: [-span, span] is used for all variables. If type is numeric and length = 2: [span[1], span[2]] is used for all variables. If a list object: [span[[1]]_i, span[[2]]_i] is used for each variable i.
}
  \item{oldspan}{
if NULL, old domains are computed from dat inspection. Otherwise, is structured as span and replaces inspected values for rescaling.
}
}
\value{
scaled data matrix.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- setDomain(irisdata, span=15)
}

