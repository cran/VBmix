\name{appendToMppca}
\alias{appendToMppca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
appendToMppca
}
\description{
appends mppca2 to mppca1.
}
\usage{
appendToMppca(mppca1, mppca2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mppca1}{
MPPCA model to be appended to.
}
  \item{mppca2}{
MPPCA to append to mod1.
}
}
\value{
appended models.
}
\author{
Pierrick Bruneau
}
\seealso{
appendToGmm appendToList
}
\examples{
temp <- appendToMppca(pcapen[[1]], pcapen[[2]])
}
