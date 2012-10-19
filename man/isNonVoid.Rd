\name{isNonVoid}
\alias{isNonVoid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
isNonVoid
}
\description{
checks if loadings contains only void columns.
}
\usage{
isNonVoid(loadings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loadings}{
matrix from which we check the columns.
}
}
\value{
TRUE if at least 1 column is not void.
}
\author{
Pierrick Bruneau
}
\seealso{
mppca newMppca
}
\examples{
isNonVoid(pcapen[[1]]$wmean[[2]])
#[1] TRUE
}

