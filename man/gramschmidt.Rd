\name{gramschmidt}
\alias{gramschmidt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gramschmidt
}
\description{
performs Gram-Schmidt orthogonalization on mat.
}
\usage{
gramschmidt(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
matrix object to orthogonalize.
}
}
\value{
orthogonalized matrix.
}
\author{
Pierrick Bruneau
}
\seealso{
mppca newMppca
}
\examples{
temp <- gramschmidt(pcapen[[3]]$wmean[[1]])
}

