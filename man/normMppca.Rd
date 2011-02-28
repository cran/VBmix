\name{normMppca}
\alias{normMppca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
normMppca
}
\description{
adjusts a MPPCA model to ensure that all factor matrices have same rank (q).
}
\usage{
normMppca(mppca1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mppca1}{
MPPCA model to be adjusted.
}
}
\value{
adjusted MPPCA.
}
\author{
Pierrick Bruneau
}
\seealso{
newMppca mppca
}
\examples{
temp <- newMppca()
for(i in 1:5) temp <- appendToMppca(temp, pcapen[[i]])
temp <- normMppca(temp)
}

