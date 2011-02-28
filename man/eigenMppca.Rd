\name{eigenMppca}
\alias{eigenMppca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
eigenMppca
}
\description{
uses eigen decompositions to align factor matrices to principal bases (see references). NB: mppca and mmppca already perform this operation during their post-processing.
}
\usage{
eigenMppca(mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
MPPCA model which components have to be aligned.
}
}
\value{
adjusted MPPCA.
}
\references{
Tipping, M. E. and Bishop, C. M. (1999) _Probabilistic principal component analysis_ ,Journal of the Royal
Statistical Society - B Series, Volume 61, Number 3, Pages 611-622.
}
\author{
Pierrick Bruneau
}
\seealso{
mppca newMppca
}
\examples{
temp <- eigenMppca(pcapen[[2]])
}

