\name{RinvDct2D}
\alias{RinvDct2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RinvDct2D
}
\description{
performs inverse 2D DCT on a real matrix.
}
\usage{
RinvDct2D(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
matrix of real values.
}
}
\value{
matrix of inverse DCT transformed values.
}
\author{
Pierrick Bruneau
}
\seealso{
Rdct Rdct2D
}
\examples{
temp <- RinvDct2D(Rdct2D(irisdata))
}

