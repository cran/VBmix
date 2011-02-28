\name{Rdct2D}
\alias{Rdct2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rdct2D
}
\description{
performs 2D DCT on a real matrix.
}
\usage{
Rdct2D(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
matrix of real values.
}
}
\value{
matrix of DCT transformed values.
}
\author{
Pierrick Bruneau
}
\seealso{
Rdct RinvDct2D
}
\examples{
temp <- Rdct2D(irisdata)
}

