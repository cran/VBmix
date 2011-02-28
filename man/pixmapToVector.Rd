\name{pixmapToVector}
\alias{pixmapToVector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pixmapToVector
}
\description{
converts a pixmapGrey object to a numeric vector. The pixel matrix is casted to a vector by appending successive columns.
}
\usage{
pixmapToVector(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
pixmapGrey object.
}
}
\value{
numeric vector containing pixel intensities.
}
\author{
Pierrick Bruneau
}
\seealso{
pixmapGrey reBuild readPixmapFile
}
\examples{
# use with path to actual train-... file
#temp <- readPixmapFile("data/train-images-idx3-ubyte")
#temp2 <- pixmapToVector(temp[[3]])
}

