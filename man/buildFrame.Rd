\name{buildFrame}
\alias{buildFrame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
buildFrame
}
\description{
builds a data frame from a matrix of elements and a vector of numeric labels.
}
\usage{
buildFrame(datamatrix, labels, dims = 1:2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamatrix}{
matrix of row-elements.
}
  \item{labels}{
vector of numeric labels.
}
  \item{dims}{
subset of variables extracted from datamatrix.
}
}
\value{
built data frame.
}
\author{
Pierrick Bruneau
}
\examples{
irisdata[c(1,7,35,56,131),]
# returns:
#         Sepal.Length Sepal.Width Petal.Length Petal.Width
#[1,]          5.1         3.5          1.4         0.2
#[2,]          4.6         3.4          1.4         0.3
#[3,]          4.9         3.1          1.5         0.2
#[4,]          5.7         2.8          4.5         1.3
#[5,]          7.4         2.8          6.1         1.9
irislabels[c(1,7,35,56,131)]
# returns:
#[1] 1 1 1 2 3
temp <- buildFrame(irisdata, irislabels, dims=1:4)
}

