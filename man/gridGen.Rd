\name{gridGen}
\alias{gridGen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gridGen
}
\description{
generates a matrix valued with a regular grid of 2D coordinates.
}
\usage{
gridGen(xlim = c(-10, 10), ylim = c(-10, 10), step = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xlim}{
x bounds.
}
  \item{ylim}{
y bounds.
}
  \item{step}{
size of the square matrix.
}
}
\value{
'grid' matrix
}
\author{
Pierrick Bruneau
}
\examples{
temp <- gridGen()
}

