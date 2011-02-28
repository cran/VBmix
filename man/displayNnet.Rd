\name{displayNnet}
\alias{displayNnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
displayNnet
}
\description{
displays the colored decision regions of a neural network model. Data symbols are also optionally displayed. Data and model should be 2D.
}
\usage{
displayNnet(nnet.model, datamatrix, datalabels, subset = NULL, displayPoints = TRUE, steps = 100, alpha = 0.4, lwd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nnet.model}{
a neural network model, as returned by nnet (nnet library)
}
  \item{datamatrix}{
a matrix of row-elements.
}
  \item{datalabels}{
matrix of binary indicator variables for labels (as used by nnet).
}
  \item{subset}{
vector of indexes of a data subset to be displayed. If NULL, all points are displayed.
}
  \item{displayPoints}{
if FALSE, only decision regions are displayed.
}
  \item{steps}{
influences the resolution of the decision regions. Low values will provoke aliasing, high values are slower to be displayed.
}
  \item{alpha}{
alpha blending parameter between decision regions and data symbols.
}
  \item{lwd}{
magnification factor for the stroke width used to plot symbols.
}
}
\value{
a new plotting window displaying decision regions associated to the parametrized neural network.
}
\author{
Pierrick Bruneau
}
\seealso{
nnet
}
\examples{
temp <- class.ind(irislabels)
temp2 <- setDomain(irisdata[,1:2], 10)
temp3 <- nnet(temp2, temp, size=10)
displayNnet(temp3, temp2, temp)
}

