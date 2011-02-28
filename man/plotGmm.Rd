\name{plotGmm}
\alias{plotGmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plotGmm
}
\description{
3D density plot of a 2D GMM.
}
\usage{
plotGmm(mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
GMM object to plot
}
}
\value{
a new plotting window with the 3D density plot.
}
\author{
Pierrick Bruneau
}
\seealso{
displayScatter
}
\examples{
plotGmm(randomGmm())
}

