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
plotGmm(mod, steps=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
GMM object to plot
}
  \item{steps}{
  specifies the hoizontal and vertical amount of vertices used to build the wireframe plot.
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
# a larger number of steps (eg 200) should be used for a visually effective 3D plot.
plotGmm(randomGmm(), steps=20)
}

