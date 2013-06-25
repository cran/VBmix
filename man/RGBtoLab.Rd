\name{RGBtoLab}
\alias{RGBtoLab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RGBtoLab
}
\description{
transform a .ppm file into a matrix of (L,a,b) pixel intensities (1 row-element per pixel).
}
\usage{
RGBtoLab(filename, filterWhite = FALSE, addCoords = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
path to a .ppm file. Alternatively, if needed, R file path manipulating routines are documented in document r-lang.pdf, section 7.1)
}
  \item{filterWhite}{
if TRUE, filter white points from result to return.
}
  \item{addCoords}{
if TRUE, append 2 normalized (x,y) coordinates for each pixel.
}
}
\note{
In order to save space, images associated to names in imgnames were not provided in this bundle.
Caltech-256 should be retrieved first, converted to .ppm (e.g. with imageMagick), and then values in imgnames associated to relevants file paths, before using RGBtoLab.
}
\value{
matrix of pixel row-elements.
}
\author{
Pierrick Bruneau
}
\examples{
# image collections are large, thus not provided. 
# The following commented example relates to a member of this image collection.
#temp <- RGBtoLab(imgnames[[2]], filterWhite=TRUE)
}

