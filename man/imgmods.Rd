\name{imgmods}
\alias{imgmods}
\docType{data}
\title{
imgmods
}
\description{
list of 200 3D GMM, sampled from the 1243 images in the 10 first categories of the Caltech-256 image collection. Built using RGBtoLab and varbayes. See reference for information about this image collection.
}
\format{
  The format is:
List of 200 GMM
}
\references{
Griffin, G., Holub, A. and Perona, P. (2007) _Caltech-256 object category dataset_, http://authors.library.caltech.edu/7694, Technical Report 7694, California Institute of Technology.
}
\examples{
temp <- gmmgen(imgmods[[10]], 1000)
}
\keyword{datasets}

