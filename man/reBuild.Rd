\name{reBuild}
\alias{reBuild}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reBuild
}
\description{
re-build a pixmapGrey object from a vector of pixel intensities. As some pixels may be irrelevant over a collection of images (e.g. pixel always white in handwritten digits), some variables may have been filtered or transformed before performing some machine learning process. These transforms are indicated as parameters, and give clues to recover objects in the original image space. NB: assumes that v is scaled in [-10,10]. Additional transformations may thus be performed as appropriate before using this function.
}
\usage{
reBuild(v, voids, nonvoids, domains, placeholder = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
vector to be converted to a pixmapGrey object.
}
  \item{voids}{
vector of position indices in the original signal (i.e. 2D matrix with its columns casted in a vector) that did not carry any information. Replaced by a placeholder in recovered image.
}
  \item{nonvoids}{
vector of positions to which v should be associated in the recovered image.
}
  \item{domains}{
original data domains of pixel intensities prior to being transformed to v's domain. Permit appropriate reconstruction in the domain of pixel intensities used by pixmap (i.e. subset of [0,1]). Formatted similarly to what is required in setDomain.
}
  \item{placeholder}{
placeholder value for pixel positions present in voids.
}
}
\value{
pixmapGrey reconstructed object.
}
\author{
Pierrick Bruneau
}
\seealso{
pixmapGrey pixmapToVector
}
\examples{
temp <- reBuild(handdat[123,], handvoid, handnonvoid, handdomains)
}

