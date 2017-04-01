\name{handdat}
\alias{handdat}
\docType{data}
\title{
handdat
}
\description{
matrix 300 x 717 of real row-elements. See reference. May be loaded into R with readDataFile. handdat was built using pixmapToVector and filtering variables with zero entropy.
}
\format{
  The format is:
 num [1:300, 1:717] 10 10 10 10 10 10 10 10 10 10 ...
}
\source{
http://yann.lecun.com/exdb/mnist/
}
\references{
LeCun, Y., Bottou, L., Bengio, Y., and Haffner, P. (1998) _Gradient-based learning applied to document
recognition_, Proceedings of the IEEE, Volume 86, Number 11, Pages 2278-2324.
}
\examples{
temp <- reBuild(handdat[123,], handvoid, handnonvoid, handdomains)
}
\keyword{datasets}
