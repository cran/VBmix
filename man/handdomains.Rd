\name{handdomains}
\alias{handdomains}
\docType{data}
\title{
handdomains
}
\description{
original domains of non-void pixels in the handwritten digits collection, to be used along with reBuild.
}
\format{
  The format is:
List of 2
 $ : num [1:717] 0.816 0.251 0.278 0.161 0.412 ...
 $ : num [1:717] 1 1 1 1 1 1 1 1 1 1 ...
}
\examples{
temp <- reBuild(handdat[123,], handvoid, handnonvoid, handdomains)
}
\keyword{datasets}

