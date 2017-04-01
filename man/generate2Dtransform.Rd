\name{generate2Dtransform}
\alias{generate2Dtransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate2Dtransform
}
\description{
generate a random matrix to transform a 2D signal to higher dimensional spaces.
}
\usage{
generate2Dtransform(dims = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dims}{
dimensionality of the target space.
}
}
\value{
a dims x 2 matrix defining the transform.
}
\author{
Pierrick Bruneau
}
\seealso{
dat1sample dat3sample
}
\examples{
temp <- generate2Dtransform()
}

