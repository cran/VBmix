\name{getElapsed}
\alias{getElapsed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getElapsed
}
\description{
get elapsed time since the moment defined by the parameter.
}
\usage{
getElapsed(stamp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stamp}{
timestamp of a specific moment, obtained with getTimestamp.
}
}
\value{
numeric vector with seconds and milliseconds of elapsed time.
}
\author{
Pierrick Bruneau
}
\seealso{
getTimestamp
}
\examples{
temp <- getTimestamp()
temp2 <- getElapsed(temp)
}

