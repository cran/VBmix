\name{readLabelFile}
\alias{readLabelFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
readLabelFile
}
\description{
reads the vector of numeric labels contained in a binary file. Labels are associated to handwritten digits, thus in [0-9].
}
\usage{
readLabelFile(name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
path to the file.
}
}
\value{
vector of digit labels.
}
\references{
LeCun, Y., Bottou, L., Bengio, Y., and Haffner, P. (1998) _Gradient-based learning applied to document
recognition_, Proceedings of the IEEE, Volume 86, Number 11, Pages 2278-2324.

http://yann.lecun.com/exdb/mnist/, repository of the labels file. 
}
\author{
Pierrick Bruneau
}
\seealso{
readPixmapFile
}
\examples{
# use with path to actual train-...
# temp <- readLabelFile("data/train-labels-idx1-ubyte")
}

