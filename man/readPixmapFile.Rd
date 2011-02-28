\name{readPixmapFile}
\alias{readPixmapFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
readPixmapFile
}
\description{
extracts a list of pixmap objects from the handwritten digits file format provided in references.
}
\usage{
readPixmapFile(name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
path to the file.
}
}
\value{
a list of pixmapGrey objects.
}
\references{
LeCun, Y., Bottou, L., Bengio, Y., and Haffner, P. (1998) _Gradient-based learning applied to document
recognition_, Proceedings of the IEEE, Volume 86, Number 11, Pages 2278-2324.

http://yann.lecun.com/exdb/mnist/, repository of the handwritten digits file. 
}
\author{
Pierrick Bruneau
}
\seealso{
pixmapGrey
}
\examples{
# use with path to actual train-... file
#temp <- readPixmapFile("data/train-images-idx3-ubyte")
}

