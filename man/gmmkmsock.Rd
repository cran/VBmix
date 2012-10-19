\name{gmmkmsock}
\alias{gmmkmsock}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gmmkmsock
}
\description{
perform k-means specifically designed for a set of GMM (see references). At each iteration, sends information about current prototypes to a server via a socket connection (see references) for info about protocol.
}
\usage{
gmmkmsock(models, names, ngroups, rho = new.env(), host = "127.0.0.1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
list of GMM objects.
}
  \item{names}{
character vector with respective names of the GMM objects.
}
  \item{ngroups}{
(maximal) number of clusters.
}
  \item{rho}{
R environment object, used for calls to R functions within C code.
}
  \item{host}{
IP address of the server for the socket (port 1979).
}
}
\value{
a set of GMM prototypes, and inferred labels (i.e. associated to the input objects).
}
\references{
Bruneau, P. , Picarougne, F. and Gelgon, M. (2010) _Interactive unsupervised classication and visu-
alization for browsing an image collection_, Pattern Recognition, Volume 43, Number 2, Pages 485-493.
}
\author{
Pierrick Bruneau
}
\note{
gmmkmsock includes a socket client that sends formatted data to a server. Detailed information about this protocol may be found in the source package (inst/doc/old_manual.pdf). Simple standalone client and server are also provided (socket/socketclient.cpp and socketserver.cpp). These can be build by running make in the source folder.
}
\examples{
temp1 <- sample(1:200, 150)
temp2 <- list()
for(i in 1:length(temp1)) temp2 <- appendToList(temp2, imgmods[[temp1[i]]])
temp3 <- imgnames[temp1]
# next command may be executed only if a server is running on 127.0.0.1:1979.
# temp4 <- gmmkmsock(temp2, temp3, 5)
}

