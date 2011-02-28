\name{mixKnn}
\alias{mixKnn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mixKnn
}
\description{
performs k-nearest neighbors over a collection of GMM. It uses jsmc to compute distances. Each elements in data is classified against all the others, and inferred class is compared to the true one (leave-one-out).
}
\usage{
mixKnn(data, labels, n = 2, KLparam = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
list of GMM.
}
  \item{labels}{
vector of numeric labels associated to data.
}
  \item{n}{
k of the algorithm.
}
  \item{KLparam}{
number of samples for jsmc.
}
}
\value{
classification error ratio in [0,1].
}
\author{
Pierrick Bruneau
}
\seealso{
mergeClassif constrClassif sampleClassif
}
\examples{
temp1 <- sample(1:1243, 150)
temp2 <- list()
for(i in temp1) temp2 <- appendToList(temp2, imgmods[[i]])
temp3 <- imglabels[temp1]
# de-activated because this process is very long...
#temp4 <- mixKnn(temp2, temp3)
}
