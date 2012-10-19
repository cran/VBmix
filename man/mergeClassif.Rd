\name{mergeClassif}
\alias{mergeClassif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mergeClassif
}
\description{
performs task analogous to mixKnn (i.e. leave-one-out classification), but uses synthetic representatives to infer labels, instead of k-NN. Each representative is obtained by concatenating all GMM (i.e. elements) of a specific label value, and applying vbcomp on this redundant mixture.
}
\usage{
mergeClassif(data, labels, KLparam = 500, rho = new.env())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
list of GMM.
}
  \item{labels}{
vector of numeric labels associated to data.
}
  \item{KLparam}{
number of samples for jsmc.
}
  \item{rho}{
R environment object. Used to issue R commands within the C routine.
}
}
\value{
classification error ratio in [0,1].
}
\author{
Pierrick Bruneau
}
\seealso{
mixKnn vbcomp
}
\examples{
temp1 <- sample(1:200, 150)
temp2 <- list()
for(i in temp1) temp2 <- appendToList(temp2, imgmods[[i]])
temp3 <- imglabels[temp1]
# de-activated because this process is very long...
#temp4 <- mergeClassif(temp2, temp3)
}

