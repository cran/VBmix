\name{appendToList}
\alias{appendToList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
appendToList
}
\description{
appends 1 list object to another.
}
\usage{
appendToList(lst, obj, appendList = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{
list object to which we append an object.
}
  \item{obj}{
object to append.
}
  \item{appendList}{
if TRUE, obj should be a list object, which elements are appended. if FALSE, obj is simply added to lst.
}
}
\value{
list object with obj appended to lst.
}
\author{
Pierrick Bruneau
}
\seealso{
appendToGmm appendToMppca
}
\examples{
temp <- list()
temp <- appendToList(temp, pcapen[[1]]$mumean, appendList=TRUE)
temp <- appendToList(temp, pcapen[[2]]$mumean, appendList=TRUE)
}

