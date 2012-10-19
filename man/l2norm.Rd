\name{l2norm}
\alias{l2norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
l2norm
}
\description{
computes Euclidian norm of vec.
}
\usage{
l2norm(vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
numeric vector.
}
}
\value{
norm value.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- l2norm(gmmpen[[2]]$mean[[1]])
}

