\name{sort_index}
\alias{sort_index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sort_index
}
\description{
returns indexes associated to the sorted values of the parameter vector.
}
\usage{
sort_index(vec, order = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
vector to be sorted.
}
  \item{order}{
if 0, ascending order, if 1, descending order.
}
}
\value{
indexes associated to the sorted input vector.
}
\author{
Pierrick Bruneau
}
\examples{
temp <- rnorm(10)
temp2 <- sort_index(temp)
}

