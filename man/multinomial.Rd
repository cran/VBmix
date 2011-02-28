\name{multinomial}
\alias{multinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
multinomial
}
\description{
samples from a k-multinomial.
}
\usage{
multinomial(weights, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weights}{
numeric vector with the weights of the multinomial. Sum to 1.
}
  \item{k}{
size of the weight vector.
}
}
\value{
an integer value in [1,k], coded as a 1-of-k variable (see reference).
}
\references{
Bishop, C. M. (2006) _Pattern Recognition and Machine Learning_ Chap. 10, Springer.
}
\author{
Pierrick Bruneau
}
\examples{
weights <- c(0.3, 0.5, 0.2)
multinomial(weights, 3)
#[1] 0 1 0
}

