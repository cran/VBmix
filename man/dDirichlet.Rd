\name{dDirichlet}
\alias{dDirichlet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
dDirichlet
}
\description{
get density of a sample w.r.t Dirichlet distribution (3D only).
}
\usage{
dDirichlet(alpha = 0.1, x1, x2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
alpha parameter of the distribution (i.e. alpha repeated 3 times).
}
  \item{x1}{
1st dimension of the sample.
}
  \item{x2}{
2nd dimension of the sample.
}
}
\value{
density value.
}
\author{
Pierrick Bruneau
}
\seealso{
rDirichlet
}
\examples{
temp <- dDirichlet(x1=0.4, x2=0.2)
# 3rd dimension is 1-x1-x2 = 0.2
}

