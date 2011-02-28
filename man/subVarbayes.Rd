\name{subVarbayes}
\alias{subVarbayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
subVarbayes
}
\description{
filters a variational posterior GMM, keeping only components with sufficient support.
}
\usage{
subVarbayes(model, thres = 2.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
variational posterior GMM.
}
  \item{thres}{
minimal support for component selection.
}
}
\value{
filtered variational posterior GMM.
}
\author{
Pierrick Bruneau
}
\seealso{
varbayes extractSimpleModel
}
\examples{
temp <- varbayes(irisdata, 20)
temp2 <- subVarbayes(temp)
}

