\name{extractSimpleModel}
\alias{extractSimpleModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
extractSimpleModel
}
\description{
extracts a GMM from a posterior variational distribution. Only relevant components (i.e. associated to a significant population) are extracted.
}
\usage{
extractSimpleModel(model = model, labels = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
variational posterior.
}
  \item{labels}{
boolean indicating wether to extract a label vector. If TRUE, model, a list object, should also contain a data attribute, used to build label vector.
}
}
\value{
GMM object.
}
\author{
Pierrick Bruneau
}
\seealso{
varbayes subVarbayes
}
\examples{
temp <- varbayes(irisdata, 20)
temp2 <- extractSimpleModel(temp)
}

