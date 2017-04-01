\name{subMppca}
\alias{subMppca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
subMppca
}
\description{
removes unused components and factor columns from model.
}
\usage{
subMppca(model, prune = FALSE, thres = 2.001, quick = FALSE, noxmean = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
MPPCA model to be shrinked.
}
  \item{prune}{
if TRUE, void factor columns are removed.
}
  \item{thres}{
threshold for component selection. A components is selected iif alpha > thres.
}
  \item{quick}{
influences method for void factor columns detection. If FALSE, a KL-based criterion is employed (more accurate). If TRUE, column norms are used (useful for very high dimensional data sets).
}
  \item{noxmean}{
should always be set to TRUE.
}
}
\value{
shrinked MPPCA model.
}
\author{
Pierrick Bruneau
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
mppca newMppca
}
\examples{
# use a subsample of pendat, for runtime (packaging) needs.
temp <- mppca(pendat[sample(1:2000,150),], 15, qmax=8, maxit=20)
temp2 <- subMppca(temp, prune=TRUE, quick=TRUE)
}

