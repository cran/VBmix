\name{randomGmm}
\alias{randomGmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
randomGmm
}
\description{
sample randomly a GMM. Number of components is sampled from a Poisson law, means uniformly from [-domain, domain], and covariance matrices using covgen function.
}
\usage{
randomGmm(domain = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{domain}{
determines the domain from which means are sampled.
}
}
\value{
randomly sampled GMM.
}
\author{
Pierrick Bruneau
}
\seealso{
covgen newGmm
}
\examples{
temp <- randomGmm()
}

