\name{displayGraph}
\alias{displayGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
displayGraph
}
\description{
displays a curve (vect, measure), and associated deviations. Typically used to present experimental results.
}
\usage{
displayGraph(measure, dev, vect, xlab = "K", ylab = "measure", main = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{measure}{
y-axis for the curve.
}
  \item{dev}{
deviations for the y-axis measures.
}
  \item{vect}{
x-axis for the curve.
}
  \item{xlab}{
label for x-axis.
}
  \item{ylab}{
label for y-axis.
}
  \item{main}{
 main label for the plotting window.
}
}
\value{
a new plotting window displaying the curve.
}
\author{
Pierrick Bruneau
}
\examples{
displayGraph(rnorm(10, mean=4, sd=3), rnorm(10, mean=0, sd=0.5), 1:10)
}

