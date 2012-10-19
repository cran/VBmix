\name{displaySVM}
\alias{displaySVM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
displaySVM
}
\description{
displays the colored decision regions of a SVM model. Data symbols are also optionally displayed. Data and model should be 2D.
}
\usage{
displaySVM(svm.model, dataframe, displayPoints = TRUE, subset = NULL, steps = 100, alpha = 0.4, lwd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{svm.model}{
a SVM model, as returned by svm (e1071 library)
}
  \item{dataframe}{
data.frame object, containing row-elements, and associated labels in the last variable.
}
  \item{displayPoints}{
if FALSE, only decision regions are displayed.
}
  \item{subset}{
vector of indexes of a data subset to be displayed. If NULL, all points are displayed.
}
  \item{steps}{
influences the resolution of the decision regions. Low values will provoke aliasing, high values are slower to be displayed.
}
  \item{alpha}{
alpha blending parameter between decision regions and data symbols.
}
  \item{lwd}{
magnification factor for the stroke width used to plot symbols.
}
}
\value{
a new plotting window displaying SVM decision regions.
}
\author{
Pierrick Bruneau
}
\seealso{
svm
}
\examples{
# extract 2 first variables and build data.frame
temp <- buildFrame(irisdata, irislabels)
iris.model <- svm(labels ~ ., data=temp, cost=100, gamma=1)
displaySVM(iris.model, temp)
}

