\name{displayScatter}
\alias{displayScatter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
displayScatter
}
\description{
general plotting function for data sets (matrix of row-elements), optionally associated to labels and a GMM. Labels influence the color and symbols of plotted data points. Gaussian envelopes of the components in the GMM are drawn. NB: data set and GMM arguments cannot be both NULL.
}
\usage{
displayScatter(data = NULL, model = NULL, labels = NULL, datasizes = NULL, 
compcolors = NULL, complabels = NULL, compstrokes = "solid", space = 1:2, 
xlim = NULL, ylim = NULL, main = "", xlab = "", ylab = "", smooth = FALSE, 
alphacol = 0.8, alphanocol = 0.5, cex.lab = 1, lwd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
matrix of row-elements. If NULL, the GMM is plotted alone.
}
  \item{model}{
GMM object.
}
  \item{labels}{
vector of numeric labels. May alternatively be present as a member of model, labels.
}
  \item{datasizes}{
vector of integer magnification factors for data symbols. If length=1, same coefficient applies to all points.
}
  \item{compcolors}{
vector of integer color indexes. These indexes are internally associated to one color among a set of appropriately chosen ones. If length=1, all GMM components are colored the same way. If length=k, each component is associated to its own color index. This k-length vector may contain NA values: associated components will be white-colored.
}
  \item{complabels}{
character vector containing text strings to be printed over Gaussian envelopes.
}
  \item{compstrokes}{
this character vector may be used to specify non default strokes for envelopes.
}
  \item{space}{
this function prints a 2D scatterplot. If data and model have higher dimensionality, this argument specifies the axes to be printed.
}
  \item{xlim}{
bounds for the first variable. If NULL, will be inferred from available data.
}
  \item{ylim}{
bounds for the second variable. If NULL, will be inferred from available data.
}
  \item{main}{
main label for the plotting window.
}
  \item{xlab}{
label for the x-axis.
}
  \item{ylab}{
label for the y-axis.
}
  \item{smooth}{
if TRUE, display the response to a kernel density function, instead of symbols for data elements.
}
  \item{alphacol}{
alpha blending parameter when a component is non-white colored.
}
  \item{alphanocol}{
alpha blending parameter when a component is white colored.
}
  \item{cex.lab}{
magnification factor for all text in the plotting window.
}
  \item{lwd}{
width of the stroke used for data symbols.
}
}
\value{
a new plotting window displaying the data set and associated model.
}
\author{
Pierrick Bruneau
}
\seealso{
plotGmm
}
\examples{
displayScatter(irisdata, NULL, irislabels)
}

