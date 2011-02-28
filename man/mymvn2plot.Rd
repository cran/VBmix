\name{mymvn2plot}
\alias{mymvn2plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mymvn2plot
}
\description{
displays mvn envelopes. For internal usage in graphical functions.
}
\usage{
mymvn2plot(w, mu, sigma, k = 15, alone = FALSE, col = NA, alphacol = 0.8, alphanocol = 0.5, lty = "solid")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
weight of the component.
}
  \item{mu}{
mean of the component.
}
  \item{sigma}{
covariance matrix of the component.
}
  \item{k}{
resolution used for drawaing the elliptic envelope.
}
  \item{alone}{
if TRUE, the component is to be plotted alone in its own window.
}
  \item{col}{
optional background color for the component.
}
  \item{alphacol}{
alpha coefficient for a component with a color.
}
  \item{alphanocol}{
alpha coefficient for a component with no color.
}
  \item{lty}{
line type for the ellipsis.
}
}

