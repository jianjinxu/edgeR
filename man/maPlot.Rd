\name{maPlot}
\Rdversion{1.1}
\alias{maPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots Log-Fold Change versus Log-Concentration (or, M versus A) for Count Data
}
\description{
To represent counts that were low (e.g. zero in 1 library and non-zero in the other) in one of the two conditions, a 'smear' of points at low A value is presented.
}
\usage{
maPlot(x, y, normalize=FALSE, smearWidth = 1, col = NULL, allCol =
"black", lowCol = "orange", deCol="red", de.tags=NULL, smooth.scatter=
FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of counts or concentrations (group 1)}
  \item{y}{vector of counts or concentrations (group 2)}
  \item{normalize}{logical, whether to divide \code{x} and \code{y} vectors by their sum}
  \item{smearWidth}{scalar, width of the smear}
  \item{col}{vector of colours for the points (if \code{NULL}, uses \code{allCol} and \code{lowCol})}
  \item{allCol}{colour of the non-smeared points}
  \item{lowCol}{colour of the smeared points}
  \item{deCol}{colour of the DE (differentially expressed) points}
  \item{de.tags}{indices for tags identified as being differentially
  expressed; use \code{exactTest} to identify DE genes}
  \item{smooth.scatter}{logical, whether to produce a 'smooth scatter'
  plot using the KernSmooth::smoothScatter function or just a regular
  scatter plot; default is \code{FALSE}, i.e. produce a regular
  scatter plot}
  \item{\dots}{further arguments passed on to \code{plot}}
}
\details{
The points to be smeared are identified as being equal to the minimum in one of the two groups.  The smear is created by using random uniform numbers of width \code{smearWidth} to the left of the minimum A value.
}
\value{a plot to the current device}

\seealso{
\code{\link{plotSmear}}
}

\author{Mark Robinson}


\examples{
y <- matrix(rnbinom(10000,mu=5,size=2),ncol=4)
maPlot(y[,1], y[,2])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
