\name{plotSmear}
\Rdversion{1.1}
\alias{plotSmear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plots log fold change versus log concentration (or, M versus A) for count data
}
\description{
To represent counts that were low (e.g. zero in 1 library and non-zero in the other) in one of the two conditions, a 'smear' of points at low A value is presented.
}
\usage{
plotSmear(object, pair = NULL, de.tags=NULL, xlab = "logConc", ylab = "logFC", pch = 19, cex = 0.2, smearWidth = 0.5, panel.first=grid(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object to plot data from (uses $conc element)}
  \item{pair}{pair of experimental conditions to plot (if \code{NULL}, the first two conditions are used)}
  \item{de.tags}{rownames for tags identified as being differentially expressed; use \code{exactTest} to identify DE genes}
  \item{xlab}{x-label of plot}
  \item{ylab}{y-label of plot}
  \item{pch}{plotting character}
  \item{cex}{character expansion factor}
  \item{smearWidth}{width of the smear}
  \item{panel.first}{an expression to be evaluated after the plot axes are set up but before any plotting takes place; the default \code{grid()} draws a background grid to aid interpretation of the plot}
  \item{\dots}{further arguments passed on to \code{plot}}
}
\details{
The points to be smeared are identified as being equal to the minimum estimated concentration in one of the two groups.  The smear is created by using random uniform numbers of width \code{smearWidth} to the left of the minimum A.
}
\value{a plot to the current device}
\author{Mark Robinson}

\seealso{
\code{\link{maPlot}}
}
\examples{
y <- matrix(rnbinom(10000,mu=5,size=2),ncol=4)
d <- DGEList(counts=y, group=rep(1:2,each=2), lib.size=colSums(y))
d <- estimateCommonDisp(d)
plotSmear(d)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
