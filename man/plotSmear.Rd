\name{plotSmear}
\alias{plotSmear}
\title{
Plots log-Fold Change versus log-Concentration (or, M versus A) for Count Data
}
\description{
Both of these functions plot the log-fold change (i.e. the log of the ratio of expression levels for each tag between two experimential groups) against the log-concentration (i.e. the overall average expression level for each tag across the two groups). To represent counts that were low (e.g. zero in 1 library and non-zero in the other) in one of the two conditions, a 'smear' of points at low A value is presented in \code{plotSmear}.
}
\usage{
plotSmear(object, pair=NULL, de.tags=NULL, xlab="Average logCPM", ylab="logFC", pch=19,
     cex=0.2, smearWidth=0.5, panel.first=grid(), smooth.scatter=FALSE, lowess=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{\code{DGEList}, \code{DGEExact} or \code{DGELRT} object containing data to produce an MA-plot.}
  \item{pair}{pair of experimental conditions to plot (if \code{NULL}, the first two conditions are used)}
  \item{de.tags}{rownames for tags identified as being differentially expressed; use \code{exactTest} to identify DE genes}
  \item{xlab}{x-label of plot}
  \item{ylab}{y-label of plot}
  \item{pch}{scalar or vector giving the character(s) to be used in the plot; default value of \code{19} gives a round point.}
  \item{cex}{character expansion factor, numerical value giving the amount by which plotting text and symbols should be magnified relative to the default; default \code{cex=0.2} to make the plotted points smaller}
  \item{smearWidth}{width of the smear}
  \item{panel.first}{an expression to be evaluated after the plot axes are set up but before any plotting takes place; the default \code{grid()} draws a background grid to aid interpretation of the plot}
  \item{smooth.scatter}{logical, whether to produce a 'smooth scatter' plot using the \code{KernSmooth::smoothScatter} function or just a regular scatter plot; default is \code{FALSE}, i.e. produce a regular scatter plot}
  \item{lowess}{logical, indicating whether or not to add a lowess curve to the MA-plot to give an indication of any trend in the log-fold change with log-concentration}
  \item{\dots}{further arguments passed on to \code{plot}}
}

\value{A plot to the current device}

\details{
\code{plotSmear} is a more sophisticated and superior way to produce an 'MA plot'. \code{plotSmear} resolves the problem of plotting tags that have a total count of zero for one of the groups by adding the 'smear' of points at low A value. The points to be smeared are identified as being equal to the minimum estimated concentration in one of the two groups.  The smear is created by using random uniform numbers of width \code{smearWidth} to the left of the minimum A. \code{plotSmear} also allows easy highlighting of differentially expressed (DE) tags.
}

\author{Mark Robinson, Davis McCarthy}

\seealso{
\code{\link{maPlot}}
}

\examples{
y <- matrix(rnbinom(10000,mu=5,size=2),ncol=4)
d <- DGEList(counts=y, group=rep(1:2,each=2), lib.size=colSums(y))
rownames(d$counts) <- paste("tag",1:nrow(d$counts),sep=".")
d <- estimateCommonDisp(d)
plotSmear(d)

# find differential expression
de <- exactTest(d)

# highlighting the top 500 most DE tags
de.tags <- rownames(topTags(de, n=500)$table)
plotSmear(d, de.tags=de.tags)
}
