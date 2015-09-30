\title{Plot exons of a differentially spliced gene}
\name{plotSpliceDGE}
\alias{plotSpliceDGE}
\description{
Plot the exon-level log-fold changes for a differentially spliced gene.
}
\usage{
plotSpliceDGE(lrt, geneid=NULL, rank=1L, FDR = 0.05)
}
\arguments{
  \item{lrt}{\code{GLMLRT} object produced by \code{diffSpliceDGE}.}
  \item{geneid}{character string, ID of the gene to plot.}
  \item{rank}{integer, if \code{geneid=NULL} then this ranked gene will be plotted.}
  \item{FDR}{numeric, mark exons with false discovery rate less than this cutoff.}
}

\details{
Plots interaction log-fold-change by exon for the specified gene.
}

\value{A plot is created on the current graphics device.}
\author{Yunshun Chen, Yifang Hu and Gordon Smyth}
\seealso{
\code{\link{diffSpliceDGE}.}
}
