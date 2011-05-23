\title{Multidimensional scaling plot of digital gene expression profiles}
\name{plotMDS.DGEList}
\alias{plotMDS.DGEList}
\description{
Calculate distances between RNA-seq or DGE libraries, then produce a multidimensional scaling plot.
Distances on the plot represent coefficient of variation of expression between samples
for the top genes that best distinguish the samples.
}
\usage{
\method{plotMDS}{DGEList}(x, top=500, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1, 2), ndim=max(dim.plot), xlab=paste("Dimension",dim.plot[1]), ylab=paste("Dimension",dim.plot[2]), ...)
}
\arguments{
  \item{x}{any matrix or \code{DGEList} object.}
  \item{top}{number of top genes used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. If \code{x} has no column names, then defaults the index of the samples.}
  \item{col}{numeric or character vector of colors for the plotting characters.}
  \item{cex}{numeric vector of plot symbol expansions.}
  \item{dim.plot}{which two dimensions should be plotted, numeric vector of length two.}
  \item{ndim}{number of dimensions in which data is to be represented}
  \item{xlab}{title for the x-axis}
  \item{ylab}{title for the y-axis}
  \item{...}{any other arguments are passed to \code{plot}.}
}

\details{
This function is a variation on the usual multdimensional scaling (or principle coordinate) plot, in that a distance measure particularly appropriate for the digital gene expression (DGE) context is used.
A set of top genes are chosen that have largest biological variation between the libraries
(those with largest tagwise dispersion treating all libraries as one group).
Then the distance between each pair of libraries (columns) is the biological coefficient of variation (square root of the common dispersion) between those two libraries alone, using
the top genes.
See \code{\link[graphics]{text}} for possible values for \code{col} and \code{cex}.

This function can be slow when there are many libraries.
}

\value{
A plot is created on the current graphics device.

An object of class \code{"MDS"} is invisibly returned.
This is a list containing the following components:
\item{distance.matrix}{numeric matrix of pairwise distances between columns of \code{x}}
\item{cmdscale.out}{output from the function \code{cmdscale} given the distance matrix}
\item{dim.plot}{dimensions plotted}
\item{x}{x-xordinates of plotted points}
\item{y}{y-cordinates of plotted points}
}

\author{Yunshun Chen and Gordon Smyth}

\seealso{
\code{\link{cmdscale}}, \code{\link{as.dist}}, \code{\link[limma]{plotMDS}}
}

\examples{
# Simulate DGE data for 1000 genes(tags) and 6 samples.
# Samples are in two groups
# First 300 genes are differentially expressed in second group

x <- 10*runif(1000)
counts <- rnbinom(6000, size = 5, mu = x)
m <- matrix(counts, 1000, 6)
rownames(m) <- paste("Gene",1:1000)
m[1:300,4:6] <- m[1:300,4:6] + 10
m <- DGEList(counts=m)
# without labels, indexes of samples are plotted.
mds <- plotMDS(m,  col=c(rep("black",3), rep("red",3)) )
# or labels can be provided, here group indicators:
plotMDS(mds,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))
}

\keyword{hplot}
