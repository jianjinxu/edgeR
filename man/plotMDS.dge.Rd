\title{Multidimensional scaling plot of SAGE data}
\name{plotMDS.dge}
\alias{plotMDS.dge}
\description{
Plot the sample relations based on Multidimensional Scaling.
}
\usage{
plotMDS.dge(x, top=500, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1,2), ndim=max(dim.plot),...)
}
\arguments{
  \item{x}{any matrix or \code{DGEList} object.}
  \item{top}{number of top genes used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. If \code{x} has no column names, then defaults the index of the samples.}
  \item{col}{numeric or character vector of colors for the plotting characters.}
  \item{cex}{numeric vector of plot symbol expansions.}
  \item{dim.plot}{which two dimensions should be plotted, numeric vector of length two.}
  \item{ndim}{number of dimensions in which data is to be represented}
  \item{...}{any other arguments are passed to \code{plot}.}
}

\details{
This function is a variation on the usual multdimensional scaling (or principle coordinate) plot, in that a distance measure particularly appropriate for the digital gene expression (DGE) context is used. The distance between each pair of samples (columns) is the square root of the common dispersion for the top \code{top} genes which best distinguish that pair of samples. These top \code{top} genes are selected according to the tagwise dispersion of all the samples.

See \code{\link[graphics]{text}} for possible values for \code{col} and \code{cex}.
}

\value{A plot is created on the current graphics device.}

\author{Yunshun Chen and Gordon Smyth}


\examples{
# Simulate DGE data for 1000 genes(tags) and 6 samples.
# Samples are in two groups
# First 300 genes are differentially expressed in second group

x <- 10*runif(1000)
counts <- rnbinom(6000, size = 5, mu = x)
m <- matrix(counts, 1000, 6)
rownames(m) <- paste("Gene",1:1000)
m[1:300,4:6] <- m[1:300,4:6] + 10
plotMDS.dge(m)

# Indexes of samples are plotted.
plotMDS.dge(m,  col=c(rep("black",3), rep("red",3)) )
}

\keyword{hplot}
