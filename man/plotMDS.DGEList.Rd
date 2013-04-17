\title{Multidimensional scaling plot of digital gene expression profiles}
\name{plotMDS.DGEList}
\alias{plotMDS.DGEList}
\description{
Calculate distances between RNA-seq or DGE libraries, then produce a multidimensional scaling plot.
Distances on the plot represent coefficient of variation of expression between samples
for the top genes that best distinguish the samples.
}
\usage{
\method{plotMDS}{DGEList}(x, top=500, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1,2),
     ndim=max(dim.plot), xlab=NULL, ylab=NULL, method="logFC", prior.count=2, gene.selection="pairwise", ...)
}
\arguments{
  \item{x}{an \code{DGEList} object.}
  \item{top}{number of top genes used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. If \code{x} has no column names, then defaults the index of the samples.}
  \item{col}{numeric or character vector of colors for the plotting characters. See \code{\link[graphics]{text}} for possible values.}
  \item{cex}{numeric vector of plot symbol expansions. See \code{\link[graphics]{text}} for possible values.}
  \item{dim.plot}{which two dimensions should be plotted, numeric vector of length two.}
  \item{ndim}{number of dimensions in which data is to be represented}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{method}{how to compute distances.  Possible values are "logFC" or \code{"bcv"}.}
  \item{prior.count}{average prior count to be added to observation to shrink the estimated log-fold-changes towards zero. Only used when \code{method="logFC"}.}
  \item{gene.selection}{character, \code{"pairwise"} to choose the top genes separately for each pairwise comparison between the samples or \code{"common"} to select the same genes for all comparisons. Only used when \code{method="logFC"}.}
  \item{...}{any other arguments are passed to \code{plot}.}
}

\details{
This function is a variation on the usual multdimensional scaling (or principle coordinate) plot, in that a distance measure particularly appropriate for the digital gene expression (DGE) context is used.
A set of top genes are chosen that have largest biological variation between the libraries
(those with largest tagwise dispersion treating all libraries as one group).
Then the distance between each pair of libraries (columns) is the biological coefficient of variation (square root of the common dispersion) between those two libraries alone, using
the top genes.

If \code{x} is a \code{DGEList}, then the library sizes and normalization factors found in the object are used.
If \code{x} is a matrix, then library sizes are computed from the column sums, but no other normalization is done.

The number \code{top} of top genes chosen for this exercise should roughly correspond to the number of differentially expressed genes with materially large fold-changes.
The default setting of 500 genes is widely effective and suitable for routine use, but a smaller value might be chosen for when the samples are distinguished by a specific focused molecular pathway.
Very large values (greater than 1000) are not usually so effective.

This function can be slow when there are many libraries.
}

\value{
A plot is created on the current graphics device.

An object of class \code{\link[limma:plotMDS]{MDS}} is invisibly returned.
}

\author{Yunshun Chen, Mark Robinson and Gordon Smyth}

\seealso{
\code{\link[limma]{plotMDS}}, \code{\link{cmdscale}}, \code{\link{as.dist}}
}

\examples{
# Simulate DGE data for 1000 genes(tags) and 6 samples.
# Samples are in two groups
# First 200 genes are differentially expressed in second group

ngenes <- 1000
nlib <- 6
counts <- matrix(rnbinom(ngenes*nlib, size=1/10, mu=20),ngenes,nlib)
rownames(counts) <- paste("Gene",1:ngenes)
group <- gl(2,3,labels=c("Grp1","Grp2"))
counts[1:200,group=="Grp2"] <- counts[1:200,group=="Grp2"] + 10
y <- DGEList(counts,group=group)
y <- calcNormFactors(y)

# without labels, indexes of samples are plotted.
col <- as.numeric(group)
mds <- plotMDS(y, top=200, col=col)

# or labels can be provided, here group indicators:
plotMDS(mds, col=col, labels=group)
}

\keyword{hplot}
