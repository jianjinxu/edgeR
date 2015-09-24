\name{cpm}
\alias{cpm}
\alias{cpm.DGEList}
\alias{cpm.default}
\alias{rpkm}
\alias{rpkm.DGEList}
\alias{rpkm.default}

\title{Counts per Million or Reads per Kilobase per Million}

\description{Computes counts per million (CPM) or reads per kilobase per million (RPKM) values.}

\usage{
\method{cpm}{DGEList}(x, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25, \dots)
\method{cpm}{default}(x, lib.size=NULL, log=FALSE, prior.count=0.25, \dots)
\method{rpkm}{DGEList}(x, gene.length=NULL, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25, \dots)
\method{rpkm}{default}(x, gene.length, lib.size=NULL, log=FALSE, prior.count=0.25, \dots)
}
\arguments{ 
\item{x}{matrix of counts or a \code{DGEList} object}
\item{normalized.lib.sizes}{logical, use normalized library sizes?}
\item{lib.size}{library size, defaults to \code{colSums(x)}.}
\item{log}{logical, if \code{TRUE} then \code{log2} values are returned.}
\item{prior.count}{average count to be added to each observation to avoid taking log of zero. Used only if \code{log=TRUE}.}
\item{gene.length}{vector of length \code{nrow(x)} giving gene length in bases, or the name of the column \code{x$genes} containing the gene lengths.}
\item{\dots}{other arguments that are not currently used.}
}

\value{numeric matrix of CPM or RPKM values.}

\details{
CPM or RPKM values are useful descriptive measures for the expression level of a gene.
By default, the normalized library sizes are used in the computation for \code{DGEList} objects but simple column sums for matrices.

If log-values are computed, then a small count, given by \code{prior.count} but scaled to be proportional to the library size, is added to \code{x} to avoid taking the log of zero.

The \code{rpkm} method for \code{DGEList} objects will try to find the gene lengths in a column of \code{x$genes} called \code{Length} or \code{length}.
Failing that, it will look for any column name containing \code{"length"} in any capitalization.
}

\note{
\code{aveLogCPM(x)}, \code{rowMeans(cpm(x,log=TRUE))} and \code{log2(rowMeans(cpm(x))} all give slightly different results.
}

\author{Davis McCarthy, Gordon Smyth}

\examples{
y <- matrix(rnbinom(20,size=1,mu=10),5,4)
cpm(y)

d <- DGEList(counts=y, lib.size=1001:1004)
cpm(d)
cpm(d,log=TRUE)

d$genes$Length <- c(1000,2000,500,1500,3000)
rpkm(d)
}

\seealso{
\code{\link{aveLogCPM}}
}
