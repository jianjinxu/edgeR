\name{dim}
\alias{dim.DGEList}
\alias{dim.deDGEList}
\alias{dim.TopTags}
\alias{length.DGEList}
\alias{length.deDGEList}
\alias{length.TopTags}
\title{Retrieve the Dimensions of a DGEList, deDGEList or TopTags Object}
\description{
Retrieve the number of rows (transcripts) and columns (libraries) for an DGEList, deDGEList or TopTags Object.
}
\usage{
\method{dim}{DGEList}(x)
\method{length}{DGEList}(x)
}
\arguments{
  \item{x}{an object of class \code{DGEList}, \code{deDGEList} or \code{TopTags}}
}
\details{
Digital gene expression data objects share many analogies with ordinary matrices in which the rows correspond to transcripts or genes and the columns to arrays.
These methods allow one to extract the size of microarray data objects in the same way that one would do for ordinary matrices.

A consequence is that row and column commands \code{nrow(x)}, \code{ncol(x)} and so on also work.
}
\value{
Numeric vector of length 2.
The first element is the number of rows (genes) and the second is the number of columns (arrays).
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{dim}} in the base package.
  
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}
\examples{
M <- A <- matrix(11:14,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("A1","A2")
MA <- new("MAList",list(M=M,A=A))
dim(M)
ncol(M)
nrow(M)
length(M)
}
\keyword{array}
