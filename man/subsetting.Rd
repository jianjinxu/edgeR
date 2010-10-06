\name{subsetting}
\alias{subsetting}
\alias{[.DGEList}
\alias{[.DGEExact}
\title{Subset DGEList and DGEExact Objects}
\description{
Extract a subset of a \code{DGEList} or a \code{DGEExact} object.
}
\usage{
\method{[}{DGEList}(object, i, j, \ldots)
\method{[}{DGEExact}(object, i, j, \ldots)
}
\arguments{
  \item{object}{object of class \code{DGEList} of \code{DGEExact}}
  \item{i,j}{elements to extract. \code{i} subsets the tags or genes while \code{j} subsets the libraries. \code{j} is not used for \code{DGEExact} objects}
  \item{\ldots}{not used}
}
\details{
If \code{object} is a \code{DGEList}, then \code{i,j} may take any values acceptable for the matrix components of \code{object}. If \code{object} is of class \code{DGEExact} object, then \code{i} may take any values acceptable for the rows of the table component of \code{object} and \code{j} should not be used (i.e.~left missing).
See the \link{Extract} help entry for more details on subsetting matrices.
}
\value{
An object of class \code{DGEList} or \code{DGEExact}, as appropriate, holding data from the specified subset of tags/genes and libraries.
}
\author{Davis McCarthy, Gordon Smyth}
\seealso{
  \code{\link{Extract}} in the base package.
}
\examples{
d <- matrix(rnbinom(8,size=1,mu=10),4,2)
rownames(d) <- c("a","b","c","d")
colnames(d) <- c("A","B")
d <- new("DGEList",list(counts=d,group=factor(c("A","B"))))
d[1:2,]
d[1:2,2]
d[,2]
}
\keyword{manip}
