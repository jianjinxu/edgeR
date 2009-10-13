\name{subsetting}
\alias{subsetting}
\alias{[.DGEList}
\title{Subset DGEist Objects}
\description{
Extract a subset of an \code{DGEList} object.
}
\usage{
\method{[}{DGEList}(object, i, j, \ldots)
}
\arguments{
  \item{object}{object of class \code{DGEList}}
  \item{i,j}{elements to extract. \code{i} subsets the tags or genes while \code{j} subsets the libraries}
  \item{\ldots}{not used}
}
\details{
\code{i,j} may take any values acceptable for the matrix components of \code{object}.
See the \link{Extract} help entry for more details on subsetting matrices.
}
\value{
An object of class \code{DGEList} holding data from the specified subset of tags/genes and libraries.
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
