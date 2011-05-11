\name{subsetting}
\alias{subsetting}
\alias{[.DGEList}
\alias{[.DGEExact}
\alias{[.DGELRT}
\alias{[.DGEGLM}
\title{Subset DGEList, DGEGLM, DGEExact and DGELRT Objects}
\description{
Extract a subset of a \code{DGEList}, \code{DGEGLM}, \code{DGEExact} or \code{DGELRT} object.
}
\usage{
\method{[}{DGEList}(object, i, j, \ldots)
\method{[}{DGEGLM}(object, i, j, \ldots)
\method{[}{DGEExact}(object, i, j, \ldots)
\method{[}{DGELRT}(object, i, j, \ldots)
}
\arguments{
 \item{object}{object of class \code{DGEList}, \code{DGEGLM}, \code{DGEExact} or \code{DGELRT}, respectively}
  \item{i,j}{elements to extract. \code{i} subsets the tags or genes while \code{j} subsets the libraries. Note, columns of \code{DGEGLM}, \code{DGEExact} and \code{DGELRT} objects cannot be subsetted.}
  \item{\ldots}{not used}
}
\details{
\code{i,j} may take any values acceptable for the matrix components of \code{object} of class \code{DGEList}.
See the \link{Extract} help entry for more details on subsetting matrices. For \code{DGEGLM}, \code{DGEExact} and \code{DGELRT} objects, only rows (i.e. \code{i}) may be subsetted.
}
\value{
An object of class \code{DGEList}, \code{DGEGLM}, \code{DGEExact} or \code{DGELRT} as appropriate, holding data from the specified subset of tags/genes and libraries.
}
\author{Davis McCarthy, Gordon Smyth}
\seealso{
  \code{\link{Extract}} in the base package.
}
\examples{
d <- matrix(rnbinom(16,size=1,mu=10),4,4)
rownames(d) <- c("a","b","c","d")
colnames(d) <- c("A1","A2","B1","B2")
d <- DGEList(counts=d,group=factor(c("A","A","B","B")))
d[1:2,]
d[1:2,2]
d[,2]
d <- estimateCommonDisp(d)
results <- exactTest(d)
results[1:2,]
# NB: cannot subset columns for DGEExact objects
}
\keyword{manip}
