\name{DGEList}
\alias{DGEList}

\title{
DGEList Constructor
}

\description{
A function to create a \code{DGEList} object from a table of counts (rows=features, columns=samples), group indicator for each column, library size (optional) and a table of annotation (optional)
}
\usage{
DGEList(counts = matrix(0, 0, 0), lib.size = NULL, group = factor(), genes = NULL, remove.zeros = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{numeric matrix containing the read counts.}
  \item{lib.size}{numeric vector containing the total to normalize against for each sample (optional)}
  \item{group}{vector giving the experimental group/condition for each sample/library}
  \item{genes}{data frame containing annotation information for the tags/transcripts/genes for which we have count data (optional).}
  \item{remove.zeros}{whether to remove rows that have 0 total count; default is FALSE so as to retain all information in the dataset}
}
\details{
If no \code{lib.size} argument is passed to the constructor, the column totals are used.

The optional \code{genes} argument is meant to be an annotation data.frame, with rows matching those in the \code{counts} argument.
}
\value{a \code{DGEList} object}

\author{Mark Robinson}

\seealso{\code{\link{DGEList}}}

\examples{
y <- matrix(rnbinom(10000,mu=5,size=2),ncol=4)
d <- DGEList(counts=y, group=rep(1:2,each=2), lib.size=colSums(y))
}

