\name{dimnames}
\alias{dimnames.DGEList}
\alias{dimnames<-.DGEList}
\alias{dimnames.DGEExact}
\alias{dimnames<-.DGEExact}
\alias{dimnames.DGELRT}
\alias{dimnames<-.DGELRT}
\alias{dimnames.TopTags}
\alias{dimnames<-.TopTags}
\title{Retrieve the Dimension Names of a DGEList Object}
\description{
Retrieve the dimension names of a digital gene expression data object.
}
\usage{
\method{dimnames}{DGEList}(x)
\method{dimnames}{DGEList}(x) <- value
}
\arguments{
  \item{x}{an object of class \code{DGEList}, \code{DGEExact}, \code{DGELRT} or \code{TopTags}}
  \item{value}{a possible value for \code{dimnames(x)}: see \code{\link{dimnames}}}
}
\details{
The dimension names of a microarray object are the same as those of the most important matrix component of that object.

A consequence is that \code{rownames} and \code{colnames} will work as expected.
}
\value{
Either \code{NULL} or a list of length 2.
If a list, its components are either \code{NULL} or a character vector the length of the appropriate dimension of \code{x}.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{dimnames}} in the base package.
  
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}
\keyword{array}
