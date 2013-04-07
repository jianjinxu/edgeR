\name{DGEList-class}
\alias{DGEList-class}
\docType{class}
\title{Digital Gene Expression data - class}

\description{
A list-based S4 class for storing read counts and associated information from digital gene expression or sequencing technologies.
}

\section{List Components}{
Objects of this class contain the following essential list components:
\tabular{ll}{
  \code{counts } \tab numeric matrix of read counts, one row for each gene (or genomic feature) and one column for each sample (library).\cr
  
  \code{samples } \tab data.frame with a row for each sample and columns \code{group}, \code{lib.size} and \code{norm.factors} containing the group labels, library sizes and normalization factors.
  Other columns can be optionally added to give more detailed sample information.
}
Optional components include:
\tabular{ll}{
  \code{genes } \tab data.frame giving annotation information for each gene.  Same number of rows as \code{counts}.\cr

  \code{AveLogCPM } \tab numeric vector giving average log2 counts per million for each gene.\cr
  
  \code{common.dispersion } \tab numeric scalar giving the overall dispersion estimate.\cr
  
  \code{trended.dispersion } \tab numeric vector giving trended dispersion estimates for each gene.\cr
    
  \code{tagwise.dispersion } \tab numeric vector giving tagwise dispersion estimates for each gene.\cr
  
  \code{offset } \tab numeric matrix of same size as \code{counts} giving offsets for use in log-linear models.
}
}

\seealso{ \code{\link{DGEList}} }

\section{Methods}{
This class inherits directly from class \code{list}, so \code{DGEList} objects can be manipulated as if they were ordinary lists.
However they can also be treated as if they were matrices for the purposes of subsetting.

The dimensions, row names and column names of a \code{DGEList} object are defined by those of \code{counts}, see \code{\link{dim.DGEList}} or \code{\link{dimnames.DGEList}}.

\code{DGEList} objects can be subsetted, see \code{\link{subsetting}}.

\code{DGEList} objects also have a \code{show} method so that printing produces a compact summary of their contents.
}

\author{edgeR team. First created by Mark Robinson.}

\keyword{classes}

