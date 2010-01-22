\name{DGEList-class}
\alias{DGEList-class}
\docType{class}
\title{Digital Gene Expression data - class}

\description{
A simple list-based class for storing read counts from digital gene expression technologies and other important information for the analysis of DGE data.
}

\section{Slots/List Components}{
Objects of this class contain (at least) the following list components:
  
  \code{counts}: numeric matrix containing the read counts.
  
  \code{samples}: data.frame containing the library size and group labels.

}

\seealso{ \code{\link{DGEList}} }

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{DGEList} objects also have a \code{show} method.
}

\author{Mark Robinson}

\keyword{classes}

