\name{DGEList-class}
\docType{class}
\alias{DGEList-class}
\alias{show,DGEList-method}
\title{Digital Gene Expression data - class}

\description{
A simple list-based class for storing read counts from digital gene expression technologies.
}

\section{Slots/List Components}{
Objects of this class contain the following list components:
\tabular{ll}{
  \code{data}:\tab numeric matrix containing the read counts.\cr
  \code{lib.size}:\tab numeric vector containing the total number of reads for each library (column of \code{code}).\cr
  \code{group}:\tab vector giving the experimental group/condiction.\cr
}
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{DGEList} objects also have a \code{show} method.
}

\author{Mark Robinson}

\keyword{classes}

