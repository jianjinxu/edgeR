\name{DGEList-class}
\docType{class}
\alias{DGEList-class}
\alias{DGEList}
\alias{show,DGEList-method}
\title{Digital Gene Expression data - class}

\description{
A simple list-based class for storing read counts from digital gene expression technologies.
}

\section{Slots/List Components}{
Objects of this class contain the following list components:
  \code{counts}: numeric matrix containing the read counts.
  \code{lib.size}: numeric vector containing the total number of reads for each library (column of \code{code}).
  \code{group}: vector giving the experimental group/condition.
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{DGEList} objects also have a \code{show} method.
}

\author{Mark Robinson}

\keyword{classes}

