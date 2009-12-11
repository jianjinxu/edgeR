\name{DGEList-class}
\docType{class}
\alias{DGEList-class}
\alias{DGEList}
\alias{show,DGEList-method}
\title{Digital Gene Expression data - class}

\description{
A simple list-based class for storing read counts from digital gene expression technologies and other important information for the analysis of DGE data.
}

\section{Slots/List Components}{
Objects of this class contain (at least) the following list components:
  
  \code{counts}: numeric matrix containing the read counts.
  
  \code{lib.size}: numeric vector containing the total number of reads for each library (i.e. each column of \code{counts}).
  
  \code{group}: vector giving the experimental group/condition for each sample/library.
  
  \code{genes}: data frame (optional) containing annotation information for the tags/transcripts/genes for which we have count data.
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{DGEList} objects also have a \code{show} method.
}

\author{Mark Robinson}

\keyword{classes}

