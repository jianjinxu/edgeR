\name{DGEExact-class}
\docType{class}
\alias{DGEExact-class}
\alias{show,DGEExact-method}
\title{differential expression of Digital Gene Expression data - class}

\description{
A simple list-based class for storing results of differential expression analysis for DGE data
}

\section{Slots/List Components}{
Objects of this class contain the following list components:

  \code{table}: data frame containing the log-concentration (i.e. expression level), the log-fold change in expression between the two groups/conditions and the exact p-value for differential expression, for each tag.
  
  \code{comparison}: vector giving the two experimental groups/conditions being compared.
  
  \code{genes}: a data frame containing information about each transcript (can be \code{NULL}).
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{DGEExact} objects also have a \code{show} method.
}

\author{Mark Robinson, Davis McCarthy}

\keyword{classes}

