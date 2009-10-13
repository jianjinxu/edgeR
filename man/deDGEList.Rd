\name{deDGEList-class}
\docType{class}
\alias{deDGEList-class}
\alias{show,deDGEList-method}
\title{differential expression of Digital Gene Expression data - class}

\description{
A simple list-based class for storing results of differential expression analysis for DGE data
}

\section{Slots/List Components}{
Objects of this class contain the following list components:
  \code{ps}: list containing estimates of \code{p} parameter.
  \code{r}: numeric vector of size parameter (\code{1/phi}) where \code{phi} is negative binomial dispersion).
  \code{pseudo}: numeric matrix with the pseudo-counts.
  \code{group}: vector giving the experimental group/condition.
  \code{M}: numeric scalar with the library size that pseudo counts are mapped to.
  
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{deDGEList} objects also have a \code{show} method.
}

\author{Mark Robinson, Davis McCarthy}

\keyword{classes}

