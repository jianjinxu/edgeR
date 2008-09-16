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
\tabular{ll}{
  \code{ps}:\tab list containing estimates of \code{p} parameter.\cr
  \code{r}:\tab numeric vector of size parameter (\code{1/phi}) where \code{phi} is negative binomial dispersion.\cr
  \code{pseudo}:\tab numeric matrix with the pseudo-counts.\cr
  \code{M}:\tab numeric scalar with the library size that pseudo counts are mapped to.\cr
  \code{exact}:\tab numeric vector of exact p-values (unadjusted).\cr
}
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{deDGEList} objects also have a \code{show} method.
}

\author{Mark Robinson}

\keyword{classes}

