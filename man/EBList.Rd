\name{EBList-class}
\docType{class}
\alias{EBList-class}
\alias{show,EBList-method}

\title{Differential Expression of Digital Gene Expression data - class}

\description{
A simple list-based class for storing results of the approximate empirical Bayes rule parameters
}

\section{Slots/List Components}{
Objects of this class contain the following list components:
  \code{sigma2.0.est}: numeric scale \code{sigma_0^2} estimate.
  \code{alpha}: numeric scalar alpha estimate.
  \code{scores}: numeric scalar (likelihood) score.
  \code{infos}: numeric vector containing the (likelihood) information for each tag.
  \code{quantileAdjusted}: list from output of \code{quantileAdjust}).
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{EBList} objects also have a \code{show} method.
}

\author{Mark Robinson, Davis McCarthy}

\keyword{classes}

