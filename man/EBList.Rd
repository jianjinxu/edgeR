\name{EBList-class}
\docType{class}
\alias{EBList-class}
\alias{show,EBList-method}

\title{differential expression of Digital Gene Expression data - class}

\description{
A simple list-based class for storing results of the approximate empirical Bayes rule parameters
}

\section{Slots/List Components}{
Objects of this class contain the following list components:
\tabular{ll}{
  \code{sigma2.0.est}:\tab numeric scale \code{sigma_0^2} estimate.\cr
  \code{alpha}:\tab numeric scalar alpha estimate.\cr
  \code{score}:\tab numeric scalar (likelihood) score.\cr
  \code{inos}:\tab numeric vector containing the (likelihood) information.\cr
  \code{exact}:\tab list containing info with respect to quantile adjustment.\cr
}
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{EBList} objects also have a \code{show} method.
}

\author{Mark Robinson}

\keyword{classes}

